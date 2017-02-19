from grako.exceptions import GrakoException

from gnomic.models import *
from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.utils import genotype_to_text, genotype_to_string, change_to_string

__all__ = (
    'DEFAULT_TYPES',
    'DEFAULT_ORGANISMS',
    'Genotype',
    'Mutation',
    'FeatureTree',
    'Fusion',
    'Plasmid',
    'Feature',
    'Organism',
    'Accession',
    'Type',
    'encoder'
)

DEFAULT_TYPES = (
    Type('promoter', ('P',)),
    Type('terminator', ('T',)),
    Type('gene')
)

DEFAULT_ORGANISMS = (
    Organism('Escherichia coli', aliases=['Ec', 'E.coli']),
    Organism('Saccaromices cerevisiae', aliases=['Sc', 'S.cerevisiae']),
)


class Match:
    NONE = 0
    PARTIAL = 1
    COMPLETE = 2


class AmbiguityError(Exception):
    pass


class Genotype(object):
    """


    .. attribute:: raw

        A :class:`tuple` containing the changes introduced in this genotype.

        A change can be one of:

        - A :class:`Mutation`
        - A :class:`Plasmid` (non-integrated plasmid)
        - A :class:`Feature` of type 'phene' with a variant (a "phenotype")

    """

    FUSION_MATCH_WHOLE = 'match-whole-only'
    FUSION_UPDATE_ON_CHANGE = 'update-on-change'
    FUSION_SPLIT_ON_CHANGE = 'split-on-change'
    FUSION_EXPLODE_ON_CHANGE = 'explode-on-change'

    MARKERS_KEEP_IN_MUTATION_TIL_CHANGE = 'keep-in-mutation-til-change'
    MARKERS_KEEP_AS_PRESENCE = 'keep-as-presence'

    def __init__(self, changes, parent=None, fusion_strategy=FUSION_MATCH_WHOLE, unambiguous_mode=False):
        self.parent = parent
        self._changes = tuple(changes)
        self.mutations = parent.mutations if parent else []
        self.presences = parent.presences if parent else []

        print "Changes:", changes
        def get_match(obj, element):
            if isinstance(obj, Feature):
                return Match.COMPLETE if isinstance(element, Feature) and element.match(obj) else Match.NONE
            elif isinstance(obj, FeatureSet):
                if isinstance(element, FeatureSet):
                    return Match.COMPLETE if element.match(obj) else Match.NONE
                elif isinstance(element, (Fusion, Feature)):
                    for part in obj:
                        result = element.match(part)
                        if result != Match.NONE:
                            return Match.PARTIAL
                        return Match.NONE
                else:
                    return Match.None
            elif isinstance(obj, Fusion):
                if isinstance(element, Fusion) and element.match(obj):
                    return Match.COMPLETE
                elif isinstance(element, (Feature, Fusion, FeatureSet)):
                    return Match.PARTIAL if obj.contains(element) else Match.NONE
                else:
                    return Match.NONE

        def substitute(obj, old, new):
            count = 0
            if isinstance(obj, Feature):
                return new if obj == old else obj, 1 if obj == old else 0
            elif isinstance(obj, FeatureSet):
                new_contents = []
                for element in obj:
                    new_element, sub_count = substitute(element, old, new)
                    count += sub_count
                    if new_element:
                        new_contents.append(new_element)

                return FeatureSet(*new_contents) if new_contents else None, count
            elif isinstance(obj, Fusion):
                new_contents = []
                for element in obj:
                    # Element can be a Feature or a FeatureSet. If element is a feature, don't process it now so that
                    # e.g. +A:B B>C:D does not parse to Fusion(A, Fusion(C, D)) - update only FeatureSets
                    new_element, sub_count = substitute(element, old, new) if isinstance(element, FeatureSet) \
                        else (element, 0)
                    count += sub_count
                    if new_element:
                        new_contents.append(new_element)

                new_fusion = Fusion(*new_contents)
                # locate range to be removed or substituted
                target_range = new_fusion.part_range(old)
                while target_range:
                    count += 1
                    pos = target_range.start
                    if new and not isinstance(new, Fusion):
                        new = [new]
                    new_fusion[target_range] = new or []
                    pos += len(new) if isinstance(new, list) else 0
                    target_range = new_fusion.part_range(old, pos)

                if len(new_fusion) > 1:
                    return new_fusion, count
                elif len(new_fusion) == 1:
                    return new_fusion[0], count
                else:
                    return None, count

        def update_mutations(change):
            new_mutations = []
            count = 0
            append_change = True
            for mutation in self.mutations:
                if isinstance(mutation, Mutation):
                    updated_mutation, new_count, new_append_change = apply_change_to_mutation(mutation, change, unambiguous_mode)
                    count += new_count
                else:
                    updated_mutation, new_append_change = apply_change_to_presence(mutation, change)
                    print "Return:", updated_mutation, new_append_change
                append_change &= new_append_change
                if updated_mutation is not None:
                    new_mutations.append(updated_mutation)

            if count > 1:
                if unambiguous_mode:
                    raise AmbiguityError("Change {} matches multiple mutations.".format(change_to_string(change)))
                else:
                    raise NotImplementedError()

            if append_change:
                new_mutations.append(change)

            return new_mutations

        def apply_change_to_mutation(mutation, change, unambiguous_mode=False):
            # presences can only influence the displayable markers within mutation
            if isinstance(change, Presence):
                if mutation.markers:
                    mutation.set_markers(*[marker for marker in mutation.markers if change.element.match(marker)])
                return mutation, 0, False  # TODO: 0, False or something else?

            # Loci does not match
            # +A@L1 -A / +A -A@L1
            # if mutation.locus != change.locus:
            #     return mutation, 0, True

            # Change is the opposite of a mutation -> Delete that mutation, no additional changes
            # +A -A / -A +A / A>B B>A / B>A A>B

            if change.is_opposite_to(mutation):
                return None, 1, False

            # Any change that is repeated after it is established is ambiguous. It has no additional effect.
            if change.is_identical_to(mutation):
                if unambiguous_mode:
                    raise AmbiguityError("Mutation {} is already in genotype.".format(change_to_string(change)))
                return mutation, 0, False

            # Before part matches the whole before part of exactly one mutation that has no after part(deletion)
            # -> Replace mutation, Ambiguous
            if change.before is not None and mutation.before is not None and get_match(mutation.before, change.before) == Match.COMPLETE:
                if unambiguous_mode:
                    raise AmbiguityError()
                if change.after is None:  # A>B -A
                    return mutation, 1, False
                else:  # -A A>B
                    return None, 1, True  # TODO: Should it be 1?

            # TODO: A>B -B = -A but +A A>B = A>B
            if change.before is not None and mutation.after is not None:
                match = get_match(mutation.after, change.before)
                print match
                # if no match or if partial match but fusion strategy is match whole, don't change mutation
                if match == Match.NONE or (match == Match.PARTIAL and fusion_strategy == Genotype.FUSION_MATCH_WHOLE):
                    return mutation, 0, True  # and add change
                # if match is complete but mutation was insertion and change was substitution (+A A>B), keep only change
                elif match == Match.COMPLETE and mutation.before is None and change.after is not None:
                        return None, 1, True
                else:
                    new_after, count = substitute(mutation.after, change.before, change.after)
                    print "Count: ", count
                    return Mutation(mutation.before, new_after), count, False

            # raise AmbiguityError("Case not considererd")
            return mutation, 0, True

# None (MATCH WHOLE/ UPDATE), PARTIAL (MATCH WHOLE), COMPLETE ()


            # if mutation.after.contains(change.before):
            #     new_after, count = substitute(mutation.after, change.before, change.after)
            #
            #     if new_after is None:
            #         # return only the markers (as presences)
            #         return None, count
            #
            #     if KEEP_NEW_MARKERS:
            #         markers = change.markers
            #     else:
            #         markers = None
            #
            #     return Mutation(mutation.before, mutation.after,
            #                     locus=mutation.locus,
            #                     multiple=mutation.multiple,
            #                     markers=change.markers), count


        def apply_change_to_presence(presence, change):
            print "Presence: ", presence, change
            if isinstance(change, Presence):
                kwargs = {"match_variant": False} if isinstance(change.element, Feature) else {}
                print change.match(presence, **kwargs), change, presence, kwargs
                if change.match(presence, **kwargs):
                    return (change, False) if change.present == presence.present else (None, False)
                else:
                    return presence, True

            # "M-, A>B::M+" becomes "A>B::M+"
            if change.markers:

                    if marker.match(presence):
                        return None
            # if any(marker.match(presence) for marker in change.markers):
            #     return  # new presence, (last match counts)

            if change.before.match(presence.element):
                return None

        for change in changes:
            print "Change:", change
            self.mutations = update_mutations(change)

        print "RESULT:", self.mutations

    @property
    def raw(self):
        return self._changes

    @classmethod
    def _parse_string(cls, string, *args, **kwargs):
        parser = GnomicParser()
        semantics = DefaultSemantics(*args, **kwargs)
        return parser.parse(string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='start')

    @classmethod
    def is_valid(cls, string):
        """
        Tests whether a gnomic genotype definition can be parsed.

        :param str string: The genotype definition to validate.
        :return: :class:`bool`
        """
        try:
            cls._parse_string(string)
            return True
        except GrakoException:
            return False

    @classmethod
    def parse(cls,
              string,
              parent=None,
              organisms=None,
              types=None,
              **kwargs):
        """

        :param str string: The genotype definition to parse
        :param Genotype parent:
        :param organisms: A list or tuple of :class:`Organism` objects used to resolve organism designations in
            features. If empty, uses: :attr:`DEFAULT_ORGANISMS`.
        :param types: A list or tuple of :class:`Type` objects used to resolve type designations in features.
            If empty, uses: :attr:`DEFAULT_TYPES`.
        :param kwargs: Any additional attributes, such as `fusion_strategy`, used with :class:`Genotype`.
        :return: A :class:`Genotype` instance with the parsed gnomic genotype.
        """
        changes = cls._parse_string(string,
                                    organisms or DEFAULT_ORGANISMS,
                                    types or DEFAULT_TYPES)
        return Genotype(changes, parent=parent, **kwargs)

    @staticmethod
    def chain_parse(definitions, **kwargs):
        genotype = Genotype.parse(definitions[0], **kwargs)
        for definition in definitions[1:]:
            genotype = Genotype.parse(definition, parent=genotype, **kwargs)
        return genotype

    def _iter_changes(self, fusions=True):
        if fusions:
            for mutation in self.mutations:
                yield mutation
        else:
            for mutation in self.mutations:
                if isinstance(mutation, Mutation):
                    if mutation.before is not None:
                        for feature in mutation.before.features():
                            yield Del(feature)
                    if mutation.after is not None:
                        for feature in mutation.after.features():
                            yield Ins(feature)
                else:
                    yield mutation

        # for presence in self.presences:
        #     yield presence

    # TODO some getter for accessing the original changes (self._changes)

    def changes(self, fusions=False):
        """
        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :return:
        """
        return set(self._iter_changes(fusions=fusions))

    def format(self, fusions=True, output='text'):
        """
        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :param str output: Output format; one of ``'text'`` and ``'string'``
        :return:
        """
        if output == 'text':
            return genotype_to_text(self, fusions=fusions)
        elif output == 'string':
            return genotype_to_string(self, fusions=fusions)
        else:
            raise NotImplementedError('Unknown output format: {}'.format(output))

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr({'parent': self.parent, 'changes': self._changes}))