from grako.exceptions import GrakoException

from gnomic.models import *
from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.utils import genotype_to_text, genotype_to_string

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

    def __init__(self, changes, parent=None, fusion_strategy=FUSION_MATCH_WHOLE):
        self.parent = parent
        self._changes = tuple(changes)
        self.mutation_list = parent.mutation_list if parent else []
        self.additions = parent.additions if parent else []

        def replace_mutation(old_mutation, new_mutation):
            # TODO: Add 'if' when locus is implemented
            # if old_mutation.locus == new_mutation.locus:
            #     return old_mutation

            is_insertion = new_mutation.new and not new_mutation.old
            target_feature = new_mutation.new if is_insertion else new_mutation.old
            substitute = new_mutation.new if new_mutation.new and new_mutation.old else None
            target_tree = old_mutation.old if is_insertion else old_mutation.new
            other_tree = old_mutation.new if is_insertion else old_mutation.old

            if not target_tree:
                return old_mutation

            if fusion_strategy == Genotype.FUSION_MATCH_WHOLE:
                new_feature = substitute if target_tree == target_feature else target_tree
            elif fusion_strategy == Genotype.FUSION_UPDATE_ON_CHANGE:
                new_feature = target_tree.updated_copy(target_feature, substitute)

            after = other_tree if is_insertion else new_feature
            before = new_feature if is_insertion else other_tree

            return Mutation(before, after) if before or after else None

        def replace_addition(old_addition, new_addition):
                kwargs = {"match_variant": False} if isinstance(new_addition.element, Feature) else {}
                if new_addition.element.match(old_addition.element, **kwargs):
                    return new_addition if new_addition.present == old_addition.present else None
                else:
                    return old_addition

        def update_component_list(component_list, change_list, replace_function):
            for change in change_list:
                new_change_list = [replace_function(component, change) for component in component_list]
                new_change_list = [new_change for new_change in new_change_list if new_change]
                if new_change_list == component_list:
                    component_list.append(change)
                else:
                    component_list = new_change_list
            return component_list

        # extract markers from changes and place them in the right positions
        updated_changes = []
        for change in changes:
            updated_changes.append(change)
            if change.markers:
                for marker in change.markers:
                    updated_changes.append(Present(marker))
                change.markers = None
        changes = updated_changes

        mutation_change_list = filter(lambda change: isinstance(change, Mutation), changes)
        self.mutation_list = update_component_list(self.mutation_list, mutation_change_list, replace_mutation)

        addition_change_list = filter(lambda change: isinstance(change, Presence), changes)
        self.additions = update_component_list(self.additions, addition_change_list, replace_addition)

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
            for mutation in self.mutation_list:
                yield mutation
        else:
            for mutation in self.mutation_list:
                if mutation.old is not None:
                    for feature in mutation.old.features():
                        yield Del(feature)
                if mutation.new is not None:
                    for feature in mutation.new.features():
                        yield Ins(feature)

        for addition in self.additions:
            yield addition

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