from grako.exceptions import GrakoException

from gnomic.models import *
from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.utils import genotype_to_text, genotype_to_string, genotype_to_html

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

        # TODO FIXME renoval strategy: do not add to removed_features list if there was a match in added_features list
        # TODO same for plasmids!

        sites = set(parent.sites if parent else ())
        markers = set(parent.markers if parent else ())
        phenotypes = set(parent.phenotypes if parent else ())
        added_plasmids = set(parent.added_plasmids if parent else ())
        removed_plasmids = set(parent.removed_plasmids if parent else ())
        added_features = set(parent.added_features if parent else ())
        removed_features = set(parent.removed_features if parent else ())
        added_fusion_features = set(parent.added_fusion_features if parent else ())
        removed_fusion_features = set(parent.removed_fusion_features if parent else ())

        def remove(features, remove, **kwargs):
            return {feature for feature in features if not remove.match(feature, **kwargs)}

        def upsert(features, addition, **kwargs):
            return remove(features, addition, **kwargs) | {addition}

        # removes a feature, but only adds it to removed features if it isn't in added features.
        def remove_or_exclude(added_features, removed_features, exclude):
            for feature in added_features:
                if exclude == feature:
                    return remove(added_features, exclude), removed_features
            return remove(added_features, exclude), upsert(removed_features, exclude)

        for change in changes:
            if isinstance(change, Plasmid):
                # add a plasmid (not integrated)
                added_plasmids |= {change}
                removed_plasmids -= {change}
            elif isinstance(change, Feature):
                # define a certain phenotype or variant of a gene:
                # - rid the genotype of other variants of this feature.
                # - replace the feature within any fusions that contain it or a variant of it.
                phenotypes = upsert(phenotypes, change, match_variant=False)

                added_features = upsert(added_features, change, match_variant=False)
                removed_features = remove(removed_features, change, match_variant=False)

                # fusion-sensitive implementation:
                added_fusion_features = upsert(added_fusion_features, change, match_variant=False)
                removed_fusion_features = remove(removed_fusion_features, change, match_variant=False)
            else:
                # mutation:
                if isinstance(change.old, Plasmid):
                    # deletion of a plasmid
                    added_plasmids, removed_plasmids = remove_or_exclude(added_plasmids,
                                                                         removed_plasmids,
                                                                         change.old)
                elif change.old:
                    # deletion of one (or more) features or fusions
                    for feature in change.old.features():
                        added_features, removed_features = remove_or_exclude(added_features,
                                                                             removed_features,
                                                                             feature)

                    # TODO fusion-sensitive implementation
                    # fusion replace/delete strategies:
                    # FUSION_MATCH_WHOLE_ONLY       | A:B:C D:E  B>X -D:E    A:B:C -B +X
                    # FUSION_UPDATE_ON_CHANGE       | A:B:C B>X              A:X:C     (uses SPLIT on delete)
                    # FUSION_SPLIT_ON_CHANGE        | A:B:C:D C>X            A:B X D
                    # FUSION_EXPLODE_ON_CHANGE      | A:B:C:D C>X            A B X D
                    if fusion_strategy not in (Genotype.FUSION_MATCH_WHOLE,):
                        raise NotImplementedError('Unsupported fusion strategy: {}'.format(fusion_strategy))

                    # TODO: support: fusion split, fusion whole.

                    if fusion_strategy == Genotype.FUSION_MATCH_WHOLE:
                        for feature_or_fusion in change.old:
                            added_fusion_features, removed_fusion_features = remove_or_exclude(added_fusion_features,
                                                                                               removed_fusion_features,
                                                                                               feature_or_fusion)
                    elif fusion_strategy == Genotype.FUSION_SPLIT_ON_CHANGE:
                        pass

                if change.new:
                    # if change.new is a plasmid, the features in the plasmid are integrated
                    if isinstance(change.new, Plasmid):
                        pass  # any further record of the integrated plasmid would go here

                    # insertion of one (or more) features or fusions
                    for feature in change.new.features():
                        removed_features, added_features = remove_or_exclude(removed_features,
                                                                             added_features,
                                                                             feature)

                        # added_features = upsert(added_features, feature)
                        # removed_features = remove(removed_features, feature)

                    # fusion-sensitive implementation:
                    for feature_or_fusion in change.new:
                        added_fusion_features = upsert(added_fusion_features, feature_or_fusion)
                        removed_fusion_features = remove(removed_fusion_features, feature_or_fusion)

                if change.old and change.new:
                    # in a replacement, the removed part must be a single feature
                    sites = upsert(sites, change.old.contents[0])

                if change.markers:
                    for marker in change.markers:
                        markers = upsert(markers, marker, match_variant=False)
                        added_features = upsert(added_features, marker, match_variant=False)
                        removed_features = remove(removed_features, marker, match_variant=False)
                        added_fusion_features = upsert(added_fusion_features, marker, match_variant=False)
                        removed_fusion_features = remove(removed_fusion_features, marker, match_variant=False)

        self.sites = tuple(sites)
        self.markers = tuple(markers)
        self.phenotypes = tuple(phenotypes)
        self.added_plasmids = tuple(added_plasmids)
        self.removed_plasmids = tuple(removed_plasmids)
        self.added_features = tuple(added_features)
        self.removed_features = tuple(removed_features)
        self.added_fusion_features = tuple(added_fusion_features)
        self.removed_fusion_features = tuple(removed_fusion_features)

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
            for feature in self.added_fusion_features:
                yield Ins(feature)

            for feature in self.removed_fusion_features:
                yield Del(feature)
        else:
            for feature in self.added_features:
                yield Ins(feature)

            for feature in self.removed_features:
                yield Del(feature)

        for plasmid in self.added_plasmids:
            yield plasmid

        for plasmid in self.removed_plasmids:
            yield Del(plasmid)

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
        :param str output: Output format; one of ``'text'``, ``'string'`` and ``'html'``
        :return:
        """
        if output == 'text':
            return genotype_to_text(self, fusions=fusions)
        elif output == 'string':
            return genotype_to_string(self, fusions=fusions)
        elif output == 'html':
            return genotype_to_html(self, fusions=fusions)
        else:
            raise NotImplementedError('Unknown output format: {}'.format(output))

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr({'parent': self.parent, 'changes': self._changes}))