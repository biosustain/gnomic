from grako.exceptions import GrakoException

from gnomic.models import *
from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics

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
    Type('promoter', 'P'),
    Type('terminator', 'T'),
    Type('gene')
)

DEFAULT_ORGANISMS = (
    Organism('Escherichia coli', aliases=['Ec', 'E.coli']),
    Organism('Saccaromices cerevisiae', aliases=['Sc', 'S.cerevisiae']),
)


class Genotype(object):
    """



    """

    FUSION_MATCH_WHOLE_ONLY = 'match-whole-only'
    FUSION_UPDATE_ON_CHANGE = 'update-on-change'
    FUSION_BREAK_ON_CHANGE = 'break-on-change'
    FUSION_EXPLODE_ON_CHANGE = 'explode-on-change'

    def __init__(self, changes, parent=None, fusion_strategy=FUSION_MATCH_WHOLE_ONLY):
        self.parent = parent
        self._changes = changes

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

        def remove(features, exclude):
            return {feature for feature in features if not exclude.match(feature)}

        def upsert(features, addition):
            return remove(features, addition) | {addition}

        for change in changes:
            if isinstance(change, Plasmid):
                # add a plasmid (not integrated)
                added_plasmids |= {change}
                removed_plasmids -= {change}
            elif isinstance(change, Feature):
                # define a certain phenotype or variant of a gene:
                # - rid the genotype of other variants of this feature.
                # - replace the feature within any fusions that contain it or a variant of it.
                phenotypes = remove(phenotypes, change)

                added_features = upsert(added_features, change)
                removed_features = remove(removed_features, change)

                # fusion-sensitive implementation:
                added_fusion_features = upsert(added_fusion_features, change)
                removed_fusion_features = remove(removed_fusion_features, change)
            else:
                # mutation:
                if isinstance(change.old, Plasmid):
                    # deletion of a plasmid
                    added_plasmids = remove(added_plasmids, change.old)
                    removed_plasmids = upsert(removed_plasmids, change.old)
                elif change.old:
                    # deletion of one (or more) features or fusions
                    for feature in change.old.features():
                        removed_features = upsert(removed_features, feature)
                        added_features = remove(added_features, feature)

                    # TODO fusion-sensitive implementation
                    # fusion replace/delete strategies:
                    # FUSION_MATCH_WHOLE_ONLY       | A:B:C D:E  B>X -D:E    A:B:C -B +X
                    # FUSION_UPDATE_ON_CHANGE       | A:B:C B>X              A:X:C     (uses SPLIT on delete)
                    # FUSION_SPLIT_ON_CHANGE        | A:B:C:D C>X            A:B X D
                    # FUSION_EXPLODE_ON_CHANGE      | A:B:C:D C>X            A B X D
                    if fusion_strategy != Genotype.FUSION_MATCH_WHOLE_ONLY:
                        raise NotImplementedError('Unsupported fusion strategy: {}'.format(fusion_strategy))

                if change.new:
                    # if change.new is a plasmid, the features in the plasmid are integrated
                    if isinstance(change.new, Plasmid):
                        pass  # any further record of the integrated plasmid would go here

                    # insertion of one (or more) features or fusions
                    for feature in change.new.features():
                        added_features = upsert(added_features, feature)
                        removed_features = remove(removed_features, feature)

                    # fusion-sensitive implementation:
                    for feature_or_fusion in change.new:
                        added_fusion_features = upsert(added_fusion_features, feature_or_fusion)
                        removed_fusion_features = remove(removed_fusion_features, feature_or_fusion)

                if change.old and change.new:
                    # in a replacement, the removed part must be a single feature
                    upsert(sites, change.old.contents[0])

                if change.marker:
                    upsert(markers, change.marker)
                    upsert(added_features, change.marker)
                    remove(removed_features, change.marker)

        self.sites = tuple(sites)
        self.markers = tuple(markers)
        self.phenotypes = tuple(phenotypes)
        self.added_plasmids = tuple(added_plasmids)
        self.removed_plasmids = tuple(removed_plasmids)
        self.added_features = tuple(added_features)
        self.removed_features = tuple(removed_features)
        self.added_fusion_features = tuple(added_fusion_features)
        self.removed_fusion_features = tuple(removed_fusion_features)

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
        try:
            cls._parse_string(string)
            return True
        except GrakoException:
            return False

    @classmethod
    def parse(cls,
              string,
              parent=None,
              organisms=DEFAULT_ORGANISMS,
              types=DEFAULT_TYPES,
              **kwargs):
        changes = cls._parse_string(string, organisms, types)
        return Genotype(changes, parent=parent, **kwargs)

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

    def format(self, fusions=True):
        """
        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :return:
        """
        pass

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr({'parent': self.parent, 'changes': self._changes}))