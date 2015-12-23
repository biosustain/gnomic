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

    def __init__(self, changes, parent=None):
        self.parent = parent
        self._changes = changes

        sites = set()
        markers = set()
        phenotypes = set()
        added_plasmids = set()
        removed_plasmids = set()
        added_features = set()
        removed_features = set()
        added_fusion_features = set()
        removed_fusion_features = set()

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
                    # MATCH_WHOLE_ONLY       | A:B:C D:E  B>X -D:E    A:B:C -B +X
                    # UPDATE                 | A:B:C B>X              A:X:C     (uses SPLIT on delete)
                    # SPLIT                  | A:B:C:D C>X            A:B X D
                    # DECOMPOSE              | A:B:C:D C>X            A B X D

                if isinstance(change.new, Plasmid):
                    # insertion of a plasmid
                    added_plasmids = upsert(added_plasmids, change.new)
                    removed_plasmids = remove(removed_plasmids, change.new)
                elif change.new:
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
              types=DEFAULT_TYPES):
        changes = cls._parse_string(string, organisms, types)
        return Genotype(changes, parent=parent)

    def changes(self, inclusive=True, fusions=True):
        """
        :param inclusive: Whether to include changes from all parents.
        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :return:
        """


        if inclusive and self.parent:
            for change in self.parent.changes(inclusive=True, fusions=fusions):
                pass


    def format(self, fusions=True):
        """
        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :return:
        """
        pass

    def __repr__(self):
        return '{}({})'.format(self.__class__.__name__, repr({'parent': self.parent, 'changes': self._changes}))