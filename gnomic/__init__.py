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
        self.changes = changes

    @classmethod
    def parse(self, value,
              parent=None,
              organisms=DEFAULT_ORGANISMS,
              types=DEFAULT_TYPES):
        parser = GnomicParser()
        semantics = DefaultSemantics(organisms, types)
        contents = parser.parse(value,
                                whitespace='',
                                semantics=semantics,
                                rule_name='start')

        return contents

    def changes(self, inclusive=True, fusions=True):
        """

        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :return:
        """
        pass

    def format(self, fusions=True):
        """

        :param bool fusions: Keeps fusions together if ``True``, otherwise includes them as individual features.
        :return:
        """
        pass

    def __str__(self):
        pass