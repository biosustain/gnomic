from unittest import TestCase

from gnomic import Feature, Organism


class ParsingTestCase(TestCase):

    def test_parse_feature(self):
        self.assertEqual(Feature(name='geneA'), Feature.parse('geneA'))
        self.assertEqual(Feature(name='A', type='gene'), Feature.parse('gene.A'))
        self.assertNotEqual(Feature(name='geneA'), Feature.parse('geneA(x)'))
        self.assertEqual(Feature(name='geneA', variant='x'), Feature.parse('geneA(x)'))
        self.assertEqual(Feature(name='geneA', variant='wild-type', organism=Organism('foo')), Feature.parse('foo/geneA+'))
        self.assertEqual(Feature(name='geneA', variant='mutant', organism=Organism('foo')), Feature.parse('foo/geneA-'))
