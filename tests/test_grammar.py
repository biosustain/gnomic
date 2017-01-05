from unittest import TestCase
from gnomic import Feature, Organism, Genotype, Ins, Accession, Type, DEFAULT_ORGANISMS, DEFAULT_TYPES, Del, Plasmid, \
    Sub, Fusion, Mutation, FeatureTree, FeatureSet


def parse(string):
    return Genotype._parse_string(string,
                                  organisms=DEFAULT_ORGANISMS,
                                  types=DEFAULT_TYPES)


class GrammarTestCase(TestCase):
    def test_empty_genotype(self):
        self.assertEqual([], parse(''))

    def test_parse_feature(self):
        self.assertEqual(Feature(name='geneA'), Feature.parse('geneA'))
        self.assertEqual(Feature(name='A', type=Type('gene')), Feature.parse('gene.A'))
        self.assertNotEqual(Feature(name='geneA'), Feature.parse('geneA(x)'))
        self.assertEqual(Feature(name='geneA', variant='x'), Feature.parse('geneA(x)'))
        self.assertEqual(Feature(name='geneA', variant='wild-type', organism=Organism('foo')),
                         Feature.parse('foo/geneA+'))
        self.assertEqual(Feature(name='geneA', variant='mutant', organism=Organism('foo')), Feature.parse('foo/geneA-'))

    def test_parse_simple_insertions(self):
        self.assertEqual([
            Ins(Feature(name='fooF'))
        ], parse('+fooF'))

        self.assertEqual([
            Ins(Feature(name='fooF', accession=Accession(identifier=123, database='FOO')))
        ], parse('+fooF#FOO:123'))

        self.assertEqual([
            Ins(Feature(accession=Accession(identifier=123, database='FOO')))
        ], parse('+#FOO:123'))

        self.assertEqual([
            Ins(Feature(accession=Accession(identifier='BAR', database='FOO')))
        ], parse('+#FOO:BAR'))

        self.assertEqual([
            Ins(Feature(accession=Accession(identifier=123)))
        ], parse('+#123'))


    def test_parse_variants(self):
        self.assertEqual([
            Feature(type=Type('phene'), name='A', variant='wild-type')
        ], parse('A+'))

        self.assertEqual([
            Feature(type=Type('phene'), name='A', variant='mutant')
        ], parse('A-'))

        self.assertEqual([
            Feature(type=Type('phene'), name='A', variant='custom')
        ], parse('A(custom)'))

        self.assertEqual([
            Feature(type=Type('phene'), name='A', variant='first, second')
        ], parse('A(first, second)'))

        self.assertEqual([
            Feature(type=Type('phene'),
                    name='A',
                    organism=Organism('Escherichia coli'),
                    variant='wild-type')
        ], parse('Ec/A+'))

        self.assertEqual([
            Feature(type=Type('phene'),
                    name='A',
                    organism=Organism('Escherichia coli'),
                    variant='wild-type')
        ], parse('Ec/A+'))

        self.assertEqual([
            Feature(type=Type('phene'),
                    accession=Accession(identifier=123, database='FOO'),
                    variant='wild-type')
        ], parse('#FOO:123+'))

    def test_parse_simple_deletions(self):
        self.assertEqual([
            Del(Feature(name='fooF'))
        ], parse('-fooF'))

        self.assertEqual([
            Del(Plasmid('pFoo', None))
        ], parse('-pFoo{}'))

        result = parse('-pFoo{fooF}')
        self.assertEqual([
            Del(Plasmid('pFoo', None))
        ], result)

        self.assertEqual((
            Feature(name='fooF'),
        ), result[0].old.contents)

    def test_parse_simple_replacements(self):
        self.assertEqual([
            Sub(Feature(name='fooF'), Feature(name='barB'))
        ], parse('fooF>barB'))

        self.assertEqual([
            Sub(Feature(name='fooF'), Feature(name='barB'), multiple=True)
        ], parse('fooF>>barB'))

        self.assertEqual([
            Sub(Feature(name='fooF'), Feature(name='barB'),
                marker=Feature(name='marker', variant='wild-type', type=Type('phene')),
                multiple=True)
        ], parse('fooF>>barB::marker+'))

    def test_parse_insert_groups(self):
        self.assertEqual([
            Ins(FeatureSet(Feature(name='fooF')))
        ], parse('+{fooF}'))

        self.assertEqual([
            Ins(FeatureSet(Feature(name='fooF'),
                            Feature(name='barB'),
                            Feature(name='bazZ')))
        ], parse('+{fooF barB bazZ}'))

        self.assertEqual([
            Ins(FeatureSet(Feature(name='fooF'),
                           Fusion(Feature(name='barB'), Feature(name='bazZ')),
                           Feature('lolO')))
        ], parse('+{fooF barB:bazZ lolO}'))

    def test_parse_plasmid_deletions(self):
        result = parse('-pFoo{barB, bazZ}')
        self.assertEqual([
            Del(Plasmid('pFoo', None))
        ], result)

        self.assertEqual((
            Feature(name='barB'),
            Feature(name='bazZ'),
        ), result[0].old.contents)

    def test_parse_fusion_featureset(self):
        self.assertEqual([
            Mutation(new=Fusion(Feature(name='featB'),  FeatureSet(Feature(name='featC'), Feature(name='featD'))),
                     old=Feature(name='featA'))
        ], parse('featA>featB:{featC featD}'))

        self.assertEqual([
            Mutation(new=Fusion(Feature(name='featB'), FeatureSet(Feature(name='featC')), Feature(name='featD')),
                     old=Feature(name='featA'))
        ], parse('featA>featB:{featC}:featD'))

    def test_parse_fusion_replacement(self):
        self.assertEqual([
            Mutation(old=Fusion(Feature(name='featA'), FeatureSet(Feature(name='featB'))),
                     new=Feature(name='featC'),
                     multiple=True,
                     marker=Feature(variant='wild-type', type=Type('phene'), name='marker'))
        ], parse('featA:{featB}>>featC::marker+'))
