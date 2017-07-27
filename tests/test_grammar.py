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

    def test_parse_plasmid_presence(self):
        self.assertEqual([Plasmid('pA', None)], parse('(pA)'))
        self.assertEqual([Plasmid('pA', None)], parse('pA{}'))

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

        self.assertEqual([
            Del(Plasmid('pBar', None))
        ], parse('-(pBar)'))

        changes = parse('-pFoo{fooF}')
        self.assertEqual([
            Del(Plasmid('pFoo', None))
        ], changes)

        self.assertEqual((
            Feature(name='fooF'),
        ), changes[0].old.contents)

        changes = parse('-(pBar barB)')
        self.assertEqual([
            Del(Plasmid('pBar', None))
        ], changes)

        self.assertEqual((
            Feature(name='barB'),
        ), changes[0].old.contents)

    def test_parse_simple_replacements(self):
        self.assertEqual([
            Sub(Feature(name='fooF'), Feature(name='barB'))
        ], parse('fooF>barB'))

        self.assertEqual([
            Sub(Feature(name='fooF'), Feature(name='barB'), multiple=True)
        ], parse('fooF>>barB'))

        self.assertEqual([
            Sub(Feature(name='fooF'), Feature(name='barB'),
                markers=[Feature(name='marker', variant='wild-type', type=Type('phene'))],
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
                     markers=[Feature(variant='wild-type', type=Type('phene'), name='marker')])
        ], parse('featA:{featB}>>featC::marker+'))

    def test_parse_markers(self):
        self.assertEqual([
            Ins(Feature(name='geneA'),
                markers=[Feature(variant='wild-type', type=Type('phene'), name='markerA')])
        ], parse('+geneA::markerA+'))

        self.assertEqual([
            Del(Feature(name='geneA'),
                markers=[Feature(variant='wild-type', type=Type('phene'), name='markerA'),
                Feature(variant='mutant', type=Type('phene'), name='markerB')])
        ], parse('-geneA::{markerA+ markerB-}'))

    def test_sequence_variants(self):
        self.assertEqual([
            Ins(Feature(name='geneA', variant='c.123A>G, p.M12N, p.Gln5*, p.X5del, foo'))
        ], parse('+geneA(c.123A>G, p.M12N, p.Gln5*, p.X5del, foo)'))

    def test_sequence_variants_dna_substitutions(self):
        substitutions = [
            'g.123C>A',
            'c.93+1G>T',
            'c.79_80delinsTT',
            'c.[12C>T;13C>G]',
            'c.12G>H',
            'c.123=',
            'c.85=/T>C',
            'c.85=//T>C',
        ]
        variants = ', '.join(substitutions)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_sequence_variants_dna_deletions(self):
        deletions = [
            'g.19del',
            'g.10_21del',
            'c.183_186+48del',
            'c.1234del',
            'c.1234+1del',
            'c.4072-1234_5155-246del',
            'c.(?_-245)_(31+1_32-1)del',
            'c.(?_-1)_(*1_?)del',
            'g.19_21=/del',
            'g.19_21del=//del'
        ]
        variants = ', '.join(deletions)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_sequence_variants_dna_duplication(self):
        duplication = [
            'g.7dup',
            'g.6_8dup',
            'c.120_123+48dup',
            'c.123dup',
            'c.4072-1234_5146-246dup',
            'c.(4071+1_4072-1)_(5145+1_5146-1)dup',
            'c.(4071+1_4072-1)_(5145+1_5146-1)[3]',
            'c.(?_-30)_(12+1_13-1)dup',
            'c.(?_-1)_(*1_?)dup',
        ]
        variants = ', '.join(duplication)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_sequence_variants_dna_insertion(self):
        insertion = [
            'g.4426_4426insA',
            'g.5756_5757insAGG',
            'g.123_124insL1234.1:23_361',
            'g.122_123ins123_234inv',
            # 'g.122_123ins213_234invinsAins123_211inv',
            'g.122_123ins212_234inv123_199inv',
            'c.(67_70)insG(p.Gly23fs)',
            'g.123_124ins(100)',
        ]
        variants = ', '.join(insertion)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_sequence_variants_dna_inversion(self):
        inversion = [
            'g.1234_1236inv',
            'c.77_80inv',
            'g.122_123ins123_234inv',
        ]
        variants = ', '.join(inversion)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_sequence_variants_dna_conversion(self):
        conversion = [
            'g.333_590con1844_2101',
            # other possibilities with reference to a file or a gene
        ]
        variants = ', '.join(conversion)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_dna_delins(self):
        delins = [
            'g.6775delinsGA',
            'g.6775_6777delinsC',
            'g.9002_9009delinsTTT',
        ]
        variants = ', '.join(delins)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_dna_repeated(self):
        repeated = [
            'g.123_124[14]',
            'g.123_124[14];[18]',
            'c.-128_-126[79]',
            'c.-128_-126[(600_800)]',
            'c.54GCA[21]',
            'c.54GCA[21]ACA[1]GCC[2]'
        ]
        variants = ', '.join(repeated)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_protein_substitution(self):
        substitution = [
            'p.Trp24Cys',
            'p.(Trp24Cys)',
            'p.Trp24Ter',
            'p.Trp24*',
            'p.Cys188=',
            'p.0',
            'p.?',
            'p.Met1?',
            'p.(Tyr4*)'
        ]
        variants = ', '.join(substitution)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_protein_deletion(self):
        deletion = [
            'p.Ala3del',
            'p.(Ala3del)',
            'p.Ala3_Ser5del',
        ]
        variants = ', '.join(deletion)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_protein_duplication(self):
        duplication = [
            'p.Ala3dup',
            'p.(Ala3dup)',
            'p.Ala3_Ser5dup',
        ]
        variants = ', '.join(duplication)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_protein_insertion(self):
        insertion = [
            'p.His4_Gln5insAla',
            'p.Lys2_Gly3insGlnSerLys',
            'p.(Met3_His4insGlyTer)',
            'p.Arg78_Gly79ins23',
            'p.Cys28delinsTrpVal',
            'p.Cys28_Lys29delinsTrp',
            'p.(Pro578_Lys579delinsLeuTer)',
            'p.[Ser44Arg;Trp46Arg]',
        ]
        variants = ', '.join(insertion)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_protein_repeated(self):
        repeated = [
            'p.Ala2[10]',
            'p.Ala2[10];[11]',
            'p.Gln18[23]',
            'p.(Gln18)[(70_80)]',
        ]
        variants = ', '.join(repeated)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_sequence_variants_protein_frameshift(self):
        frameshift = [
            'p.Arg97ProfsTer23',
            'p.Arg97fs',
            'p.Ile327Argfs*?',
            'p.Gln151Thrfs*9',
        ]
        variants = ', '.join(frameshift)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))
        
    def test_sequence_variants_protein_extension(self):
        extension = [
            'p.Met1ext-5',
            'p.Met1Valext-12',
            'p.Ter110Clnext*17',
            'p.(Ter315TyrextAsnLysGlyThrTer)',
            'p.Ter327Argext*?',
            'p.*327Argext*',
        ]
        variants = ', '.join(extension)
        self.assertEqual([
            Ins(Feature(name='geneA', variant=variants))
        ], parse('+geneA({})'.format(variants)))

    def test_variable_variants(self):
        self.assertEqual([
            Ins(Feature(name='geneA', variant='c.123A>G, org.foo.fooBar=42.0, fooString="bar"'))
        ], parse('+geneA(c.123A>G, org.foo.fooBar=42.0, fooString="bar")'))
