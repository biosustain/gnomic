import pytest

from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.types import Feature, Plasmid, Change, Accession, Fusion, CompositeAnnotation

# http://varnomen.hgvs.org/recommendations/


@pytest.fixture
def parse():
    def parse_(gnomic_string, *args, **kwargs):
        parser = GnomicParser()
        semantics = DefaultSemantics(*args, **kwargs)
        return parser.parse(gnomic_string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='FEATURE')
    return parse_


def test_sequence_variants_dna_substitutions(parse):
    variants = [
        'g.123C>A',
        'c.93+1G>T',
        'c.79_80delinsTT',
        'c.[12C>T;13C>G]',
        'c.12G>H',
        'c.123=',
        'c.85=/T>C',
        'c.85=//T>C',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_deletions(parse):
    variants = [
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
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_duplication(parse):
    variants = [
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
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_insertion(parse):
    variants = [
        'g.4426_4426insA',
        'g.5756_5757insAGG',
        'g.123_124insL1234.1:23_361',
        'g.122_123ins123_234inv',
        # 'g.122_123ins213_234invinsAins123_211inv',
        'g.122_123ins212_234inv123_199inv',
        'c.(67_70)insG(p.Gly23fs)',
        'g.123_124ins(100)',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_inversion(parse):
    variants = [
        'g.1234_1236inv',
        'c.77_80inv',
        'g.122_123ins123_234inv',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_conversion(parse):
    variants = [
        'g.333_590con1844_2101',
        # other possibilities with reference to a file or a gene
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_delins(parse):
    variants = [
        'g.6775delinsGA',
        'g.6775_6777delinsC',
        'g.9002_9009delinsTTT',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_dna_repeated(parse):
    variants = [
        'g.123_124[14]',
        'g.123_124[14];[18]',
        'c.-128_-126[79]',
        'c.-128_-126[(600_800)]',
        'c.54GCA[21]',
        'c.54GCA[21]ACA[1]GCC[2]'
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_substitution(parse):
    variants = [
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
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_deletion(parse):
    variants = [
        'p.Ala3del',
        'p.(Ala3del)',
        'p.Ala3_Ser5del',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_duplication(parse):
    variants = [
        'p.Ala3dup',
        'p.(Ala3dup)',
        'p.Ala3_Ser5dup',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_insertion(parse):
    variants = [
        'p.His4_Gln5insAla',
        'p.Lys2_Gly3insGlnSerLys',
        'p.(Met3_His4insGlyTer)',
        'p.Arg78_Gly79ins23',
        'p.Cys28delinsTrpVal',
        'p.Cys28_Lys29delinsTrp',
        'p.(Pro578_Lys579delinsLeuTer)',
        'p.[Ser44Arg;Trp46Arg]',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_repeated(parse):
    variants = [
        'p.Ala2[10]',
        'p.Ala2[10];[11]',
        'p.Gln18[23]',
        'p.(Gln18)[(70_80)]',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_frameshift(parse):
    variants = [
        'p.Arg97ProfsTer23',
        'p.Arg97fs',
        'p.Ile327Argfs*?',
        'p.Gln151Thrfs*9',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))


def test_sequence_variants_protein_extension(parse):
    variants = [
        'p.Met1ext-5',
        'p.Met1Valext-12',
        'p.Ter110Clnext*17',
        'p.(Ter315TyrextAsnLysGlyThrTer)',
        'p.Ter327Argext*?',
        'p.*327Argext*',
    ]
    assert parse('geneA({})'.format('; '.join(variants))) == Feature('geneA', variant=tuple(variants))
