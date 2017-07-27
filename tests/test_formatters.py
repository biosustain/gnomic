import pytest

from gnomic.formatters import GnomicFormatter, TextFormatter
from gnomic.types import Feature, Change, Fusion, Plasmid, AtLocus
from gnomic import Genotype


@pytest.fixture
def gnomic_formatter():
    return GnomicFormatter()


@pytest.fixture
def text_formatter():
    return TextFormatter()


def test_feature_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_annotation(Feature(name='foo')) == 'foo'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f)')) == 'foo(f)'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f, g)')) == 'foo(f; g)'
    assert gnomic_formatter.format_annotation(Feature.parse('organism/type.gene#db:123(f, g)')) \
        == 'organism/type.gene#db:123(f; g)'
    assert gnomic_formatter.format_annotation(Feature.parse('#db:123')) == '#db:123'
    assert gnomic_formatter.format_annotation(Feature.parse('gene#123')) == 'gene#123'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(wild-type)')) == 'foo+'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(mutant)')) == 'foo-'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(mutant; variant)')) == 'foo-(variant)'


def test_change_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_change(Change(before=Feature('foo'))) == '-foo'
    assert gnomic_formatter.format_change(Change(after=Feature('foo'))) == '+foo'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('bar'))) == 'foo>bar'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('bar'), multiple=True)) == 'foo>>bar'


def test_fusion_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_fusion(Fusion(Feature('foo'), Feature('bar'))) == 'foo:bar'
    assert gnomic_formatter.format_fusion(Fusion(Feature('foo'), Feature('bar'), Plasmid('p'))) == 'foo:bar:(p)'


def test_plasmid_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_plasmid(Plasmid('foo')) == '(foo)'
    assert gnomic_formatter.format_plasmid(Plasmid('foo', annotations=(Feature('bar'),))) == '(foo bar)'


def test_at_locus_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_at_locus(AtLocus(Feature('foo'), Feature('bar'))) == 'foo@bar'


def test_genotype_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB site>feature')) \
           == '+geneA -geneB site>feature'
    # assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB site>>feature')) \
    #        == '+geneA -geneB site>>feature'
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB -geneA')) == '-geneB'


def test_genotype_text_format(text_formatter):
    assert text_formatter.format_genotype(Genotype.parse('+geneA')) == 'geneA'
    assert text_formatter.format_genotype(Genotype.parse('-geneA')) == u'\u0394geneA'
    #assert text_formatter.format_genotype(Genotype.parse('siteA>(pA)')) == u'\u0394siteA'
    assert text_formatter.format_genotype(Genotype.parse('foo>bar')) == u'\u0394foo'
    assert text_formatter.format_genotype(Genotype.parse('-(pA)')) == u'\u0394(pA)'
    assert text_formatter.format_genotype(Genotype.parse('+geneA(x)')) == 'geneA(x)'
    assert text_formatter.format_genotype(Genotype.parse('+geneA(wild-type, var)')) == 'geneA+(var)'
