import pytest

from gnomic.formatters import GnomicFormatter
from gnomic.types import Feature, Change
from gnomic import Genotype


@pytest.fixture
def gnomic_formatter():
    return GnomicFormatter()


def test_feature_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_annotation(Feature(name='foo')) == 'foo'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f)')) == 'foo(f)'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f, g)')) == 'foo(f; g)'


def test_change_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_change(Change(before=Feature('foo'))) == '-foo'
    assert gnomic_formatter.format_change(Change(after=Feature('foo'))) == '+foo'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('bar'))) == 'foo>bar'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('bar'), multiple=True)) == 'foo>>bar'


def test_genotype_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB site>feature')) \
           == '+geneA -geneB site>feature'
    # assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB site>>feature')) \
    #        == '+geneA -geneB site>>feature'
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB -geneA')) == '-geneB'
