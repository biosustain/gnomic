import pytest

from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.types import Feature, Plasmid, Change


@pytest.fixture
def parse():
    def parse_(gnomic_string, *args, **kwargs):
        parser = GnomicParser()
        semantics = DefaultSemantics(*args, **kwargs)
        return parser.parse(gnomic_string,
                            whitespace='',
                            semantics=semantics,
                            rule_name='start')
    return parse_


def test_parse_feature():
    assert Feature('foo') == Feature.parse('foo')
    assert Feature('foo', type='gene') == Feature.parse('gene.foo')


def test_genotype_parse_plasmid_insertion(parse):
    assert [Change(after=Plasmid('foo'))] == parse('(foo)')
    assert [Change(after=Plasmid('foo', [Feature.parse('gene.A')]))] == parse('(foo gene.A)')
