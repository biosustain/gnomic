import pytest

from gnomic.grammar import GnomicParser
from gnomic.semantics import DefaultSemantics
from gnomic.types import Feature, Plasmid, Change, Accession, Fusion, CompositeAnnotation


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


def test_parse_deletions(parse):
    assert [Change(before=Feature('foo'))] == parse('-foo')


def test_parse_insertions(parse):
    assert [Change(after=Feature('foo'))] == parse('+foo')
    assert [Change(after=Feature('foo', variant=('variant', )))] == parse('foo(variant)')


def test_parse_replacements(parse):
    assert [Change(before=Feature('foo'), after=Feature('bar'))] == parse('foo>bar')
    assert [Change(before=Feature('foo'), after=Feature('bar'), multiple=True)] == parse('foo>>bar')


def test_parse_feature():
    assert Feature('foo') == Feature.parse('foo')
    assert Feature('foo', type='gene') == Feature.parse('gene.foo')
    assert Feature('foo', organism='organism') == Feature.parse('organism/foo')
    assert Feature('foo', accession=Accession('id', 'database')) == Feature.parse('foo#database:id')
    assert Feature('foo', accession=Accession('id123')) == Feature.parse('foo#id123')
    assert Feature('foo', type='gene', organism='organism', accession=Accession('id', 'database'),
                   variant=tuple('variant')) == Feature.parse('organism/gene.foo#database:id(variant)')


def test_genotype_parse_plasmid_insertion(parse):
    assert [Change(after=Plasmid('foo'))] == parse('(foo)')
    assert [Change(after=Plasmid('foo', [Feature.parse('gene.A')]))] == parse('(foo gene.A)')


def test_parse_fusion(parse):
    assert [Change(after=Fusion(Feature('foo'), Feature('bar')))] == parse('+foo:bar')
    assert [Change(before=Fusion(Feature('foo'), Feature('bar')), after=Feature('foo'))] == parse('foo:bar>foo')


def test_parse_composite_annotation(parse):
    assert [Change(after=CompositeAnnotation(Feature('geneA'), Feature('geneB')))] == parse('+{geneA, geneB}')
    assert [Change(before=CompositeAnnotation(Feature('geneA'), Feature('geneB')))] == parse('-{geneA, geneB}')
