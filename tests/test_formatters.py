# -*- coding: utf-8 -*-
from __future__ import unicode_literals
import pytest

from gnomic.formatters import GnomicFormatter, TextFormatter, HTMLFormatter
from gnomic.types import Feature, Change, Fusion, Plasmid, AtLocus, Accession
from gnomic import Genotype


@pytest.fixture
def gnomic_formatter():
    return GnomicFormatter()


@pytest.fixture
def text_formatter():
    return TextFormatter()


@pytest.fixture
def html_formatter():
    return HTMLFormatter()


def test_feature_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_annotation(Feature(name='foo')) == 'foo'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f)')) == 'foo(f)'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f, g)')) == 'foo(f; g)'
    assert gnomic_formatter.format_annotation(Feature.parse('organism/type.gene#db:123(f, g)')) \
        == 'organism/type.gene#db:123(f; g)'
    assert gnomic_formatter.format_annotation(Feature.parse('#db:123')) == '#db:123'
    assert gnomic_formatter.format_annotation(Feature.parse('gene#123')) == 'gene#123'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(mutant; variant)')) == 'foo(mutant; variant)'


def test_change_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_change(Change(before=Feature('foo'))) == '-foo'
    assert gnomic_formatter.format_change(Change(after=Feature('foo'))) == '+foo'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('bar'))) == 'foo>bar'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('bar'), multiple=True)) == 'foo>>bar'
    assert gnomic_formatter.format_change(Change(after=Plasmid('foo'))) == '(foo)'
    assert gnomic_formatter.format_change(Change(before=Feature('foo'),
                                                 after=Feature('foo', variant=['x']))) == 'foo(x)'
    assert gnomic_formatter.format_change(Change(before=Feature('foo', variant=['x']),
                                                 after=Feature('foo', variant=['y']))) == 'foo(y)'
    assert gnomic_formatter.format_change(Change(before=AtLocus(Feature('foo', variant=['x']),
                                                                Feature('f')),
                                                 after=Feature('foo', variant=['y']))) == 'foo(x)@f>foo(y)'


def test_fusion_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_fusion(Fusion(Feature('foo'), Feature('bar'))) == 'foo:bar'
    assert gnomic_formatter.format_fusion(Fusion(Feature('foo'), Feature('bar'), Plasmid('p'))) == 'foo:bar:(p)'


def test_plasmid_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_plasmid(Plasmid('foo')) == '(foo)'
    assert gnomic_formatter.format_plasmid(Plasmid('foo', annotations=(Feature('bar'),))) == '(foo bar)'


def test_at_locus_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_at_locus(AtLocus(Feature('foo'), Feature('bar'))) == 'foo@bar'
    assert gnomic_formatter.format_genotype(Genotype.parse('foo@bar>new')) == 'foo@bar>new'


def test_genotype_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB site>feature')) \
        == '+geneA -geneB site>feature'
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB site>>feature')) \
        == '+geneA -geneB site>feature'
    assert gnomic_formatter.format_genotype(Genotype.parse('+geneA -geneB -geneA')) == '-geneB'
    assert gnomic_formatter.format_genotype(Genotype.parse('+{geneA, geneB} -{geneC, geneD}')) \
        == '+{geneA, geneB} -{geneC, geneD}'


def test_genotype_text_format(text_formatter):
    assert text_formatter.format_genotype(Genotype.parse('+geneA')) == 'geneA'
    assert text_formatter.format_genotype(Genotype.parse('-geneA')) == '\u0394geneA'
    assert text_formatter.format_genotype(Genotype.parse('siteA>(pA)')) == '\u0394siteA\u2192(pA)'
    assert text_formatter.format_genotype(Genotype.parse('foo>bar')) == '\u0394foo\u2192bar'
    assert text_formatter.format_genotype(Genotype.parse('foo>>bar')) == '\u0394foo\u2192bar'
    assert text_formatter.format_genotype(Genotype.parse('-(pA)')) == '\u0394(pA)'
    assert text_formatter.format_genotype(Genotype.parse('+geneA(x)')) == 'geneA(x)'
    assert text_formatter.format_genotype(Genotype.parse('+geneA(var1, var2)')) == 'geneA(var1; var2)'
    assert text_formatter.format_genotype(Genotype.parse('+{geneA, geneB}')) == '{geneA, geneB}'


def test_genotype_html_format(html_formatter):
    assert html_formatter.format_genotype(Genotype.parse('+geneA')) == '<span class="gnomic-feature">geneA</span>'
    assert html_formatter.format_genotype(Genotype.parse('-geneB')) \
        == '\u0394<span class="gnomic-feature">geneB</span>'
    assert html_formatter.format_genotype(Genotype.parse('foo>>bar')) \
        == '\u0394<span class="gnomic-feature">foo</span>\u2192<span class="gnomic-feature">bar</span>'
    assert html_formatter.format_genotype(Genotype.parse('(pA)')) \
        == '<span class="gnomic-plasmid">(<span class="gnomic-plasmid-name">pA</span>)</span>'
    assert html_formatter.format_genotype(Genotype.parse('+geneA(x)')) \
        == '<span class="gnomic-feature">geneA<sup>x</sup></span>'
    assert html_formatter.format_genotype(Genotype.parse('+geneA(x; wild-type; mutant; y)')) \
        == '<span class="gnomic-feature">geneA<sup>x; wild-type; mutant; y</sup></span>'
    assert html_formatter.format_genotype(Genotype.parse('+foo:bar')) \
        == '<span class="gnomic-fusion"><span class="gnomic-feature">foo</span>:' \
           '<span class="gnomic-feature">bar</span></span>'


def test_feature_text_format(text_formatter):
    assert text_formatter.format_feature(Feature('foo')) == 'foo'
    assert text_formatter.format_feature(Feature('foo', organism='org', type='gene')) == 'org/gene.foo'
    assert text_formatter.format_feature(Feature('foo', accession=Accession('123', 'db'))) == 'foo#db:123'


def test_plasmid_text_format(text_formatter):
    assert text_formatter.format_plasmid(Plasmid('foo')) == '(foo)'
    assert text_formatter.format_plasmid(Plasmid('foo', annotations=(Feature('bar'), Feature('x')))) == '(foo bar x)'


def test_feature_html_format_escape_html(html_formatter):
    assert html_formatter.format_feature(Feature('geneB&')) == '<span class="gnomic-feature">geneB&amp;</span>'
    assert html_formatter.format_feature(Feature('geneB<')) == '<span class="gnomic-feature">geneB&lt;</span>'
    assert html_formatter.format_feature(Feature('geneB>')) == '<span class="gnomic-feature">geneB&gt;</span>'
