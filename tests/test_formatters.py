import pytest

from gnomic.formatters import GnomicFormatter
from gnomic.types import Feature


@pytest.fixture
def gnomic_formatter():
    return GnomicFormatter()


def test_feature_gnomic_format(gnomic_formatter):
    assert gnomic_formatter.format_annotation(Feature(name='foo')) == 'foo'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f)')) == 'foo(f)'
    assert gnomic_formatter.format_annotation(Feature.parse('foo(f, g)')) == 'foo(f; g)'
