import pytest

from gnomic.types import Change, Feature, AtLocus, Fusion


def test_remove_feature():
    assert Change(before=Feature('foo')) == -Feature('foo')


def test_insert_feature():
    assert Change(after=Feature('foo')) == +Feature('foo')


def test_replace_feature():
    assert Change(before=Feature('a'), after=Feature('b')) == (Feature('a') > Feature('b'))
    assert Change(before=Feature('a'), after=Feature('b'), multiple=True) == (Feature('a') >> Feature('b'))


def test_feature_at_locus():
    assert AtLocus(Feature('foo'), Feature('site')) == Feature('foo') % Feature('site')
    assert AtLocus(Feature('foo'), Feature('site')) == Feature('foo').__matmul__(Feature('site'))


def test_replace_feature_at_locus():
    assert Change(before=AtLocus(Feature('a'), Feature('site')),
                  after=Feature('b')) == (Feature('a') % Feature('site') > Feature('b'))


def test_remove_feature_at_locus():
    assert Change(before=AtLocus(Feature('foo'), Feature('site'))) == -Feature('foo') % Feature('site')

    with pytest.raises(ValueError):
        -Feature('foo') % Feature('a') % Feature('b')


def test_fusion():
    assert Fusion(Feature('a'), Feature('b')) == Feature('a') ** Feature('b')
    assert Fusion(Feature('a'), Feature('b'), Feature('c')) == Feature('a') ** Feature('b') ** Feature('c')


def test_insert_fusion():
    assert Change(after=Fusion(Feature('a'), Feature('b'))) == +Feature('a') ** Feature('b')


def test_remove_fusion():
    assert Change(before=Fusion(Feature('a'), Feature('b'))) == -Feature('a') ** Feature('b')


def test_replace_fusion():
    assert Change(before=Feature('a'),
                  after=Fusion(Feature('b'), Feature('c'))) == (Feature('a') > Feature('b') ** Feature('c'))

    assert Change(before=Fusion(Feature('a'), Feature('b')),
                  after=Feature('c')) == (Feature('a') ** Feature('b') > Feature('c'))

