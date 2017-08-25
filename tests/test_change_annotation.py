from gnomic.genotype import change_annotation
from gnomic.types import Feature as F, CompositeAnnotation


def test_replace_feature_in_composite_annotation():
    assert change_annotation(CompositeAnnotation(F('a'), F('b')), F('a'), None) == CompositeAnnotation(F('b'))

    assert change_annotation(CompositeAnnotation(F('b')), F('a'), None) == CompositeAnnotation(F('b'))

    assert change_annotation(CompositeAnnotation(F('a') ** F('b'),
                                                 F('c')),
                             F('a') ** F('b'),
                             F('d')) == CompositeAnnotation(F('d'), F('c'))

    assert change_annotation(CompositeAnnotation(F('a') ** F('b'),
                                                 F('c')),
                             F('a') ** F('b'),
                             F('c')) == CompositeAnnotation(F('c'))


def test_replace_feature_in_fusion():
    assert change_annotation(F('a') ** F('b'), F('a')) == F('b')
    assert change_annotation(F('a') ** F('b') ** F('c'), F('b')) == F('a') ** F('c')
    assert change_annotation(F('a') ** F('b') ** F('c'), F('b'), F('x')) == F('a') ** F('x') ** F('c')
    assert change_annotation(F('a') ** F('b') ** F('c'), F('b'),
                             F('x') ** F('y')) == F('a') ** F('x') ** F('y') ** F('c')


def test_remove_fusion_in_fusion():
    assert change_annotation(F('a') ** F('b') ** F('c') ** F('d'), F('b') ** F('c')) == F('a') ** F('d')
    assert change_annotation(F('a') ** F('b') ** F('c'), F('a') ** F('b')) == F('c')


def test_replace_fusion_in_fusion():
    assert change_annotation(F('a') ** F('b') ** F('c') ** F('d'),
                             F('b') ** F('c'), F('x')) == F('a') ** F('x') ** F('d')
    assert change_annotation(F('a') ** F('b') ** F('c'), F('a') ** F('b'),
                             F('x') ** F('y')) == F('x') ** F('y') ** F('c')


def test_insert_feature_in_fusion():
    assert change_annotation(F('a') ** F('b'), None, F('x')) == CompositeAnnotation(F('a') ** F('b'), F('x'))
