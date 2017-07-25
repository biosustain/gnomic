
from gnomic.types import Change, Feature as F, CompositeAnnotation, Fusion


def test_fusion_contains():
    assert Fusion(F('a'), F('b'), F('c'), F('d')).contains(Fusion(F('b'), F('c'))) is True
    assert Fusion(F('a'), F('b'), F('c'), F('d')).contains(Fusion(F('a'), F('c'))) is False
    assert Fusion(F('a'), F('b'), F('c'), F('d')).contains(F('a')) is True
    assert Fusion(F('a'), F('b'), F('c'), F('d')).contains(F('x')) is False

