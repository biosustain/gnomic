import pytest

from gnomic.genotype import GenotypeState
from gnomic.types import Change, Feature as F, Present, AtLocus


@pytest.fixture
def state():
    return GenotypeState()


def test_insertion(state):
    state.change(Change(None, F('foo')))

    assert state.changes == (
        Change(None, F('foo')),
    )


def test_repeat_insertion(state):
    state.change(Change(None, F('foo')))
    state.change(Change(None, F('foo')))

    assert state.changes == (
        Change(None, F('foo')),
    )


def test_deletion(state):
    state.change(Change(F('foo'), None))

    assert state.changes == (
        Change(F('foo'), None),
    )


def test_repeat_deletion(state):
    state.change(Change(F('foo'), None))
    state.change(Change(F('foo'), None))

    assert state.changes == (
        Change(F('foo'), None),
    )

def test_deletion_followed_by_insertion(state):
    state.change(Change(F('foo'), None))
    assert state.changes == (
        Change(F('foo'), None),
    )

    state.change(Change(None, F('foo')))
    assert state.changes == ()


def test_insertion_followed_by_deletion(state):
    state.change(Change(None, F('foo')))

    assert state.changes == (
        Change(None, F('foo')),
    )

    state.change(Change(F('foo'), None))
    assert state.changes == ()


def test_replacement(state):
    state.change(Change(F('a'), F('b')))

    assert state.changes == (
        Change(F('a'), F('b')),
    )


def test_repeat_replacement(state):
    state.change(Change(F('a'), F('b')))
    state.change(Change(F('a'), F('b')))

    assert state.changes == (
        Change(F('a'), F('b')),
    )


def test_replacement_of_variants(state):
    state.change(Change(F.parse('a(x)'), F('b')))

    assert state.changes == (
        Change(F.parse('a(x)'), F('b')),
    )

    state.change(Change(F.parse('a(y)'), F('b')))
    assert state.changes == (
        Change(F.parse('a(x)'), F('b')),
        Change(F.parse('a(y)'), F('b')),
    )

    state.change(Change(F.parse('a(y)'), F('c')))
    assert state.changes == (
        Change(F.parse('a(x)'), F('b')),
        Change(F.parse('a(y)'), F('c')),
    )


def test_replacement_replace_variant(state):
    state.change(Change(F('a'), F.parse('a(x)')))
    assert state.changes == (
        Change(F('a'), F.parse('a(x)')),
    )

    state.change(Change(F('a'), F.parse('a(y)')))
    assert state.changes == (
        Change(F('a'), F.parse('a(y)')),
    )


def test_replacement_at_unknown_locus(state):
    state.change(Change(F('a') % F('locus'), F('b')))

    assert state.changes == (
        Change(F('a') % F('locus'), F('b')),
    )


def test_replacement_in_fusion_at_locus(state):
    state.change(Change(F('x'), F('a') ** F('b')))
    state.change(Change(F('y'), F('a') ** F('b')))

    assert state.changes == (
        Change(F('x'), F('a') ** F('b')),
        Change(F('y'), F('a') ** F('b')),
    )

    state.change(Change(F('b') % F('y'), F('b') ** F('c')))
    assert state.changes == (
        Change(F('x'), F('a') ** F('b')),
        Change(F('y'), F('a') ** F('b') ** F('c')),
    )


def test_change_at_locus_that_cancels_out(state):
    state.change(Change(F('a') % F('locus'), F('b')))

    assert state.changes == (
        Change(F('a') % F('locus'), F('b')),
    )

    state.change(Change(F('b') % F('locus'), F('a'), ))
    assert state.changes == ()


def test_change_that_cancels_out(state):
    state.change(Change(F('a'), F('b')))

    assert state.changes == (
        Change(F('a'), F('b')),
    )

    state.change(Change(F('b'), F('a')))
    assert state.changes == ()


def test_standalone_change_at_replaced_locus(state):
    state.change(Change(F('a'), F('b')))

    assert state.changes == (
        Change(F('a'), F('b')),
    )

    state.change(Change(AtLocus(F('b'), F('a')), F('c')))
    assert state.changes == (
        Change(F('a'), F('c')),
    )


def test_replacement_at_locus_of_locus(state):
    # XXX consider not simplifying (foo@foo>... and -foo@foo) if there is a logical use
    state.change(Change(F('a') % F('a'), F('b')))

    assert state.changes == (
        Change(F('a'), F('b')),
    )


def test_change_at_locus_when_before_is_at_locus(state):
    state.change(Change(F('a') % F('a'), F('b')))

    assert state.changes == (
        Change(F('a'), F('b')),
    )

    state.change(Change(F('c') % F('c'), None))
    assert state.changes == (
        Change(F('a'), F('b')),
        Change(F('c'), None),
    )


def test_change_presence(state):
    state.change(Present(F.parse('gene.foo(A)')))

    assert state.changes == (
        Change(None, F.parse('gene.foo(A)')),
    )

    state.change(Present(F.parse('gene.foo(B)')))

    assert state.changes == (
        Change(None, F.parse('gene.foo(B)')),
    )
