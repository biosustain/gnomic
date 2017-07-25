import pytest

from gnomic.genotype import GenotypeState
from gnomic.types import Change, Feature as F


@pytest.fixture
def state():
    return GenotypeState()


def test_scenario_1a_disambiguate_site_using_loci(state):
    state.change(Change(F.parse('gene.A'), F.parse('promoter.X') ** F.parse('gene.B')))
    state.change(Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D')))

    assert state.changes == (
        Change(F.parse('gene.A'), F.parse('promoter.X') ** F.parse('gene.B')),
        Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D'))
    )

    state.change(Change(F.parse('promoter.X') % F.parse('gene.A'), F.parse('promoter.Y')))

    assert state.changes == (
        Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D')),
        Change(F.parse('gene.A'), F.parse('promoter.Y') ** F.parse('gene.B'))
    )


def test_scenario_1b_disambiguate_site_using_fusion(state):
    state.change(Change(F.parse('gene.A'), F.parse('promoter.X') ** F.parse('gene.B')))
    state.change(Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D')))

    state.change(Change(F.parse('promoter.X') ** F.parse('gene.B'), F.parse('promoter.Y') ** F.parse('gene.B')))
    assert state.changes == (
        Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D')),
        Change(F.parse('gene.A'), F.parse('promoter.Y') ** F.parse('gene.B'))
    )


def test_scenario_1c_disambiguate_site_using_repeat_site(state):
    state.change(Change(F.parse('gene.A'), F.parse('promoter.X') ** F.parse('gene.B')))
    state.change(Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D')))

    state.change(Change(F.parse('gene.A'), F.parse('promoter.Y') ** F.parse('gene.B')))
    assert state.changes == (
        Change(F.parse('gene.C'), F.parse('promoter.X') ** F.parse('gene.D')),
        Change(F.parse('gene.A'), F.parse('promoter.Y') ** F.parse('gene.B'))
    )


def test_scenario_2_change_at_unknown_locus(state):
    state.change(Change(F.parse('promoter.A') % F.parse('gene.foo'), F.parse('promoter.B')))

    assert state.changes == (
        Change(F.parse('promoter.A') % F.parse('gene.foo'), F.parse('promoter.B')),
    )

    state.change(Change(None, F.parse('gene.foo')))
    assert state.changes == (
        Change(F.parse('promoter.A') % F.parse('gene.foo'), F.parse('promoter.B')),
        Change(None, F.parse('gene.foo')),
    )

    state.change(Change(F.parse('promoter.B') % F.parse('gene.foo'), F.parse('promoter.A')))
    assert state.changes == (
        Change(None, F.parse('gene.foo')),
    )
