import pytest

from gnomic import Genotype
from gnomic.types import Plasmid


def test_added_plasmids():
    genotype = Genotype.parse('+gene.B (pA) (pB) -(pB) -(pC)')
    assert genotype.added_plasmids == {Plasmid('pA')}


def test_removed_plasmids():
    genotype = Genotype.parse('+gene.A (pA) (pB) -(pB) -(pC)')
    assert genotype.removed_plasmids == {Plasmid('pC')}


# TODO added ...