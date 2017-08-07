import pytest

from gnomic import Genotype
from gnomic.types import Plasmid, Fusion, Feature


def test_added_plasmids():
    genotype = Genotype.parse('+gene.B (pA) (pB) -(pB) -(pC)')
    assert genotype.added_plasmids == {Plasmid('pA')}


def test_removed_plasmids():
    genotype = Genotype.parse('+gene.A (pA) (pB) -(pB) -(pC)')
    assert genotype.removed_plasmids == {Plasmid('pC')}


def test_added_fusions():
    genotype = Genotype.parse('+A:B +B:C -B:C -C:D')
    assert genotype.added_fusions == {Fusion(Feature('A'), Feature('B'))}


def test_removed_fusions():
    genotype = Genotype.parse('+A:B +B:C -B:C -C:D')
    assert genotype.removed_fusions == {Fusion(Feature('C'), Feature('D'))}