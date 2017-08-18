from gnomic import Genotype
from gnomic.types import Plasmid, Fusion, Feature, CompositeAnnotation


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


def test_added_features():
    genotype = Genotype.parse('+geneA -geneB +geneB -geneC')
    assert genotype.added_features == {Feature('geneA')}
    genotype = Genotype.parse('+{geneA, geneB:geneC}')
    assert genotype.added_features == {Feature('geneA'), Feature('geneB'), Feature('geneC')}


def test_removed_features():
    genotype = Genotype.parse('+geneA -geneB +geneB -geneC')
    assert genotype.removed_features == {Feature('geneC')}


def test_added_fusion_features():
    genotype = Genotype.parse('+geneA -geneB:geneC +geneB:geneC +{geneA, geneB}')
    assert genotype.added_fusion_features == {Feature('geneA'), CompositeAnnotation(Feature('geneA'), Feature('geneB'))}


def test_removed_fusion_features():
    genotype = Genotype.parse('+geneA -geneB:geneC -geneA +{geneA, geneB}')
    assert genotype.removed_fusion_features == {Fusion(Feature('geneB'), Feature('geneC'))}
