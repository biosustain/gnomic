import pytest

from gnomic.genotype import GenotypeState
from gnomic.types import Feature as F, Plasmid, CompositeAnnotation


@pytest.fixture
def state():
    return GenotypeState()


def test_match_whole_fusion_on_delete(state):
    state.change(+F.parse('gene.A') ** F.parse('gene.B'))
    state.change(+F.parse('gene.C'))
    state.change(-F('gene.A'))

    assert state.changes == (
        +F.parse('gene.A') ** F.parse('gene.B'),
        +F.parse('gene.C'),
        -F('gene.A'),
    )

    state.change(-F.parse('gene.C'))
    assert state.changes == (
        +F.parse('gene.A') ** F.parse('gene.B'),
        -F('gene.A'),
    )

    state.change(-F.parse('gene.A') ** F.parse('gene.B'))
    assert state.changes == (
        -F('gene.A'),
    )


def test_match_whole_fusion_on_delete_with_variant(state):
    state.change(+F.parse('gene.A') ** F.parse('gene.B'))
    state.change(+F.parse('gene.C'))
    state.change(-F.parse('gene.A') ** F.parse('gene.B(x)'))

    assert state.changes == (
        +F.parse('gene.A') ** F.parse('gene.B'),
        +F.parse('gene.C'),
        -F.parse('gene.A') ** F.parse('gene.B(x)'),
    )


def test_integrate_plasmid_with_fusion(state):
    state.change(F('site') > Plasmid('pA', [F.parse('gene.A') ** F.parse('gene.B')]))

    assert state.changes == (
        F('site') > Plasmid('pA', [F.parse('gene.A') ** F.parse('gene.B')]),
    )

    state.change(F.parse('gene.A') ** F.parse('gene.B') > F.parse('gene.C'))
    assert state.changes == (
        F('site') > CompositeAnnotation(F.parse('gene.C')),
    )


# NOTE for now fusions have to be updated explictly (entirely) through replacement or using @locus
def test_update_fusion_in_deletion(state):
    state.change(+F.parse('gene.A') ** F.parse('gene.B'))
    state.change(+F.parse('gene.C'))
    state.change(-F.parse('gene.A'))

    assert state.changes == (
        +F.parse('gene.A') ** F.parse('gene.B'),
        +F.parse('gene.C'),
        -F.parse('gene.A'),
    )

    state.change(-F.parse('gene.A') ** F.parse('gene.B'))
    assert state.changes == (
        +F.parse('gene.C'),
        -F.parse('gene.A'),
    )


def test_update_fusion_in_deletion_at_locus(state):
    state.change(F.parse('gene.A') > F.parse('gene.B') ** F.parse('gene.C'))
    state.change(-F.parse('gene.C') % F.parse('gene.A'))

    assert state.changes == (
        F.parse('gene.A') > F.parse('gene.B'),
    )

# TODO for now fusions have to be updated explictly (entirely)
#     def test_fusion_delete_update_on_change(self):
#         self.assertEqual({
#             Ins(Feature(name='geneB')),
#             Ins(Feature(name='geneC'))
#         }, self.chain('+geneA:geneB +geneC', '-geneA',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Feature(name='geneC')),
#         }, self.chain('+geneA:geneB +geneC', '-geneA:geneB',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
#             Del(Fusion(Feature(name='geneA'), Feature(name='geneB', variant='x'))),
#             Ins(Feature(name='geneC')),
#         }, self.chain('+geneA:geneB +geneC', '-geneA:geneB(x)',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
#         }, self.chain('+geneA:geneB:geneC', '-geneB',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Feature(name='geneA')),
#         }, self.chain('+geneA:geneB:geneC', '-geneB:geneC',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'), Feature(name='geneC'))),
#             Del(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
#         }, self.chain('+geneA:geneB:geneC', '-geneA:geneC',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         with self.assertRaises(AmbiguityError):
#             self.chain('+geneA:geneB:geneC +geneA', '-geneA',
#                        fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True)
#
#         self.assertEqual({
#             Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneB')))),
#         }, self.chain('+geneA:{geneB geneC:geneD}', '-geneC:geneD',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneB'), Feature(name='geneD')))),
#         }, self.chain('+geneA:{geneB geneC:geneD}', '-geneC',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Fusion(Feature(name='geneA'), FeatureSet(Fusion(Feature(name='geneC'), Feature(name='geneD'))))),
#         }, self.chain('+geneA:{geneB geneC:geneD}', '-geneB',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         self.assertEqual({
#             Ins(Feature(name='geneA')),
#         }, self.chain('+geneA:{geneB geneC:geneD}', '-geneB', '-geneC:geneD',
#                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
#
#         with self.assertRaises(AmbiguityError):
#             self.chain('+geneA:geneB:geneC:geneA:geneB', '-geneA:geneB',
#                        fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True)


def test_replace_gene_with_fusion(state):
    state.change(+F.parse('gene.A'))
    state.change(F.parse('gene.A') > F.parse('gene.A') ** F.parse('gene.B'))

    assert state.changes == (
        +F.parse('gene.A') ** F.parse('gene.B'),
    )


def test_update_fusion_in_replacement(state):
    state.change(+F.parse('gene.A') ** F.parse('gene.B'))
    state.change(F.parse('gene.A') ** F.parse('gene.B') > F.parse('gene.C') ** F.parse('gene.D'))

    assert state.changes == (
        +F.parse('gene.C') ** F.parse('gene.D'),
    )

    state.change(F.parse('gene.C') ** F.parse('gene.D') > F.parse('gene.E'))
    assert state.changes == (
        +F.parse('gene.E'),
    )


def test_update_part_of_fusion_in_replacement(state):
    state.change(+F.parse('gene.A') ** F.parse('gene.B'))
    state.change(F.parse('gene.A') > F.parse('gene.C'))

    assert state.changes == (
        +F.parse('gene.C') ** F.parse('gene.B'),
    )

    state.change(F.parse('gene.C') ** F.parse('gene.B') > F.parse('gene.A') ** F.parse('gene.B') ** F.parse('gene.C'))
    assert state.changes == (
        +F.parse('gene.A') ** F.parse('gene.B') ** F.parse('gene.C'),
    )

    state.change(F.parse('gene.A') ** F.parse('gene.B') > F.parse('gene.D') ** F.parse('gene.E'))
    assert state.changes == (
        +F.parse('gene.D') ** F.parse('gene.E') ** F.parse('gene.C'),
    )


def test_update_fusion_with_composite_annotation_in_replacement(state):
    state.change(+F.parse('gene.A') ** F.parse('gene.B'))
    state.change(F.parse('gene.B') > CompositeAnnotation(F.parse('gene.C'), F.parse('gene.D')))

    assert state.changes == (
        +F.parse('gene.A') ** CompositeAnnotation(F.parse('gene.C'), F.parse('gene.D')),
    )
