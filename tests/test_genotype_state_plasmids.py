import pytest

from gnomic.genotype import GenotypeState
from gnomic.types import Plasmid


@pytest.fixture
def state():
    return GenotypeState()


def test_insert_plasmid(state):
    state.change(+Plasmid(name='pA'))

    assert state.changes == (
        +Plasmid(name='pA'),
    )


def test_plasmid_insertion_followed_by_removal(state):
    state.change(+Plasmid(name='pA'))
    state.change(-Plasmid(name='pA'))

    assert state.changes == ()

    state.change(+Plasmid(name='pB'))
    assert state.changes == (
        +Plasmid(name='pB'),
    )


def test_plasmid_removal_followed_by_insertion(state):
    state.change(-Plasmid(name='pA'))
    state.change(+Plasmid(name='pA'))

    assert state.changes == ()

# def test_plasmids(self):
#     self.assertEqual({
#         Present(Plasmid('p3', []))
#     }, self.chain('-p1{} p2{} p3{}', '-p2{} p1{}').changes())
#
#     self.assertEqual({
#         Present(Plasmid('p2', [])),
#         Absent(Plasmid('p1', []))
#     }, self.chain('-p1{} p2{}').changes())
#
#     self.assertEqual({
#         Absent(Plasmid('p1', [])),
#         Present(Plasmid('p2', []))
#     }, self.chain('-p1{} -p1{} p2{} p2{}', unambiguous_mode=False).changes())
#
#     with self.assertRaises(AmbiguityError):
#         self.chain('-p1{} -p1{}')
#
#     with self.assertRaises(AmbiguityError):
#         self.chain('p1{} p1{}')


# def test_integrated_plasmid_vector(self):
#     self.assertEqual({
#         Del(Feature(name='siteA')),
#     }, self.chain('siteA>pA{}').changes())
#
#     self.assertEqual({
#         Sub(Feature(name='siteA'), FeatureSet(Feature(name='geneA'), Feature(name='geneB'))),
#     }, self.chain('siteA>pA{geneA geneB}').changes(True))
#
#     self.assertEqual({
#         Del(Feature(name='siteA')),
#         Ins(Feature(name='geneA')),
#         Ins(Feature(name='geneB')),
#     }, self.chain('siteA>pA{geneA geneB}').changes())
#
#     self.assertEqual({
#         Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneC'), Feature(name='geneD')))),
#     }, self.chain('+geneA:geneB', 'geneB>pA{geneC geneD}',
#                   fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))
#
#     self.assertEqual({
#         Ins(Fusion(Feature(name='geneA'), FeatureSet(), Feature(name='geneC'))),
#     }, self.chain('+geneA:geneB:geneC', 'geneB>pA{}',
#                   fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))
#
#     self.assertEqual({
#         Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
#     }, self.chain('+geneA:geneB:geneC', 'geneC>pA{}',
#                   fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))
