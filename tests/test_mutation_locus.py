import pytest

from gnomic.types import Feature, Change, Fusion
from gnomic.utils import chain


@pytest.mark.skip
def test_imply_locus_with_substitution():

    assert chain('gene.A>promoter.X:gene.B',
                 'gene.C>promoter.X:gene.D',
                 'promoter.X@gene.A>promoter.Y').changes() == (
        Change(Feature.parse('gene.C'), Fusion(Feature.parse('promoter.X'), Feature.parse('gene.D'))),
        Change(Feature.parse('gene.A'), Fusion(Feature.parse('promoter.Y'), Feature.parse('gene.B'))),
    )


