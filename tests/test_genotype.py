from unittest import TestCase

from gnomic import Genotype, Feature, Ins, Del

class BaseTestCase(TestCase):

    def chain(self, *definitions, **kwargs):
        genotype = Genotype.parse(definitions[0], **kwargs)
        for definition in definitions[1:]:
            genotype = Genotype.parse(definition, parent=genotype, **kwargs)
        return genotype

class GenotypeTestCase(BaseTestCase):

    def test_chain_propagate_added_features(self):
        genotype = self.chain('+geneA', '+geneB')

        self.assertEqual({
            Ins(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
        }, genotype.changes())

    def test_chain_propagate_removed_features(self):
        self.assertEqual({
            Del(Feature(name='geneA')),
            Del(Feature(name='geneB')),
        }, self.chain('-geneA', '-geneB').changes())

        self.assertEqual({
            Del(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
            Ins(Feature(name='geneC')),
        }, self.chain('-geneA -geneB', '+geneB', '+geneC').changes())

    def test_integrated_plasmid_vector(self):
        self.assertEqual({
            Del(Feature(name='siteA')),
        }, self.chain('siteA>pA{}').changes())

        self.assertEqual({
            Del(Feature(name='siteA')),
            Ins(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
        }, self.chain('siteA>pA{geneA geneB}').changes())

    def test_variants(self):
        self.assertEqual({
            Del(Feature(name='geneA')),
        }, self.chain('+geneA(x)', '-geneA').changes())

        self.assertEqual({
            Ins(Feature(name='geneA')),
            Del(Feature(name='geneA', variant='x')),
        }, self.chain('+geneA', '-geneA(x)').changes())

        self.assertEqual({
            Del(Feature(name='geneA', variant='x')),
            Ins(Feature(name='geneA', variant='y')),
        }, self.chain('-geneA(x)', '+geneA(y)').changes())


class GenotypeFusionsTestCase(BaseTestCase):

    def test_break_fusions_on_deletion(self):
        genotype = self.chain('+P.promoterA:geneB:T.terminatorC +geneD',
                              '-geneB')


        print(genotype)

        # should have P.promoterA, T.terminatorC, geneD

        self.fail()

        genotype = self.chain('+geneA:geneB',
                              '+geneB',
                              '-geneA:geneB')

        # should have neither geneA, nor geneB

    def test_update_fusions_where_possible(self):
        genotype = self.chain('+P.promoterA:geneB:T.terminatorC +geneD',
                              '-geneB')


        print(genotype)

        # should have P.promoterA, T.terminatorC, geneD

        self.fail()

        genotype = self.chain('+geneA:geneB',
                              '+geneB',
                              '-geneA:geneB')

        # should have neither geneA, nor geneB

    def omit_changes_already_in_fusion(self):
        # XXX should they be omitted or not?

        genotype = self.chain('+P.promoterA:geneA +geneA')

        # should only return the fusion?
        # maybe there should be a flag.
