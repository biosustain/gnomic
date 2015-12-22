from unittest import TestCase

from gnomic import Genotype


class GenotypeTestCase(TestCase):

    def chain(self, *definitions, **kwargs):
        genotype = Genotype.parse(definitions[0], **kwargs)

        for definition in definitions[1:]:
            genotype = Genotype.parse(definition, parent=genotype, **kwargs)

        return genotype

    def test_keep_fusions_where_possible(self):
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
