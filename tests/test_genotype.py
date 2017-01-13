from unittest import TestCase, SkipTest

from gnomic import Genotype, Feature, Ins, Del, Fusion, Sub, Type, Range, FeatureSet


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
        # self.assertEqual({
        #     Del(Feature(name='siteA')),
        # }, self.chain('siteA>pA{}').changes())

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

    def test_phenotypes_replace_variants(self):

        # when variants are used (default case):
        self.assertEqual({
            Ins(Feature(name='geneA', variant='x')),
            Ins(Feature(name='geneA', variant='y')),
            Ins(Feature(name='geneA', variant='z')),
        }, self.chain('+geneA(x) +geneA(y)', '+geneA(z)').changes())

        # when phenotypes are used:
        self.assertEqual({
            Ins(Feature(name='pheneA', type=Type('phene'), variant='mutant')),
        }, self.chain('pheneA+', 'pheneA-').changes())

        # when variants are mixed:
        self.assertEqual({
            Ins(Feature(name='geneA', type=Type('phene'), variant='z')),
        }, self.chain('+geneA(x) +geneA(y)', 'geneA(z)').changes())

        self.assertEqual({
            Ins(Feature(name='geneA', variant='x')),
            Ins(Feature(name='geneA', type=Type('phene'), variant='z')),
        }, self.chain('+geneA(x) +geneA(y)', 'geneA(z)', '+geneA(x)').changes())

    def test_markers_as_phenotypes(self):
        self.assertEqual({
            Ins(Feature(name='geneA')),
            Ins(Feature(name='pheneA', type=Type('phene'), variant='mutant')),
        }, self.chain('+geneA::pheneA+', 'pheneA-').changes())

        self.assertEqual({
            Ins(Feature(name='pheneA', type=Type('phene'), variant='wild-type')),
        }, self.chain('+geneA::pheneA+', 'pheneA-', '-geneA::pheneA+').changes())

        self.assertEqual({
            Ins(Feature(name='pheneA', type=Type('phene'), variant='mutant')),
            Ins(Feature(name='pheneB', type=Type('phene'), variant='mutant')),
        }, self.chain('+geneA::{pheneA+ pheneB+}', 'pheneA-', '-geneA::pheneB-').changes())

    def test_multiple_insertion(self):
        self.assertEqual(Sub(Feature(name='siteA'), Feature(name='geneA'), multiple=True),
                         self.chain('siteA>>geneA').raw[0])

        self.assertNotEqual(Sub(Feature(name='siteA'), Feature(name='geneA'), multiple=True),
                            self.chain('siteA>geneA').raw[0])

    def test_no_delete_if_present(self):
        # geneA(x) is not marked as deleted as the removal was an exact match
        self.assertEqual({
            Ins(Feature(name='geneB')),
        }, self.chain('+geneA(x) +geneB', '-geneA(x)').changes())

        # geneA(x) is now marked as deleted as it was not present
        self.assertEqual({
            Ins(Feature(name='geneB')),
            Del(Feature(name='geneA', variant='x')),
        }, self.chain('+geneB', '-geneA(x)').changes())

        # geneA is marked as deleted because the match was not exact
        self.assertEqual({
            Del(Feature(name='geneA')),
        }, self.chain('+geneA(x)', '-geneA').changes())


class GenotypeRangeTestCase(BaseTestCase):

    def test_delete_range_basic(self):
        self.assertEqual({
            Del(Feature(name='geneA', range=Range('coding', 5, 10))),
        }, self.chain('-geneA[c.5_10]').changes())

        self.assertEqual({
            Del(Feature(name='geneA', range=Range('protein', 5, 5))),
        }, self.chain('-geneA[p.5]').changes())

    def test_delete_insert(self):
        self.assertEqual({
            Ins(Feature(name='geneA')),
        }, self.chain('-geneA[c.5_10]', '+geneA').changes())

    @SkipTest
    def test_delete_multiple_ranges(self):
        # TODO in the current implementation, only the most recently deleted range is accounted for.
        # TODO this implementation may change

        self.assertEqual({
        #    Del(Feature(name='geneA', range=Range('coding', 5, 10))),
            Del(Feature(name='geneA', range=Range('coding', 11, 12))),
        }, self.chain('-geneA[c.5_10]', '-geneA[c.11_12]').changes())

        self.assertEqual({
            Del(Feature(name='geneA'))
        }, self.chain('-geneA[c.5_10]', '-geneA').changes())


    # TODO detailed tracking of different bits & pieces of features.

class GenotypeFusionsTestCase(BaseTestCase):

    @SkipTest
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

    @SkipTest
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

    @SkipTest
    def omit_changes_already_in_fusion(self):
        # XXX should they be omitted or not?

        genotype = self.chain('+P.promoterA:geneA +geneA')

        # should only return the fusion?
        # maybe there should be a flag.

    def test_replace_fusion_feature_set(self):
        self.assertEqual({
            Ins(Feature(name='promoterA')),
            Ins(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
        }, self.chain('+promoterA:{geneA geneB}').changes())

        self.assertEqual({
            Ins(Fusion(Feature(name='promoterA'), FeatureSet(Feature(name='geneA'), Feature(name='geneB')))),
        }, self.chain('+promoterA:{geneA geneB}').changes(True))

        self.assertEqual({
            Ins(Feature(name='promoterA')),
            Ins(Feature(name='geneB'))
        }, self.chain('+promoterA:{geneA geneB}', '-geneA').changes())

        self.assertEqual({
            Ins(Fusion(Feature(name='promoterA'), FeatureSet(Feature(name='geneA'), Feature(name='geneB')))),
            Del(Feature(name='geneA'))
        }, self.chain('+promoterA:{geneA geneB}', '-geneA',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        # TODO fusion strategy FUSION_SPLIT_ON_CHANGE:
        # self.assertEqual({
        #     Ins(Fusion(Feature(name='promoterA'), FeatureSet(Feature(name='geneB')))),
        # }, self.chain('+promoterA:{geneA geneB}', '-geneA',
        #               fusion_strategy=Genotype.FUSION_SPLIT_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='promoterA'), FeatureSet(Feature(name='geneA'), Feature(name='geneB')))),
            Del(Feature(name='geneA'))
        }, self.chain('+promoterA:{geneA geneB}', '-geneA',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

    def test_fusion_delete_match_whole(self):
        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Ins(Feature(name='geneC')),
            Del(Feature(name='geneA')),
        }, self.chain('+geneA:geneB +geneC', '-geneA',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        # self.assertEqual({
        #     Ins(Fusion(Feature(name='geneA'), Feature(name='geneB')))
        # }, self.chain('+geneA:geneB +geneC', '-geneC',
        #               fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))
        #
        # self.assertEqual({
        #     Ins(Feature(name='geneC')),
        # }, self.chain('+geneA:geneB +geneC', '-geneA:geneB',
        #               fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))
        #
        # self.assertEqual({
        #     Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
        #     Del(Fusion(Feature(name='geneA'), Feature(name='geneB', variant='x'))),
        #     Ins(Feature(name='geneC')),
        # }, self.chain('+geneA:geneB +geneC', '-geneA:geneB(x)',
        #               fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

    def test_fusion_replace_match_whole(self):
        self.assertEqual({
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB', 'geneA:geneB>geneC',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
            Ins(Feature(name='geneD')),
        }, self.chain('+geneA:geneB', 'geneA:geneC>geneD',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB', variant='x'))),
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB', 'geneA:geneB(x)>geneC',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

    def test_integrated_plasmid_vector_fusion(self):
        # self.assertEqual({
        #     Del(Feature(name='siteA')),
        #     Ins(Feature(name='geneA')),
        #     Ins(Feature(name='geneB')),
        # }, self.chain('siteA>pA{geneA:geneB}').changes())

        self.assertEqual({
            Del(Feature(name='siteA')),
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
        }, self.chain('siteA>pA{geneA:geneB}').changes(fusions=True))


class GenotypeFusionsUpdateOnChangeTestCase(BaseTestCase):

    def test_fusion_insert_update_on_change(self):
        self.assertEqual({
            Ins(Feature(name='geneB')),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneC')))
        }, self.chain('-geneA:geneB:geneC', '+geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'), Feature(name='geneC')))
        }, self.chain('-geneA:geneB:geneC', '+geneA:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneB'), Feature(name='geneC'))),
            Del(Feature(name='geneA'))
        }, self.chain('-geneA:geneB:geneC', '+geneB:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneB')),
            Del(Feature(name='geneA'))
        }, self.chain('-geneA:geneB', '+geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB')))
        }, self.chain('-geneA:geneB', '+geneA:geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Ins(Fusion(Feature(name='geneB'), Feature(name='geneA')))
        }, self.chain('-geneA:geneB', '+geneB:geneA',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'), Feature(name='geneC'))),
            Ins(Feature(name='geneD'))
        }, self.chain('-geneA:geneB:{geneC geneD}', '+geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Ins(Feature(name='geneC')),
            Ins(Feature(name='geneD'))
        }, self.chain('-geneA:geneB:{geneC geneD}', '+geneC', '+geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Ins(FeatureSet(Feature(name='geneC'), Feature(name='geneD')))
        }, self.chain('-geneA:geneB:{geneC geneD}', '+{geneC geneD}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneC'), Feature(name='geneD')))),
            Ins(Feature(name='geneB'))
        }, self.chain('-geneA:geneB:{geneC geneB:geneD}', '+geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'),
                       FeatureSet(Feature(name='geneC'), Fusion(Feature(name='geneB'), Feature(name='geneE'))))),
            Ins(Feature(name='geneD'))
        }, self.chain('-geneA:{geneC geneB:geneD:geneE}', '+geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

    def test_fusion_delete_update_on_change(self):
        self.assertEqual({
            Ins(Feature(name='geneB')),
            Ins(Feature(name='geneC'))
        }, self.chain('+geneA:geneB +geneC', '-geneA',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB +geneC', '-geneA:geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB', variant='x'))),
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB +geneC', '-geneA:geneB(x)',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
        }, self.chain('+geneA:geneB:geneC', '-geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneA')),
        }, self.chain('+geneA:geneB:geneC', '-geneB:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'), Feature(name='geneC'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
        }, self.chain('+geneA:geneB:geneC', '-geneA:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneB'), Feature(name='geneC'))),
        }, self.chain('+geneA:geneB:geneC +geneA', '-geneA',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneC:geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneB'), Feature(name='geneD')))),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneA')),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneB', '-geneC:geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

    def test_fusion_replace_update_on_change(self):
        self.assertEqual({
            Ins(Fusion(Feature(name='geneB'), Feature(name='geneC'))),
        }, self.chain('+geneA', 'geneA>geneB:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB', 'geneA:geneB>geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
            Ins(Feature(name='geneD')),
        }, self.chain('+geneA:geneB', 'geneA:geneC>geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneC'), Feature(name='geneB'))),
        }, self.chain('+geneA:geneB', 'geneA>geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneD'), Feature(name='geneC'))),
        }, self.chain('+geneA:geneB:geneC', 'geneA:geneB>geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneC'), Feature(name='geneD'))),
        }, self.chain('+geneA:geneB', 'geneB>geneC:geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(FeatureSet(Feature(name='geneB'), Feature(name='geneC')), Feature(name='geneD'))),
        }, self.chain('+geneA:geneD', 'geneA>{geneB geneC}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(FeatureSet(Feature(name='geneA'), Feature(name='geneC'))),
        }, self.chain('+{geneA geneB}', 'geneB>geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneB'), Feature(name='geneE')))),
        }, self.chain('+geneA:{geneB geneC:geneD}', 'geneC:geneD>geneE',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'),
                       FeatureSet(Feature(name='geneC'), Feature(name='geneF')),
                       Feature(name='geneF'))),
        }, self.chain('+geneA:geneB:{geneC geneD:geneE}:geneD:geneE', 'geneD:geneE>geneF',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneD'), Feature(name='geneE')))),
        }, self.chain('+geneA:geneB:geneC', 'geneB:geneC>{geneD geneE}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneE'), Feature(name='geneD'))),
        }, self.chain('+geneA:{geneB geneC}:geneD', 'geneA:{geneB geneC}>geneE',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))


