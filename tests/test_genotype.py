from unittest import TestCase, SkipTest

from gnomic import Genotype, Feature, Ins, Del, Fusion, Sub, Type, Range, Plasmid, FeatureTree, Organism, FeatureSet, \
    Present, Absent, AmbiguityError
from gnomic.utils import genotype_to_text, feature_to_text, genotype_to_string, change_to_string


class BaseTestCase(TestCase):
    def chain(self, *definitions, **kwargs):
        return Genotype.chain_parse(list(definitions), **kwargs)


class GenotypeTestCase(BaseTestCase):
    def test_chain_propagate_added_features(self):
        self.assertEqual({
            Ins(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
        }, self.chain('+geneA', '+geneB').changes())

    def test_chain_propagate_removed_features(self):
        self.assertEqual({
            Del(Feature(name='geneA')),
            Del(Feature(name='geneB')),
        }, self.chain('-geneA', '-geneB').changes())

        self.assertEqual({
            Del(Feature(name='geneA')),
            Ins(Feature(name='geneC')),
        }, self.chain('-geneA -geneB', '+geneB', '+geneC').changes())

    def test_opposite_mutations(self):
        self.assertEqual(set(), self.chain('+geneA', '-geneA').changes())

        self.assertEqual(set(), self.chain('-geneA', '+geneA').changes())

        self.assertEqual(set(), self.chain('geneA>geneB', 'geneB>geneA').changes())

        self.assertEqual(set(), self.chain('geneA(x)>geneB', 'geneB>geneA').changes())
        #
        self.assertEqual({
            Sub(Feature(name='geneA'), Feature(name='geneB')),
            Sub(Feature(name='geneB', variant='x'), Feature(name='geneA'))
        }, self.chain('geneA>geneB', 'geneB(x)>geneA').changes(True))

    def test_identical_mutations(self):
        self.assertEqual({
            Ins(Feature(name='geneA')),
            Del(Feature(name='geneB')),
            Sub(Feature(name='geneC'), Feature(name='geneD'))
        }, self.chain('+geneA', '+geneA', '-geneB', '-geneB', 'geneC>geneD', 'geneC>geneD',
                      unambiguous_mode=False).changes(True))

        with self.assertRaises(AmbiguityError):
            self.chain('+geneA', '+geneA')

        with self.assertRaises(AmbiguityError):
            self.chain('-geneA', '-geneA')

        with self.assertRaises(AmbiguityError):
            self.chain('geneA>geneB', 'geneA>geneB')

    def test_simple_deletion_and_substitution(self):
        self.assertEqual({
            Sub(Feature(name='geneA'), Feature(name='geneB'))
        }, self.chain('-geneA', 'geneA>geneB', unambiguous_mode=False).changes(True))

        self.assertEqual({
            Sub(Feature(name='geneA'), Feature(name='geneB'))
        }, self.chain('geneA>geneB', '-geneA', unambiguous_mode=False).changes(True))

        with self.assertRaises(AmbiguityError):
            self.chain('geneA>geneB', '-geneA')

        with self.assertRaises(AmbiguityError):
            self.chain('-geneA', 'geneA>geneB')

        self.assertEqual({
            Del(Feature(name='geneA')),
            Sub(Feature(name='geneB'), Feature(name='geneA'))
        }, self.chain('-geneA', 'geneB>geneA').changes(True))

        self.assertEqual({
            Del(Feature(name='geneB'))
        }, self.chain('geneB>geneA', '-geneA').changes(True))

    def test_simple_insertion_and_substitution(self):
        self.assertEqual({
            Sub(Feature(name='geneA'), Feature(name='geneB'))
        }, self.chain('+geneA', 'geneA>geneB').changes(True))

        self.assertEqual({
            Ins(Feature(name='geneA')),
            Sub(Feature(name='geneB'), Feature(name='geneA'))
        }, self.chain('+geneA', 'geneB>geneA').changes(True))

        self.assertEqual({
            Ins(Feature(name='geneA')),
            Sub(Feature(name='geneA'), Feature(name='geneB'))
        }, self.chain('geneA>geneB', '+geneA').changes(True))

        self.assertEqual({
            Ins(Feature(name='geneA')),
            Sub(Feature(name='geneB'), Feature(name='geneA'))
        }, self.chain('geneB>geneA', '+geneA').changes(True))

    def test_mutations_with_fusions(self):
        self.assertEqual({
            Sub(Feature(name='geneB'), Fusion(Feature(name='geneA'), Feature(name='geneB')))
        }, self.chain('geneB>geneA', 'geneA>geneA:geneB').changes(True))

        self.assertEqual({
            Del(Feature(name='geneB')),
        }, self.chain('geneB>geneA:geneB', '-geneA:geneB').changes(True))

        self.assertEqual({
            Sub(Feature(name='geneB'), Fusion(Feature(name='geneD'), Feature(name='geneC')))
        }, self.chain('geneB>geneA:geneB:geneC', 'geneA:geneB>geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneA'))),
            Sub(Feature(name='geneA'), Feature(name='geneB')),
        }, self.chain('+geneA:geneA', 'geneA>geneB', unambiguous_mode=False,
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

        # +A:A A>B normally throws AmbiguityError for FUSION_UPDATE_ON_CHANGE strategy as there are two matches
        with self.assertRaises(AmbiguityError):
            self.chain('+geneA:geneA', 'geneA>geneB', fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE)

        # but for FUSION_MATCH_WHOLE +A:A A>B gives no matches so the test should pass
        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneA'))),
            Sub(Feature(name='geneA'), Feature(name='geneB')),
        }, self.chain('+geneA:geneA', 'geneA>geneB', fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(True))

        self.assertEqual({
            Sub(Fusion(Feature(name='geneA'), Feature(name='geneA')),
                Fusion(Feature(name='geneB'), Feature(name='geneB')))
        }, self.chain('+geneA:geneA', 'geneA:geneA>geneB:geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneB'), Feature(name='geneC')))
        }, self.chain('+geneA:geneC', 'geneA>geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneD'))),
            Sub(Feature(name='geneA'), Feature(name='geneB')),
        }, self.chain('+geneA:geneC +geneA:geneD', 'geneA>geneB', unambiguous_mode=False,
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))
        #
        with self.assertRaises(AmbiguityError):
            self.chain('+geneA:geneC +geneA:geneD', 'geneA>geneB', fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE)

    def test_loci(self):
        self.assertEqual({
            Ins(Feature(name='geneA'), locus='L1'),
            Del(Feature(name='geneA'), locus='L2'),
        }, self.chain('+geneA@L1', '-geneA@L2').changes())

        self.assertEqual({
            Ins(Feature(name='geneA'), locus='L1'),
            Sub(Feature(name='geneA'), Feature(name='geneB')),
        }, self.chain('+geneA@L1', 'geneA>geneB').changes(True))

        self.assertEqual({
            Ins(Feature(name='geneA'), locus='L1'),
            Ins(Feature(name='geneA')),
        }, self.chain('+geneA', '+geneA@L1').changes())

        self.assertEqual({
            Sub(Feature(name='geneA'), Feature(name='geneB'), locus='L1'),
            Del(Feature(name='geneB')),
        }, self.chain('geneA@L1>geneB', '-geneB').changes(True))

    def test_plasmids(self):
        self.assertEqual({
            Present(Plasmid('p3', []))
        }, self.chain('-p1{} p2{} p3{}', '-p2{} p1{}').changes())

        self.assertEqual({
            Present(Plasmid('p2', [])),
            Absent(Plasmid('p1', []))
        }, self.chain('-p1{} p2{}').changes())

        self.assertEqual({
            Absent(Plasmid('p1', [])),
            Present(Plasmid('p2', []))
        }, self.chain('-p1{} -p1{} p2{} p2{}', unambiguous_mode=False).changes())

        with self.assertRaises(AmbiguityError):
            self.chain('-p1{} -p1{}')

        with self.assertRaises(AmbiguityError):
            self.chain('p1{} p1{}')

    def test_phenotypes(self):
        self.assertEqual({
            Present(Feature(name='geneA', type=Type('phene'), variant='y')),
        }, self.chain('geneA(x)', 'geneA(y)').changes())

        self.assertEqual({
            Present(Feature(name='geneB', type=Type('phene'), variant='x')),
            Del(Feature(name='geneA', variant='x'))
        }, self.chain('geneA(x) geneB(x)', '-geneA(x)').changes())

        self.assertEqual({
            Present(Feature(name='geneA', type=Type('phene'), variant='y')),
            Ins(Feature(name='geneA', variant='x'))
        }, self.chain('+geneA(x)', 'geneA(y)').changes())

        self.assertEqual({
            Present(Feature(name='geneB', type=Type('phene'), variant='x')),
            Del(Feature(name='geneA'))
        }, self.chain('geneA(x) geneB(x)', '-geneA').changes())

        self.assertEqual({
            Present(Feature(name='geneA', type=Type('phene'), variant='x')),
        }, self.chain('geneA(x) geneA(x)', unambiguous_mode=False).changes())

        with self.assertRaises(AmbiguityError):
            self.chain('geneA(x) geneA(x)')

        self.assertEqual({
            Present(Feature(name='geneA', type=Type('phene'), variant='x')),
            Ins(Feature(name='geneA', variant='x'))
        }, self.chain('geneA(x)', '+geneA(x)').changes())

        self.assertEqual(set(), self.chain('geneA(x)', '+geneA(x)', '-geneA(x)').changes())

        self.assertEqual({
            Sub(Feature(name='geneA', variant='x'), Feature(name='geneB'))
        }, self.chain('geneA(x)', '+geneA(x)', 'geneA(x)>geneB').changes(True))

    def test_marker_presences(self):
        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='wild-type')),
            Ins(Feature(name='geneA'), markers=[Feature(name='M', type=Type('phene'), variant='wild-type')])
        }, self.chain('+geneA::M+').changes(True))

        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='wild-type')),
            Ins(Feature(name='geneA'), markers=[Feature(name='M', type=Type('phene'), variant='wild-type')])
        }, self.chain('M-', '+geneA::M+').changes(True))

        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='mutant')),
            Ins(Feature(name='geneA'))
        }, self.chain('+geneA::M+', 'M-').changes(True))

        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='mutant')),
            Present(Feature(name='N', type=Type('phene'), variant='wild-type')),
            Ins(Feature(name='geneA'), markers=[Feature(name='N', type=Type('phene'), variant='wild-type')])
        }, self.chain('+geneA::{M+ N+}', 'M-').changes(True))

        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='wild-type')),
            Sub(Feature(name='geneA'), Feature(name='geneC'))
        }, self.chain('geneA>geneB::M+', 'geneB>geneC').changes(True))

        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='x')),
            Present(Feature(name='N', type=Type('phene'), variant='y')),
            Sub(Feature(name='geneA'), Feature(name='geneC'),
                markers=[Feature(name='N', type=Type('phene'), variant='y')])
        }, self.chain('geneA>geneB::M(x)', 'geneB>geneC::N(y)').changes(True))

        self.assertEqual({
            Present(Feature(name='M', type=Type('phene'), variant='y')),
            Ins(Feature(name='geneA'), markers=[Feature(name='M', type=Type('phene'), variant='y')])
        }, self.chain('+geneA::{M(x) M(y)}', unambiguous_mode=False).changes(True))

        with self.assertRaises(AmbiguityError):
            self.chain('+geneA::{M(x) M(y)}')

    def test_integrated_plasmid_vector(self):
        self.assertEqual({
            Del(Feature(name='siteA')),
        }, self.chain('siteA>pA{}').changes())

        self.assertEqual({
            Sub(Feature(name='siteA'), FeatureSet(Feature(name='geneA'), Feature(name='geneB'))),
        }, self.chain('siteA>pA{geneA geneB}').changes(True))

        self.assertEqual({
            Del(Feature(name='siteA')),
            Ins(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
        }, self.chain('siteA>pA{geneA geneB}').changes())

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneC'), Feature(name='geneD')))),
        }, self.chain('+geneA:geneB', 'geneB>pA{geneC geneD}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(), Feature(name='geneC'))),
        }, self.chain('+geneA:geneB:geneC', 'geneB>pA{}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
        }, self.chain('+geneA:geneB:geneC', 'geneC>pA{}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(True))

    # def test_multiple_insertion(self):
    #     self.assertEqual(Sub(Feature(name='siteA'), Feature(name='geneA'), multiple=True),
    #                      self.chain('siteA>>geneA').raw[0])
    #
    #     self.assertNotEqual(Sub(Feature(name='siteA'), Feature(name='geneA'), multiple=True),
    #                         self.chain('siteA>geneA').raw[0])


class GenotypeRangeTestCase(BaseTestCase):
    def test_delete_range_basic(self):
        self.assertEqual({
            Del(Feature(name='geneA', range=Range('coding', 5, 10))),
        }, self.chain('-geneA[c.5_10]').changes())

        self.assertEqual({
            Del(Feature(name='geneA', range=Range('protein', 5, 5))),
        }, self.chain('-geneA[p.5]').changes())

    @SkipTest
    def test_delete_insert(self):
        self.assertEqual({
            Ins(Feature(name='geneA')),  # XXX needs rethinking.
        }, self.chain('-geneA[c.5_10]', '+geneA').changes())

        self.assertEqual({
            Ins(Feature(name='geneA')),
            Del(Feature(name='geneA', range=Range('coding', 5, 10))),
        }, self.chain('+geneA', '-geneA[c.5_10]').changes())

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

    def test_fusion_delete_match_whole(self):
        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Ins(Feature(name='geneC')),
            Del(Feature(name='geneA')),
        }, self.chain('+geneA:geneB +geneC', '-geneA',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB')))
        }, self.chain('+geneA:geneB +geneC', '-geneC',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB +geneC', '-geneA:geneB',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB', variant='x'))),
            Ins(Feature(name='geneC')),
        }, self.chain('+geneA:geneB +geneC', '-geneA:geneB(x)',
                      fusion_strategy=Genotype.FUSION_MATCH_WHOLE).changes(fusions=True))

    def test_integrated_plasmid_vector_fusion(self):
        self.assertEqual({
            Del(Feature(name='siteA')),
            Ins(Feature(name='geneA')),
            Ins(Feature(name='geneB')),
        }, self.chain('siteA>pA{geneA:geneB}').changes())

        self.assertEqual({
            Sub(Feature(name='siteA'), FeatureSet(Fusion(Feature(name='geneA'), Feature(name='geneB')))),
        }, self.chain('siteA>pA{geneA:geneB}').changes(fusions=True))


class GenotypeToTextTestCase(BaseTestCase):
    def test_added_features(self):
        self.assertEqual(
            genotype_to_text(self.chain('+geneA')),
            'geneA')

    def test_removed_features(self):
        self.assertEqual(genotype_to_text(self.chain('-geneA')),
                         u"\u0394geneA")

    def test_added_and_removed_features(self):
        self.assertEqual(genotype_to_text(self.chain('-geneB')),
                         u"\u0394geneB")

    def test_plasmid(self):
        self.assertEqual(genotype_to_text(self.chain('siteA>pA{}')),
                         u"\u0394siteA")

        self.assertEqual(genotype_to_text(self.chain('-pA{}')),
                         u"\u0394(pA)")

    def test_variants(self):
        self.assertEqual(genotype_to_text(self.chain('-geneA(x)')),
                         u"\u0394geneA(x)")

        self.assertEqual(genotype_to_text(self.chain('+geneA(x)')),
                         u"geneA(x)")

        self.assertEqual(genotype_to_text(self.chain('+geneA(wild-type)')),
                         u"geneA\u207A")

        self.assertEqual(genotype_to_text(self.chain('+geneA(mutant)')),
                         u"geneA\u207B")


class FeatureToTextTestCase(BaseTestCase):
    def test_plasmid_without_contents(self):
        feature = Plasmid("foo", contents=None)
        self.assertEqual(feature_to_text(feature), "foo")
        self.assertEqual(feature_to_text(feature, integrated=False), "(foo)")

    def test_plasmid_with_contents(self):
        feature = Plasmid("foo", contents=[Feature("bar")])
        self.assertEqual(feature_to_text(feature), "foo(bar)")
        self.assertEqual(feature_to_text(feature, integrated=False), "(foo bar)")

    def test_fusion(self):
        feature = Fusion(Feature(name="foo"), Feature(name="bar"))
        self.assertEqual(feature_to_text(feature), "foo:bar")

    def test_featuretree(self):
        feature = FeatureTree(Feature(name="foo"), Feature(name="bar"))
        self.assertEqual(feature_to_text(feature), "foo bar")

    def test_feature(self):
        feature = Feature(name="foo")
        self.assertEqual(feature_to_text(feature), "foo")

    def test_feature_with_organism(self):
        feature = Feature(name="foo", organism=Organism("bar"))
        self.assertEqual(feature_to_text(feature), "bar/foo")

    def test_feature_with_variant(self):
        feature_wild = Feature(name="foo", variant="wild-type")
        feature_mutant = Feature(name="foo", variant="mutant")
        feature_other = Feature(name="foo", variant="x")

        self.assertEqual(feature_to_text(feature_wild), u"foo\u207A")
        self.assertEqual(feature_to_text(feature_mutant), u"foo\u207B")
        self.assertEqual(feature_to_text(feature_other), "foo(x)")

    def test_feature_is_maker(self):
        feature = Feature(name="foo")
        self.assertEqual(feature_to_text(feature, is_marker=True), "::foo")


class GenotypeToStringTestCase(BaseTestCase):

    def test_genotype_to_string(self):
        changes = [
            '-e.coli/geneA',
            '+geneB(a)',
            'vectorC{geneC}',
            '-vectorD{}',
            '+geneE::markerF+',
            '+gene.G::{markerH+ markerI+}'
        ]
        gnomic = genotype_to_string(self.chain(*changes))

        self.assertEqual(gnomic, " ".join(changes))

        self.assertIsInstance(Genotype.parse(gnomic), Genotype)

    def test_change_to_string(self):
        self.assertEqual('-geneX', change_to_string(Del(Feature('geneX'))))
        self.assertEqual('-geneA:geneB', change_to_string(Del(Fusion(Feature('geneA'), Feature('geneB')))))

        self.assertEqual('-geneX::marker+',
                         change_to_string(Del(Feature('geneX'),
                                              markers=[Feature('marker', variant='wild-type')])))
        self.assertEqual('-geneX::{markerA+ markerB+}',
                         change_to_string(Del(Feature('geneX'),
                                              markers=[Feature('markerA', variant='wild-type'),
                                                       Feature('markerB', variant='wild-type')])))

        self.assertEqual('+geneX',
                         change_to_string(Ins(Feature('geneX'))))

        self.assertEqual('+{geneA geneB}',
                         change_to_string(Ins(FeatureSet(Feature('geneA'), Feature('geneB')))))

        self.assertEqual('plasmid{}',
                         change_to_string(Present(Plasmid('plasmid', []))))

        self.assertEqual('-plasmid{}',
                         change_to_string(Del(Plasmid('plasmid', []))))

        self.assertEqual('siteX>geneY',
                         change_to_string(Sub(Feature('siteX'), Feature('geneY'))))

        self.assertEqual('siteX>>geneY',
                         change_to_string(Sub(Feature('siteX'), Feature('geneY'), multiple=True)))


class GenotypeFusionsUpdateOnChangeTestCase(BaseTestCase):

    @SkipTest
    def test_fusion_insert_update_on_change(self):
        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneC')))
        }, self.chain('-geneA:geneB:geneC', '+geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneC'))),
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'), Feature(name='geneC')))
        }, self.chain('-geneA:geneB:geneC', '+geneA:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Feature(name='geneA'))
        }, self.chain('-geneA:geneB:geneC', '+geneB:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Feature(name='geneA'))
        }, self.chain('-geneA:geneB', '+geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual(set(), self.chain('-geneA:geneB', '+geneA:geneB',
                                           fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Ins(Fusion(Feature(name='geneB'), Feature(name='geneA')))
        }, self.chain('-geneA:geneB', '+geneB:geneA',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB'), FeatureSet(Feature(name='geneC'))))
        }, self.chain('-geneA:geneB:{geneC geneD}', '+geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB')))
        }, self.chain('-geneA:geneB:{geneC geneD}', '+geneC', '+geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), Feature(name='geneB')))
        }, self.chain('-geneA:geneB:{geneC geneD}', '+{geneC geneD}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneC'), Feature(name='geneD'))))
        }, self.chain('-geneA:geneB:{geneC geneB:geneD}', '+geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Del(Fusion(Feature(name='geneA'),
                       FeatureSet(Feature(name='geneC'), Fusion(Feature(name='geneB'), Feature(name='geneE')))))
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

        with self.assertRaises(AmbiguityError):
            self.chain('+geneA:geneB:geneC +geneA', '-geneA',
                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True)

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneB')))),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneC:geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneB'), Feature(name='geneD')))),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Fusion(Feature(name='geneC'), Feature(name='geneD'))))),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneB',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Feature(name='geneA')),
        }, self.chain('+geneA:{geneB geneC:geneD}', '-geneB', '-geneC:geneD',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        with self.assertRaises(AmbiguityError):
            self.chain('+geneA:geneB:geneC:geneA:geneB', '-geneA:geneB',
                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True)

    def test_fusion_replace_update_on_change(self):
        self.assertEqual({
            Sub(Feature(name='geneA'), Fusion(Feature(name='geneB'), Feature(name='geneC'))),
        }, self.chain('+geneA', 'geneA>geneB:geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Sub(Fusion(Feature(name='geneA'), Feature(name='geneB')), Feature(name='geneC')),
        }, self.chain('+geneA:geneB', 'geneA:geneB>geneC',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
            Sub(Fusion(Feature(name='geneA'), Feature(name='geneC')), Feature(name='geneD'))
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

        with self.assertRaises(AmbiguityError):
            self.chain('+geneA:geneB:{geneC geneD:geneE}:geneD:geneE', 'geneD:geneE>geneF',
                       fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True)

        self.assertEqual({
            Ins(Fusion(Feature(name='geneA'), FeatureSet(Feature(name='geneD'), Feature(name='geneE')))),
        }, self.chain('+geneA:geneB:geneC', 'geneB:geneC>{geneD geneE}',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))

        self.assertEqual({
            Ins(Fusion(Feature(name='geneE'), Feature(name='geneD'))),
        }, self.chain('+geneA:{geneB geneC}:geneD', 'geneA:{geneB geneC}>geneE',
                      fusion_strategy=Genotype.FUSION_UPDATE_ON_CHANGE).changes(fusions=True))
