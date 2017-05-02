# -*- coding: utf-8 -*-
from __future__ import unicode_literals
from unittest import TestCase, SkipTest

from gnomic import Genotype, Feature, Ins, Del, Fusion, Sub, Type, Range, Plasmid, FeatureTree, Organism, FeatureSet, \
    Accession
from gnomic.utils import genotype_to_text, feature_to_text, genotype_to_string, change_to_string


class BaseTestCase(TestCase):
    def chain(self, *definitions, **kwargs):
        return Genotype.chain_parse(list(definitions), **kwargs)


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

    def test_deletion_of_plasmid_vector(self):
        self.assertEqual({
            Plasmid('pA', [])
        }, self.chain('pA{}').changes())

        self.assertEqual({
            Del(Plasmid('pA', []))
        }, self.chain('-pA{}').changes())

        self.assertEqual(set(), self.chain('pA{}', '-pA{}').changes())
        self.assertEqual(set(), self.chain('pA{foo}', '-pA{}').changes())

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

    def test_insertion_followed_by_deletion(self):
        self.assertEqual({
            Ins(Feature(name='geneX')),
        }, self.chain('+geneX +geneA', '-geneA').changes())

        self.assertEqual({
            Ins(Feature(name='geneX')),
            Ins(Feature(name='geneY')),
        }, self.chain('+geneX +geneA', '+geneY', '-geneA').changes())

    def test_deletion_followed_by_insertion(self):
        self.assertEqual({
            Ins(Feature(name='geneX')),
        }, self.chain('+geneX -geneA', '+geneA').changes())

        self.assertEqual({
            Ins(Feature(name='geneX')),
            Ins(Feature(name='geneY')),
        }, self.chain('+geneX -geneA', '+geneY', '+geneA').changes())

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
            Del(Feature(name='siteA')),
            Ins(Fusion(Feature(name='geneA'), Feature(name='geneB'))),
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

    def test_markers(self):
        self.assertEqual(genotype_to_text(self.chain('(p)::m+')), '(p)::m⁺')
        self.assertEqual(genotype_to_text(self.chain('(p)::{m+ n+}')), '(p)::{m⁺ n⁺}')

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

    def test_feature_tree(self):
        feature = FeatureTree(Feature(name="foo"), Feature(name="bar"))
        self.assertEqual(feature_to_text(feature), "{foo bar}")

    def test_feature_set(self):
        feature = FeatureSet(Feature(name="foo"), Feature(name="bar"))
        self.assertEqual(feature_to_text(feature), "{foo bar}")

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

    def test_feature_is_marker(self):
        feature = Feature(name="foo")
        self.assertEqual(feature_to_text(feature, is_marker=True), "::foo")

    def test_feature_with_accession(self):
        feature = Feature(name="foo", accession=Accession(identifier='bar', database='database'), type=Type('type'))
        self.assertEqual(feature_to_text(feature), "type.foo#database:bar")


class GenotypeToStringTestCase(BaseTestCase):

    def test_genotype_to_string(self):
        gnomic = genotype_to_string(self.chain('-e.coli/geneA',
                                               '+geneB(a)',
                                               '(vectorC geneC)',
                                               '-(vectorD)',
                                               '+geneE::markerF+',
                                               '+gene.G::{markerH+, markerI+}',
                                               '-Ecoli/reaction.PGK#bigg:PGK'))
        self.assertEqual(
            set(gnomic.split()),
            set('+geneE '
                '+markerH+ '
                '(vectorC geneC) '
                '-e.coli/geneA '
                '+markerF+ '
                '+geneB(a) '
                '+markerI+ '
                '+gene.G '
                '-(vectorD) '
                '-Ecoli/reaction.PGK#bigg:PGK'.split())
        )

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
                         change_to_string(Ins([Feature('geneA'), Feature('geneB')])))

        self.assertEqual('(plasmid)',
                         change_to_string(Plasmid('plasmid', [])))

        self.assertEqual('-(plasmid)',
                         change_to_string(Del(Plasmid('plasmid', []))))

        self.assertEqual('siteX>geneY',
                         change_to_string(Sub(Feature('siteX'), Feature('geneY'))))

        self.assertEqual('siteX>>geneY',
                         change_to_string(Sub(Feature('siteX'), Feature('geneY'), multiple=True)))
