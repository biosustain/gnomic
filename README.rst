Gnomic
======

*(A Python implementation of gnomic, currently in development.)*

Gnomic is a human– and computer–readable representation of microbial genotypes and phenotypes. The ``gnomic``
Python package contains a parser for the Gnomic grammar able to interpret changes over multiple generations.

The first formal guidelines for microbial genetic nomenclature were drawn up in the 1960s and 70s. They are often hopelessly
ambiguous and not useful for modern computer-assisted genome engineering. The Gnomic grammar is an improvement over existing nomenclatures
designed to be clear, unambiguous and computer–readable and describe genotypes at various levels of details.

A JavaScript (Node) version of the package is available on NPM as `gnomic-grammar <https://www.npmjs.com/package/gnomic-grammar>`_.

Example usage
-------------

In this example, we parse *"geneA ΔsiteA::promoterB:geneB ΔgeneC"* and *"ΔgeneA"* in *gnomic* syntax:

::

   >>> from gnomic import *
   >>> g1 = Genotype.parse('+E.coli/geneA siteA>P.promoterB:E.coli/geneB -geneC')
   >>> g1.added_features
   (Feature(organism=Organism('Escherichia coli'), name='geneB'),
    Feature(organism=Organism('Escherichia coli'), name='geneA'),
    Feature(type=Type('promoter'), name='promoterB'))
   >>> g1.removed_features
   (Feature(name='geneC'),
    Feature(name='siteA'))
   >>> g2 = Genotype.parse('-geneA', parent=g1)
   >>> g2.added_features
   (Feature(type=Type('promoter'), name='promoterB'),
    Feature(name='geneB', organism=Organism('Escherichia coli')))
   >>> g2.removed_features
   (Feature(name='siteA'),
    Feature(name='geneC'),
    Feature(name='geneA'))

Language grammar
----------------

The grammar consists of a list of genotype or phenotype designations, separated by
spaces and/or commas. The designations are described using the following nomenclature:

============================================================= ==================================
 Designation                               Grammar expression
============================================================= ==================================
 ``feature`` deleted                                          ``-feature``
 ``feature`` inserted                                         ``+feature``
 ``feature`` inserted at ``site``                             ``+site::feature``
 ``site`` replaced with ``feature``                           ``site>feature``
 ``site`` (multiple integration) replaced with ``feature``     ``site>>feature``
 ``feature`` of ``organism``                                   ``organism/feature``
 ``feature`` with variant                                      ``feature(variant)``
 ``feature`` with accession number                             ``feature#GB:123456``
 ``feature`` by accession number                              ``#GB:123456``
 fusion of ``feature1`` and ``feature2``                      ``feature1:feature2``
 insertion of two fused features                              ``+feature1:feature2``
 insertion of a list of features or fusions                    ``+{..insertables}``
 phenotype: wild-type                                         ``phene+`` or ``phene(wild-type)``
 phenotype: mutant                                            ``phene-`` or ``phene(mutant)``
 selection marker: used (wild-type)                           ``marker+``
 selection marker: available (mutantt)                        ``marker-``
 a non-integrated plasmid                                     ``plasmid{}`` or ``plasmid{..insertables}``
 plasmid with selection marker                                ``plasmid{..insertables}::marker+``
 integrated plasmid vector with required insertion site       ``site>vector{..insertables}``
 nucleotide range of a ``feature``                            ``feature[startBase_endBase]``
 coding nucleotide range of a ``gene``                        ``gene[c.startBase_endBase]``
 protein amino-acid range of a ``gene``                       ``gene[p.startAA_endAA]``
 protein amino-acid of a ``gene``                             ``gene[p.AA]``
============================================================= ==================================

Development
-----------

To rebuild the gnomic parser using `grako`, run:

::

    grako genotype.enbf -o gnomic/grammar.py -m Gnomic