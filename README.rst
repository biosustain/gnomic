Gnomic
======

.. image:: https://travis-ci.org/biosustain/gnomic-python.svg?branch=master
    :target: https://travis-ci.org/biosustain/gnomic-python

.. image:: https://zenodo.org/badge/47830031.svg
   :target: https://zenodo.org/badge/latestdoi/47830031

Gnomic is a human– and computer–readable representation of microbial genotypes and phenotypes. The ``gnomic``
Python package contains a parser for the Gnomic grammar able to interpret changes over multiple generations.

The first formal guidelines for microbial genetic nomenclature were drawn up in the 1960s. These traditional nomenclatures are too
ambiguous to be useful for modern computer-assisted genome engineering. The *gnomic* grammar is an improvement over existing nomenclatures, designed to be clear, unambiguous, computer–readable and describe genotypes at various levels of detail.

A JavaScript (Node) version of the package is available on NPM as `gnomic-grammar <https://www.npmjs.com/package/gnomic-grammar>`_.

Installation
------------

.. code-block:: bash

    pip install gnomic

Example usage
-------------

In this example, we parse `"EcGeneA ΔsiteA::promoterB:EcGeneB ΔgeneC"` and `"ΔgeneA"` in *gnomic* syntax:

.. code-block:: python

   >>> from gnomic import *
   >>> g1 = Genotype.parse('+Ec/geneA siteA>P.promoterB:Ec/geneB -geneC')
   >>> g1.added_features
   (Feature(organism=Organism('Escherichia coli'), name='geneB'),
    Feature(organism=Organism('Escherichia coli'), name='geneA'),
    Feature(type=Type('promoter'), name='promoterB'))
   >>> g1.removed_features
   (Feature(name='geneC'),
    Feature(name='siteA'))
   >>> g1.raw
   (Mutation(new=FeatureTree(Feature(organism=Organism('Escherichia coli'), name='geneA'))),
    Mutation(old=FeatureTree(Feature(name='siteA')),
             new=FeatureTree(Fusion(Feature(type=Type('promoter'), name='promoterB'),
                                    Feature(organism=Organism('Escherichia coli'), name='geneB')))),
    Mutation(old=FeatureTree(Feature(name='geneC'))))
   >>>
   >>> g2 = Genotype.parse('-geneA', parent=g1)
   >>> g2.added_features
   (Feature(type=Type('promoter'), name='promoterB'),
    Feature(name='geneB', organism=Organism('Escherichia coli')))
   >>> g2.removed_features
   (Feature(name='siteA'),
    Feature(name='geneC'),
    Feature(name='geneA'))
    >>> g2.changes()
    {Mutation(old=FeatureTree(Feature(name='siteA'))), 
     Mutation(new=FeatureTree(Feature(name='promoterB', type=Type('promoter')))), 
     Mutation(new=FeatureTree(Feature(organism=Organism('Escherichia coli'), name='geneB'))), 
     Mutation(old=FeatureTree(Feature(name='geneC'))),
     Mutation(old=FeatureTree(Feature(name='geneA')))}
    >>> g2.changes(fusions=True)
    {Mutation(old=FeatureTree(Feature(name='siteA'))), 
     Mutation(new=FeatureTree(Fusion(Feature(name='promoterB', type=Type('promoter')), 
                                     Feature(organism=Organism('Escherichia coli'), name='geneB')))),
     Mutation(old=FeatureTree(Feature(name='geneC'))),
     Mutation(old=FeatureTree(Feature(name='geneA')))}


Language grammar
----------------

The grammar consists of a list of genotype or phenotype designations, separated by
spaces and/or commas. The designations are described using the following nomenclature:

============================================================= ==================================
Designation                                                   Grammar expression
============================================================= ==================================
``feature`` deleted                                           ``-feature``
``feature`` inserted                                          ``+feature``
``site`` replaced with ``feature``                            ``site>feature``
``site`` (multiple integration) replaced with ``feature``     ``site>>feature``
``feature`` of ``organism``                                   ``organism/feature``
``feature`` with variant                                      ``feature(variant)``
``feature`` with accession number                             ``feature#GB:123456``
``feature`` by accession number                               ``#GB:123456``
fusion of ``feature1`` and ``feature2``                       ``feature1:feature2``
insertion of two fused features                               ``+feature1:feature2``
insertion of a list of features or fusions                    ``+{..insertables}``
phenotype: wild-type                                          ``phene+`` or ``phene(wild-type)``
phenotype: mutant                                             ``phene-`` or ``phene(mutant)``
selection marker: used (wild-type)                            ``marker+``
selection marker: available (mutant)                          ``marker-``
a non-integrated plasmid                                      ``plasmid{}`` or ``plasmid{..insertables}``
plasmid with single selection marker                          ``plasmid{..insertables}::marker+``
plasmid with multiple selection markers                       ``plasmid{..insertables}::{markerA+ markerB+}``
integrated plasmid vector with required insertion site        ``site>vector{..insertables}``
nucleotide range of a ``feature``                             ``feature[startBase_endBase]``
coding nucleotide range of a ``gene``                         ``gene[c.startBase_endBase]``
protein amino-acid range of a ``gene``                        ``gene[p.startAA_endAA]``
protein amino-acid of a ``gene``                              ``gene[p.AA]``
============================================================= ==================================

Development
-----------

To rebuild the gnomic parser using `grako`, run:

::

    grako genotype.enbf -o gnomic/grammar.py -m Gnomic
