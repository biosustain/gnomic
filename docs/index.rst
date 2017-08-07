.. Gnomic documentation master file, created by
   sphinx-quickstart on Mon Jan  4 11:40:21 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Gnomic's documentation!
==================================

Gnomic is a human-- and computer--readable representation of microbial genotypes and phenotypes. The ``gnomic``
Python package contains a parser for the Gnomic grammar able to interpret changes over multiple generations.

The first formal guidelines for microbial genetic nomenclature were drawn up in the 1960s. These traditional nomenclatures are too
ambiguous to be useful for modern computer-assisted genome engineering. The *gnomic* grammar is an improvement over existing nomenclatures
designed to be clear, unambiguous and computer–readable and describe genotypes at various levels of detail.

A JavaScript (Node) version of the package is available on NPM as `gnomic-grammar <https://www.npmjs.com/package/gnomic-grammar>`_.

Example usage
-------------

In this example, we parse `"EcGeneA ΔsiteA::promoterB:EcGeneB ΔgeneC"` and `"ΔgeneA"` in *gnomic* syntax:

.. code-block:: python

   >>> from gnomic import Genotype
   >>> g1 = Genotype.parse('+Ec/geneA(variant) siteA>P.promoterB:Ec/geneB -geneC')
   >>> g1.added_features
   {Feature(organism=Organism('Ec'), name='geneA', variant=('variant',)),
    Feature(organism=Organism('Ec'), name='geneB'),
    Feature(type='P', name='promoterB')}
   >>> g1.removed_features
   {Feature(name='geneC'),
    Feature(name='siteA')}

   >>> g2 = Genotype.parse('-geneA', parent=g1)
   >>> g2.added_features
   {Feature(type='P', name='promoterB'),
    Feature(name='geneB', organism='Ec')}
   >>> g2.removed_features
   {Feature(name='siteA'),
    Feature(name='geneC')}
    >>> g2.changes()
    (Change(multiple=False,
            after=Fusion(annotations=(Feature(type='P', name='promoterB'), Feature(organism='Ec', name='geneB'))),
            before=Feature(name='siteA')),
     Change(multiple=False, before=Feature(name='geneC')))

    >>> g2.format()
    'ΔsiteA P.promoterB:Ec/geneB ΔgeneC'


User's guide
------------

.. toctree::
   :maxdepth: 2

   grammar
   models

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
