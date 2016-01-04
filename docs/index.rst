.. Gnomic documentation master file, created by
   sphinx-quickstart on Mon Jan  4 11:40:21 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Gnomic's documentation!
==================================

Gnomic is a human-- and computer--readable representation of microbial genotypes and phenotypes. The ``gnomic``
Python package contains a parser for the Gnomic grammar able to interpret changes over multiple generations.

The first formal guidelines for microbial genetic nomenclature were drawn up in the 1960s and 70s. They are often hopelessly
ambiguous and not useful for modern computer-assisted genome engineering. The Gnomic grammar is an improvement over existing nomenclatures
designed to be clear, unambiguous and computer-readable and describe genotypes at various levels of details.

A JavaScript (Node) version of the package is available on NPM as `gnomic-grammar <https://www.npmjs.com/package/gnomic-grammar>`_.

Example usage
-------------

In this example, we parse `"geneA ΔsiteA::promoterB:geneB ΔgeneC"` and `"ΔgeneA"` in *gnomic* syntax:


.. code-block:: python

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


User's guide
------------

.. toctree::
   :maxdepth: 2

   grammar

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
