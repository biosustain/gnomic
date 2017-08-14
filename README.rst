Gnomic
======

.. image:: https://travis-ci.org/biosustain/gnomic.svg?branch=master
    :target: https://travis-ci.org/biosustain/gnomic

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

Language grammar
----------------

The grammar consists of a list of genotype or phenotype designations, separated by
spaces and/or commas. The designations are described using the following nomenclature:

============================================================= ==================================
Designation                                                   Grammar expression
============================================================= ==================================
``feature`` deleted                                           ``-feature``
``feature`` at ``locus`` deleted                              ``-feature@locus``
``feature`` inserted                                          ``+feature``
``site`` replaced with ``feature``                            ``site>feature``
``site`` (multiple integration) replaced with ``feature``     ``site>>feature``
``site`` at ``locus`` replaced with ``feature``               ``site@locus>feature``
``feature`` of ``organism``                                   ``organism/feature``
``feature`` with ``type``                                     ``type.feature``
``feature`` with variant                                      ``feature(variant)``
``feature`` with list of variants                             ``feature(var1, var2)`` or ``feature(var1; var2)``
``feature`` with accession number                             ``feature#GB:123456``
``feature`` by accession number                               ``#GB:123456``
accession number                                              ``#database:id`` or ``#id``
fusion of ``feature1`` and ``feature2``                       ``feature1:feature2``
insertion of two fused features                               ``+feature1:feature2``
insertion of a list of features or fusions                    ``+{..insertables}``
fusion of a list and a feature                                ``{..insertables}:feature``
a non-integrated plasmid                                      ``(plasmid)`` or ``(plasmid ...insertables)``
integrated plasmid vector with required insertion site        ``site>(vector ..insertables)``
============================================================= ==================================


Feature variants
^^^^^^^^^^^^^^^^

Features may have one or more variants, separated by colon ";" or comma ",".

For example: ``geneX(cold-resistant; heat-resistant)``

Variants can either be identifiers (using the characters a-z, 0-9, "-" and "_") or be sequence variants following
the HGVS `Sequence Variant Nomenclature <http://www.hgvs.org/varnomen>`_.

For example: ``geneY(c.123G>T)``


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


Development
-----------

To rebuild the gnomic parser using `grako` (version 3.18.1), run:

::

    grako gnomic-grammar/genotype.enbf -o gnomic/grammar.py -m Gnomic
    
References
-----------

- `Wikipedia — Bacterial genetic nomenclature <http://en.wikipedia.org/wiki/Bacterial_genetic_nomenclature>`_
- `Journal of Bacteriology — Instructions to Authors <http://jb.asm.org/site/misc/journal-ita_nom.xhtml#03>`_
- `Human Genome Variation Society — Recommendations for the description of sequence variants <http://www.hgvs.org/mutnomen/recs.html>`_
- `Databases cross-referenced in UniProtKB <http://www.uniprot.org/docs/dbxref>`_


