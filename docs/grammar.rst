
Language grammar
----------------

The grammar consists of a list of genotype or phenotype designations, separated by
spaces and/or commas. The designations are described using the following nomenclature:

============================================================= =======================================================
 Designation                                                  Grammar expression
============================================================= =======================================================
 ``feature`` deleted                                          ``-feature``
 ``feature`` inserted                                         ``+feature``
 ``feature`` inserted at ``site``                             ``+site::feature``
 ``site`` replaced with ``feature``                           ``site>feature``
 ``site`` (multiple integration) replaced with ``feature``    ``site>>feature``
 ``feature`` of ``organism``                                  ``organism/feature``
 ``feature`` with variant                                     ``feature(variant)``
 ``feature`` with accession number                            ``feature#GB:123456``
 ``feature`` by accession number                              ``#GB:123456``
 fusion of ``feature1`` and ``feature2``                      ``feature1:feature2``
 insertion of two fused features                              ``+feature1:feature2``
 insertion of a list of features or fusions                   ``+{..insertables}``
 phenotype: wild-type                                         ``phene+`` or ``phene(wild-type)``
 phenotype: mutant                                            ``phene-`` or ``phene(mutant)``
 selection marker: used (wild-type)                           ``marker+``
 selection marker: available (mutant)                         ``marker-``
 a non-integrated plasmid                                     ``plasmid{}`` or ``plasmid{..insertables}``
 plasmid with single selection marker                         ``plasmid{..insertables}::marker+``
 plasmid with multiple selection markers                      ``plasmid{..insertables}::{markerA+ markerB+}``
 integrated plasmid vector with required insertion site       ``site>vector{..insertables}``
 nucleotide range of a ``feature``                            ``feature[startBase_endBase]``
 coding nucleotide range of a ``gene``                        ``gene[c.startBase_endBase]``
 protein amino-acid range of a ``gene``                       ``gene[p.startAA_endAA]``
 protein amino-acid of a ``gene``                             ``gene[p.AA]``
============================================================= =======================================================