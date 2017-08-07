
Language grammar
----------------

The grammar consists of a list of genotype or phenotype designations, separated by
spaces and/or commas. The designations are described using the following nomenclature:

============================================================= =======================================================
 Designation                                                  Grammar expression
============================================================= =======================================================
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
fusion of a list and a feature                                ``{..insertables}:feature
a non-integrated plasmid                                      ``(plasmid)`` or ``(plasmid ...insertables)``
integrated plasmid vector with required insertion site        ``site>(vector ..insertables)``
============================================================= =======================================================