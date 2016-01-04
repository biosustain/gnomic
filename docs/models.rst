Gnomic API
==========

Genotype
--------

.. module:: gnomic


.. autoclass:: Genotype
    :members:


.. py:attribute:: DEFAULT_ORGANISMS

    The default set of known organisms used when parsing.

.. py:attribute:: DEFAULT_TYPES

    The default set of known :class:`Feature` type modifiers used when parsing.


Mutations
---------

.. autoclass:: Mutation
    :members:

.. py:function:: Ins(insert, **kwargs)

    A shortcut for `Mutation(None, insert, **kwargs)`

.. py:function:: Sub(before, after, **kwargs)

    A shortcut for `Mutation(before, after, **kwargs)`

.. py:function:: Del(delete, **kwargs)

    A shortcut for `Mutation(delete, None, **kwargs)`

Features
--------

.. autoclass:: Feature
    :members:

.. autoclass:: FeatureTree
    :members:

.. autoclass:: Plasmid
    :members:

.. autoclass:: Fusion
    :members:

.. autoclass:: Organism
    :members:

.. autoclass:: Type
    :members:
