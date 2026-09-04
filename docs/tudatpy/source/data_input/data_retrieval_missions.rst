.. _data_retrieval_missions:

``missions``
============

.. toctree::
   :maxdepth: 2
   :caption: Modules

   data_retrieval_missions/cassini
   data_retrieval_missions/grail
   data_retrieval_missions/juice
   data_retrieval_missions/mex
   data_retrieval_missions/mro
   data_retrieval_missions/rosetta

.. automodule:: tudatpy.data_input.data_retrieval.missions

This submodule contains mission-archive download and discovery utilities. The
:class:`LoadPDS` class is the main mission-data retrieval interface; mission
specific support is implemented in the per-mission submodules listed above.

Classes
-------
.. currentmodule:: tudatpy.data_input.data_retrieval.missions

.. autosummary::

   LoadPDS

.. autoclass:: LoadPDS
   :members:
