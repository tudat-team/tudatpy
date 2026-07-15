.. _environment_data_coma:

``coma``
========

This modules contains functionality to load files to be used for the definition of comet coma density and wind modelling, using the
specific implementation of TODO-AI (Reichel et al. 2026, to be submitted, other details pending, add to bib). TODO-AI: add links to specific
environment_setup functions where these are to be used.

.. automodule:: tudatpy.data_input.environment_data.coma

Functions
---------
.. currentmodule:: tudatpy.data_input.environment_data.coma

.. autosummary::

   coma_model_file_processor
   coma_wind_file_processor

.. autofunction:: coma_model_file_processor
.. autofunction:: coma_wind_file_processor

Classes
-------

.. autosummary::

   ComaModelFileProcessor
   ComaWindModelFileProcessor
   ComaPolyDataset
   ComaStokesDataset
   ComaWindDatasetCollection

.. autoclass:: ComaModelFileProcessor
   :members:

.. autoclass:: ComaWindModelFileProcessor
   :members:

.. autoclass:: ComaPolyDataset
   :members:

.. autoclass:: ComaStokesDataset
   :members:

.. autoclass:: ComaWindDatasetCollection
   :members:
