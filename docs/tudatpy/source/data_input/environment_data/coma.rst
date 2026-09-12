.. _environment_data_coma:

``coma``
========

This module contains functionality for loading comet-coma density and wind data
from files. The loaded data are used to define coma atmosphere and wind model
settings through :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_from_poly_data`,
:func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_model_from_stokes_data`,
and :func:`~tudatpy.dynamics.environment_setup.atmosphere.coma_wind_model`.

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
