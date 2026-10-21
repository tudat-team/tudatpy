.. _coma_model_data:

``coma_model``
==============
This module exposes data containers and file-processing utilities for coma
atmosphere and coma wind models. The settings objects that consume these data
containers are defined in :ref:`atmosphere`.

Functions
---------
.. currentmodule:: tudatpy.data.coma_model

.. autosummary::

   coma_model_file_processor
   coma_wind_file_processor

.. autofunction:: tudatpy.data.coma_model.coma_model_file_processor

.. autofunction:: tudatpy.data.coma_model.coma_wind_file_processor


Classes
-------
.. currentmodule:: tudatpy.data.coma_model

.. autosummary::

   ComaPolyDataset
   ComaStokesDataset
   ComaWindDatasetCollection
   ComaModelFileProcessor
   ComaWindModelFileProcessor

.. autoclass:: tudatpy.data.coma_model.ComaPolyDataset
   :members:

.. autoclass:: tudatpy.data.coma_model.ComaStokesDataset
   :members:

.. autoclass:: tudatpy.data.coma_model.ComaWindDatasetCollection
   :members:

.. autoclass:: tudatpy.data.coma_model.ComaModelFileProcessor
   :members:

.. autoclass:: tudatpy.data.coma_model.ComaWindModelFileProcessor
   :members:
