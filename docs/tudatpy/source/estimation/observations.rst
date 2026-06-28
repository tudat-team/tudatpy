.. _observations:

``observations``
=====================

This module contains the objects in Tudat that are used to store and process observations. The :class:`~tudatpy.estimation.observations.ObservationDataset` is the backend storage object for observations, residuals, weights, dependent variables, metadata and row-level selection state. Dealing with observations in Tudat is discussed extensively on the `user guide <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation.html>`_.




.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observations/observations_geometry


Functions
---------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   compute_residuals_and_dependent_variables_for_dataset

   observation_simulation_settings_from_dataset


.. autofunction:: tudatpy.estimation.observations.compute_residuals_and_dependent_variables_for_dataset

.. autofunction:: tudatpy.estimation.observations.observation_simulation_settings_from_dataset


Classes
-------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   ObservationDataset

   ObservationDatasetViewer

   ObservationCondition

   FlattenedObservationData

   ObservationSetMetadata

   ObservationDatasetRow

   ObservationScalarComponentRow

   ObservationWeightBlock

.. autoclass:: tudatpy.estimation.observations.ObservationDataset
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationDatasetViewer
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationCondition
   :members:

.. autoclass:: tudatpy.estimation.observations.FlattenedObservationData
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationSetMetadata
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationDatasetRow
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationScalarComponentRow
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationWeightBlock
   :members:
