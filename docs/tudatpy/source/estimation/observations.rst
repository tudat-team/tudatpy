.. _observations:

``observations``
================

This module provides the objects and functions used to store, inspect, select,
and process observations for estimation. The central container is
:class:`~tudatpy.estimation.observations.ObservationDataset`, which stores
observation events together with their link definitions, times, residuals,
dependent variables, and weights. Observation events with common metadata are
grouped into logical sets, while vector observables remain single events with
multiple scalar components.

Selections are expressed with
:data:`~tudatpy.estimation.observations.observation_query` and may be used to
inspect data, assign weights, reject or restore observations, and create reduced
datasets. A detailed description of the observation-data workflow will be
provided on the `ObservationDataset user guide page <https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-simulation/observation-dataset.html>`_.

Functions
---------

.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   create_observation_dataset_from_tracking_data
   create_observation_dataset_from_arrays
   create_single_type_observation_dataset_from_arrays
   create_pseudo_observation_dataset_and_models
   create_pseudo_observation_dataset_and_models_from_observation_times
   simulate_observation_dataset
   create_compressed_doppler_dataset
   observation_simulation_settings_from_dataset
   compute_residuals_and_dependent_variables
   set_tracking_supplementary_data_in_bodies

.. autofunction:: create_observation_dataset_from_tracking_data
.. autofunction:: create_observation_dataset_from_arrays
.. autofunction:: create_single_type_observation_dataset_from_arrays
.. autofunction:: create_pseudo_observation_dataset_and_models
.. autofunction:: create_pseudo_observation_dataset_and_models_from_observation_times
.. autofunction:: simulate_observation_dataset
.. autofunction:: create_compressed_doppler_dataset
.. autofunction:: observation_simulation_settings_from_dataset
.. autofunction:: compute_residuals_and_dependent_variables
.. autofunction:: set_tracking_supplementary_data_in_bodies

Classes
-------

.. autosummary::

   ObservationDataset
   ObservationSelectionCondition
   ObservationSelectionConditionType
   ObservationVectorData
   ObservationSetMetadata
   ObservationDatasetRow
   ObservationScalarComponentRow
   ObservationWeightSettings

.. autoclass:: ObservationDataset
   :members:
   :special-members: __init__

.. autoclass:: ObservationSelectionCondition
   :members:

.. autoclass:: ObservationSelectionConditionType
   :members:

.. autoclass:: ObservationVectorData
   :members:

.. autoclass:: ObservationSetMetadata
   :members:

.. autoclass:: ObservationDatasetRow
   :members:

.. autoclass:: ObservationScalarComponentRow
   :members:

.. autoclass:: ObservationWeightSettings
   :members:

Observation query
-----------------

Use :data:`observation_query` to build row-level
:class:`ObservationSelectionCondition` objects. Conditions may be combined with
``&`` and ``|`` and negated with ``~``. Parenthesize comparisons, and do not use
Python's ``and``, ``or``, or ``not`` operators with conditions.

.. data:: observation_query

   Query builder exposing selectors for observable type, link definition, link
   ends, set ID, time, active/rejected status, observations, residuals, and
   dependent variables.

Legacy compatibility
--------------------

:class:`SingleObservationSet` and :class:`ObservationCollection` remain
available as compatibility facades. New code should use
:class:`ObservationDataset` and condition-based selection.

.. autoclass:: SingleObservationSet
   :members:

.. autoclass:: ObservationCollection
   :members:

Submodules
----------

.. toctree::
   :maxdepth: 1

   /estimation/observations/observations_geometry
