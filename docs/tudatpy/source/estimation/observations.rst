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

User-defined numerical metadata is stored separately as
``dataset.numerical_observation_metadata[key][observation_id]``. Use
``add_numerical_observation_metadata(key, observation_ids, values)`` to attach
values and ``get_numerical_observation_metadata(key, observation_ids)`` to
retrieve a vector in an explicit observation-ID order. Sorting and rejection
retain these values; copying and transfers preserve them, remapping IDs when
needed. Physical removal deletes the corresponding metadata entries. Metadata
does not change weights or observation models.

``dataset.reject_observations(condition, sigma_limit=5.0, reason="...")`` rejects
an active selected observation when any residual component exceeds five times
its measurement uncertainty. Residuals and weights remain unchanged. For
correlated observations, the uncertainty includes the full inverse set weight
matrix. Omitting ``sigma_limit`` retains selection-only rejection.

``ObservationVectorData.inverse_weight_matrix_for_observation(observation_id)``
returns a measurement covariance block using the complete correlated set weights,
including rejected observations. This uses the existing covariance cache.

Translational ephemerides created or reset from tracking supplementary data use
linear interpolation with silent extrapolation outside the supplied history.
Their original frame origin is preserved, including Earth's origin for Gaia.

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
   /estimation/observations/observation_corrections
