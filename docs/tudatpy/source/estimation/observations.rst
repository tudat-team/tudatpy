.. _observations:

``observations``
=====================

This module contains the objects in Tudat that are used to store and process observations. The :class:`~tudatpy.estimation.observations.ObservationDataset` is the object in which a complete set of observations (e.g. all observations to be used in a single estimation) is stored, together with their times, link ends, residuals, weights and dependent variables. Each observation is a measurement at a single time, which may have several components (e.g. right ascension and declination). Observations are stored in observation sets: all observations of a set have the same observable type, link ends, etc. (several sets may share these properties, since each simulation settings object, tracking data object or call of :meth:`~tudatpy.estimation.observations.ObservationDataset.add_observation_set` creates its own set). When using simulated data, an observation dataset is created with :func:`~tudatpy.estimation.observations.simulate_observation_dataset`; real data is converted from :ref:`tracking_data` objects with :func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data`. The residuals are computed during an estimation (see :ref:`estimation_analysis`). Dependent variables are computed when simulating observations, or for any dataset with :func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`, which also recomputes the residuals.

Most functions of the :class:`~tudatpy.estimation.observations.ObservationDataset` act on the observations selected by a condition, which is created with :data:`~tudatpy.estimation.observations.observation_query` (see :ref:`below <observation_query>`). Such conditions are used to retrieve data (e.g. observation times, residuals or weights), to set weights, to reject observations (rejected observations remain in the dataset, but are not used in the estimation), or to create a new dataset with a subset of the observations. The weights may be diagonal, or may contain correlations between the components of an observation, between the observations in a set, or between observations in different sets.

The :class:`~tudatpy.estimation.observations.ObservationCollection` and :class:`~tudatpy.estimation.observations.SingleObservationSet` classes, and the functions that use them, are deprecated, but remain available (see :ref:`below <observation_collection_deprecated>`). Objects of these classes can be converted to and from an :class:`~tudatpy.estimation.observations.ObservationDataset`. Note that the conversion to an observation collection keeps only the diagonal of the weight matrix, and is not possible for a dataset that contains rejected observations.




.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observations/observations_geometry
   /estimation/observations/observations_processing


Functions
---------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::




   compute_residuals_and_dependent_variables





















Classes
-------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   ObservationDataset

   ObservationSelectionCondition

   ObservationSelectionConditionType

   ObservationVectorData

   ObservationSetMetadata

   ObservationDatasetRow

   ObservationScalarComponentRow

.. autoclass:: tudatpy.estimation.observations.ObservationDataset
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.estimation.observations.ObservationSelectionCondition
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.estimation.observations.ObservationSelectionConditionType
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationVectorData
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationSetMetadata
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationDatasetRow
   :members:

.. autoclass:: tudatpy.estimation.observations.ObservationScalarComponentRow
   :members:




   create_single_observation_set

   create_filtered_observation_set

   split_observation_set

   create_filtered_observation_collection

   merge_observation_collections

   split_observation_collection

   create_new_observation_collection

   create_observation_collection_from_tracking_data

   create_count_time_compressed_doppler_collection

   create_observation_collection_from_arrays

   create_single_type_observation_collection_from_arrays

   simulate_observations

   simulate_pseudo_observations

   set_tracking_supplementary_data_in_bodies
   create_compressed_doppler_collection



.. autofunction:: tudatpy.estimation.observations.compute_residuals_and_dependent_variables



.. autofunction:: tudatpy.estimation.observations.create_single_observation_set

.. autofunction:: tudatpy.estimation.observations.create_filtered_observation_set

.. autofunction:: tudatpy.estimation.observations.split_observation_set

.. autofunction:: tudatpy.estimation.observations.merge_observation_collections

.. autofunction:: tudatpy.estimation.observations.create_filtered_observation_collection

.. autofunction:: tudatpy.estimation.observations.split_observation_collection

.. autofunction:: tudatpy.estimation.observations.create_new_observation_collection

.. autofunction:: tudatpy.estimation.observations.create_observation_collection_from_tracking_data

.. autofunction:: tudatpy.estimation.observations.create_observation_collection_from_arrays

.. autofunction:: tudatpy.estimation.observations.create_single_type_observation_collection_from_arrays

.. autofunction:: tudatpy.estimation.observations.simulate_observations

.. autofunction:: tudatpy.estimation.observations.simulate_pseudo_observations

.. autofunction:: tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies
.. autofunction:: tudatpy.estimation.observations.create_compressed_doppler_collection

.. autofunction:: tudatpy.estimation.observations.create_count_time_compressed_doppler_collection


Classes
^^^^^^^
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   SingleObservationSet

   ObservationCollection

.. autoclass:: tudatpy.estimation.observations.SingleObservationSet
   :members:
   :special-members: __init__

.. autoclass:: tudatpy.estimation.observations.ObservationCollection
   :members:
   :special-members: __init__
