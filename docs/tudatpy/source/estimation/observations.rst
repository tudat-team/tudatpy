.. _observations:

``observations``
=====================

This module contains the objects in Tudat that are used to store and process observations. The :class:`~tudatpy.estimation.observations.ObservationDataset` is the object in which a complete set of observations (e.g. all observations to be used in a single estimation) is stored, together with their times, link ends, residuals, weights and dependent variables. Each observation is a measurement at a single time, which may have several components (e.g. right ascension and declination). Observations with a common observable type, link ends, etc. are grouped into observation sets. When using simulated data, an observation dataset is created with :func:`~tudatpy.estimation.observations.simulate_observation_dataset`; real data is converted from :ref:`tracking_data` objects with :func:`~tudatpy.estimation.observations.create_observation_dataset_from_tracking_data`. The residuals and dependent variables are normally computed during an estimation (see :ref:`estimation_analysis`), but can also be computed outside an estimation loop using the function :func:`~tudatpy.estimation.observations.compute_residuals_and_dependent_variables`.

Most functions of the :class:`~tudatpy.estimation.observations.ObservationDataset` act on the observations selected by a condition, which is created with :data:`~tudatpy.estimation.observations.observation_query` (see :ref:`below <observation_query>`). Such conditions are used to retrieve data (e.g. observation times, residuals or weights), to set weights, to reject observations (rejected observations remain in the dataset, but are not used in the estimation), or to create a new dataset with a subset of the observations. The weights may be diagonal, or may contain correlations between the components of an observation, between the observations in a set, or between observations in different sets.

The :class:`~tudatpy.estimation.observations.ObservationCollection` and :class:`~tudatpy.estimation.observations.SingleObservationSet` classes, and the functions that use them, are deprecated, but remain available (see :ref:`below <observation_collection_deprecated>`). Objects of these classes can be converted to and from an :class:`~tudatpy.estimation.observations.ObservationDataset`. Note that the conversion to an observation collection keeps only the diagonal of the weight matrix, and is not possible for a dataset that contains rejected observations.




.. toctree::
   :maxdepth: 2
   :caption: Modules

   /estimation/observations/observations_geometry
   /estimation/observations/observations_processing
   /estimation/observations/observation_corrections


Functions
---------
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   simulate_observation_dataset

   create_observation_dataset_from_tracking_data

   set_tracking_supplementary_data_in_bodies

   compute_residuals_and_dependent_variables

   observation_simulation_settings_from_dataset

   create_observation_dataset_from_arrays

   create_single_type_observation_dataset_from_arrays

   create_pseudo_observation_dataset_and_models

   create_pseudo_observation_dataset_and_models_from_observation_times

   create_compressed_doppler_dataset

   create_count_time_compressed_doppler_dataset


.. autofunction:: tudatpy.estimation.observations.simulate_observation_dataset

.. autofunction:: tudatpy.estimation.observations.create_observation_dataset_from_tracking_data

.. autofunction:: tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies

.. autofunction:: tudatpy.estimation.observations.compute_residuals_and_dependent_variables

.. autofunction:: tudatpy.estimation.observations.observation_simulation_settings_from_dataset

.. autofunction:: tudatpy.estimation.observations.create_observation_dataset_from_arrays

.. autofunction:: tudatpy.estimation.observations.create_single_type_observation_dataset_from_arrays

.. autofunction:: tudatpy.estimation.observations.create_pseudo_observation_dataset_and_models

.. autofunction:: tudatpy.estimation.observations.create_pseudo_observation_dataset_and_models_from_observation_times

.. autofunction:: tudatpy.estimation.observations.create_compressed_doppler_dataset

.. autofunction:: tudatpy.estimation.observations.create_count_time_compressed_doppler_dataset


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


.. _observation_query:

Observation query
-----------------

Observations in an :class:`~tudatpy.estimation.observations.ObservationDataset` are selected with conditions, created by comparing an attribute of ``observation_query`` to a value. For instance, the following retrieves the residuals of the observations received at the DSS-14 ground station that are not rejected:

.. code-block:: python

    from tudatpy.estimation.observable_models_setup import links
    from tudatpy.estimation.observations import observation_query

    station = links.body_reference_point_link_end_id("Earth", "DSS-14")
    condition = (observation_query.receiver == station) & observation_query.active
    residuals = observation_dataset.get_residuals(condition, ordering="estimation")

.. autodata:: tudatpy.estimation.observations.observation_query
   :annotation:


.. _observation_collection_deprecated:

Observation collections (deprecated)
------------------------------------

The :class:`~tudatpy.estimation.observations.ObservationCollection` (composed of a list of :class:`~tudatpy.estimation.observations.SingleObservationSet` objects, each of which stores observations of a common type, link ends, etc.) is deprecated, and replaced by the :class:`~tudatpy.estimation.observations.ObservationDataset`. The classes and functions below remain available for existing code. Observation collections are converted to and from an observation dataset with :func:`~tudatpy.estimation.observations.create_observation_dataset_from_collection` and :func:`~tudatpy.estimation.observations.create_observation_collection_from_dataset`. The observation sets in a collection are selected, filtered and split with the functionality in :ref:`observations_processing`.

Functions
^^^^^^^^^
.. currentmodule:: tudatpy.estimation.observations

.. autosummary::

   create_observation_dataset_from_collection

   create_observation_dataset_from_single_observation_set

   create_observation_collection_from_dataset

   create_single_observation_set_from_dataset

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

   create_compressed_doppler_collection

   create_count_time_compressed_doppler_collection


.. autofunction:: tudatpy.estimation.observations.create_observation_dataset_from_collection

.. autofunction:: tudatpy.estimation.observations.create_observation_dataset_from_single_observation_set

.. autofunction:: tudatpy.estimation.observations.create_observation_collection_from_dataset

.. autofunction:: tudatpy.estimation.observations.create_single_observation_set_from_dataset

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
