.. _observations:

``observations``
================

This module stores and manipulates observation data for estimation. The central
object is :class:`~tudatpy.estimation.observations.ObservationDataset`.
Selections are built with :data:`~tudatpy.estimation.observations.observation_query` and applied
through condition-based methods.

Minimal example
---------------

.. code-block:: python

   from tudatpy.estimation import observations
   observation_query = observations.observation_query

   viewer = dataset.create_viewer(
       condition=(
           (observation_query.observable_type == observations.one_way_doppler)
           & observation_query.time.between(t0, t1)
           & ~observation_query.rejected
       )
   )

Functions
---------

.. autofunction:: tudatpy.estimation.observations.create_observation_dataset_from_tracking_data
.. autofunction:: tudatpy.estimation.observations.create_observation_dataset_from_arrays
.. autofunction:: tudatpy.estimation.observations.create_single_type_observation_dataset_from_arrays
.. autofunction:: tudatpy.estimation.observations.create_pseudo_observation_dataset_and_models
.. autofunction:: tudatpy.estimation.observations.create_pseudo_observation_dataset_and_models_from_observation_times
.. autofunction:: tudatpy.estimation.observations.simulate_observation_dataset
.. autofunction:: tudatpy.estimation.observations.create_compressed_doppler_dataset
.. autofunction:: tudatpy.estimation.observations.set_tracking_supplementary_data_in_bodies
.. autofunction:: tudatpy.estimation.observations.compute_residuals_and_dependent_variables_for_dataset
.. autofunction:: tudatpy.estimation.observations.observation_simulation_settings_from_dataset

Selecting observations
----------------------

.. currentmodule:: tudatpy.estimation.observations

Use :data:`observation_query` to build :class:`ObservationSelectionCondition` objects. Conditions select
observation rows/events, not individual scalar components.

.. data:: observation_query

   Query builder for creating :class:`ObservationSelectionCondition` objects.

.. code-block:: python

   condition = (
       (observation_query.observable_type == observations.one_way_doppler)
       & (observation_query.receiver == observations.LinkEndId("Earth", "DSS-63"))
       & observation_query.time.between(t0, t1)
       & ~observation_query.rejected
   )

Combine conditions with ``&`` and ``|`` and negate them with ``~``. Do not use
Python ``and``, ``or`` or ``not``; those operators ask Python to convert the
condition to ``bool``. Parentheses are required around comparisons because of
Python operator precedence.

Common selectors
~~~~~~~~~~~~~~~~

.. code-block:: python

   observation_query.observable_type == observations.one_way_doppler
   observation_query.set_id == 3
   observation_query.link_definition == link_definition
   observation_query.receiver == observations.LinkEndId("Earth", "DSS-63")
   observation_query.link_end(observations.transmitter) == observations.LinkEndId("Earth", "DSS-25")
   observation_query.time.between(t0, t1)
   observation_query.time >= t0
   observation_query.time < t1
   observation_query.active
   observation_query.rejected
   ~observation_query.rejected

Value thresholds
~~~~~~~~~~~~~~~~

.. code-block:: python

   observation_query.residual.abs_greater_than(3.0 * sigma)
   observation_query.observation.abs_greater_than(1.0e9)
   observation_query.dependent_variable(elevation_settings).greater_than(15.0 * np.pi / 180.0)

For vector observables, scalar limits are broadcast to all components. A
threshold condition selects a row if any component exceeds the supplied limit.

.. autoclass:: tudatpy.estimation.observations.ObservationSelectionCondition
.. autoclass:: tudatpy.estimation.observations.ObservationSelectionConditionType

ObservationDataset
------------------

.. autoclass:: tudatpy.estimation.observations.ObservationDataset

Creating datasets
~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.add_observation_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.add_observation_set_with_weights
.. automethod:: tudatpy.estimation.observations.ObservationDataset.add_observation_set_from_dataset

Inspecting datasets
~~~~~~~~~~~~~~~~~~~

.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.number_of_observation_sets
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.number_of_observations
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.total_scalar_size
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_observation_set_metadata
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_link_end_reference_point
.. automethod:: tudatpy.estimation.observations.ObservationDataset.ancillary_settings_for_set

Creating views and reduced datasets
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.create_viewer
.. automethod:: tudatpy.estimation.observations.ObservationDataset.create_new_and_keep
.. automethod:: tudatpy.estimation.observations.ObservationDataset.create_new_and_drop

Rejecting, restoring, and removing observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.reject_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.restore_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.remove_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.remove_rejected_observations

Flattening data
~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.estimation_flattened_observation_data
.. automethod:: tudatpy.estimation.observations.ObservationDataset.computation_flattened_observation_data
.. automethod:: tudatpy.estimation.observations.ObservationDataset.ordered_flattened_observation_data

Diagnostics
~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.rms_residuals_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.mean_residuals_for_set

Observation weights
-------------------

Most weights are scalar per observation or ``N x N`` blocks per observation. A
scalar weight for a vector observable stays compact and expands to ``w I_N``
only when flattened matrix data are materialized. Diagonal vectors remain
diagonal storage. Cross-observation and cross-set correlations should use
sparse/block machinery rather than dense global matrices.

Common cases
~~~~~~~~~~~~

.. code-block:: python

   dataset.set_constant_single_observation_scalar_weight(
       condition=(observation_query.observable_type == observations.angular_position),
       weight=1.0e-9,
   )

   dataset.set_constant_single_observation_diagonal_weight(
       condition=(observation_query.set_id == set_id),
       weight=weight_vector,
   )

   dataset.set_weight_vector_for_set(set_id, weight_vector)
   dataset.set_weight_matrix_for_observation(observation_id, weight_matrix)
   dataset.set_weight_matrix_for_set(set_id, weight_matrix)

When a set-level matrix, per-observation weights, and extra scalar-component
blocks overlap, the effective matrix is the one returned by flattened data, for
example ``dataset.estimation_flattened_observation_data().weight_matrix``.
``weight_matrix_for_set`` returns the stored set-level block when one is present,
otherwise it materializes the set's compact per-observation weights.

Off-diagonal blocks
~~~~~~~~~~~~~~~~~~~

.. code-block:: python

   ids_a = dataset.observation_ids_matching_condition(observation_query.receiver == dss63)
   ids_b = dataset.observation_ids_matching_condition(observation_query.receiver == dss25)

   dataset.set_weight_block(
       row_observation_ids=ids_a,
       column_observation_ids=ids_b,
       weight_block=block,
   )

Blocks are symmetrized automatically; the transposed block is stored for the
opposite row/column selection.

Weight API
~~~~~~~~~~

.. autoclass:: tudatpy.estimation.observations.ObservationWeightSettings

.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_constant_single_observation_scalar_weight
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_constant_single_observation_diagonal_weight
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_constant_single_observation_matrix_weight
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_weight_vector_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_weight_matrix_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.has_weight_matrix_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_weight_matrix_for_observation
.. automethod:: tudatpy.estimation.observations.ObservationDataset.has_weight_matrix_for_observation
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_weight_block
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.extra_weight_blocks
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.has_extra_weight_blocks

.. autoclass:: tudatpy.estimation.observations.ObservationWeightBlock

Supporting dataset objects
--------------------------

ObservationDatasetViewer
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.estimation.observations.ObservationDatasetViewer

.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetViewer.number_of_observations
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetViewer.observation_ids
.. automethod:: tudatpy.estimation.observations.ObservationDatasetViewer.observation_row
.. automethod:: tudatpy.estimation.observations.ObservationDatasetViewer.observation_value
.. automethod:: tudatpy.estimation.observations.ObservationDatasetViewer.observation_time
.. automethod:: tudatpy.estimation.observations.ObservationDatasetViewer.create_viewer
.. automethod:: tudatpy.estimation.observations.ObservationDatasetViewer.estimation_flattened_observation_data
.. automethod:: tudatpy.estimation.observations.ObservationDatasetViewer.ordered_flattened_observation_data

FlattenedObservationData
~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.estimation.observations.FlattenedObservationData

.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.observation_vector
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.residual_vector
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.weight_vector
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.weight_matrix
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.sparse_weight_matrix
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.is_diagonal_weight_only
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.has_off_diagonal_weights
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.times
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.observation_ids
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.set_ids
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.scalar_component_ids
.. autoattribute:: tudatpy.estimation.observations.FlattenedObservationData.set_ids_in_row_order
.. automethod:: tudatpy.estimation.observations.FlattenedObservationData.unique_observation_ids_for_set
.. automethod:: tudatpy.estimation.observations.FlattenedObservationData.flattened_row

ObservationSetMetadata
~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.estimation.observations.ObservationSetMetadata

.. autoattribute:: tudatpy.estimation.observations.ObservationSetMetadata.observable_type
.. autoattribute:: tudatpy.estimation.observations.ObservationSetMetadata.link_definition_id
.. autoattribute:: tudatpy.estimation.observations.ObservationSetMetadata.reference_link_end
.. autoattribute:: tudatpy.estimation.observations.ObservationSetMetadata.observable_size
.. autoattribute:: tudatpy.estimation.observations.ObservationSetMetadata.ancillary_settings_id
.. autoattribute:: tudatpy.estimation.observations.ObservationSetMetadata.dependent_variable_layout_id

ObservationDatasetRow
~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.estimation.observations.ObservationDatasetRow

.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.time
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.set_id
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.first_scalar_component
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.scalar_size
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.index_in_set
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.is_active
.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.rejection_reason

ObservationScalarComponentRow
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. autoclass:: tudatpy.estimation.observations.ObservationScalarComponentRow

.. autoattribute:: tudatpy.estimation.observations.ObservationScalarComponentRow.observation_id
.. autoattribute:: tudatpy.estimation.observations.ObservationScalarComponentRow.component_index

Advanced dataset API
--------------------

These methods expose row ids, scalar-component rows, set-level vectors and
metadata registries. They are useful for diagnostics, ingestion tools and
advanced weight blocks, but they are not the primary observation-selection
workflow.

Rows and identifiers
~~~~~~~~~~~~~~~~~~~~

.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.observation_rows
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_row
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.scalar_component_rows
.. automethod:: tudatpy.estimation.observations.ObservationDataset.scalar_component_row
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_ids_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_ids_matching_condition

Per-set and per-observation data
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.observations_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_vector_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.computed_observations_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.computed_observation_vector_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_value
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_times_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.observation_time
.. automethod:: tudatpy.estimation.observations.ObservationDataset.weights_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.weight_vector_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.weight_value
.. automethod:: tudatpy.estimation.observations.ObservationDataset.weight_matrix_for_observation
.. automethod:: tudatpy.estimation.observations.ObservationDataset.weight_matrix_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.residuals_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.residual_vector_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.residual_value

Dependent variables
~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.dependent_variables_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.dependent_variables
.. automethod:: tudatpy.estimation.observations.ObservationDataset.single_dependent_variable_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.single_dependent_variable_for_set_by_index
.. automethod:: tudatpy.estimation.observations.ObservationDataset.compatible_dependent_variable_settings_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.all_compatible_dependent_variables_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_dependent_variables_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.clear_dependent_variables_for_set

Mutation helpers
~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_observations_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_residuals_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.add_observations_to_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.remove_observations_from_set

Registries and set metadata
~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.time_bounds_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.number_of_observations_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.total_scalar_size_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.link_definition
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.number_of_link_definitions
.. automethod:: tudatpy.estimation.observations.ObservationDataset.ancillary_settings
.. automethod:: tudatpy.estimation.observations.ObservationDataset.dependent_variable_bookkeeping

Advanced weight conveniences
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_constant_single_observation_scalar_weight_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_constant_single_observation_diagonal_weight_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.set_constant_single_observation_matrix_weight_for_set

Submodules
----------

.. toctree::
   :maxdepth: 1

   /estimation/observations/observations_geometry
