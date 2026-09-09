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

   data = dataset.get_data(
       condition=(
           (observation_query.observable_type == observations.one_way_doppler)
           & observation_query.time.between(t0, t1)
           & ~observation_query.rejected
       ),
       fields=("times", "observations", "residuals"),
       ordering="estimation",
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

Ownership, identities and projections
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The dataset owns observations, residuals, dependent-variable values and weights.
Each event has a stable ``observation_id`` within its dataset. Removing another
row, sorting a set, or appending data never reuses or renumbers that identity.
``set_id`` identifies metadata grouping. Scalar storage positions may change;
use a fresh projection to obtain indices for numerical work.
In C++, ``projection.getFlattenedRow(observation_id, component_index)`` gives
the scalar offset in that projection. A set can span noncontiguous storage
after appends, so set counts alone do not define scalar offsets.

``create_estimation_projection()`` selects active rows in estimation order: observable type, link ends, set, event within the set, then component.
Estimation and covariance use this same route. Computation projections include
rejected rows by default, so residual diagnostics can inspect them. Restoring a
row retains its last rejection reason as well as its identity and weights.

A projection is a snapshot. Its observations, times, residuals, weights,
dependent variables and link metadata describe one selection in one order.
Residual writeback checks the originating dataset, structure, observed values,
selection and ancillary settings. Rebuild the projection after adding, removing,
regrouping, rejecting or restoring rows, changing observed values, or replacing
link/layout metadata. Updating residuals or weights does not invalidate an iteration
mapping, but does not refresh the projection's captured values. Estimation creates
its input projection once per call and computes fresh residuals and design matrices
on each iteration. Covariance also prepares its projection once per call. Set
observations, selection and weights before starting these numerical operations.

Independent dataset copies also clone mutable
ancillary settings and dependent-variable settings. Filtered copies preserve
metadata identifiers, including groups left empty by the selection. Custom C++
dependent-variable setting subclasses must implement ``clone()`` to participate
in metadata snapshots and independent dataset copies; unknown derived types fail explicitly instead
of being sliced or shared silently.

Inspecting observation data
~~~~~~~~~~~~~~~~~~~~~~~~~~

Use the dataset directly to obtain independent value snapshots:

.. code-block:: python

   times = dataset.get_times(condition, ordering="internal")
   residuals = dataset.get_residuals(condition, ordering="estimation")
   data = dataset.get_data(
       condition=condition,
       fields=("times", "observations", "residuals", "observation_ids"),
       ordering="estimation",
   )

Every new getter accepts the same ``condition`` (default: all observations) and
``ordering`` (default: ``"internal"``). Internal order preserves the selected
rows' relative dataset storage order. Estimation order uses the same
observable/link/set/event/component ordering as numerical estimation, including
existing set order and equal-time ties. Getters never reorder the dataset.
Changing ordering never changes membership: **rejected observations are included
by default in both orders**. To reproduce estimator membership as well, combine
your condition with ``observation_query.active``.

Snapshots remain valid after any dataset mutation or destruction. Updating
values, residuals, correlations, metadata, rejection status or event membership
cannot change a saved result. Conversely, editing a returned array, nested list,
dictionary or mutable setting cannot change the dataset. Snapshots can themselves
be modified; they need no revision or invalidation checks. Separate getter calls
capture separate dataset states. ``get_data`` resolves one selection and ordering
and gathers its requested fields once, within the existing supported mutation
model. Its default fields are ``("times", "observations")``.

.. list-table:: Fields returned by ``get_data`` and the corresponding ``get_<field>`` methods
   :header-rows: 1
   :widths: 27 18 55

   * - Field
     - Alignment
     - Owning return value
   * - ``times``
     - Event
     - List of ``Time`` objects; no conversion to Python float.
   * - ``observations``, ``residuals``
     - Event
     - List of arrays, one vector per event; mixed observable dimensions need no rectangular matrix.
   * - ``observation_ids``, ``set_ids``
     - Event
     - Lists of stable IDs, without renumbering; set IDs may repeat.
   * - ``dependent_variables``
     - Event
     - List of arrays; empty arrays represent missing values. Layouts are in ``metadata``.
   * - ``rows``
     - Event
     - List of dictionaries with IDs, precise time, component size, original storage/set positions, dependent values, active status and rejection reason.
   * - ``scalar_components``
     - Scalar component
     - List of ``(observation_id, component_index)`` tuples in event/component order.
   * - ``weight_diagonal``
     - Scalar component
     - NumPy vector containing effective diagonal weights.
   * - ``weight_matrix``
     - Scalar component on both axes
     - SciPy sparse principal matrix, including all selected correlations. Both axes follow ``scalar_components``.
   * - ``metadata``
     - Selected metadata groups
     - Dictionary keyed by set ID, containing set attributes, copied ``link_definition``, cloned ``ancillary_settings`` and ``dependent_variable_layout``.

``dependent_variable_layout`` maps ``(start, size)`` to detached settings for slices
of the corresponding event's dependent-variable vector. Missing ancillary settings
are ``None``; missing layouts are empty dictionaries. Only sets represented in the
selection appear. Metadata IDs still refer to the original dataset registries;
the returned dictionary contains the corresponding copied values for interpretation.
``rows.first_scalar_component`` describes the source storage at extraction time;
use ``scalar_components`` to interpret extracted scalar arrays and matrices.

Omitted residuals retain the existing dataset convention of stored zeros; there
is no separate flag distinguishing uncomputed residuals. Stored NaNs, if present,
are copied unchanged. Empty selections return empty lists/vectors/dictionaries
and a 0-by-0 sparse weight matrix. Invalid ordering, unknown fields and duplicate
fields raise ``ValueError``.

Extraction preserves the represented observation scalar and time precision. C++
getters retain ``ObservationScalarType`` and ``TimeType``; the standard Python
build uses double observation arrays and precise ``Time`` objects. Explicitly
converting a time to ``float`` can lose precision. A flattened inspection array
can be assembled with ``np.concatenate(data["observations"])``; its component
association is ``dataset.get_scalar_components`` for the same selection and order.
Prefer requesting both fields in one ``get_data`` call when alignment matters.

Copying costs memory proportional to requested values and selected metadata.
Single-field getters do not construct a ``FlattenedObservationData``, and asking
for times does not copy residuals, weights or dependent-variable values. Full
weight access remains sparse; call ``.toarray()`` only when a dense quadratic-size
matrix is explicitly needed. ``get_weight_diagonal`` does not assemble correlations.

C++ provides typed ``getTimes``, ``getObservations``, ``getResiduals``,
``getObservationIds``, ``getSetIds``, ``getRows``, ``getDependentVariableValues``,
``getScalarComponents``, ``getWeightDiagonal``, ``getWeightMatrix`` and
``getMetadata`` methods with ``ObservationOrdering::internal`` or
``ObservationOrdering::estimation``. Return values are owning STL/Eigen values.
``InspectionMetadata`` is an ordinary map of tuples containing set metadata,
link definition, cloned ancillary settings and the cloned dependent-variable
layout, in that order. Existing compatibility accessors retain their defaults
and ordering conventions.

The removed ``ObservationDatasetViewer`` and ``create_viewer`` API are replaced
by these direct getters. Narrow a condition with ``&`` or reuse selected stable
IDs through the existing condition mechanism. ``FlattenedObservationData`` and
its existing builders remain numerical inputs for estimation, covariance and
residual computation, with their source identity and safe-writeback checks.

.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_times
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_residuals
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_observation_ids
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_set_ids
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_rows
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_dependent_variables
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_scalar_components
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_weight_diagonal
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_weight_matrix
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_metadata
.. automethod:: tudatpy.estimation.observations.ObservationDataset.get_data

Legacy ownership and conversion
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

``SingleObservationSet`` is a facade over one dataset metadata group. Legacy
collections group these facades and preserve shared-set behavior: putting the
same set into another collection never migrates its data or changes its owner.
Converting such a collection to a dataset creates an independent snapshot.
Conversion preserves correlations between selected sets from the same source
dataset. Repeating correlated source rows in a conversion is ambiguous and
raises an error; repeated diagonal-only legacy sets remain supported.

A collection explicitly created from a dataset is a live facade over that
backend. Editing its set membership makes it a legacy grouping of the retained
sets; it does not delete data from the source dataset. Dataset-native rejection
and removal should be used to change the backend itself. Estimation and covariance
inputs constructed from a legacy collection prepare a fresh snapshot when an
operation starts, and estimation writes computed residuals back to those sets.
Binary files saved with the base-branch ``SingleObservationSet`` and
``ObservationCollection`` layouts remain readable through a compatibility
reader. New files store the shared dataset, preserving row status, identities
and correlations across groups.
Covariance and estimation result archives from the base branch also retain
their diagonal-weight interpretation when loaded.

An empty dataset or metadata group produces empty observation vectors. Simulating
an empty tabulated group preserves its metadata without requiring an observation
model. Legacy estimation and covariance inputs accept ``None`` for configuring
settings; operations that need observations then report a missing-source error.
Estimation
and covariance analysis report an explicit error when no active observations
remain, including after rejecting every row. Residuals
not supplied at creation default to zero, matching legacy behavior; zero alone
does not indicate that residual computation has run. Dependent-variable values
may be absent for an entire set, but a partially populated set is rejected.

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

.. automethod:: tudatpy.estimation.observations.ObservationDataset.create_new_and_keep
.. automethod:: tudatpy.estimation.observations.ObservationDataset.create_new_and_drop

Rejecting, restoring, and removing observations
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.reject_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.restore_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.remove_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.delete_rejected_observations
.. automethod:: tudatpy.estimation.observations.ObservationDataset.remove_rejected_observations

Flattening data
~~~~~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.create_estimation_projection
.. automethod:: tudatpy.estimation.observations.ObservationDataset.estimation_flattened_observation_data
.. automethod:: tudatpy.estimation.observations.ObservationDataset.computation_flattened_observation_data
.. automethod:: tudatpy.estimation.observations.ObservationDataset.ordered_flattened_observation_data

Diagnostics
~~~~~~~~~~~

.. automethod:: tudatpy.estimation.observations.ObservationDataset.rms_residuals_for_set
.. automethod:: tudatpy.estimation.observations.ObservationDataset.mean_residuals_for_set

Observation weights
-------------------

The dataset stores one effective symmetric weight matrix. Its diagonal is a
compact vector; only nonzero off-diagonal coefficients need additional storage.
Small per-observation blocks and cross-observation correlations update this same
matrix, without allocating a dense global matrix.

Assignments replace the addressed entries. A full-set diagonal or matrix replaces
that set's principal block. A per-observation diagonal or matrix replaces only
that observation's principal block and preserves correlations with other rows.
All weight getters report the current effective coefficients. Diagonal values
must be finite and nonnegative; blocks must be finite and consistent with a
symmetric matrix, including overlapping or permuted selections.

Selection, rejection and removal restrict weights to the corresponding principal
submatrix. Restoring rejected observations restores their original correlations.

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

Use observation identities to address cross-observation blocks:

.. code-block:: python

   ids_a = dataset.observation_ids_matching_condition(observation_query.receiver == dss63)
   ids_b = dataset.observation_ids_matching_condition(observation_query.receiver == dss25)

   dataset.set_weight_block(
       row_observation_ids=ids_a,
       column_observation_ids=ids_b,
       weight_block=block,
   )

The transpose is assigned automatically. If row and column selections overlap,
entries that address the same symmetric pair must agree; inconsistent requests
raise an exception before changing any weights.

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
.. autoattribute:: tudatpy.estimation.observations.ObservationDataset.has_extra_weight_blocks


Supporting dataset objects
--------------------------

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

.. autoattribute:: tudatpy.estimation.observations.ObservationDatasetRow.observation_id
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
