# PR 905 findings that may also affect PR 1000

This file records findings from the PR 905 validation that may apply to the
already completed PR 1000 work. These entries are records only: no PR 1000
change is to be made from them during this work.

## P1000-1: observation metadata becomes model ancillary settings

- **Found while running:**
  `examples/tudatpy/estimation/retrieving_mpc_observation_data.py` and
  `examples/tudatpy/estimation/estimation_with_mpc.py` with the PR 905 kernel.
- **Problem:** `create_observation_dataset_from_tracking_data` attempted to
  convert MPC per-observation metadata, such as the object `number`, into
  observation-model ancillary settings. This first raised an unknown-key error.
  After metadata was skipped, it still produced an empty ancillary-settings
  object, which angular-position observation models reject.
- **PR 905 correction:** observation metadata is skipped for every ancillary
  value category, and the conversion returns no ancillary settings when no
  model setting was transferred.
- **PR 1000 status:** possible shared issue. Record only; do not change or
  retest PR 1000 until the PR 905 work is complete.

## P1000-2: routine residual plotting needs scalar-aligned metadata

- **Found while running:**
  `examples/tudatpy/estimation/estimation_with_mpc.py` with the PR 905 kernel.
- **Problem:** `ObservationDataset.get_times()` and `get_set_ids()` return one
  value per observation event. Estimator residual histories are scalar vectors,
  so two-component angular-position observations require a scalar-aligned time
  and set-ID sequence for ordinary residual plots.
- **PR 905 correction:** the general scalar-aligned representation is now named
  `ObservationVectorData` and is returned by `observation_vector_data()`. The
  MPC estimation example uses its `times` and `set_ids` directly, without
  expanding values in example code.
- **PR 1000 status:** possible shared design issue. Record only; do not change
  or retest PR 1000 until the PR 905 work is complete.

## P1000-3: full daily TNF input makes the residual test unnecessarily slow

- **Found while running:**
  `tests/test_tudatpy/test_data_input_tracking_residuals.py` with the PR 905
  kernel.
- **Problem:** the MRO test loaded a 78,995,114-byte daily TNF file, then removed
  thousands of observations individually to retain a one-hour interval. It
  exceeded 125 s locally. Interrupted downloads were also accepted whenever a
  nonempty partial file remained.
- **PR 905 correction:** a 5,109,334-byte fixture retains all earlier frequency
  ramps needed at light-time-shifted transmit epochs and only the bounded
  Doppler record section. The same 59-sample and residual assertions pass in
  1.42 s. Test downloads are written to a temporary file and renamed only when
  complete.
- **PR 1000 status:** likely shared test issue. Record only; apply or retest it
  on PR 1000 after the PR 905 work is complete.

## P1000-4: migrated examples need direct dataset operations

- **Found while converting:** the GRAIL, MEX, and MRO examples to the PR 905
  `ObservationDataset` interface.
- **Problem:** the dataset initially lacked direct operations for precise global
  time bounds, observable-specific link definitions, selected dependent
  variables, reference-point switch histories, and scalar-aligned link IDs.
  Residual computation also required a differently named dataset function.
  Recreating these operations in each example would make the examples longer
  and easier to misuse.
- **PR 905 correction:** these operations are available on the dataset, and
  `compute_residuals_and_dependent_variables` accepts either representation.
  The compatibility overload for the earlier collection continues to emit its
  deprecation warning.
- **PR 1000 status:** possible shared interface issue. Record only; do not
  change or retest PR 1000 until the PR 905 work is complete.

## P1000-5: `time_bounds` is already occupied by a compatibility property

- **Found while validating:** the first name chosen for precise dataset-wide
  time bounds resolved to the older `time_bounds` compatibility property and
  emitted a deprecation warning.
- **PR 905 correction:** the modern property is named
  `observation_time_bounds`, leaving the compatibility property and its warning
  intact.
- **PR 1000 status:** naming collision to check after PR 905 work is complete.

## P1000-6: legacy collection changes could detach live dataset state

- **Found while running:** the PR 1000 and current-master MEX residual examples
  with the PR 905 kernel, followed by the GRAIL residual and estimation
  examples.
- **Problem:** legacy collection removal, splitting, replacement, and
  reference-point operations rebuilt wrapper lists without consistently
  retaining a live `ObservationDataset`. Later residual updates could therefore
  be absent from the legacy collection. Splitting also assumed dataset set ID
  zero and failed for a populated set following an empty set.
- **PR 905 correction:** no-op removal preserves the existing dataset;
  membership-changing operations rebuild one live dataset and rebind existing
  wrappers; reference-point overloads retain or adopt that dataset; empty sets
  are ignored during time-varying splitting; and splitting retains the source
  set ID. Four focused compatibility tests cover these cases.
- **PR 1000 status:** possible shared compatibility issue. The PR 1000 branch
  has not been modified during this PR 905 stage.
