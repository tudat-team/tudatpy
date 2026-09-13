# Follow-up Report for AI Review

## Context

This report describes the changes made after the static Claude review of the
`feature/data-refactor` branch. The review was based on the branch state at:

- Reviewed head: `4e2672635098`
- Current base commit when this report was generated: `1bb8c430d`
- Local working-tree edits after that commit are included in the diff.
- Incremental diff file: `claude_review_followup_vs_reviewed_head.diff`

The incremental diff was generated with:

```bash
git diff --no-ext-diff 4e2672635098 -- . ':!examples/tudatpy' > claude_review_followup_vs_reviewed_head.diff
```

The `examples/tudatpy` submodule was excluded because it has unrelated local
state and the MRO example was explicitly not to be touched.

## Summary of Changes Made

### PERF-1 / BLK-1: repeated time-converter construction

The review found direct construction of
`earth_orientation::TerrestrialTimeScaleConverter` in new observation-data paths,
bypassing the cached default converter. This has been addressed.

Changes:

- Replaced direct converter construction in the reviewed data-flow paths with
  `earth_orientation::createDefaultTimeConverter()`.
- The default behavior of `createDefaultTimeConverter` now reuses the cached
  converter when no custom EOP reader is supplied.
- The TrackingData-to-ObservationCollection path only obtains the converter when
  the input time scale is not already TDB.
- Added a C++ test comparing default-constructed and factory-created converters
  for UTC/TDB conversion consistency.

### CORR-1: observation corrections dropped by public collection path

The review found that correction data stored in `TrackingData` could not be
applied through the public collection creation entry point.

Changes:

- Added an `applyCorrections` argument to
  `createObservationCollection`.
- Forwarded the flag into `createSingleObservationSetFromTrackingData`.
- Added Python coverage checking that corrections are ignored by default and
  applied when requested.

### CORR-2: concatenated epochs for vector observables

The review found that `TrackingData::getObservationEpochsVector()` incorrectly
assigned scalar epochs into vector segments.

Changes:

- Fixed the implementation to repeat each observation epoch for every component
  of a vector-valued observable.
- Added a dedicated C++ `TrackingData` test target.
- The test file now covers:
  - concatenated vector-valued observations,
  - concatenated epochs for vector-valued observations,
  - concatenated vector-valued weights,
  - validation of inconsistent weight/correction sizes.

### COMPAT-5 / BLK-2: removed `observations_wrapper` path

The review found that the old
`tudatpy.estimation.observations_setup.observations_wrapper` import path had
been removed without a compatibility shim.

Changes:

- Reintroduced the old module path as a compatibility wrapper.
- It lazily forwards to `tudatpy.estimation.observations`.
- It emits `DeprecationWarning` only when the old interface is used.
- Added Python deprecation/compatibility coverage.

### `tudatpy.data` backwards compatibility

After review, the compatibility requirement was clarified: everything that used
to work through `tudatpy.data` and its submodules must keep working during the
deprecation period.

Changes:

- Re-enabled the old `tudatpy.data` compatibility surface.
- Restored old raw-reader behavior for:
  - `read_crd_file`,
  - `read_tracking_txt_file`,
  - `read_ifms_file`,
  - `read_fdets_file`.
- Restored old weather-data functions:
  - `set_dsn_weather_data_in_ground_stations`,
  - `set_estrack_weather_data_in_ground_stations`.
- These functions emit deprecation warnings and delegate to active
  compatibility bindings in `add_data_to_kernel.cpp` that preserve the old
  signatures/return behavior.
- The old `src/tudatpy/data/expose_data.*` and
  `src/tudatpy/data/coma_model/expose_coma_model.*` C++ files were removed
  after verifying they are not part of the kernel build and have no remaining
  references. This does not remove the public compatibility API: old
  `tudatpy.data` raw-reader/weather functions are served by
  `add_data_to_kernel.cpp` and old `tudatpy.data.coma_model` names are served by
  the Python shim forwarding to `tudatpy.data_input.environment_data.coma`.

### `get_weights_VFCC17`

The review initially treated this function as intentionally removed. The
compatibility requirement was later clarified: it must remain available, even
though there is no one-to-one replacement in the new workflow.

Changes:

- Restored `tudatpy.data.mpc.get_weights_VFCC17`.
- Restored access through `tudatpy.data.mpc.mpc.get_weights_VFCC17`.
- Added a deprecation warning explaining:
  - there is no one-to-one equivalent in the new workflow,
  - the old function remains supported during the deprecation period,
  - the new workflow should request VFCC17 weighting with `add_weights=True`,
  - direct weight arrays/intermediate tables are no longer the new public
    interface.
- Added a link to the TudatPy migration guide in the warning.
- Added Python tests for both old access paths.
- Deferred the `astroquery`/`astropy` imports into the function body so importing
  `tudatpy.data.mpc` does not eagerly require those heavier dependencies.

### DEPR-1: warning stack levels

The review found that deprecation warnings pointed at TudatPy internals instead
of user call sites.

Changes:

- Updated stack levels in `data/_compat.py`.
- Added tests that check warning source locations for representative deprecated
  access.

### DOC-1: committed generated `docs/tudatpy/build.zip`

The review flagged a committed generated binary.

Changes:

- Deleted `docs/tudatpy/build.zip`.
- No `.gitignore` change was retained.

Note: deleting the file from this branch does not remove the blob from existing
git history. No history rewrite was attempted, and no `.gitignore` change was
retained because the maintainer requested that `.gitignore` be left unchanged.

### API-2: TNF frequency bands used raw enum-number duplicates

The review flagged Python-side raw frequency-band numbers such as `1.0` for
X-band. The maintainer later clarified that the data-input layer should use
strings, not higher-level Tudat enums.

Changes:

- Removed the TNF `frequencyBandIds` numeric map.
- TNF Doppler and range converters now store:
  - `"frequency bands"` as string vectors, e.g. `["X-band", "X-band"]`,
  - `"DSN reference frequency band at reception"` as a string, e.g. `"X-band"`.
- The central C++ observation-collection conversion remains responsible for
  translating these strings into observation-model frequency-band values.
- Updated affected Python tests to use the string path.

### API-3: `TrackingData` const-correctness and copies

The review flagged read-only `TrackingData` getters as non-`const` and
copy-returning.

Changes:

- Marked read-only getters `const`.
- Returned stored strings, vectors, maps, weights, and corrections by `const&`
  where appropriate.
- Left flattened vector getters returning by value because they construct new
  concatenated vectors by design.

## Points Not Addressed, With Maintainer Approval

### API-1: string ancillary setting keys

The review flagged string keys such as `"frequency bands"` and
`"DSN Doppler reference frequency"` as fragile.

Maintainer decision:

- This is intentional.
- The `data_input` layer should not depend on higher-level observation-model
  enums.
- Reviewers should not request enum use in `data_input` for these keys.
- Reviewers may still flag actual typos, missing validation, or inconsistent
  key spelling.

### NEEDS-RUN-1: full old-vs-new numerical regression campaign

The review asked for broad numerical equivalence testing against `develop`
across representative ODF/TNF/IFMS/FDETS/PSF/MPC workflows.

Maintainer decision:

- This remains not fully performed.
- Targeted tests and residual checks exist, but there is no exhaustive
  branch-to-branch golden-output comparison in this follow-up.
- This residual risk is approved for now.

### COMPAT-3: dynamic extraction of old `tudatpy.data` public symbols

The review suggested mechanically extracting the full old public symbol set from
`develop` and comparing it against the compatibility layer.

Maintainer decision:

- This was not done dynamically.
- Instead, compatibility is covered with an explicit old-symbol inventory and
  direct tests of representative deprecated functions/classes.
- This limitation is approved for now.

### API-3 remaining copy behavior for flattened vectors

The review suggested reducing copies in `TrackingData` getters.

Maintainer decision:

- Stored containers now return by `const&`.
- Flattened vector getters still return by value because they assemble new
  vectors.
- This remaining by-value behavior is intentional and approved.

## Verification Already Run

The following focused checks passed after the follow-up changes:

```bash
conda run -n tudatpy-dev cmake --build cmake-build-release --target kernel -j6
conda run -n tudatpy-dev cmake --build cmake-build-release --target test_io_TrackingData -j6
conda run -n tudatpy-dev ./cmake-build-release/tests/test_io_TrackingData
PYTHONPATH=cmake-build-release/src:src conda run -n tudatpy-dev python -m pytest tests/test_tudatpy/test_data_deprecation_compatibility.py -q
PYTHONPATH=cmake-build-release/src:src conda run -n tudatpy-dev python -m pytest tests/test_tudatpy/test_tnf.py -q
PYTHONPATH=cmake-build-release/src:src conda run -n tudatpy-dev python -m pytest tests/test_tudatpy/test_data_input_tracking_residuals.py::test_fdets_juice_short_arc_residual_scatter_is_millihertz_level -q
git diff --check
```

Observed results:

- `test_io_TrackingData`: passed, 4 C++ cases.
- `test_data_deprecation_compatibility.py`: passed, 161 Python cases.
- `test_tnf.py`: passed, 13 Python cases.
- FDETS residual test: passed.
- `git diff --check`: clean.

## Suggested Focus for the Next AI Review

Please review the incremental diff with these priorities:

1. Confirm that the old `tudatpy.data` compatibility functions preserve old
   signatures and behavior while emitting appropriate `DeprecationWarning`s.
2. Check that `get_weights_VFCC17` faithfully preserves the old public behavior
   and warning message.
3. Verify that TNF band strings are consistently converted only at the central
   TrackingData-to-observation-collection boundary.
4. Check that returning `const&` from `TrackingData` getters is safe for existing
   C++ and pybind call sites.
5. Do not re-raise the approved design choices listed above unless there is a
   concrete implementation bug.
