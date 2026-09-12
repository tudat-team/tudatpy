# Deprecated data API: compatibility and removal

This package preserves the pre-refactor `tudatpy.data` API during deprecation.
Keep new functionality in `tudatpy.data_input`, `tudatpy.dynamics`, or
`tudatpy.estimation`. Dependencies must point from this package to the current
API; current Python modules must not import this package or `tudatpy.kernel.data`.
`test_data_compatibility_boundaries.py` checks that source dependency boundary.

## Warning contract

- Public aliases warn with `DeprecationWarning` when retrieved through attribute
  access or `from ... import ...`. They return the actual target class/function,
  preserving identity; reusing an already imported alias need not warn again.
- Historical nested modules have their own lookup functions so the message names
  the public import path actually used.
- Functions with retained signatures, such as raw-file readers and
  `get_weights_VFCC17`, warn when called. Private `_legacy` modules are implementation
  details, not additional public entry points.
- `_compat.py` owns warning emission and maps retained implementations to public
  migration targets. Warnings point past compatibility/import machinery to the
  caller. Preserve Python's normal warning filters; do not enable warnings globally.
- Calls through current API paths must not emit compatibility warnings. Do not
  import a deprecated alias internally just to reuse its implementation.

The MPC, TNF, and OMMUtils behavior retained by the review lives in `_legacy.py`
and `_legacy_converters.py` files beneath this package. The standalone legacy
weighting implementation is `mpc/_vfcc17.py`. These may depend on current helpers;
current helpers must not depend on them.

## Removal when the deprecation period ends

1. Delete `src/tudatpy/data/`, including its alias tables, warning helper, legacy
   implementations, and this note. Retire the compatibility-only tests
   `test_data_deprecation_compatibility.py`, `test_data_legacy_calls.py`, and
   `test_data_compatibility_boundaries.py` in `tests/test_tudatpy/` together.
2. Remove `src/tudatpy/add_data_to_kernel.cpp` and its source entry in
   `src/tudatpy/CMakeLists.txt`. In `src/tudatpy/kernel.cpp`, remove the
   `add_data_to_kernel` declaration/call and the deprecated `data` submodule
   creation. Preserve `data_input` registration and shared C++ readers, which
   remain part of the current API.
3. Remove the old API documentation at `docs/tudatpy/source/data.rst` and
   `docs/tudatpy/source/data/`, and update any documentation references to them.
4. Audit remaining references to `tudatpy.data` and `tudatpy.kernel.data`, migrate
   examples/callers, and run the current data-input and environment tests. Keep
   the optical, frequency, metadata, and other correctness regressions introduced
   by the review: those test retained functionality, not the deprecated API.

The deprecated `HorizonsQuery`/`HorizonsBatch` wrappers in
`dynamics/environment_setup/ephemeris/horizons_wrapper.py` also serve historical
imports outside `tudatpy.data`. Retire those classes and their exports separately
when that deprecation ends. Keep `jpl_horizons`, `jpl_horizons_from_query`, and
`add_horizons_batch_ephemerides`: they implement the supported environment workflow.
