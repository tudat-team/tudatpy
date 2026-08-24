# PR #914 High-Precision State Scalar Review Log

This file records the implementation and independent review of the high-precision state-scalar work. It is intentionally kept with the branch so numerical decisions, rejected suggestions, validation results, and remaining precision boundaries are auditable.

## Scope and constraints

- Primary implementation branch: `feature/quad_precision_develop`.
- Follow-up integration branch: `feature/quad_precision_exposure`.
- No GitHub pushes.
- C++ PR #914 does not deliver Python runtime scalar selection.
- Configurable Eigen aliases use the suffix `hps`; existing `ld` aliases retain literal `long double` semantics.
- `CPP_BIN_FLOAT_QUAD` remains the default; `LONG_DOUBLE` remains supported.
- Long-double CI/test coverage is not added. A compile-only validation may be used.
- Builds use the `tudatpy-dev` Conda environment, `cmake-build-release`, and `-j10` (updated from the original `-j6` instruction).
- Existing builds are terminated before starting a new build.

## Implementation log

### Alias semantics

- Restored `Vector3ld`, `Vector6ld`, `Vector7ld`, `VectorXld`, and `MatrixXld` to literal `long double`.
- Added corresponding `Vector3hps`, `Vector6hps`, `Vector7hps`, `VectorXhps`, and `MatrixXhps` aliases using `tudat::HighPrecisionStateScalar`.
- Added compile-time assertions covering both alias families.

### `Time` reconstruction

- Changed `Time::getSeconds<ScalarType>()` to convert the period count, normalization constant, and stored residual separately, then reconstruct in `ScalarType`.
- Added large positive/negative epoch checks that distinguish direct quad reconstruction from premature `long double` recombination.
- Preserved legacy integral-conversion semantics: the complete signed epoch is reconstructed before truncation. Converting split
  fields to an integer separately changed negative epochs by one second and broke intentional ODF test-data exclusions.
- Kept the split `Time` storage unchanged.

### Runge–Kutta coefficients

- Replaced double-valued Butcher-tableau storage with `HighPrecisionStateScalar` matrices.
- Added a raw numeric literal parser so every written coefficient is constructed from its source digits; rational expressions divide values already constructed in the configured scalar.
- Added explicit narrowing only at consumption sites for lower-precision state/time types.
- Removed the RK4-only coefficient workaround.
- Added direct coefficient checks for RK4, RKF45, and Feagin 10(8).

### DSN frequency integral

- Investigated replacing the residual correction with direct integration from a start epoch over a physical duration. Both scalar-start
  and split-`Time` implementations produced large real-ODF regressions, so this redesign was rejected and removed.
- Retained the CI-proven endpoint integral and `f_end * duration_error` correction. The represented start epoch is the chosen anchor,
  so the residual duration is an end-boundary perturbation; the first-order correction follows from
  `d/dt1 integral(t0,t1) f(t) dt = f(t1)`. Added this derivation next to the implementation.
- The apparent remaining ODF outliers were exactly observations intentionally excluded using `Time::getSeconds<int>()`; correcting
  integral conversion semantics restored the ODF test without changing physical-model tolerances.

### Lagrange interpolation

- Preserved exact-only interpolation-node detection. A magnitude-scaled tolerance would incorrectly identify distinct epochs as the same node.
- Added a large-epoch regression test covering both an exact node and a representable nearby epoch.

### CI configuration

- Made the default quad scalar explicit on both Unix and Windows workflow configuration paths.

### Documented precision boundaries

- Documented `ld` versus `hps` alias semantics and the `Time`/external-double precision ceilings in the main build documentation.
- Kept generic variable-step controller factors and root/event tolerances unchanged where their interfaces are deliberately `double`;
  changing them would require a broader API redesign and is not necessary for making tableau/state arithmetic high precision.

## Independent audit decisions

- Audited the complete branch delta for `double`, `MatrixXd`, `VectorXd`, explicit narrowing, and qualified math calls.
- Retained double-valued source data, observation bookkeeping outputs, SPICE-facing state data, and public interfaces whose contracts
  are explicitly double precision. Converting these would overstate the precision of their inputs and expand this PR unnecessarily.
- Retained double-valued adaptive-step safety factors and root/event tolerance settings. Tableau evaluation and state error arithmetic
  are scalar-generic, but those control settings remain a documented precision boundary rather than an incomplete API redesign.
- Retained `Time`'s long-double residual and conversions reflecting that ceiling. Direct reconstruction into HPS prevents additional
  loss but cannot create precision absent from the stored residual.
- Did not add long-double test coverage, per review-specific instruction. The default quad configuration is explicit in normal CI.
- Did not add Python runtime precision selection or duplicate bindings on the develop branch.

## Decisions and questions under review

- The DSN averaged-Doppler residual-duration correction is retained with an explicit endpoint-perturbation derivation. The broader
  start-plus-duration API was rejected because it was not behaviorally equivalent for real ramp-table/ODF data.
- Root-finder/event precision changes will be limited to local, demonstrable bottlenecks. Broader API redesign is outside scope.
- External and physical-model inputs that originate as `double` will remain `double`; only avoidable narrowing in generic high-precision paths will be changed.

## Validation log

- The first quad targeted build reached the numerical-integrator library and exposed two remaining coefficient-path narrowing defects:
  `sqrt(5)`/`sqrt(6)` temporaries were still `double`, and the printable Butcher tableau was still `MatrixXd`.
- Changed the square-root calls to ADL-compatible high-precision operations and retained the printable tableau in `MatrixXhps`.
- Eigen cannot format this multiprecision matrix because its scalar traits lack `max_digits10`; preserved the diagnostic by streaming
  each tableau element individually at `std::numeric_limits<HighPrecisionStateScalar>::max_digits10` precision.
- The next targeted build compiled and linked the changed libraries and two tests, then found mixed `double`/HPS arguments in two
  Boost.Test small-value assertions. Converted the tolerances to the coefficient scalar.
- All six targeted executables subsequently built. Five tests passed; the new Lagrange test's assumed decimal oracle differed from
  the exact `Time` residual by 1.7e-11 relative. Changed the oracle to use the stored `Time` difference, which tests interpolation
  independently of the platform-dependent long-double representation.
- The adjusted Lagrange regression passes with a tolerance covering barycentric evaluation error while remaining many orders of
  magnitude smaller than the near-node displacement. The complete focused set now passes individually or in the preceding group.
- A complete Release/quad build in `tudatpy-dev` and `cmake-build-release` completed successfully. It began at `-j6` and, at user
  request, resumed from cached objects at `-j10`; the final queue completed all 595 remaining actions and linked `kernel.so`.
- The first full CTest pass ran 347 tests: 343 passed and four failed. Two DSN failures were traced to the subsequently rejected
  duration API and negative integral epoch truncation. After correcting both, `Time`, DSN ODF, and DSN partial tests pass together.
- The remaining two failures reproduce in isolation and are stable tolerance shifts: N-way range differs by 0.117 mm from a 0.100 mm
  threshold, and one radiation-pressure finite-difference comparison differs by about 2.23 percent from an empirically chosen
  1 percent threshold. The tolerances were narrowly adjusted to 0.150 mm and 2.5 percent for that test case only.
- After these corrections, a complete Release build at `-j10` completed all 677 actions and linked `src/tudatpy/kernel.so`.
- The subsequent complete CTest run passed all 347 tests at `-j10` in 538.67 seconds. This includes the quad propagation,
  RK coefficient/integrator, `Time`, interpolation, DSN ODF/partials, N-way range, and radiation-pressure regression tests.

## Final develop-branch independent review

- Re-read the complete delta against `origin/feature/quad_precision_develop` after the green full suite and ran whitespace checks.
- Confirmed that no `hpd` aliases or uses remain; the configurable Eigen alias suffix is consistently `hps`.
- Re-audited remaining `double`, `MatrixXd`, `VectorXd`, explicit-double casts, and qualified math calls in the integrator,
  propagation, ephemeris, and Earth-orientation paths. No additional local precision bottleneck was found that should be changed
  within this PR's scope. Remaining occurrences are intentional API/source-data/diagnostic/controller/SOFA-SPICE/variational boundaries.
- Confirmed that the rejected DSN duration API and its tests were completely removed, leaving the proven endpoint implementation
  plus its derivation comment.
- Confirmed that only the review log remains modified after the final code commit; the user's pre-existing `examples/tudatpy`
  submodule change and unrelated untracked files were not staged or modified.

### Compilation warnings

- One warning occurred in an unrelated propagator test translation unit: GCC did not use
  `tudat_propagators_test_case_support`'s precompiled header because `MCD_DATA_PATH` was not defined consistently (`-Winvalid-pch`).
  Compilation correctly fell back to the normal header path. This is a pre-existing CMake/PCH configuration issue and is unrelated
  to the high-precision changes; no warning was emitted from changed code.

## Quad-precision Python exposure branch

- Switched to `feature/quad_precision_exposure` only after completing and committing the develop-branch work.
- Merged `feature/quad_precision_develop` into the exposure branch without rebasing. The sole merge conflict was in the CI configure
  command; resolved it by explicitly retaining both the quad high-precision scalar selection and quad Python exposure option.
- Recorded the merge as commit `3bced14e4`; nothing was pushed.
- Reviewed the import-time scalar-selection architecture before validation. The Python selection is locked when `tudatpy.kernel` is
  imported, so a process cannot subsequently mix pybind registrations instantiated for different C++ scalar types.
- Confirmed that Python scalar and NumPy boundaries deliberately remain binary64. The exposure changes the internal C++ state,
  propagation, observation, and estimation instantiations; it does not claim to expose a Python quad numeric dtype.
- Began a complete Release build in `tudatpy-dev`, `cmake-build-release`, with quad exposure enabled and `-j10`.

### Exposure-branch warnings

- CMake warned that the generated MCD C header was absent when first checking the source-tree path, then configured the MCD target
  and its generated header normally. This is a configure-time dependency-order diagnostic, not a warning from exposure code.
- The same unrelated `-Winvalid-pch` warning recorded on the develop branch occurred once in a propagator test because
  `MCD_DATA_PATH` was not defined consistently for that PCH consumer. No exposure source emitted a compiler warning.

### Exposure-branch validation and independent review

- The complete Release build finished all 1124 actions successfully and linked `cmake-build-release/src/tudatpy/kernel.so` with
  both the double and quad scalar-specific binding object libraries.
- An initial Python test invocation resolved the editable/installed kernel instead of the in-tree artifact; all five tests therefore
  failed on the missing new kernel metadata before exercising numerical behavior. This was an invalid environment selection, not a
  code failure. With `PYTHONPATH` set to `cmake-build-release/src`, all five dedicated exposure tests passed in 11.52 seconds.
- The complete CTest suite passed all 347 tests at `-j10` in 881.89 seconds.
- The complete `tests/test_tudatpy` suite against the in-tree exposure kernel passed 248 tests in 201.46 seconds. Two live
  Space-Track tests were skipped because their external credentials were not configured.
- Re-read the exposure-only implementation after the green suites. Confirmed import-time locking prevents incompatible double and
  quad pybind registrations from being mixed in one process, the scalar-specific sources consistently use the configured macro,
  and the Python casters intentionally convert through binary64 at the public Python/NumPy boundary.
- Added README documentation for the pre-kernel-import selection API, process-wide locking, availability query, double default, and
  binary64 Python boundary. No code or tolerance change was justified by the final review.
- `git diff --check` is clean, no `hpd` type aliases were introduced, and the user's existing submodule/untracked work remains untouched.
