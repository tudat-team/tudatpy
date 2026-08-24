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
- Builds use the `tudatpy-dev` Conda environment, `cmake-build-release`, and `-j6`.
- Existing builds are terminated before starting a new build.

## Implementation log

### Alias semantics

- Restored `Vector3ld`, `Vector6ld`, `Vector7ld`, `VectorXld`, and `MatrixXld` to literal `long double`.
- Added corresponding `Vector3hps`, `Vector6hps`, `Vector7hps`, `VectorXhps`, and `MatrixXhps` aliases using `tudat::HighPrecisionStateScalar`.
- Added compile-time assertions covering both alias families.

### `Time` reconstruction

- Changed `Time::getSeconds<ScalarType>()` to convert the period count, normalization constant, and stored residual separately, then reconstruct in `ScalarType`.
- Added large positive/negative epoch checks that distinguish direct quad reconstruction from premature `long double` recombination.
- Kept the split `Time` storage unchanged.

### Runge–Kutta coefficients (in progress)

- Replaced double-valued Butcher-tableau storage with `HighPrecisionStateScalar` matrices.
- Added a raw numeric literal parser so every written coefficient is constructed from its source digits; rational expressions divide values already constructed in the configured scalar.
- Added explicit narrowing only at consumption sites for lower-precision state/time types.
- Removed the RK4-only coefficient workaround.
- Added direct coefficient checks for RK4, RKF45, and Feagin 10(8).

### DSN frequency integral

- Replaced the unexplained `f_end * duration_error` correction with integration from the represented transmission-start epoch over the independently derived physical UTC duration.
- Implemented this directly for constant and piecewise-linear ramp models, including negative durations, invalid gaps, ramp crossings, and multiple segments.
- Added analytic constant, linear-ramp, independently perturbed-boundary, single-boundary, and multiple-segment tests.

### Lagrange interpolation

- Preserved exact-only interpolation-node detection. A magnitude-scaled tolerance would incorrectly identify distinct epochs as the same node.
- Added a large-epoch regression test covering both an exact node and a representable nearby epoch.

### CI configuration

- Made the default quad scalar explicit on both Unix and Windows workflow configuration paths.

## Decisions and questions under review

- The DSN averaged-Doppler residual-duration correction is not accepted merely because it improves a regression test. Its endpoint-error model and frequency-ramp behavior must be derived before retention or replacement.
- Root-finder/event precision changes will be limited to local, demonstrable bottlenecks. Broader API redesign is outside scope.
- External and physical-model inputs that originate as `double` will remain `double`; only avoidable narrowing in generic high-precision paths will be changed.

## Validation log

No post-change build or test has been run yet.
