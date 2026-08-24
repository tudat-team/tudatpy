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

### Documented precision boundaries

- Documented `ld` versus `hps` alias semantics and the `Time`/external-double precision ceilings in the main build documentation.
- Kept generic variable-step controller factors and root/event tolerances unchanged where their interfaces are deliberately `double`;
  changing them would require a broader API redesign and is not necessary for making tableau/state arithmetic high precision.

## Decisions and questions under review

- The DSN averaged-Doppler residual-duration correction is not accepted merely because it improves a regression test. Its endpoint-error model and frequency-ramp behavior must be derived before retention or replacement.
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
