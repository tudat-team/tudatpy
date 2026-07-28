# Remaining stub-generation issues

This report records the remaining issues found during a source-level audit of
the changes intended to address the original `pybind11-stubgen` errors.
No compilation was performed, and the existing `kernel.so` predates the
current source changes.

## Addressed categories

- All 23 enum classes reported by `pybind11-stubgen` have corresponding
  `--enum-class-locations` mappings.
- All 102 unstable object or container defaults from the original output now
  use stable `py::arg_v` representations.
- Every raw C++ type or signature reported in the original output has a
  corresponding exposure-order change, registration, or wrapper in the
  current source.

These findings are based on source inspection. A build and a new stubgen run
are still required to verify the result at runtime.

## Remaining issues

### Defaults that do not yet use their documented factories

- The 47 `ObservationCollectionParser` defaults still use
  `std::make_shared<ObservationCollectionParser>()`, while their documentation
  refers to `observations_processing.observation_parser()`.
- The 24 `LightTimeConvergenceCriteria` defaults still use
  `std::make_shared<LightTimeConvergenceCriteria>()`, while their documentation
  refers to `light_time_convergence_settings()`.
- The `EstimationConvergenceChecker` default still uses
  `std::make_shared<EstimationConvergenceChecker>()`, while its documentation
  refers to `estimation_convergence_checker()`.

The factory functions are exact equivalents for their respective
default-constructed objects and should be used consistently in both the C++
defaults and documentation.

### Lagrange interpolator defaults

Three defaults construct `LagrangeInterpolatorSettings` directly but document
the default as `lagrange_interpolation(order)`. These are not equivalent:

- The constructor defaults to `extrapolate_at_boundary`.
- `lagrange_interpolation(order)` defaults to
  `extrapolate_at_boundary_with_warning`.

The public factory can reproduce the existing behavior only if
`boundary_interpolation=extrapolate_at_boundary` is specified explicitly.

### Missing `set_vmf_troposphere_data` parameter documentation

`set_vmf_troposphere_data` now uses `py::arg_v` for its interpolator default,
but its docstring does not contain a `Parameters` section. Consequently, the
actual default factory is not documented there.

### Quaternion callback behavior

`CustomThrustOrientationSettings.thrust_orientation_function` originally
exposed the stored callback directly. That callback returns an
`Eigen::Quaternion<double>`.

The current wrapper instead returns a new callback that converts every
quaternion result to an `Eigen::Matrix3d`. This avoids the raw
`Eigen::Quaternion<double>` type in the generated signature, but it changes the
Python API:

- Previously, the property represented the original quaternion-returning
  callback.
- With the wrapper, it represents a different matrix-returning callback.
- User code expecting quaternion semantics, quaternion components, or the
  original callable identity can therefore behave differently.

The stubgen warning is addressed, but by changing behavior rather than only
changing type registration or representation. This should be reconsidered so
that fixing the stubs does not silently alter the public API.

### Defaults without an existing public factory

Two descriptions were deliberately left unresolved:

- The default DSN station-name map is produced by
  `getDefaultDsnStationNamesPerComplex()`, which is not currently exposed as a
  public Python factory.
- The empty ancillary-settings default uses
  `ObservationAncillarySimulationSettings()`. A public class constructor
  exists, but no dedicated factory function exists.

A decision is required on whether to document the existing literal/constructor
or add public factory functions.

## Checks performed without compiling

- `git diff --check` passes.
- `python -m py_compile build.py` passes.
- The inspected files show no IDE code-model errors.
- The current built extension cannot validate these source changes because it
  predates them.
