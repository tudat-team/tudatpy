# Remaining stub-generation issues

This report records the resolution of issues found during a source-level audit
of the changes intended to address the original `pybind11-stubgen` errors.

## Addressed categories

- All 23 enum classes reported by `pybind11-stubgen` have corresponding
  `--enum-class-locations` mappings.
- All unstable object or container defaults from the original output now
  use stable `py::arg_v` representations.
- Every raw C++ type or signature reported in the original output has a
  corresponding exposure-order change, registration, or wrapper in the
  current source.
- The `ObservationCollectionParser`, `LightTimeConvergenceCriteria`, and
  `EstimationConvergenceChecker` defaults now use their documented factories.
- The three Lagrange defaults document the explicit
  `boundary_interpolation=extrapolate_at_boundary` argument needed to preserve
  their existing behavior.
- `set_vmf_troposphere_data` documents all inputs and its interpolator default.
- Public factories are exposed for the default DSN station-name map and empty
  ancillary settings, and the corresponding defaults refer to them.

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

The wrapper is intentionally left unchanged.
