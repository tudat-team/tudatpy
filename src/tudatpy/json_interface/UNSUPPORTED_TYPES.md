# Unsupported input types

This list contains types for which no JSON representation or settings-based
alternative exists. Raw acceleration-, torque-, and mass-rate-model maps remain
unavailable as JSON values, but their corresponding `*_from_*_settings`
factories provide the supported JSON path.

## Required inputs

| Needed type | Contract factory inputs |
|---|---|
| `dict[float, tuple[float[2], float]]` | `environment_setup.rotation_model.iau_rotation_model.pole_periodic_terms` |
| `tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` | `propagation_setup.dependent_variable.acceleration_derivative_partial_wrt_parameter.parameter_settings`; `total_acceleration_derivative_partial_wrt_parameter.parameter_settings` |

## Optional inputs using Tudat defaults

| Needed type | Contract factory inputs |
|---|---|
| `tudatpy.math.interpolators.InterpolatorGenerationSettings` | `environment_setup.rotation_model.gcrs_to_itrs.cio_interpolation_settings`; `tdb_to_tt_interpolation_settings`; `short_term_eop_interpolation_settings` |
