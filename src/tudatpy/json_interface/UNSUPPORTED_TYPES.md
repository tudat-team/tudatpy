# Unsupported input types

This list contains types for which no JSON representation or settings-based
alternative exists. Raw acceleration-, torque-, and mass-rate-model maps remain
unavailable as JSON values, but their corresponding `*_from_*_settings`
factories provide the supported JSON path.

## Required inputs

| Needed type | Contract factory inputs |
|---|---|
| `dict[float, tuple[float[2], float]]` | `dynamics.environment_setup.rotation_model.iau_rotation_model.pole_periodic_terms` |
| `dict[float, numpy.ndarray[7]]` | `dynamics.environment_setup.rotation_model.tabulated_rotation.rotational_state_history` |
| `dict[int, list[tuple[int, int]]]` | `dynamics.parameters_setup.polynomial_gravity_field_variation_amplitudes.cosine_indices_per_power`; `sine_indices_per_power`; `periodic_gravity_field_variation_amplitudes.cosine_indices_per_period`; `sine_indices_per_period` |
| `dict[tuple[int, int], list[tuple[int, int]]]` | `dynamics.parameters_setup.mode_coupled_k_love_numbers.love_number_indices` |
| `dict[float, numpy.ndarray]` / `dict[float, float]` | `estimation.observable_models_setup.biases.arcwise_absolute_bias_per_time.bias_values_per_start_time`; `arcwise_relative_bias_per_time.bias_values_per_start_time`; `arc_wise_time_drift_bias_per_time.bias_value_per_start_time`; `arcwise_time_bias.time_bias_per_arc_start_time` |
| `dict[str, dict[float, numpy.ndarray]]` | `estimation.observable_models_setup.light_time_corrections.measured_tropospheric_light_time_correction.zenith_delays_per_station` |

## Optional inputs using Tudat defaults

| Needed type | Contract factory inputs |
|---|---|
| `tudatpy.math.interpolators.InterpolatorGenerationSettings` | `dynamics.environment_setup.rotation_model.gcrs_to_itrs.cio_interpolation_settings`; `tdb_to_tt_interpolation_settings`; `short_term_eop_interpolation_settings` |
| `dict[int, str]` | `estimation.observable_models_setup.light_time_corrections.dsn_tabulated_ionospheric_light_time_correction.spacecraft_name_per_id`; `quasar_name_per_id` |

## Inputs supplied from runtime context

These inputs are live objects that are supplied from the runtime context when
the settings are combined into models, not from JSON.

| Runtime input | Contract factory inputs |
|---|---|
| `tudatpy.dynamics.environment.SystemOfBodies` (required) | `dynamics.parameters_setup.initial_states.bodies` |
| `tudatpy.dynamics.environment.SystemOfBodies` (optional) | `dynamics.propagation_setup.propagator.translational_from_acceleration_settings.bodies`; `rotational_from_torque_settings.bodies`; `mass_from_mass_rate_settings.bodies` |

## Factories without a contract

Settings factories with a required Python callable (the `custom_*` factories,
`scaled_by_function`, `scaled_by_vector_function`, `time_variable_luminosity`,
`irradiance_based_time_variable_luminosity`, `time_varying_panel_geometry`,
`custom_inertial_direction_based`,
`step_size_control_custom_blockwise_*_tolerance`), factories taking runtime
data objects (`sp3_from_file_contents`, `observation_settings_from_collection`),
the Coma factories, the build-dependent Mars Climate Database factories, and
the deprecated integrator factories have no contract.
