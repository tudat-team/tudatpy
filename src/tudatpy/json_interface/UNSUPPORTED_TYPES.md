# Unsupported input types

The contract `unsupported` fields are authoritative for this list. Required
unsupported inputs prevent the associated factory from being loaded. Optional
unsupported inputs are omitted so that the Tudat factory uses its own default.

## Required inputs

| Needed type | Contract factory inputs |
|---|---|
| `tudatpy.data.coma_model.ComaPolyDataset` | `environment_setup.atmosphere.coma_model_from_poly_data.poly_data` |
| `tudatpy.data.coma_model.ComaStokesDataset` | `environment_setup.atmosphere.coma_model_from_stokes_data.stokes_data` |
| `tudatpy.data.coma_model.ComaWindDatasetCollection` | `environment_setup.atmosphere.coma_wind_model.dataset_collection` |
| `tudatpy.math.interpolators.InterpolatorSettings` | `environment_setup.gravity_field_variation.tabulated.interpolation_settings` |
| `dict[float, tuple[float[2], float]]` | `environment_setup.rotation_model.iau_rotation_model.pole_periodic_terms` |
| `tudatpy.dynamics.parameters_setup.EstimatableParameterSettings` | `propagation_setup.dependent_variable.acceleration_derivative_partial_wrt_parameter.parameter_settings`; `total_acceleration_derivative_partial_wrt_parameter.parameter_settings` |
| Acceleration model map | `propagation_setup.propagator.translational.acceleration_models` |
| Torque model map | `propagation_setup.propagator.rotational.torque_models` |
| Mass-rate model map | `propagation_setup.propagator.mass.mass_rate_models` |

The acceleration model map is intentionally replaced for JSON use by
`propagation_setup.propagator.translational_from_acceleration_settings`, which
accepts acceleration settings and creates the models using the supplied bodies.

## Optional inputs using Tudat defaults

| Needed type | Contract factory inputs |
|---|---|
| `tudatpy.math.interpolators.InterpolatorSettings` | `environment_setup.aerodynamic_coefficients.tabulated.interpolator_settings`; `tabulated_force_only.interpolator_settings`; `tabulated_from_files.interpolator_settings`; `tabulated_force_only_from_files.interpolator_settings`; `environment_setup.ephemeris.interpolated_spice.interpolator_settings`; `tabulated_from_existing.interpolator_settings` |
| `tudatpy.math.interpolators.InterpolatorGenerationSettings` | `environment_setup.rotation_model.gcrs_to_itrs.cio_interpolation_settings`; `tdb_to_tt_interpolation_settings`; `short_term_eop_interpolation_settings` |
| `tudatpy.astro.element_conversion.PositionElementTypes` | `environment_setup.ground_station.basic_station.station_position_element_type` |
| `tudatpy.dynamics.propagation_setup.propagator.HybridArcPropagatorProcessingSettings` | `propagation_setup.propagator.hybrid_arc.processing_settings` |
| `tudatpy.dynamics.propagation_setup.propagator.MultiArcPropagatorProcessingSettings` | `propagation_setup.propagator.multi_arc.processing_settings` |
| `tudatpy.dynamics.propagation_setup.propagator.SingleArcPropagatorProcessingSettings` | `propagation_setup.propagator.translational.processing_settings`; `translational_from_acceleration_settings.processing_settings`; `rotational.processing_settings`; `mass.processing_settings`; `multitype.processing_settings`; `direct_relativistic_time_settings.output_settings` |
| `tudatpy.math.root_finders.RootFinderSettings` | `propagation_setup.propagator.dependent_variable_termination.termination_root_finder_settings` |
