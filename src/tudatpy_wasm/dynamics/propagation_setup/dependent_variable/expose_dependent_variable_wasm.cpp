/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../../wasm_module.h"
#include "../../../stl_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/propagation_setup/propagationOutputSettings.h>

namespace tp = tudat::propagators;

WASM_MODULE_PATH("dynamics_propagation_setup_dependent_variable")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_dependent_variable) {
    using namespace emscripten;

    // PropagationDependentVariables enum (subset of most common ones)
    enum_<tp::PropagationDependentVariables>("dynamics_propagation_setup_dependent_variable_PropagationDependentVariables")
        .value("mach_number_dependent_variable", tp::mach_number_dependent_variable)
        .value("altitude_dependent_variable", tp::altitude_dependent_variable)
        .value("airspeed_dependent_variable", tp::airspeed_dependent_variable)
        .value("local_density_dependent_variable", tp::local_density_dependent_variable)
        .value("relative_speed_dependent_variable", tp::relative_speed_dependent_variable)
        .value("relative_position_dependent_variable", tp::relative_position_dependent_variable)
        .value("relative_distance_dependent_variable", tp::relative_distance_dependent_variable)
        .value("relative_velocity_dependent_variable", tp::relative_velocity_dependent_variable)
        .value("radiation_pressure_dependent_variable", tp::radiation_pressure_dependent_variable)
        .value("total_acceleration_norm_dependent_variable", tp::total_acceleration_norm_dependent_variable)
        .value("single_acceleration_norm_dependent_variable", tp::single_acceleration_norm_dependent_variable)
        .value("total_acceleration_dependent_variable", tp::total_acceleration_dependent_variable)
        .value("single_acceleration_dependent_variable", tp::single_acceleration_dependent_variable)
        .value("aerodynamic_force_coefficients_dependent_variable", tp::aerodynamic_force_coefficients_dependent_variable)
        .value("aerodynamic_moment_coefficients_dependent_variable", tp::aerodynamic_moment_coefficients_dependent_variable)
        .value("stagnation_point_heat_flux_dependent_variable", tp::stagnation_point_heat_flux_dependent_variable)
        .value("local_temperature_dependent_variable", tp::local_temperature_dependent_variable)
        .value("geodetic_latitude_dependent_variable", tp::geodetic_latitude_dependent_variable)
        .value("control_surface_deflection_dependent_variable", tp::control_surface_deflection_dependent_variable)
        .value("total_mass_rate_dependent_variables", tp::total_mass_rate_dependent_variables)
        .value("tnw_to_inertial_frame_rotation_dependent_variable", tp::tnw_to_inertial_frame_rotation_dependent_variable)
        .value("periapsis_altitude_dependent_variable", tp::periapsis_altitude_dependent_variable)
        .value("apoapsis_altitude_dependent_variable", tp::apoapsis_altitude_dependent_variable)
        .value("keplerian_state_dependent_variable", tp::keplerian_state_dependent_variable)
        .value("modified_equinocial_state_dependent_variable", tp::modified_equinocial_state_dependent_variable)
        .value("spherical_harmonic_acceleration_terms_dependent_variable", tp::spherical_harmonic_acceleration_terms_dependent_variable)
        .value("spherical_harmonic_acceleration_norm_terms_dependent_variable", tp::spherical_harmonic_acceleration_norm_terms_dependent_variable)
        .value("body_fixed_airspeed_based_velocity_variable", tp::body_fixed_airspeed_based_velocity_variable)
        .value("total_aerodynamic_g_load_variable", tp::total_aerodynamic_g_load_variable)
        .value("body_fixed_groundspeed_based_velocity_variable", tp::body_fixed_groundspeed_based_velocity_variable)
        .value("local_dynamic_pressure_dependent_variable", tp::local_dynamic_pressure_dependent_variable)
        .value("current_body_mass_dependent_variable", tp::current_body_mass_dependent_variable)
        .value("radiation_pressure_coefficient_dependent_variable", tp::radiation_pressure_coefficient_dependent_variable)
        .value("rsw_to_inertial_frame_rotation_dependent_variable", tp::rsw_to_inertial_frame_rotation_dependent_variable)
        .value("custom_dependent_variable", tp::custom_dependent_variable);

    // VariableSettings base class
    class_<tp::VariableSettings>("dynamics_propagation_setup_dependent_variable_VariableSettings")
        .smart_ptr<std::shared_ptr<tp::VariableSettings>>("shared_ptr_VariableSettings");

    // SingleDependentVariableSaveSettings
    class_<tp::SingleDependentVariableSaveSettings, base<tp::VariableSettings>>(
        "dynamics_propagation_setup_dependent_variable_SingleDependentVariableSaveSettings")
        .smart_ptr<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>(
            "shared_ptr_SingleDependentVariableSaveSettings");

    // Factory functions for common dependent variables
    function("dynamics_propagation_setup_dependent_variable_mach_number",
        &tp::machNumberDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_altitude",
        &tp::altitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_density",
        &tp::densityDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_position",
        &tp::relativePositionDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_velocity",
        &tp::relativeVelocityDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_distance",
        &tp::relativeDistanceDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_speed",
        &tp::relativeSpeedDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_keplerian_state",
        &tp::keplerianStateDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_latitude",
        &tp::latitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_longitude",
        &tp::longitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_heading",
        &tp::headingDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_flight_path_angle",
        &tp::flightPathAngleDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_total_acceleration",
        &tp::totalAccelerationDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_total_acceleration_norm",
        &tp::totalAccelerationNormDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_single_acceleration",
        &tp::singleAccelerationDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_single_acceleration_norm",
        &tp::singleAccelerationNormDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_aerodynamic_force_coefficients",
        &tp::aerodynamicForceCoefficientDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_current_body_mass",
        &tp::bodyMassVariable);

    function("dynamics_propagation_setup_dependent_variable_radiation_pressure",
        &tp::radiationPressureDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_periapsis_altitude",
        &tp::periapsisAltitudeVariable);
}

#endif
