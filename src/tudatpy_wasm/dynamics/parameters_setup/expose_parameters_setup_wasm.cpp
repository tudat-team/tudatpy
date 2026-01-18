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
#include "../../wasm_module.h"
#include "../../eigen_wasm.h"
#include "../../stl_wasm.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/estimation_setup/estimatableParameterSettings.h>
#include <tudat/simulation/estimation_setup/createEstimatableParameters.h>
#include <tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h>

namespace tep = tudat::estimatable_parameters;
namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_parameters_setup")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_parameters_setup) {
    using namespace emscripten;

    // EstimatebleParametersEnum (note: typo is in original Tudat API)
    enum_<tep::EstimatebleParametersEnum>("dynamics_parameters_setup_EstimatebleParametersEnum")
        .value("arc_wise_initial_body_state", tep::arc_wise_initial_body_state)
        .value("initial_body_state", tep::initial_body_state)
        .value("initial_rotational_body_state", tep::initial_rotational_body_state)
        .value("initial_mass_state", tep::initial_mass_state)
        .value("gravitational_parameter", tep::gravitational_parameter)
        .value("constant_drag_coefficient", tep::constant_drag_coefficient)
        .value("radiation_pressure_coefficient", tep::radiation_pressure_coefficient)
        .value("arc_wise_radiation_pressure_coefficient", tep::arc_wise_radiation_pressure_coefficient)
        .value("spherical_harmonics_cosine_coefficient_block", tep::spherical_harmonics_cosine_coefficient_block)
        .value("spherical_harmonics_sine_coefficient_block", tep::spherical_harmonics_sine_coefficient_block)
        .value("constant_rotation_rate", tep::constant_rotation_rate)
        .value("rotation_pole_position", tep::rotation_pole_position)
        .value("constant_additive_observation_bias", tep::constant_additive_observation_bias)
        .value("arcwise_constant_additive_observation_bias", tep::arcwise_constant_additive_observation_bias)
        .value("constant_relative_observation_bias", tep::constant_relative_observation_bias)
        .value("arcwise_constant_relative_observation_bias", tep::arcwise_constant_relative_observation_bias)
        .value("empirical_acceleration_coefficients", tep::empirical_acceleration_coefficients)
        .value("arc_wise_empirical_acceleration_coefficients", tep::arc_wise_empirical_acceleration_coefficients)
        .value("full_degree_tidal_love_number", tep::full_degree_tidal_love_number)
        .value("single_degree_variable_tidal_love_number", tep::single_degree_variable_tidal_love_number)
        .value("direct_dissipation_tidal_time_lag", tep::direct_dissipation_tidal_time_lag)
        .value("inverse_tidal_quality_factor", tep::inverse_tidal_quality_factor)
        .value("ground_station_position", tep::ground_station_position)
        .value("equivalence_principle_lpi_violation_parameter", tep::equivalence_principle_lpi_violation_parameter)
        .value("ppn_parameter_gamma", tep::ppn_parameter_gamma)
        .value("ppn_parameter_beta", tep::ppn_parameter_beta)
        .value("desaturation_delta_v_values", tep::desaturation_delta_v_values)
        .value("mean_moment_of_inertia", tep::mean_moment_of_inertia)
        .value("periodic_spin_variation", tep::periodic_spin_variation)
        .value("polar_motion_amplitude", tep::polar_motion_amplitude)
        .value("core_factor", tep::core_factor)
        .value("free_core_nutation_rate", tep::free_core_nutation_rate);

    // EstimatableParameterSettings base class
    class_<tep::EstimatableParameterSettings>("dynamics_parameters_setup_EstimatableParameterSettings")
        .smart_ptr<std::shared_ptr<tep::EstimatableParameterSettings>>(
            "shared_ptr_EstimatableParameterSettings");

    // SphericalHarmonicEstimatableParameterSettings
    class_<tep::SphericalHarmonicEstimatableParameterSettings,
           base<tep::EstimatableParameterSettings>>(
        "dynamics_parameters_setup_SphericalHarmonicEstimatableParameterSettings")
        .smart_ptr<std::shared_ptr<tep::SphericalHarmonicEstimatableParameterSettings>>(
            "shared_ptr_SphericalHarmonicEstimatableParameterSettings");

    // EstimatableParameterSet
    class_<tep::EstimatableParameterSet<double>>("dynamics_parameters_setup_EstimatableParameterSet")
        .smart_ptr<std::shared_ptr<tep::EstimatableParameterSet<double>>>(
            "shared_ptr_EstimatableParameterSet")
        .function("getParameterSetSize", &tep::EstimatableParameterSet<double>::getParameterSetSize)
        .function("getEstimatedParameterSetSize", &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize)
        .function("getInitialDynamicalStateParameterSize",
            &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize)
        .function("getNonDynamicalStateParameterSize",
            &tep::EstimatableParameterSet<double>::getNonDynamicalStateParameterSize)
        .function("getParametersDescriptions", &tep::EstimatableParameterSet<double>::getParametersDescriptions);

    // Factory functions for common parameter settings
    function("dynamics_parameters_setup_gravitational_parameter",
        &tep::gravitationalParameter);

    function("dynamics_parameters_setup_constant_drag_coefficient",
        &tep::constantDragCoefficient);

    function("dynamics_parameters_setup_radiation_pressure_coefficient",
        &tep::radiationPressureCoefficient);

    function("dynamics_parameters_setup_constant_rotation_rate",
        &tep::constantRotationRate);

    function("dynamics_parameters_setup_rotation_pole_position",
        &tep::rotationPolePosition);

    function("dynamics_parameters_setup_ground_station_position",
        &tep::groundStationPosition);

    function("dynamics_parameters_setup_ppn_parameter_gamma",
        &tep::ppnParameterGamma);

    function("dynamics_parameters_setup_ppn_parameter_beta",
        &tep::ppnParameterBeta);

    // Factory function to create parameters to estimate
    function("dynamics_parameters_setup_create_parameters_to_estimate",
        select_overload<std::shared_ptr<tep::EstimatableParameterSet<double>>(
            const std::vector<std::shared_ptr<tep::EstimatableParameterSettings>>&,
            const tss::SystemOfBodies&,
            const std::shared_ptr<tudat::propagators::PropagatorSettings<double>>,
            const std::vector<std::shared_ptr<tep::EstimatableParameterSettings>>&)>(
            &tss::createParametersToEstimate<double, double>));
}

#endif
