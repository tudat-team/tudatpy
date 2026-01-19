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
#include <tudat/astro/basic_astro/accelerationModelTypes.h>
#include <tudat/astro/basic_astro/torqueModelTypes.h>

namespace tp = tudat::propagators;
namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("dynamics_propagation_setup_dependent_variable")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_dependent_variable) {
    using namespace emscripten;

    // ========================================================================
    // PropagationDependentVariables enum - Full list
    // ========================================================================
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
        .value("inertial_to_body_fixed_rotation_matrix_variable", tp::inertial_to_body_fixed_rotation_matrix_variable)
        .value("intermediate_aerodynamic_rotation_matrix_variable", tp::intermediate_aerodynamic_rotation_matrix_variable)
        .value("relative_body_aerodynamic_orientation_angle_variable", tp::relative_body_aerodynamic_orientation_angle_variable)
        .value("body_fixed_airspeed_based_velocity_variable", tp::body_fixed_airspeed_based_velocity_variable)
        .value("total_aerodynamic_g_load_variable", tp::total_aerodynamic_g_load_variable)
        .value("stagnation_point_heat_flux_dependent_variable", tp::stagnation_point_heat_flux_dependent_variable)
        .value("local_temperature_dependent_variable", tp::local_temperature_dependent_variable)
        .value("geodetic_latitude_dependent_variable", tp::geodetic_latitude_dependent_variable)
        .value("control_surface_deflection_dependent_variable", tp::control_surface_deflection_dependent_variable)
        .value("total_mass_rate_dependent_variables", tp::total_mass_rate_dependent_variables)
        .value("tnw_to_inertial_frame_rotation_dependent_variable", tp::tnw_to_inertial_frame_rotation_dependent_variable)
        .value("rsw_to_inertial_frame_rotation_dependent_variable", tp::rsw_to_inertial_frame_rotation_dependent_variable)
        .value("periapsis_altitude_dependent_variable", tp::periapsis_altitude_dependent_variable)
        .value("apoapsis_altitude_dependent_variable", tp::apoapsis_altitude_dependent_variable)
        .value("total_torque_norm_dependent_variable", tp::total_torque_norm_dependent_variable)
        .value("single_torque_norm_dependent_variable", tp::single_torque_norm_dependent_variable)
        .value("total_torque_dependent_variable", tp::total_torque_dependent_variable)
        .value("single_torque_dependent_variable", tp::single_torque_dependent_variable)
        .value("body_fixed_groundspeed_based_velocity_variable", tp::body_fixed_groundspeed_based_velocity_variable)
        .value("keplerian_state_dependent_variable", tp::keplerian_state_dependent_variable)
        .value("modified_equinocial_state_dependent_variable", tp::modified_equinocial_state_dependent_variable)
        .value("spherical_harmonic_acceleration_terms_dependent_variable", tp::spherical_harmonic_acceleration_terms_dependent_variable)
        .value("spherical_harmonic_acceleration_norm_terms_dependent_variable", tp::spherical_harmonic_acceleration_norm_terms_dependent_variable)
        .value("body_fixed_relative_cartesian_position", tp::body_fixed_relative_cartesian_position)
        .value("body_fixed_relative_spherical_position", tp::body_fixed_relative_spherical_position)
        .value("total_gravity_field_variation_acceleration", tp::total_gravity_field_variation_acceleration)
        .value("single_gravity_field_variation_acceleration", tp::single_gravity_field_variation_acceleration)
        .value("single_gravity_field_variation_acceleration_terms", tp::single_gravity_field_variation_acceleration_terms)
        .value("acceleration_partial_wrt_body_translational_state", tp::acceleration_partial_wrt_body_translational_state)
        .value("local_dynamic_pressure_dependent_variable", tp::local_dynamic_pressure_dependent_variable)
        .value("euler_angles_to_body_fixed_313", tp::euler_angles_to_body_fixed_313)
        .value("current_body_mass_dependent_variable", tp::current_body_mass_dependent_variable)
        .value("radiation_pressure_coefficient_dependent_variable", tp::radiation_pressure_coefficient_dependent_variable)
        .value("gravity_field_potential_dependent_variable", tp::gravity_field_potential_dependent_variable)
        .value("gravity_field_laplacian_of_potential_dependent_variable", tp::gravity_field_laplacian_of_potential_dependent_variable)
        .value("custom_dependent_variable", tp::custom_dependent_variable);

    // ========================================================================
    // VariableSettings base class
    // ========================================================================
    class_<tp::VariableSettings>("dynamics_propagation_setup_dependent_variable_VariableSettings")
        .smart_ptr<std::shared_ptr<tp::VariableSettings>>("shared_ptr_VariableSettings");

    // ========================================================================
    // SingleDependentVariableSaveSettings
    // ========================================================================
    class_<tp::SingleDependentVariableSaveSettings, base<tp::VariableSettings>>(
        "dynamics_propagation_setup_dependent_variable_SingleDependentVariableSaveSettings")
        .smart_ptr<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>(
            "shared_ptr_SingleDependentVariableSaveSettings")
        .function("getDependentVariableType", &tp::SingleDependentVariableSaveSettings::getDependentVariableType)
        .function("getAssociatedBody", &tp::SingleDependentVariableSaveSettings::getAssociatedBody)
        .function("getSecondaryBody", &tp::SingleDependentVariableSaveSettings::getSecondaryBody)
        .function("getComponentIndex", &tp::SingleDependentVariableSaveSettings::getComponentIndex);

    // ========================================================================
    // SingleAccelerationDependentVariableSaveSettings
    // ========================================================================
    class_<tp::SingleAccelerationDependentVariableSaveSettings, base<tp::SingleDependentVariableSaveSettings>>(
        "dynamics_propagation_setup_dependent_variable_SingleAccelerationDependentVariableSaveSettings")
        .smart_ptr<std::shared_ptr<tp::SingleAccelerationDependentVariableSaveSettings>>(
            "shared_ptr_SingleAccelerationDependentVariableSaveSettings")
        .function("getAccelerationModelType", &tp::SingleAccelerationDependentVariableSaveSettings::getAccelerationModelType);

    // ========================================================================
    // SingleTorqueDependentVariableSaveSettings
    // ========================================================================
    class_<tp::SingleTorqueDependentVariableSaveSettings, base<tp::SingleDependentVariableSaveSettings>>(
        "dynamics_propagation_setup_dependent_variable_SingleTorqueDependentVariableSaveSettings")
        .smart_ptr<std::shared_ptr<tp::SingleTorqueDependentVariableSaveSettings>>(
            "shared_ptr_SingleTorqueDependentVariableSaveSettings")
        .property("torqueModelType_", &tp::SingleTorqueDependentVariableSaveSettings::torqueModelType_);

    // ========================================================================
    // Utility functions
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_get_dependent_variable_id",
        &tp::getDependentVariableId);

    // ========================================================================
    // Atmospheric dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_mach_number",
        &tp::machNumberDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_altitude",
        &tp::altitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_airspeed",
        &tp::airspeedDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_density",
        &tp::densityDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_temperature",
        &tp::localTemperatureDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_dynamic_pressure",
        &tp::localDynamicPressureDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_local_aerodynamic_g_load",
        &tp::totalAerodynamicGLoadDependentVariable);

    // ========================================================================
    // Relative state dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_relative_position",
        &tp::relativePositionDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_velocity",
        &tp::relativeVelocityDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_distance",
        &tp::relativeDistanceDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_relative_speed",
        &tp::relativeSpeedDependentVariable);

    // ========================================================================
    // Orbital element dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_keplerian_state",
        &tp::keplerianStateDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_modified_equinoctial_state",
        &tp::modifiedEquinoctialStateDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_periapsis_altitude",
        &tp::periapsisAltitudeVariable);

    function("dynamics_propagation_setup_dependent_variable_apoapsis_altitude",
        &tp::apoapsisAltitudeVariable);

    // ========================================================================
    // Acceleration dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_total_acceleration",
        &tp::totalAccelerationDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_total_acceleration_norm",
        &tp::totalAccelerationNormDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_single_acceleration",
        &tp::singleAccelerationDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_single_acceleration_norm",
        &tp::singleAccelerationNormDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_spherical_harmonic_terms_acceleration",
        &tp::sphericalHarmonicAccelerationTermsDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_spherical_harmonic_terms_acceleration_norm",
        &tp::sphericalHarmonicAccelerationTermsNormDependentVariable);

    // ========================================================================
    // Torque dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_total_torque",
        &tp::totalTorqueDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_total_torque_norm",
        &tp::totalTorqueNormDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_single_torque",
        &tp::singleTorqueVariable);

    function("dynamics_propagation_setup_dependent_variable_single_torque_norm",
        &tp::singleTorqueNormVariable);

    // ========================================================================
    // Aerodynamic coefficient dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_aerodynamic_force_coefficients",
        &tp::aerodynamicForceCoefficientDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_aerodynamic_moment_coefficients",
        &tp::aerodynamicMomentCoefficientDependentVariable);

    // ========================================================================
    // Angular dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_latitude",
        &tp::latitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_longitude",
        &tp::longitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_geodetic_latitude",
        &tp::geodeticLatitudeDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_heading_angle",
        &tp::headingDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_flight_path_angle",
        &tp::flightPathAngleDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_angle_of_attack",
        &tp::angleOfAttackDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_sideslip_angle",
        &tp::sideslipAngleDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_bank_angle",
        &tp::bankAngleDependentVariable);

    // ========================================================================
    // Radiation pressure dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_radiation_pressure",
        &tp::radiationPressureDependentVariable);

    // ========================================================================
    // Mass dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_current_body_mass",
        &tp::bodyMassVariable);

    function("dynamics_propagation_setup_dependent_variable_total_mass_rate",
        &tp::totalMassRateDependentVariable);

    // ========================================================================
    // Rotation matrix dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_inertial_to_body_fixed_rotation_frame",
        &tp::inertialToBodyFixedRotationMatrixVariable);

    function("dynamics_propagation_setup_dependent_variable_tnw_to_inertial_rotation_matrix",
        &tp::tnwToInertialFrameRotationMatrixVariable);

    function("dynamics_propagation_setup_dependent_variable_rsw_to_inertial_rotation_matrix",
        &tp::rswToInertialFrameRotationMatrixVariable);

    function("dynamics_propagation_setup_dependent_variable_inertial_to_body_fixed_313_euler_angles",
        &tp::eulerAnglesToBodyFixed313Variable);

    // ========================================================================
    // Body-fixed velocity dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_body_fixed_airspeed_velocity",
        &tp::bodyFixedAirspeedBasedVelocityVariable);

    function("dynamics_propagation_setup_dependent_variable_body_fixed_groundspeed_velocity",
        &tp::bodyFixedGroundspeedBasedVelocityVariable);

    // ========================================================================
    // Body-fixed position dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_central_body_fixed_cartesian_position",
        &tp::centralBodyFixedCartesianPositionVariable);

    function("dynamics_propagation_setup_dependent_variable_central_body_fixed_spherical_position",
        &tp::centralBodyFixedSphericalPositionVariable);

    // ========================================================================
    // Gravity field variation dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_total_gravity_field_variation_acceleration",
        &tp::totalGravityFieldVariationAccelerationContributionVariable);

    function("dynamics_propagation_setup_dependent_variable_single_gravity_field_variation_acceleration",
        &tp::singleGravityFieldVariationAccelerationContributionVariable);

    function("dynamics_propagation_setup_dependent_variable_single_per_term_gravity_field_variation_acceleration",
        &tp::singleGravityFieldVariationSeparateTermsAccelerationContributionVariable);

    // ========================================================================
    // Gravity field potential dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_gravity_field_potential",
        &tp::gravityFieldPotentialDependentVariable);

    function("dynamics_propagation_setup_dependent_variable_gravity_field_laplacian_of_potential",
        &tp::gravityFieldLaplacianOfPotentialDependentVariable);

    // ========================================================================
    // Spherical harmonic coefficient variation dependent variables
    // ========================================================================
    function("dynamics_propagation_setup_dependent_variable_total_spherical_harmonic_cosine_coefficient_variations",
        &tp::totalSphericalHarmonicCosineCoefficientVariation);

    function("dynamics_propagation_setup_dependent_variable_total_spherical_harmonic_sine_coefficient_variations",
        &tp::totalSphericalHarmonicSineCoefficientVariation);

    function("dynamics_propagation_setup_dependent_variable_total_spherical_harmonic_cosine_coefficient_variations_from_indices",
        &tp::totalSphericalHarmonicCosineCoefficientVariationFromIndices);

    function("dynamics_propagation_setup_dependent_variable_total_spherical_harmonic_sine_coefficient_variations_from_indices",
        &tp::totalSphericalHarmonicSineCoefficientVariationFromIndices);
}

#endif
