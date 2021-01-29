/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagation_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

#define TUDAT_NAN std::numeric_limits< double >::signaling_NaN()

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;

namespace tudatpy {


void expose_dependent_variable_setup(py::module &m) {


    py::enum_<tp::PropagationDependentVariables>(m, "PropagationDependentVariables")
            // C++ legacy variable names.
            .value("mach_number_type", tp::PropagationDependentVariables::mach_number_dependent_variable)
            .value("altitude_type", tp::PropagationDependentVariables::altitude_dependent_variable)
            .value("airspeed_type", tp::PropagationDependentVariables::airspeed_dependent_variable)
            .value("local_density_type", tp::PropagationDependentVariables::local_density_dependent_variable)
            .value("relative_speed_type", tp::PropagationDependentVariables::relative_speed_dependent_variable)
            .value("relative_position_type", tp::PropagationDependentVariables::relative_position_dependent_variable)
            .value("relative_distance_type", tp::PropagationDependentVariables::relative_distance_dependent_variable)
            .value("relative_velocity_type", tp::PropagationDependentVariables::relative_velocity_dependent_variable)
            .value("radiation_pressure_type", tp::PropagationDependentVariables::radiation_pressure_dependent_variable)
            .value("total_acceleration_norm_type", tp::PropagationDependentVariables::total_acceleration_norm_dependent_variable)
            .value("single_acceleration_norm_type", tp::PropagationDependentVariables::single_acceleration_norm_dependent_variable)
            .value("total_acceleration_type", tp::PropagationDependentVariables::total_acceleration_dependent_variable)
            .value("single_acceleration_type", tp::PropagationDependentVariables::single_acceleration_dependent_variable)
            .value("aerodynamic_force_coefficients_type", tp::PropagationDependentVariables::aerodynamic_force_coefficients_dependent_variable)
            .value("aerodynamic_moment_coefficients_type", tp::PropagationDependentVariables::aerodynamic_moment_coefficients_dependent_variable)
            .value("rotation_matrix_to_body_fixed_frame_type", tp::PropagationDependentVariables::rotation_matrix_to_body_fixed_frame_variable)
            .value("intermediate_aerodynamic_rotation_matrix_type", tp::PropagationDependentVariables::intermediate_aerodynamic_rotation_matrix_variable)
            .value("relative_body_aerodynamic_orientation_angle_type", tp::PropagationDependentVariables::relative_body_aerodynamic_orientation_angle_variable)
            .value("body_fixed_airspeed_based_velocity_type", tp::PropagationDependentVariables::body_fixed_airspeed_based_velocity_variable)
            .value("total_aerodynamic_g_load_type", tp::PropagationDependentVariables::total_aerodynamic_g_load_variable)
            .value("stagnation_point_heat_flux_type", tp::PropagationDependentVariables::stagnation_point_heat_flux_dependent_variable)
            .value("local_temperature_type", tp::PropagationDependentVariables::local_temperature_dependent_variable)
            .value("geodetic_latitude_type", tp::PropagationDependentVariables::geodetic_latitude_dependent_variable)
            .value("control_surface_deflection_type", tp::PropagationDependentVariables::control_surface_deflection_dependent_variable)
            .value("total_mass_rate_type", tp::PropagationDependentVariables::total_mass_rate_dependent_variables)
            .value("lvlh_to_inertial_frame_rotation_type", tp::PropagationDependentVariables::lvlh_to_inertial_frame_rotation_dependent_variable)
            .value("periapsis_altitude_type", tp::PropagationDependentVariables::periapsis_altitude_dependent_variable)
            .value("total_torque_norm_type", tp::PropagationDependentVariables::total_torque_norm_dependent_variable)
            .value("single_torque_norm_type", tp::PropagationDependentVariables::single_torque_norm_dependent_variable)
            .value("total_torque_type", tp::PropagationDependentVariables::total_torque_dependent_variable)
            .value("single_torque_type", tp::PropagationDependentVariables::single_torque_dependent_variable)
            .value("body_fixed_groundspeed_based_velocity_type", tp::PropagationDependentVariables::body_fixed_groundspeed_based_velocity_variable)
            .value("keplerian_state_type", tp::PropagationDependentVariables::keplerian_state_dependent_variable)
            .value("modified_equinocial_state_type", tp::PropagationDependentVariables::modified_equinocial_state_dependent_variable)
            .value("spherical_harmonic_acceleration_terms_type", tp::PropagationDependentVariables::spherical_harmonic_acceleration_terms_dependent_variable)
            .value("spherical_harmonic_acceleration_norm_terms_type", tp::PropagationDependentVariables::spherical_harmonic_acceleration_norm_terms_dependent_variable)
            .value("body_fixed_relative_cartesian_position_type", tp::PropagationDependentVariables::body_fixed_relative_cartesian_position)
            .value("body_fixed_relative_spherical_position_type", tp::PropagationDependentVariables::body_fixed_relative_spherical_position)
            .value("total_gravity_field_variation_acceleration_type", tp::PropagationDependentVariables::total_gravity_field_variation_acceleration)
            .value("single_gravity_field_variation_acceleration_type", tp::PropagationDependentVariables::single_gravity_field_variation_acceleration)
            .value("single_gravity_field_variation_acceleration_terms_type", tp::PropagationDependentVariables::single_gravity_field_variation_acceleration_terms)
            .value("acceleration_partial_wrt_body_translational_state_type", tp::PropagationDependentVariables::acceleration_partial_wrt_body_translational_state)
            .value("local_dynamic_pressure_type", tp::PropagationDependentVariables::local_dynamic_pressure_dependent_variable)
            .value("local_aerodynamic_heat_rate_type", tp::PropagationDependentVariables::local_aerodynamic_heat_rate_dependent_variable)
            .value("euler_angles_to_body_fixed_type", tp::PropagationDependentVariables::euler_angles_to_body_fixed_313)
            .value("current_body_mass_type", tp::PropagationDependentVariables::current_body_mass_dependent_variable)
            .value("radiation_pressure_coefficient_type", tp::PropagationDependentVariables::radiation_pressure_coefficient_dependent_variable)
            .export_values();


    py::class_<tp::VariableSettings,
            std::shared_ptr<tp::VariableSettings>>
            VariableSettings_(m, "VariableSettings");

    py::class_<tp::SingleDependentVariableSaveSettings,
            std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
            tp::VariableSettings>(m, "tp::SingleDependentVariableSaveSettings")
            .def(py::init<
                 const tp::PropagationDependentVariables,
                 const std::string &,
                 const std::string &,
                 const int>(),
                 py::arg("dependent_variable_type"),
                 py::arg("associated_body"),
                 py::arg("secondary_body") = "",
                 py::arg("component_idx") = -1);

    //////////////////////////////////////////////////////////////////////////////
    // astro/basic_astro/accelerationModels.h   TODO: This needs to be moved in tudat source.
    //////////////////////////////////////////////////////////////////////////////
    py::class_<
            tp::SingleAccelerationDependentVariableSaveSettings,
            std::shared_ptr<tp::SingleAccelerationDependentVariableSaveSettings>,
            tp::SingleDependentVariableSaveSettings>(m, "SingleAccelerationDependentVariableSaveSettings")
            .def(py::init<
                 const tudat::basic_astrodynamics::AvailableAcceleration,
                 const std::string &,
                 const std::string &,
                 const bool,
                 const int>(),
                 py::arg("acceleration_model_type"),
                 py::arg("body_undergoing_acceleration"),
                 py::arg("body_exerting_acceleration"),
                 py::arg("use_norm") = 0,
                 py::arg("component_index") = -1);

    m.def("create",
          &tp::createDependentVariableSaveSettings,
          py::arg("dependent_variable_list"),
          py::arg("print_variable_indices")=1);

    m.def("mach_number",
          &tp::machNumberDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("altitude",
          &tp::altitudeDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("airspeed",
          &tp::airspeedDependentVariable,
          py::arg("body"),
          py::arg("body_with_atmosphere"));

    m.def("density",
          &tp::densityDependentVariable,
          py::arg("body"),
          py::arg("body_with_atmosphere"));

    m.def("relative_speed",
          &tp::relativeSpeedDependentVariable,
          py::arg("body"),
          py::arg("relative_body"));

    m.def("relative_position",
          &tp::relativePositionDependentVariable,
          py::arg("body"),
          py::arg("relative_body"));

    m.def("relative_distance",
          &tp::relativeDistanceDependentVariable,
          py::arg("body"),
          py::arg("relative_body"));

    m.def("relative_velocity",
          &tp::relativeVelocityDependentVariable,
          py::arg("body"),
          py::arg("relative_body"));

    m.def("keplerian_state",
          &tp::keplerianStateDependentVariable,
          py::arg("associated_body"),
          py::arg("central_body"));

    m.def("single_acceleration",
          &tp::singleAccelerationDependentVariable,
          py::arg("acceleration_type"),
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"));

    m.def("single_acceleration_norm",
          &tp::singleAccelerationNormDependentVariable,
          py::arg("acceleration_type"),
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"));

    m.def("spherical_harmonic_terms_acceleration",
          &tp::sphericalHarmonicAccelerationTermsDependentVariable,
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"),
          py::arg("component_indices"));

    m.def("spherical_harmonic_terms_acceleration_norm",
          &tp::sphericalHarmonicAccelerationTermsNormDependentVariable,
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"),
          py::arg("component_indices"));

    m.def("total_acceleration",
          &tp::totalAccelerationDependentVariable,
          py::arg("body_undergoing_acceleration"));

    m.def("total_acceleration_norm",
          &tp::totalAccelerationNormDependentVariable,
          py::arg("body_undergoing_acceleration"));

    m.def("aerodynamic_force_coefficients",
          &tp::aerodynamicForceCoefficientDependentVariable,
          py::arg("body"));

    m.def("aerodynamic_moment_coefficients",
          &tp::aerodynamicMomentCoefficientDependentVariable,
          py::arg("body"));

    m.def("latitude",
          &tp::latitudeDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("longitude",
          &tp::longitudeDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("heading_angle",
          &tp::headingDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("flight_path_angle",
          &tp::flightPathAngleDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("angle_of_attack",
          &tp::angleOfAttackDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("sideslip_angle",
          &tp::sideslipAngleDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("bank_angle",
          &tp::bankAngleDependentVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("radiation_pressure",
          &tp::radiationPressureDependentVariable,
          py::arg("body"),
          py::arg("radiating_body"));


    m.def("central_body_fixed_spherical_position",
          &tp::centralBodyFixedSphericalPositionVariable,
          py::arg("associatedBody"),
          py::arg("central_body"));

    m.def("central_body_fixed_cartesian_position",
          &tp::centralBodyFixedCartesianPositionVariable,
          py::arg("associatedBody"),
          py::arg("central_body"));

    m.def("body_mass",
          &tp::bodyMassVariable,
          py::arg("associatedBody"));

    m.def("total_gravity_field_variation_acceleration",
          &tp::totalGravityFieldVariationAccelerationContributionVariable,
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"));

    m.def("single_gravity_field_variation_acceleration",
          &tp::singleGravityFieldVariationAccelerationContributionVariable,
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"),
          py::arg("deformation_type"),
          py::arg("identifier") = "" );

    m.def("single_per_terms_gravity_field_variation_acceleration",
          &tp::singleGravityFieldVariationSeparateTermsAccelerationContributionVariable,
          py::arg("body_undergoing_acceleration"),
          py::arg("body_exerting_acceleration"),
          py::arg("component_indices"),
          py::arg("deformation_type"),
          py::arg("identifier") = "" );

    m.def("rotation_matrix_to_body_fixed_frame",
          &tp::rotationMatrixToBodyFixedFrameVariable,
          py::arg("body") );

    m.def("intermediate_aerodynamic_rotation_matrix_variable",
          &tp::intermediateAerodynamicRotationMatrixVariable,
          py::arg("body"),
          py::arg("base_frame"),
          py::arg("target_frame"));

    m.def("body_fixed_airspeed_velocity",
          &tp::bodyFixedAirspeedBasedVelocityVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("body_fixed_groundspeed_velocity",
          &tp::bodyFixedGroundspeedBasedVelocityVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("lvlh_to_inertial_rotation_matrix",
          &tp::lvlhToInertialFrameRotationMatrixVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("periapsis_altitude",
          &tp::periapsisAltitudeVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("central_body_fixed_spherical_position",
          &tp::centralBodyFixedSphericalPositionVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("central_body_fixed_cartesian_position",
          &tp::centralBodyFixedCartesianPositionVariable,
          py::arg("body"),
          py::arg("central_body"));

    m.def("inertial_to_body_fixed_313_euler_angles",
          &tp::eulerAnglesToBodyFixed313Variable,
          py::arg("body"));


    m.def("body_mass",
          &tp::bodyMassVariable,
          py::arg("body"));

    m.def("radiation_pressure_coefficient",
          &tp::radiationPressureCoefficientVariable,
          py::arg("body"),
          py::arg("emitting_body") );


    //    inline std::shared_ptr< SingleDependentVariableSaveSettings > singleTorqueNormVariable(
    //            const basic_astrodynamics::AvailableTorque torqueModelType,
    //            const std::string& bodyUndergoingTorque,
    //            const std::string& bodyExertingTorque )
    //    {
    //        return std::make_shared< SingleTorqueDependentVariableSaveSettings >(
    //                    torqueModelType, bodyUndergoingTorque, bodyExertingTorque, true );
    //    }

    //    inline std::shared_ptr< SingleDependentVariableSaveSettings > singleTorqueVariable(
    //            const basic_astrodynamics::AvailableTorque torqueModelType,
    //            const std::string& bodyUndergoingTorque,
    //            const std::string& bodyExertingTorque )
    //    {
    //        return std::make_shared< SingleTorqueDependentVariableSaveSettings >(
    //                    torqueModelType, bodyUndergoingTorque, bodyExertingTorque, false );
    //    }


}

void expose_acceleration_setup(py::module &m) {

    /*
     * This contains the addition of IntegratorSettings and AvailableIntegrators
     * and AvailableAccelerations which should be relocated in the tudat source.
     */

    py::enum_<tba::AvailableAcceleration>(m, "AvailableAcceleration")
            .value("undefined_acceleration_type", tba::AvailableAcceleration::undefined_acceleration)
            .value("point_mass_gravity_type", tba::AvailableAcceleration::point_mass_gravity)
            .value("central_gravity_type", tba::AvailableAcceleration::central_gravity)
            .value("aerodynamic_type", tba::AvailableAcceleration::aerodynamic)
            .value("cannonball_radiation_pressure_type", tba::AvailableAcceleration::cannon_ball_radiation_pressure)
            .value("spherical_harmonic_gravity_type", tba::AvailableAcceleration::spherical_harmonic_gravity)
            .value("mutual_spherical_harmonic_gravity_type", tba::AvailableAcceleration::mutual_spherical_harmonic_gravity)
            .value("third_body_point_mass_gravity_type", tba::AvailableAcceleration::third_body_point_mass_gravity)
            .value("third_body_central_gravity_type", tba::AvailableAcceleration::third_body_central_gravity)
            .value("third_body_spherical_harmonic_gravity_type", tba::AvailableAcceleration::third_body_spherical_harmonic_gravity)
            .value("third_body_mutual_spherical_harmonic_gravity_type", tba::AvailableAcceleration::third_body_mutual_spherical_harmonic_gravity)
            .value("thrust_acceleration_type", tba::AvailableAcceleration::thrust_acceleration)
            .value("relativistic_correction_acceleration_type", tba::AvailableAcceleration::relativistic_correction_acceleration)
            .value("empirical_acceleration_type", tba::AvailableAcceleration::empirical_acceleration)
            .value("direct_tidal_dissipation_in_central_body_acceleration_type", tba::AvailableAcceleration::direct_tidal_dissipation_in_central_body_acceleration)
            .value("direct_tidal_dissipation_in_orbiting_body_acceleration_type", tba::AvailableAcceleration::direct_tidal_dissipation_in_orbiting_body_acceleration)
            .value("panelled_radiation_pressure_acceleration_type", tba::AvailableAcceleration::panelled_radiation_pressure_acceleration)
            .value("momentum_wheel_desaturation_acceleration_type", tba::AvailableAcceleration::momentum_wheel_desaturation_acceleration)
            .value("solar_sail_acceleration_type", tba::AvailableAcceleration::solar_sail_acceleration)
            .export_values();

    //////////////////////////////////////////////////////////////////////////////
    // accelerationSettings.h
    //////////////////////////////////////////////////////////////////////////////

    // Unified interface functions for acceleration settings
    //  m.def("acceleration", &tss::acceleration, py::arg("acceleration_type"));
    m.def("point_mass_gravity", &tss::pointMassGravityAcceleration);

    m.def("aerodynamic", &tss::aerodynamicAcceleration);

    m.def("cannonball_radiation_pressure", &tss::cannonBallRadiationPressureAcceleration);

    m.def("spherical_harmonic_gravity", &tss::sphericalHarmonicAcceleration,
          py::arg( "maximum_degree" ),
          py::arg( "maximum_order" ) );

    m.def("mutual_spherical_harmonic_gravity", &tss::mutualSphericalHarmonicAcceleration,
          py::arg( "maximum_degree_body_exerting" ),
          py::arg( "maximum_order_body_exerting" ),
          py::arg( "maximum_degree_body_undergoing" ),
          py::arg( "maximum_order_body_undergoing" ),
          py::arg( "maximum_degree_central_body" ) = 0,
          py::arg( "maximum_order_central_body" ) = 0 );

    m.def("relativistic_correction", &tss::relativisticAccelerationCorrection,
          py::arg( "use_schwarzschild" ) = 0,
          py::arg( "use_lense_thirring" ) = 0,
          py::arg( "use_de_sitter" ) = 0,
          py::arg( "de_sitter_central_body" ) = "",
          py::arg( "lense_thirring_angular_momentum" ) = Eigen::Vector3d::Zero( ) );

    m.def("empirical", &tss::empiricalAcceleration,
          py::arg( "constant_acceleration" ) = Eigen::Vector3d::Zero( ),
          py::arg( "sine_acceleration" ) = Eigen::Vector3d::Zero( ),
          py::arg( "cosine_acceleration" ) = Eigen::Vector3d::Zero( ) );

    // TODO: add overloaded methods
    m.def("thrust_acceleration", py::overload_cast<const std::shared_ptr<tss::ThrustDirectionGuidanceSettings>,
          const std::shared_ptr<tss::ThrustMagnitudeSettings>>(&tss::thrustAcceleration),
          py::arg("thrust_direction_guidance_settings"),
          py::arg("thrust_magnitude_settings"));
    m.def("thrust_acceleration", py::overload_cast<
          const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>>&,
          const std::function<double(const double)>,
          const tss::ThrustFrames,
          const std::string>(&tss::thrustAcceleration),
          py::arg("data_interpolation_settings"),
          py::arg("specific_impulse_function"),
          py::arg("thrust_frame") = tss::ThrustFrames::unspecified_thrust_frame,
          py::arg("central_body") = "");
    m.def("thrust_acceleration", py::overload_cast<
          const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>>&,
          const double,
          const tss::ThrustFrames,
          const std::string>(&tss::thrustAcceleration),
          py::arg("data_interpolation_settings"),
          py::arg("constant_specific_impulse"),
          py::arg("thrust_frame") = tss::ThrustFrames::unspecified_thrust_frame,
          py::arg("central_body") = "");
    m.def("direct_tidal_dissipation_acceleration", &tss::directTidalDissipationAcceleration);
    m.def("momentum_wheel_desaturation_acceleration", &tss::momentumWheelDesaturationAcceleration);

    py::class_<tss::AccelerationSettings,
            std::shared_ptr<tss::AccelerationSettings>>(m, "AccelerationSettings")
            .def(py::init<const tudat::basic_astrodynamics::AvailableAcceleration>(),
                 py::arg("acceleration_type"));

    py::class_<tss::SphericalHarmonicAccelerationSettings,
            std::shared_ptr<tss::SphericalHarmonicAccelerationSettings>,
            tss::AccelerationSettings>(m, "SphericalHarmonicAccelerationSettings")
            .def(py::init<const int, const int>(), py::arg("maximum_degree"),
                 py::arg("maximum_order"));

    py::class_<tss::MutualSphericalHarmonicAccelerationSettings,
            std::shared_ptr<tss::MutualSphericalHarmonicAccelerationSettings>,
            tss::AccelerationSettings>(m, "MutualSphericalHarmonicAccelerationSettings");

    py::class_<tss::EmpiricalAccelerationSettings,
            std::shared_ptr<tss::EmpiricalAccelerationSettings>,
            tss::AccelerationSettings>(m, "EmpiricalAccelerationSettings");


    py::class_<tss::ThrustAccelerationSettings,
            std::shared_ptr<tss::ThrustAccelerationSettings>,
            tss::AccelerationSettings>(m, "ThrustAccelerationSettings")
            .def(py::init<//ctor 1
                 const std::shared_ptr<tss::ThrustDirectionGuidanceSettings>,
                 const std::shared_ptr<tss::ThrustMagnitudeSettings>>(),
                 py::arg("thrust_direction_settings"),
                 py::arg("thrust_magnitude_settings"))
            .def(py::init<//ctor 2
                 const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>> &,
                 const std::function<double(const double)>,
                 const tss::ThrustFrames,
                 const std::string>(),
                 py::arg("data_interpolation_settings"),
                 py::arg("specific_impulse_function"),
                 py::arg("thrust_frame"),
                 py::arg("central_body") = "")
            .def(py::init<//ctor 3
                 const std::shared_ptr<tinterp::DataInterpolationSettings<double, Eigen::Vector3d>> &,
                 const double,
                 const tss::ThrustFrames,
                 const std::string>(),
                 py::arg("data_interpolation_settings"),
                 py::arg("constant_specific_impulse"),
                 py::arg("thrust_frame"),
                 py::arg("central_body") = "")
            .def_readwrite("direction_settings", &tss::ThrustAccelerationSettings::thrustMagnitudeSettings_ );


    //////////////////////////////////////////////////////////////////////////////
    // createThrustModelGuidance.h / createThrustModelGuidance.cpp
    //////////////////////////////////////////////////////////////////////////////
    m.def("get_combined_thrust_direction",
          &tss::getCombinedThrustDirection,
          py::arg("thrust_directions"),
          py::arg("thrust_magnitudes"));

    m.def("get_body_fixed_thrust_direction",
          &tss::getBodyFixedThrustDirection,
          py::arg("thrust_magnitude_settings"),
          py::arg("body_system"),
          py::arg("body_name"));

    m.def("create_thrust_magnitude_wrapper",
          &tss::createThrustMagnitudeWrapper,
          py::arg("thrust_magnitude_settings"),
          py::arg("body_system"),
          py::arg("name_of_body_with_guidance"),
          py::arg("magnitude_update_settings"));

    m.def("update_thrust_magnitude_and_direction",
          &tss::updateThrustMagnitudeAndDirection,
          py::arg("thrust_magnitude_wrapper"),
          py::arg("thrust_direction_guidance"),
          py::arg("current_time"));

    m.def("reset_thrust_magnitude_and_direction_time",
          &tss::resetThrustMagnitudeAndDirectionTime,
          py::arg("thrust_magnitude_wrapper"),
          py::arg("thrust_direction_guidance"),
          py::arg("current_time"));

    //////////////////////////////////////////////////////////////////////////////
    // thrustSettings.h / thrustSettings.cpp
    //////////////////////////////////////////////////////////////////////////////
    py::enum_<tss::ThrustDirectionGuidanceTypes>(m, "ThrustDirectionGuidanceTypes")
            .value("colinear_with_state_segment_thrust_direction", tss::ThrustDirectionGuidanceTypes::colinear_with_state_segment_thrust_direction)
            .value("thrust_direction_from_existing_body_orientation", tss::ThrustDirectionGuidanceTypes::thrust_direction_from_existing_body_orientation)
            .value("custom_thrust_direction", tss::ThrustDirectionGuidanceTypes::custom_thrust_direction)
            .value("custom_thrust_orientation", tss::ThrustDirectionGuidanceTypes::custom_thrust_orientation)
            .value("mee_costate_based_thrust_direction", tss::ThrustDirectionGuidanceTypes::mee_costate_based_thrust_direction);

    m.def("get_propulsion_input_variables",
          &tss::getPropulsionInputVariables,
          py::arg("body_with_guidance") = std::shared_ptr<tss::Body>(),
          py::arg("independent_variables") = std::vector<tudat::propulsion::ThrustIndependentVariables>(),
          py::arg("guidance_input_functions") = std::vector<std::function<double()>>());

    py::class_<
            tss::ThrustDirectionGuidanceSettings,
            std::shared_ptr<tss::ThrustDirectionGuidanceSettings>>(m, "ThrustDirectionGuidanceSettings")
            .def(py::init<
                 const tss::ThrustDirectionGuidanceTypes,
                 const std::string>(),
                 py::arg("thrust_direction_type"),
                 py::arg("relative_body"))
            .def_readonly("thrust_direction_type", &tss::ThrustDirectionGuidanceSettings::thrustDirectionType_)
            .def_readonly("relative_body", &tss::ThrustDirectionGuidanceSettings::relativeBody_);

    py::class_<
            tss::ThrustDirectionFromStateGuidanceSettings,
            std::shared_ptr<tss::ThrustDirectionFromStateGuidanceSettings>,
            tss::ThrustDirectionGuidanceSettings>(m, "ThrustDirectionFromStateGuidanceSettings")
            .def(py::init<const std::string &,
                 const bool,
                 const bool>(),
                 py::arg("central_body"),
                 py::arg("is_colinear_with_velocity"),
                 py::arg("direction_is_opposite_to_vector"))
            .def_readonly("is_colinear_with_velocity", &tss::ThrustDirectionFromStateGuidanceSettings::isColinearWithVelocity_)
            .def_readonly("direction_is_opposite_to_vector", &tss::ThrustDirectionFromStateGuidanceSettings::directionIsOppositeToVector_);

    py::class_<
            tss::CustomThrustDirectionSettings,
            std::shared_ptr<tss::CustomThrustDirectionSettings>,
            tss::ThrustDirectionGuidanceSettings>(m, "CustomThrustDirectionSettings")
            .def(py::init<const std::function<Eigen::Vector3d(const double)>>(),
                 py::arg("thrust_direction_function") );

    py::class_<
            tss::CustomThrustOrientationSettings,
            std::shared_ptr<tss::CustomThrustOrientationSettings>,
            tss::ThrustDirectionGuidanceSettings>(m, "CustomThrustOrientationSettings")
            .def(py::init<const std::function<Eigen::Quaterniond(const double)>>(),
                 py::arg("thrust_orientation_function"))
            .def_readonly("thrust_orientation_function", &tss::CustomThrustOrientationSettings::thrustOrientationFunction_);

    py::class_<
            tss::MeeCostateBasedThrustDirectionSettings,
            std::shared_ptr<tss::MeeCostateBasedThrustDirectionSettings>,
            tss::ThrustDirectionGuidanceSettings>(m, "MeeCostateBasedThrustDirectionSettings")
            .def(py::init<const std::string &,//ctor 1
                 const std::string &,
                 const std::function<Eigen::VectorXd(const double)>>(),
                 py::arg("vehicle_name"),
                 py::arg("central_body_name"),
                 py::arg("costate_function"))
            .def(py::init<const std::string &,//ctor 2
                 const std::string &,
                 std::shared_ptr<tinterp::OneDimensionalInterpolator<double, Eigen::VectorXd>>>(),
                 py::arg("vehicle_name"),
                 py::arg("central_body_name"),
                 py::arg("costate_interpolator"))
            .def(py::init<const std::string &,//ctor 3
                 const std::string &,
                 const Eigen::VectorXd>(),
                 py::arg("vehicle_name"),
                 py::arg("central_body_name"),
                 py::arg("constant_costates"))
            .def_readonly("vehicle_name", &tss::MeeCostateBasedThrustDirectionSettings::vehicleName_)
            .def_readonly("costate_function", &tss::MeeCostateBasedThrustDirectionSettings::costateFunction_);

    py::enum_<tss::ThrustMagnitudeTypes>(m, "ThrustMagnitudeTypes")
            .value("constant_thrust_magnitude", tss::ThrustMagnitudeTypes::constant_thrust_magnitude)
            .value("from_engine_properties_thrust_magnitude", tss::ThrustMagnitudeTypes::from_engine_properties_thrust_magnitude)
            .value("thrust_magnitude_from_time_function", tss::ThrustMagnitudeTypes::thrust_magnitude_from_time_function)
            .value("thrust_magnitude_from_dependent_variables", tss::ThrustMagnitudeTypes::thrust_magnitude_from_dependent_variables)
            .value("bang_bang_thrust_magnitude_from_mee_costates", tss::ThrustMagnitudeTypes::bang_bang_thrust_magnitude_from_mee_costates);

    py::class_<
            tss::ThrustMagnitudeSettings,
            std::shared_ptr<tss::ThrustMagnitudeSettings>>(m, "ThrustMagnitudeSettings")
            .def(py::init<
                 const tss::ThrustMagnitudeTypes,
                 const std::string &>(),
                 py::arg("thrust_magnitude_guidance_type"),
                 py::arg("thrust_origin_id"))
            .def_readonly("thrust_magnitude_guidance_type", &tss::ThrustMagnitudeSettings::thrustMagnitudeGuidanceType_)
            .def_readonly("thrust_origin_id", &tss::ThrustMagnitudeSettings::thrustOriginId_);

    py::class_<
            tss::ConstantThrustMagnitudeSettings,
            std::shared_ptr<tss::ConstantThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "ConstantThrustMagnitudeSettings")
            .def(py::init<
                 const double,
                 const double,
                 const Eigen::Vector3d>(),
                 py::arg("thrust_magnitude"),
                 py::arg("specific_impulse"),
                 py::arg("body_fixed_thrust_direction") = Eigen::Vector3d::UnitX())
            .def_readonly("thrust_magnitude", &tss::ConstantThrustMagnitudeSettings::thrustMagnitude_)
            .def_readonly("specific_impulse", &tss::ConstantThrustMagnitudeSettings::specificImpulse_)
            .def_readonly("body_fixed_thrust_direction", &tss::ConstantThrustMagnitudeSettings::bodyFixedThrustDirection_);

    py::class_<
            tss::FromBodyThrustMagnitudeSettings,
            std::shared_ptr<tss::FromBodyThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "FromBodyThrustMagnitudeSettings")
            .def(py::init<
                 const double,
                 const std::string &>(),
                 py::arg("use_all_engines"),
                 py::arg("thrust_origin"))
            .def_readonly("use_all_engines", &tss::FromBodyThrustMagnitudeSettings::useAllEngines_);

    py::class_<
            tss::FromFunctionThrustMagnitudeSettings,
            std::shared_ptr<tss::FromFunctionThrustMagnitudeSettings>,
            tss::ThrustMagnitudeSettings>(m, "FromFunctionThrustMagnitudeSettings")
            .def(py::init<
                 const std::function< double( const double ) >,
                 const std::function< double( const double ) >,
                 const std::function< bool( const double ) >,
                 const std::function< Eigen::Vector3d( ) >,
                 const std::function< void( const double ) > >(),
                 py::arg("thrust_magnitude_function"),
                 py::arg("specific_impulse_function"),
                 py::arg("is_engine_on_function" ) =
            std::function< bool( const double ) >( [ ]( const double ){ return true; } ),
                 py::arg("body_fixed_thrust_direction" ) =
            std::function< Eigen::Vector3d( ) >( [ ]( ){ return  Eigen::Vector3d::UnitX( ); } ),
                 py::arg("custom_thrust_reset_function" ) = std::function< void( const double ) >( ) );

    m.def("custom_thrust_direction", &tss::customThrustDirectionSettings,
          py::arg( "thrust_direction_function" ) );

    m.def("custom_thrust_magnitude", &tss::fromFunctionThrustMagnitudeSettings,
          py::arg("thrust_magnitude_function"),
          py::arg("specific_impulse_function"),
          py::arg("is_engine_on_function" ) =
            std::function< bool( const double ) >( [ ]( const double ){ return true; } ),
          py::arg("body_fixed_thrust_direction" ) =
            std::function< Eigen::Vector3d( ) >( [ ]( ){ return  Eigen::Vector3d::UnitX( ); } ),
          py::arg("custom_thrust_reset_function" ) = std::function< void( const double ) >( ) );
}

void expose_mass_rate_setup(py::module &m)
{
    py::enum_<tba::AvailableMassRateModels>(m, "AvailableMassRateModels")
            .value("undefined_mass_rate_type", tba::AvailableMassRateModels::undefined_mass_rate_model)
            .value("custom_mass_rate_type", tba::AvailableMassRateModels::custom_mass_rate_model)
            .value("from_thrust_mass_rate_type", tba::AvailableMassRateModels::from_thrust_mass_rate_model)
            .export_values();

    py::class_<tss::MassRateModelSettings,
            std::shared_ptr<tss::MassRateModelSettings>>(m, "MassRateModelSettings")
            .def(py::init<const tudat::basic_astrodynamics::AvailableMassRateModels>(),
                 py::arg("mass_rate_type"));

    py::class_<tss::FromThrustMassModelSettings,
            std::shared_ptr<tss::FromThrustMassModelSettings>,
            tss::MassRateModelSettings>(m, "FromThrustMassModelSettings")
            .def(py::init<const bool, const std::string&>(),
                 py::arg("use_all_thrust_models") = 1,
                 py::arg("associated_thrust_source") = "" );

    m.def("custom", &tss::customMassRate,
          py::arg( "mass_rate_function" ) );

    m.def("from_thrust", &tss::fromThrustMassRate,
          py::arg( "use_all_thrust_models" ) = 1,
          py::arg( "associated_thrust_source" ) = "" );

}

void expose_integrator_setup(py::module &m) {


    py::enum_<tni::AvailableIntegrators>(m, "AvailableIntegrators")
            .value("euler_type", tni::AvailableIntegrators::euler)
            .value("runge_kutta_4_type", tni::AvailableIntegrators::rungeKutta4)
            .value("runge_kutta_variable_step_size_type", tni::AvailableIntegrators::rungeKuttaVariableStepSize)
            .value("bulirsch_stoer_type", tni::AvailableIntegrators::bulirschStoer)
            .value("adams_bashforth_moulton_type", tni::AvailableIntegrators::adamsBashforthMoulton)
            .export_values();

    py::enum_<tni::RungeKuttaCoefficients::CoefficientSets>(m, "RKCoefficientSets")
            .value("rkf_45", tni::RungeKuttaCoefficients::rungeKuttaFehlberg45)
            .value("rkf_56", tni::RungeKuttaCoefficients::rungeKuttaFehlberg56)
            .value("rkf_78", tni::RungeKuttaCoefficients::rungeKuttaFehlberg78)
            .value("rkdp_87", tni::RungeKuttaCoefficients::rungeKutta87DormandPrince)
            .export_values();

    py::enum_<tni::ExtrapolationMethodStepSequences>(m, "ExtrapolationMethodStepSequences")
            .value("bulirsch_stoer_sequence", tni::ExtrapolationMethodStepSequences::bulirsch_stoer_sequence)
            .value("deufelhard_sequence", tni::ExtrapolationMethodStepSequences::deufelhard_sequence)
            .export_values();


    py::class_<
            tni::IntegratorSettings<double>,
            std::shared_ptr<tni::IntegratorSettings<double>>>(m, "IntegratorSettings")
            .def(py::init<
                 const tni::AvailableIntegrators,
                 const double,
                 const double,
                 const int,
                 const bool>(),
                 py::arg("integrator_type"),
                 py::arg("initial_time"),
                 py::arg("initial_time_step"),
                 py::arg("save_frequency") = 1,
                 // TODO: Discuss length of this argument: assess_propagation_termination_condition_during_integration_substeps.
                 py::arg("assess_propagation_termination_condition_during_integration_substeps") = false)
            .def_readwrite("initial_time", &tni::IntegratorSettings<double>::initialTime_ );

    m.def("euler",
          &tni::eulerSettings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false);

    m.def("runge_kutta_4",
          &tni::rungeKutta4Settings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false);

    //    m.def("runge_kutta_variable_step_size",
    //		  &tni::rungeKuttaVariableStepSettings< double >,
    //		  py::arg("initial_time"),
    //		  py::arg("initial_time_step"),
    //		  py::arg("coefficient_set"),
    //		  py::arg("minimum_step_size"),
    //		  py::arg("maximum_step_size"),
    //		  py::arg("relative_error_tolerance"),
    //		  py::arg("absolute_error_tolerance"),
    //		  py::arg("save_frequency") = 1,
    //		  py::arg("assess_termination_on_minor_steps") = false,
    //		  py::arg("safety_factor") = 0.8,
    //		  py::arg("maximum_factor_increase") = 4.0,
    //		  py::arg("minimum_factor_increase") = 0.1);

    //! Function defined twice (here with shorter name)
    m.def("runge_kutta_variable_step_size",
          &tni::rungeKuttaVariableStepSettingsScalarTolerances,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("coefficient_set"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance"),
          py::arg("absolute_error_tolerance"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("safety_factor") = 0.8,
          py::arg("maximum_factor_increase") = 4.0,
          py::arg("minimum_factor_increase") = 0.1 );


    m.def("runge_kutta_variable_step_size_scalar_tolerances",
          &tni::rungeKuttaVariableStepSettingsScalarTolerances,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("coefficient_set"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance"),
          py::arg("absolute_error_tolerance"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("safety_factor") = 0.8,
          py::arg("maximum_factor_increase") = 4.0,
          py::arg("minimum_factor_increase") = 0.1 );

    m.def("runge_kutta_variable_step_size_vector_tolerances",
          &tni::rungeKuttaVariableStepSettingsVectorTolerances,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("coefficient_set"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance"),
          py::arg("absolute_error_tolerance"),
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("safety_factor") = 0.8,
          py::arg("maximum_factor_increase") = 4.0,
          py::arg("minimum_factor_increase") = 0.1);

    m.def("bulirsch_stoer",
          &tni::bulirschStoerIntegratorSettings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("extrapolation_sequence"),
          py::arg("maximum_number_of_steps"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance") = 1.0E-12,
          py::arg("absolute_error_tolerance") = 1.0E-12,
          py::arg("save_frequency") = 1,
          py::arg("check_termination_on_minor_steps") = 0,
          py::arg("safety_factor") = 0.7,
          py::arg("maximum_factor_increase") = 10.0,
          py::arg("minimum_factor_increase") = 10.0 );

    m.def("adams_bashforth_moulton",
          &tni::adamsBashforthMoultonSettings< double >,
          py::arg("initial_time"),
          py::arg("initial_time_step"),
          py::arg("minimum_step_size"),
          py::arg("maximum_step_size"),
          py::arg("relative_error_tolerance") = 1.0E-12,
          py::arg("absolute_error_tolerance") = 1.0E-12,
          py::arg("minimum_order") = 6,
          py::arg("maximum_order") = 11,
          py::arg("save_frequency") = 1,
          py::arg("assess_termination_on_minor_steps") = false,
          py::arg("bandwidth") = 200.0 );
}

void expose_propagator_setup(py::module &m)
{

    py::enum_<tp::TranslationalPropagatorType>(m, "TranslationalPropagatorType")
            .value("undefined_translational_propagator",
                   tp::TranslationalPropagatorType::undefined_translational_propagator)
            .value("cowell",
                   tp::TranslationalPropagatorType::cowell)
            .value("encke",
                   tp::TranslationalPropagatorType::encke)
            .value("gauss_keplerian",
                   tp::TranslationalPropagatorType::gauss_keplerian)
            .value("gauss_modified_equinoctial",
                   tp::TranslationalPropagatorType::gauss_modified_equinoctial)
            .value("unified_state_model_quaternions",
                   tp::TranslationalPropagatorType::unified_state_model_quaternions)
            .value("unified_state_model_modified_rodrigues_parameters",
                   tp::TranslationalPropagatorType::unified_state_model_modified_rodrigues_parameters)
            .value("unified_state_model_exponential_map",
                   tp::unified_state_model_exponential_map)
            .export_values();

    py::class_<tp::DependentVariableSaveSettings,
            std::shared_ptr<tp::DependentVariableSaveSettings>>(m, "DependentVariableSaveSettings")
            .def(py::init<
                 const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>,
                 const bool>(),
                 py::arg("dependent_variables"),
                 py::arg("print_dependent_variable_types") = true);

    py::class_<
            tp::PropagatorSettings<double>,
            std::shared_ptr<tp::PropagatorSettings<double>>>(m, "PropagatorSettings")
            .def("reset_initial_states", &tp::PropagatorSettings<double>::resetInitialStates);

    py::class_<
            tp::SingleArcPropagatorSettings<double>,
            std::shared_ptr<tp::SingleArcPropagatorSettings<double>>,
            tp::PropagatorSettings<double>>(m, "SingleArcPropagatorSettings")
            .def_property("termination_settings",
                          &tp::SingleArcPropagatorSettings<double>::getTerminationSettings,
                          &tp::SingleArcPropagatorSettings<double>::resetTerminationSettings);
    py::class_<
            tp::TranslationalStatePropagatorSettings<double>,
            std::shared_ptr<tp::TranslationalStatePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "TranslationalStatePropagatorSettings")
            .def(// ctor 1
                 py::init<
                 const std::vector<std::string> &,
                 const tba::AccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const std::shared_ptr<tp::PropagationTerminationSettings>,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_models"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_settings"),
                 py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def(// ctor 2
                 py::init<const std::vector<std::string> &,
                 const tss::SelectedAccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const std::shared_ptr<tp::PropagationTerminationSettings>,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_settings"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_settings"),
                 py::arg("propagator") = tp::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def(// ctor 3
                 py::init<const std::vector<std::string> &,
                 const tba::AccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const double,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_models"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_time"),
                 py::arg("propagator") = tp::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def(// ctor 4
                 py::init<const std::vector<std::string> &,
                 const tss::SelectedAccelerationMap &,
                 const std::vector<std::string> &,
                 const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                 const double,
                 const tp::TranslationalPropagatorType,
                 const std::shared_ptr<tp::DependentVariableSaveSettings>,
                 const double>(),
                 py::arg("central_bodies"),
                 py::arg("acceleration_settings"),
                 py::arg("bodies_to_integrate"),
                 py::arg("initial_states"),
                 py::arg("termination_time"),
                 py::arg("propagator") = tp::cowell,
                 py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
                 py::arg("print_interval") = TUDAT_NAN)
            .def("recreate_state_derivative_models", &tp::TranslationalStatePropagatorSettings<double>::resetIntegratedStateModels,
                 py::arg("bodies") )
            .def("get_propagated_state_size", &tp::TranslationalStatePropagatorSettings<double>::getPropagatedStateSize)
            .def("reset_and_recreate_acceleration_models", &tp::TranslationalStatePropagatorSettings<double>::resetAccelerationModelsMap,
                 py::arg("new_acceleration_settings"),
                 py::arg("bodies") )
            .def_property_readonly("acceleration_settings", &tp::TranslationalStatePropagatorSettings<double>::getAccelerationSettingsMap);


    py::class_<
            tp::MultiTypePropagatorSettings<double>,
            std::shared_ptr<tp::MultiTypePropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "MultiTypePropagatorSettings")
            .def("reset_initial_states", &tp::MultiTypePropagatorSettings<double>::resetInitialStates)
            .def("recreate_state_derivative_models", &tp::MultiTypePropagatorSettings<double>::resetIntegratedStateModels,
                 py::arg("bodies") )
            .def("single_type_settings", &tp::MultiTypePropagatorSettings<double>::getSingleTypePropagatorSettings,
                 py::arg("state_type") )
            .def_property_readonly("propagator_settings_per_type", &tp::MultiTypePropagatorSettings<double>::getPropagatorSettingsMap);




    py::class_<
            tp::MassPropagatorSettings<double>,
            std::shared_ptr<tp::MassPropagatorSettings<double>>,
            tp::SingleArcPropagatorSettings<double>>(m, "MassPropagatorSettings");


    m.def("combine_initial_states",
          &tp::createCombinedInitialState<double>,
          py::arg("propagator_settings_per_type") );

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >&,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") =  std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);


    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const double,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tba::AccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const double,
          const tp::TranslationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >&,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_models"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") =  std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tss::SelectedAccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const double,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_settings"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tss::SelectedAccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::shared_ptr<tp::DependentVariableSaveSettings>,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_settings"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("translational",
          py::overload_cast<
          const std::vector< std::string >&,
          const tss::SelectedAccelerationMap&,
          const std::vector< std::string >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr<tp::PropagationTerminationSettings>,
          const tp::TranslationalPropagatorType,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >&,
          const double>(&tp::translationalStatePropagatorSettings<double>),
          py::arg("central_bodies"),
          py::arg("acceleration_settings"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::shared_ptr< tba::MassRateModel > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::vector< std::shared_ptr< tba::MassRateModel > > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const tss::SelectedMassRateModelMap&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_settings"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);




    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::shared_ptr< tba::MassRateModel > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const std::map< std::string, std::vector< std::shared_ptr< tba::MassRateModel > > >&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_models"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);

    m.def("mass",
          py::overload_cast<
          const std::vector< std::string >,
          const tss::SelectedMassRateModelMap&,
          const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >(&tp::massPropagatorSettings<double>),
          py::arg("bodies_with_mass_to_propagate"),
          py::arg("mass_rate_settings"),
          py::arg("initial_body_masses"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >(),
          py::arg("print_interval") = TUDAT_NAN);





    m.def("multitype",
          py::overload_cast<
          const std::vector< std::shared_ptr< tp::SingleArcPropagatorSettings< double > > >,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::shared_ptr< tp::DependentVariableSaveSettings >,
          const double >( &tp::multiTypePropagatorSettings<double> ),
          py::arg("propagator_settings_list"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>( ),
          py::arg("print_interval") = TUDAT_NAN );


    m.def("multitype",
          py::overload_cast<
          const std::vector< std::shared_ptr< tp::SingleArcPropagatorSettings< double > > >,
          const std::shared_ptr< tp::PropagationTerminationSettings >,
          const std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >,
          const double >( &tp::multiTypePropagatorSettings<double> ),
          py::arg("propagator_settings_list"),
          py::arg("termination_settings"),
          py::arg("output_variables") = std::vector< std::shared_ptr< tp::SingleDependentVariableSaveSettings > >( ),
          py::arg("print_interval") = TUDAT_NAN );

    py::class_<tp::PropagationTerminationSettings,
            std::shared_ptr<tp::PropagationTerminationSettings>>
            PropagationTerminationSettings_(m, "PropagationTerminationSettings");

    py::class_<
            tp::PropagationDependentVariableTerminationSettings,
            std::shared_ptr<tp::PropagationDependentVariableTerminationSettings>,
            tp::PropagationTerminationSettings>(m, "PropagationDependentVariableTerminationSettings")
            .def(py::init<
                 const std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
                 const double,
                 const bool,
                 const bool,
                 const std::shared_ptr<tudat::root_finders::RootFinderSettings>>(),
                 py::arg("dependent_variadble_settings"),
                 py::arg("limit_value"),
                 py::arg("use_as_lower_limit"),
                 py::arg("terminate_exactly_on_final_condition") = false,
                 py::arg("termination_root_finder_settings") = nullptr);

    m.def("time_termination",
          &tp::propagationTimeTerminationSettings,
          py::arg("termination_time"),
          py::arg("terminate_exactly_on_final_condition") = false);

    m.def("dependent_variable_termination",
          &tp::propagationDependentVariableTerminationSettings,
          py::arg("dependent_variable_settings"),
          py::arg("limit_value"),
          py::arg("use_as_lower_limit"),
          py::arg("terminate_exactly_on_final_condition") = false,
          py::arg("termination_root_finder_settings") = nullptr);

    m.def("hybrid_termination",
          &tp::propagationHybridTerminationSettings,
          py::arg("termination_settings"),
          py::arg("fulfill_single_condition") );
}

void expose_propagation_setup(py::module &m) {

    py::enum_<tss::ThrustFrames>(m, "ThrustFrames")
            .value("unspecified_thrust_frame", tss::ThrustFrames::unspecified_thrust_frame)
            .value("inertial_thurst_frame", tss::ThrustFrames::inertial_thurst_frame)
            .value("lvlh_thrust_frame", tss::ThrustFrames::lvlh_thrust_frame)
            .export_values();
    /*
   * propagation_setup
   *  ├── accelerationSettings.h
   *  ├── createAccelerationModels.h
   *  ├── createEnvironmentUpdater.h
   *  ├── createMassRateModels.h
   *  ├── createStateDerivativeModel.h
   *  ├── createThrustModelGuidance.h
   *  ├── createTorqueModel.h
   *  ├── dynamicsSimulator.h
   *  ├── environmentUpdater.h
   *  ├── propagationCR3BPFullProblem.h
   *  ├── propagationLambertTargeterFullProblem.h
   *  ├── propagationOutput.h
   *  ├── propagationOutputSettings.h
   *  ├── propagationPatchedConicFullProblem.h
   *  ├── propagationSettings.h
   *  ├── propagationTermination.h
   *  ├── propagationTerminationSettings.h
   *  ├── setNumericallyIntegratedStates.h
   *  ├── thrustSettings.h
   *  └── torqueSettings.h
   *
   * propagation_setup/
   *  ├── createAccelerationModels.cpp
   *  ├── createEnvironmentUpdater.cpp
   *  ├── createMassRateModels.cpp
   *  ├── createStateDerivativeModel.cpp
   *  ├── createThrustModelGuidance.cpp
   *  ├── createTorqueModel.cpp
   *  ├── dynamicsSimulator.cpp
   *  ├── environmentUpdater.cpp
   *  ├── propagationCR3BPFullProblem.cpp
   *  ├── propagationLambertTargeterFullProblem.cpp
   *  ├── propagationOutput.cpp
   *  ├── propagationOutputSettings.cpp
   *  ├── propagationPatchedConicFullProblem.cpp
   *  ├── propagationSettings.cpp
   *  ├── propagationTermination.cpp
   *  ├── setNumericallyIntegratedStates.cpp
   *  └── thrustSettings.cpp
   *
   */


    //////////////////////////////////////////////////////////////////////////////
    // propagationTerminationSettings.h
    //////////////////////////////////////////////////////////////////////////////
    //  py::enum_<tss::PropagationTerminationTypes,
    //            std::shared_ptr<>>
    //  enum PropagationTerminationTypes
    //  {
    //    time_stopping_condition = 0,
    //    cpu_time_stopping_condition = 1,
    //    dependent_variable_stopping_condition = 2,
    //    hybrid_stopping_condition = 3,
    //    custom_stopping_condition = 4
    //  };

    //  py::class_<tss::ThrustAccelerationSettings,
    //             std::shared_ptr<tss::ThrustAccelerationSettings>,
    //             tss::AccelerationSettings>(m, "ThrustAccelerationSettings")
    //      .def(py::init<//ctor 1
    //               const std::shared_ptr<tss::ThrustDirectionGuidanceSettings>,
    //               const std::shared_ptr<tss::ThrustMagnitudeSettings>>(),
    //           py::arg("thrust_direction_settings"),
    //           py::arg("thrust_magnitude_settings"));

    //////////////////////////////////////////////////////////////////////////////
    // createAccelerationModels.cpp
    //////////////////////////////////////////////////////////////////////////////
    m.def("create_acceleration_models",// overload [1/2]
          py::overload_cast<const tss::SystemOfBodies &,
          const tss::SelectedAccelerationMap &,
          const std::vector<std::string> &,
          const std::vector<std::string> &>(
              &tss::createAccelerationModelsMap),
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"));

    m.def("create_acceleration_models",// overload [2/2]
          py::overload_cast<const tss::SystemOfBodies &,
          const tss::SelectedAccelerationMap &,
          const std::map<std::string, std::string> &>(
              &tss::createAccelerationModelsMap),
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("central_bodies"));

    //////////////////////////////////////////////////////////////////////////////
    // dynamicsSimulator.h / dynamicsSimulator.cpp
    //////////////////////////////////////////////////////////////////////////////
    //  m.def("get_initial_state_of_bodies",// overload [1/2]
    //        py::overload_cast<const std::vector<std::string> &,
    //                          const std::vector<std::string> &,
    //                          const tss::SystemOfBodies &,
    //                          const double,
    //                          std::shared_ptr<te::ReferenceFrameManager>>(
    //            &tp::getInitialStatesOfBodies<>));

    py::enum_<tp::IntegratedStateType>(m, "StateType")
            .value("hybrid_type", tp::IntegratedStateType::hybrid)
            .value("translational_type", tp::IntegratedStateType::translational_state)
            .value("rotational_type", tp::IntegratedStateType::rotational_state)
            .value("mass_type", tp::IntegratedStateType::body_mass_state)
            .value("custom_type", tp::IntegratedStateType::custom_state)
            .export_values();

    m.def("get_initial_state_of_bodies",// overload [2/2]
          py::overload_cast<const std::vector<std::string> &,
          const std::vector<std::string> &,
          const tss::SystemOfBodies &,
          const double>(
              &tp::getInitialStatesOfBodies<>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"));

    py::class_<
            tp::SingleArcDynamicsSimulator<double, double>,
            std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double>>>(m, "SingleArcDynamicsSimulator")
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies &,
                 const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double>>,
                 const std::shared_ptr<tp::PropagatorSettings<double>>,
                 const bool,
                 const bool,
                 const bool,
                 const bool,
                 const std::chrono::steady_clock::time_point,
                 const std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>> &>(),
                 py::arg("body_map"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("are_equations_of_motion_to_be_integrated") = true,
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_number_of_function_evaluations") = false,
                 py::arg("initial_clock_time") = std::chrono::steady_clock::now(),
                 py::arg("state_derivative_models") =
            std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>())
            .def("integrate_equations_of_motion",
                 &tp::SingleArcDynamicsSimulator<double, double>::integrateEquationsOfMotion,
                 py::arg("initial_states"))
            .def("get_equations_of_motion_numerical_solution",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
            .def_property_readonly("state_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
            .def("get_equations_of_motion_numerical_solution_raw",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
            .def("get_dependent_variable_history",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory)
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory)
            .def("get_cumulative_computation_time_history",
                 &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistory)
            .def("get_cumulative_number_of_function_evaluations",
                 &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeNumberOfFunctionEvaluations)
            .def("get_equations_of_motion_numerical_solution_base",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionBase)
            .def("get_dependent_variable_numerical_solution_base",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableNumericalSolutionBase)
            .def("get_cumulative_computation_time_history_base",
                 &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistoryBase)
            .def("manually_set_and_process_raw_numerical_equations_of_motion_solution",
                 &tp::SingleArcDynamicsSimulator<double, double>::manuallySetAndProcessRawNumericalEquationsOfMotionSolution,
                 py::arg("equations_of_motion_numerical_solution"),
                 py::arg("dependent_variable_history"),
                 py::arg("process_solution"))
            .def("get_integrator_settings",
                 &tp::SingleArcDynamicsSimulator<double, double>::getIntegratorSettings)
            .def("get_state_derivative_function",
                 &tp::SingleArcDynamicsSimulator<double, double>::getStateDerivativeFunction)
            .def("get_double_state_derivative_function",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDoubleStateDerivativeFunction)
            .def("get_environment_updater",
                 &tp::SingleArcDynamicsSimulator<double, double>::getEnvironmentUpdater)
            .def("get_dynamics_state_derivative",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDynamicsStateDerivative)
            .def("get_propagation_termination_condition",
                 &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationCondition)
            .def("get_integrated_state_processors",
                 &tp::SingleArcDynamicsSimulator<double, double>::getIntegratedStateProcessors)
            .def("get_propagation_termination_reason",
                 &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationReason)
            .def("integration_completed_successfully",
                 &tp::SingleArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully)
            .def("get_dependent_variable_ids",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableIds)
            .def("get_initial_propagation_time",
                 &tp::SingleArcDynamicsSimulator<double, double>::getInitialPropagationTime)
            .def("reset_initial_propagation_time",
                 &tp::SingleArcDynamicsSimulator<double, double>::resetInitialPropagationTime)
            .def("get_dependent_variables_functions",
                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariablesFunctions)
            .def("reset_propagation_termination_conditions",
                 &tp::SingleArcDynamicsSimulator<double, double>::resetPropagationTerminationConditions)
            .def("process_numerical_equations_of_motion_solution",
                 &tp::SingleArcDynamicsSimulator<double, double>::processNumericalEquationsOfMotionSolution);


    //        py::enum_<tp::VariableType>(m, "VariableType")
    //                .value("independent_variable", tp::VariableType::independentVariable)
    //                .value("cpu_time_variable", tp::VariableType::cpuTimeVariable)
    //                .value("state_variable", tp::VariableType::stateVariable)
    //                .value("dependent_variable", tp::VariableType::dependentVariable)
    //                .export_values();



    auto acceleration_setup = m.def_submodule("acceleration");
    expose_acceleration_setup(acceleration_setup);

    auto integrator_setup = m.def_submodule("integrator");
    expose_integrator_setup(integrator_setup);

    auto propagator_setup = m.def_submodule("propagator");
    expose_propagator_setup(propagator_setup);

    auto mass_setup = m.def_submodule("mass");
    expose_mass_rate_setup(mass_setup);

    auto dependent_variable_setup = m.def_submodule("dependent_variable");
    expose_dependent_variable_setup(dependent_variable_setup);
}

}// namespace tudatpy
