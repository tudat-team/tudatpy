/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_dependent_variable_setup.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;
namespace tinterp = tudat::interpolators;
namespace te = tudat::ephemerides;
namespace tni = tudat::numerical_integrators;
namespace trf = tudat::reference_frames;
namespace tmrf = tudat::root_finders;

namespace tudatpy {

    void expose_dependent_variable_setup(py::module &m) {


        py::enum_<tp::PropagationDependentVariables>(m, "PropagationDependentVariables")
                // C++ legacy variable names.
                .value("mach_number_type", tp::PropagationDependentVariables::mach_number_dependent_variable)
                .value("altitude_type", tp::PropagationDependentVariables::altitude_dependent_variable)
                .value("airspeed_type", tp::PropagationDependentVariables::airspeed_dependent_variable)
                .value("local_density_type", tp::PropagationDependentVariables::local_density_dependent_variable)
                .value("relative_speed_type", tp::PropagationDependentVariables::relative_speed_dependent_variable)
                .value("relative_position_type",
                       tp::PropagationDependentVariables::relative_position_dependent_variable)
                .value("relative_distance_type",
                       tp::PropagationDependentVariables::relative_distance_dependent_variable)
                .value("relative_velocity_type",
                       tp::PropagationDependentVariables::relative_velocity_dependent_variable)
                .value("radiation_pressure_type",
                       tp::PropagationDependentVariables::radiation_pressure_dependent_variable)
                .value("total_acceleration_norm_type",
                       tp::PropagationDependentVariables::total_acceleration_norm_dependent_variable)
                .value("single_acceleration_norm_type",
                       tp::PropagationDependentVariables::single_acceleration_norm_dependent_variable)
                .value("total_acceleration_type",
                       tp::PropagationDependentVariables::total_acceleration_dependent_variable)
                .value("single_acceleration_type",
                       tp::PropagationDependentVariables::single_acceleration_dependent_variable)
                .value("aerodynamic_force_coefficients_type",
                       tp::PropagationDependentVariables::aerodynamic_force_coefficients_dependent_variable)
                .value("aerodynamic_moment_coefficients_type",
                       tp::PropagationDependentVariables::aerodynamic_moment_coefficients_dependent_variable)
                .value("rotation_matrix_to_body_fixed_frame_type",
                       tp::PropagationDependentVariables::rotation_matrix_to_body_fixed_frame_variable)
                .value("intermediate_aerodynamic_rotation_matrix_type",
                       tp::PropagationDependentVariables::intermediate_aerodynamic_rotation_matrix_variable)
                .value("relative_body_aerodynamic_orientation_angle_type",
                       tp::PropagationDependentVariables::relative_body_aerodynamic_orientation_angle_variable)
                .value("body_fixed_airspeed_based_velocity_type",
                       tp::PropagationDependentVariables::body_fixed_airspeed_based_velocity_variable)
                .value("total_aerodynamic_g_load_type",
                       tp::PropagationDependentVariables::total_aerodynamic_g_load_variable)
                .value("stagnation_point_heat_flux_type",
                       tp::PropagationDependentVariables::stagnation_point_heat_flux_dependent_variable)
                .value("local_temperature_type",
                       tp::PropagationDependentVariables::local_temperature_dependent_variable)
                .value("geodetic_latitude_type",
                       tp::PropagationDependentVariables::geodetic_latitude_dependent_variable)
                .value("control_surface_deflection_type",
                       tp::PropagationDependentVariables::control_surface_deflection_dependent_variable)
                .value("total_mass_rate_type", tp::PropagationDependentVariables::total_mass_rate_dependent_variables)
                .value("lvlh_to_inertial_frame_rotation_type",
                       tp::PropagationDependentVariables::lvlh_to_inertial_frame_rotation_dependent_variable)
                .value("periapsis_altitude_type",
                       tp::PropagationDependentVariables::periapsis_altitude_dependent_variable)
                .value("total_torque_norm_type",
                       tp::PropagationDependentVariables::total_torque_norm_dependent_variable)
                .value("single_torque_norm_type",
                       tp::PropagationDependentVariables::single_torque_norm_dependent_variable)
                .value("total_torque_type", tp::PropagationDependentVariables::total_torque_dependent_variable)
                .value("single_torque_type", tp::PropagationDependentVariables::single_torque_dependent_variable)
                .value("body_fixed_groundspeed_based_velocity_type",
                       tp::PropagationDependentVariables::body_fixed_groundspeed_based_velocity_variable)
                .value("keplerian_state_type", tp::PropagationDependentVariables::keplerian_state_dependent_variable)
                .value("modified_equinoctial_state_type",
                       tp::PropagationDependentVariables::modified_equinocial_state_dependent_variable)
                .value("spherical_harmonic_acceleration_terms_type",
                       tp::PropagationDependentVariables::spherical_harmonic_acceleration_terms_dependent_variable)
                .value("spherical_harmonic_acceleration_norm_terms_type",
                       tp::PropagationDependentVariables::spherical_harmonic_acceleration_norm_terms_dependent_variable)
                .value("body_fixed_relative_cartesian_position_type",
                       tp::PropagationDependentVariables::body_fixed_relative_cartesian_position)
                .value("body_fixed_relative_spherical_position_type",
                       tp::PropagationDependentVariables::body_fixed_relative_spherical_position)
                .value("total_gravity_field_variation_acceleration_type",
                       tp::PropagationDependentVariables::total_gravity_field_variation_acceleration)
                .value("single_gravity_field_variation_acceleration_type",
                       tp::PropagationDependentVariables::single_gravity_field_variation_acceleration)
                .value("single_gravity_field_variation_acceleration_terms_type",
                       tp::PropagationDependentVariables::single_gravity_field_variation_acceleration_terms)
                .value("acceleration_partial_wrt_body_translational_state_type",
                       tp::PropagationDependentVariables::acceleration_partial_wrt_body_translational_state)
                .value("local_dynamic_pressure_type",
                       tp::PropagationDependentVariables::local_dynamic_pressure_dependent_variable)
                .value("local_aerodynamic_heat_rate_type",
                       tp::PropagationDependentVariables::local_aerodynamic_heat_rate_dependent_variable)
                .value("euler_angles_to_body_fixed_type",
                       tp::PropagationDependentVariables::euler_angles_to_body_fixed_313)
                .value("current_body_mass_type",
                       tp::PropagationDependentVariables::current_body_mass_dependent_variable)
                .value("radiation_pressure_coefficient_type",
                       tp::PropagationDependentVariables::radiation_pressure_coefficient_dependent_variable)
                .export_values();


        py::class_<tp::VariableSettings,
                std::shared_ptr<tp::VariableSettings>>
                VariableSettings_(m, "VariableSettings");

        py::class_<tp::SingleDependentVariableSaveSettings,
                std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
                tp::VariableSettings>(m, "SingleDependentVariableSaveSettings");
//            .def(py::init<
//                 const tp::PropagationDependentVariables,
//                 const std::string &,
//                 const std::string &,
//                 const int>(),
//                 py::arg("dependent_variable_type"),
//                 py::arg("associated_body"),
//                 py::arg("secondary_body") = "",
//                 py::arg("component_idx") = -1);

        //////////////////////////////////////////////////////////////////////////////
        // astro/basic_astro/accelerationModels.h   TODO: This needs to be moved in tudat source.
        //////////////////////////////////////////////////////////////////////////////
        py::class_<
                tp::SingleAccelerationDependentVariableSaveSettings,
                std::shared_ptr<tp::SingleAccelerationDependentVariableSaveSettings>,
                tp::SingleDependentVariableSaveSettings>(m, "SingleAccelerationDependentVariableSaveSettings");
//            .def(py::init<
//                 const tudat::basic_astrodynamics::AvailableAcceleration,
//                 const std::string &,
//                 const std::string &,
//                 const bool,
//                 const int>(),
//                 py::arg("acceleration_model_type"),
//                 py::arg("body_undergoing_acceleration"),
//                 py::arg("body_exerting_acceleration"),
//                 py::arg("use_norm") = 0,
//                 py::arg("component_index") = -1);

        m.def("create",
              &tp::createDependentVariableSaveSettings,
              py::arg("dependent_variable_list"),
              py::arg("print_variable_indices") = 1);

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
              py::arg("body"),
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

        m.def("total_gravity_field_variation_acceleration",
              &tp::totalGravityFieldVariationAccelerationContributionVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"));

        m.def("single_gravity_field_variation_acceleration",
              &tp::singleGravityFieldVariationAccelerationContributionVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              py::arg("deformation_type"),
              py::arg("identifier") = "");

        m.def("single_per_term_gravity_field_variation_acceleration",
              &tp::singleGravityFieldVariationSeparateTermsAccelerationContributionVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              py::arg("component_indices"),
              py::arg("deformation_type"),
              py::arg("identifier") = "");

        m.def("body_fixed_airspeed_velocity",
              &tp::bodyFixedAirspeedBasedVelocityVariable,
              py::arg("body"),
              py::arg("central_body"));

        m.def("body_fixed_groundspeed_velocity",
              &tp::bodyFixedGroundspeedBasedVelocityVariable,
              py::arg("body"),
              py::arg("central_body"));

        // TODO: I propose to name it "inertial_to_body_fixed_rotation_frame"
        m.def("rotation_matrix_to_body_fixed_frame",
              &tp::rotationMatrixToBodyFixedFrameVariable,
              py::arg("body"));

        m.def("lvlh_to_inertial_rotation_matrix",
              &tp::lvlhToInertialFrameRotationMatrixVariable,
              py::arg("body"),
              py::arg("central_body"));

        // TODO: Inconsistent (why Euler angles?)
        m.def("inertial_to_body_fixed_313_euler_angles",
              &tp::eulerAnglesToBodyFixed313Variable,
              py::arg("body"));

        // TODO: inconsistent naming and signature
        m.def("intermediate_aerodynamic_rotation_matrix_variable",
              &tp::intermediateAerodynamicRotationMatrixVariable,
              py::arg("body"),
              py::arg("base_frame"),
              py::arg("target_frame"));

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

        m.def("body_mass",
              &tp::bodyMassVariable,
              py::arg("body"));

        m.def("radiation_pressure_coefficient",
              &tp::radiationPressureCoefficientVariable,
              py::arg("body"),
              py::arg("emitting_body"));

        m.def("dynamic_pressure",
              &tp::dynamicPressureVariable,
              py::arg("body"),
              py::arg("central_body"));

        m.def("aerodynamic_g_load",
              &tp::aerodynamicGLoadVariable,
              py::arg("body"),
              py::arg("central_body"));

        m.def("atmospheric_temperature",
              &tp::atmosphericTemperatureVariable,
              py::arg("body"),
              py::arg("central_body"));

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

} // namespace tudatpy