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

#include "tudatpy/docstrings.h"
#include <tudat/simulation/propagation_setup.h>

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
namespace numerical_simulation {
namespace propagation_setup {
namespace dependent_variable {

    void expose_dependent_variable_setup(py::module &m) {

        //////////////////////////////////////////////////////////////////////////////////////
        /// ENUMS ////////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////

        py::enum_<tp::PropagationDependentVariables>(m, "PropagationDependentVariables",
                                                     get_docstring("PropagationDependentVariables").c_str())
                // C++ legacy variable names.
                .value("mach_number_type",
                       tp::PropagationDependentVariables::mach_number_dependent_variable,
                       get_docstring("PropagationDependentVariables.mach_number_type").c_str())
                .value("altitude_type",
                       tp::PropagationDependentVariables::altitude_dependent_variable,
                       get_docstring("PropagationDependentVariables.altitude_type").c_str())
                .value("airspeed_type",
                       tp::PropagationDependentVariables::airspeed_dependent_variable,
                       get_docstring("PropagationDependentVariables.airspeed_type").c_str())
                .value("local_density_type",
                       tp::PropagationDependentVariables::local_density_dependent_variable,
                       get_docstring("PropagationDependentVariables.local_density_type").c_str())
                .value("relative_speed_type",
                       tp::PropagationDependentVariables::relative_speed_dependent_variable,
                       get_docstring("PropagationDependentVariables.relative_speed_type").c_str())
                .value("relative_position_type",
                       tp::PropagationDependentVariables::relative_position_dependent_variable,
                       get_docstring("PropagationDependentVariables.relative_position_type").c_str())
                .value("relative_distance_type",
                       tp::PropagationDependentVariables::relative_distance_dependent_variable,
                       get_docstring("PropagationDependentVariables.relative_distance_type").c_str())
                .value("relative_velocity_type",
                       tp::PropagationDependentVariables::relative_velocity_dependent_variable,
                       get_docstring("PropagationDependentVariables.relative_velocity_type").c_str())
                .value("radiation_pressure_type",
                       tp::PropagationDependentVariables::radiation_pressure_dependent_variable,
                       get_docstring("PropagationDependentVariables.radiation_pressure_type").c_str())
                .value("total_acceleration_norm_type",
                       tp::PropagationDependentVariables::total_acceleration_norm_dependent_variable,
                       get_docstring("PropagationDependentVariables.total_acceleration_norm_type").c_str())
                .value("single_acceleration_norm_type",
                       tp::PropagationDependentVariables::single_acceleration_norm_dependent_variable,
                       get_docstring("PropagationDependentVariables.single_acceleration_norm_type").c_str())
                .value("total_acceleration_type",
                       tp::PropagationDependentVariables::total_acceleration_dependent_variable,
                       get_docstring("PropagationDependentVariables.total_acceleration_type").c_str())
                .value("single_acceleration_type",
                       tp::PropagationDependentVariables::single_acceleration_dependent_variable,
                       get_docstring("PropagationDependentVariables.single_acceleration_type").c_str())
                .value("aerodynamic_force_coefficients_type",
                       tp::PropagationDependentVariables::aerodynamic_force_coefficients_dependent_variable,
                       get_docstring("PropagationDependentVariables.aerodynamic_force_coefficients_type").c_str())
                .value("aerodynamic_moment_coefficients_type",
                       tp::PropagationDependentVariables::aerodynamic_moment_coefficients_dependent_variable,
                       get_docstring("PropagationDependentVariables.aerodynamic_moment_coefficients_type").c_str())
                .value("rotation_matrix_to_body_fixed_frame_type",
                       tp::PropagationDependentVariables::inertial_to_body_fixed_rotation_matrix_variable,
                       get_docstring("PropagationDependentVariables.rotation_matrix_to_body_fixed_frame_type").c_str())
                .value("intermediate_aerodynamic_rotation_matrix_type",
                       tp::PropagationDependentVariables::intermediate_aerodynamic_rotation_matrix_variable,
                       get_docstring("PropagationDependentVariables.intermediate_aerodynamic_rotation_matrix_type").c_str())
                .value("relative_body_aerodynamic_orientation_angle_type",
                       tp::PropagationDependentVariables::relative_body_aerodynamic_orientation_angle_variable,
                       get_docstring("PropagationDependentVariables.relative_body_aerodynamic_orientation_angle_type").c_str())
                .value("body_fixed_airspeed_based_velocity_type",
                       tp::PropagationDependentVariables::body_fixed_airspeed_based_velocity_variable,
                       get_docstring("PropagationDependentVariables.body_fixed_airspeed_based_velocity_type").c_str())
                .value("total_aerodynamic_g_load_type",
                       tp::PropagationDependentVariables::total_aerodynamic_g_load_variable,
                       get_docstring("PropagationDependentVariables.total_aerodynamic_g_load_type").c_str())
                .value("stagnation_point_heat_flux_type",
                       tp::PropagationDependentVariables::stagnation_point_heat_flux_dependent_variable,
                       get_docstring("PropagationDependentVariables.stagnation_point_heat_flux_type").c_str())
                .value("local_temperature_type",
                       tp::PropagationDependentVariables::local_temperature_dependent_variable,
                       get_docstring("PropagationDependentVariables.local_temperature_type").c_str())
                .value("geodetic_latitude_type",
                       tp::PropagationDependentVariables::geodetic_latitude_dependent_variable,
                       get_docstring("PropagationDependentVariables.geodetic_latitude_type").c_str())
                .value("control_surface_deflection_type",
                       tp::PropagationDependentVariables::control_surface_deflection_dependent_variable,
                       get_docstring("PropagationDependentVariables.control_surface_deflection_type").c_str())
                .value("total_mass_rate_type",
                       tp::PropagationDependentVariables::total_mass_rate_dependent_variables,
                       get_docstring("PropagationDependentVariables.total_mass_rate_type").c_str())
                .value("tnw_to_inertial_frame_rotation_type",
                       tp::PropagationDependentVariables::tnw_to_inertial_frame_rotation_dependent_variable,
                       get_docstring("PropagationDependentVariables.tnw_to_inertial_frame_rotation_type").c_str())
                .value("rsw_to_inertial_frame_rotation_type",
                       tp::PropagationDependentVariables::rsw_to_inertial_frame_rotation_dependent_variable,
                       get_docstring("PropagationDependentVariables.tnw_to_inertial_frame_rotation_type").c_str())
                .value("periapsis_altitude_type",
                       tp::PropagationDependentVariables::periapsis_altitude_dependent_variable,
                       get_docstring("PropagationDependentVariables.periapsis_altitude_type").c_str())
                .value("total_torque_norm_type",
                       tp::PropagationDependentVariables::total_torque_norm_dependent_variable,
                       get_docstring("PropagationDependentVariables.total_torque_norm_type").c_str())
                .value("single_torque_norm_type",
                       tp::PropagationDependentVariables::single_torque_norm_dependent_variable,
                       get_docstring("PropagationDependentVariables.single_torque_norm_type").c_str())
                .value("total_torque_type",
                       tp::PropagationDependentVariables::total_torque_dependent_variable,
                       get_docstring("PropagationDependentVariables.total_torque_type").c_str())
                .value("single_torque_type",
                       tp::PropagationDependentVariables::single_torque_dependent_variable,
                       get_docstring("PropagationDependentVariables.single_torque_type").c_str())
                .value("body_fixed_groundspeed_based_velocity_type",
                       tp::PropagationDependentVariables::body_fixed_groundspeed_based_velocity_variable,
                       get_docstring("PropagationDependentVariables.body_fixed_groundspeed_based_velocity_type").c_str())
                .value("keplerian_state_type",
                       tp::PropagationDependentVariables::keplerian_state_dependent_variable,
                       get_docstring("PropagationDependentVariables.keplerian_state_type").c_str())
                .value("modified_equinoctial_state_type",
                       tp::PropagationDependentVariables::modified_equinocial_state_dependent_variable,
                       get_docstring("PropagationDependentVariables.modified_equinoctial_state_type").c_str())
                .value("spherical_harmonic_acceleration_terms_type",
                       tp::PropagationDependentVariables::spherical_harmonic_acceleration_terms_dependent_variable,
                       get_docstring("PropagationDependentVariables.spherical_harmonic_acceleration_terms_type").c_str())
                .value("spherical_harmonic_acceleration_norm_terms_type",
                       tp::PropagationDependentVariables::spherical_harmonic_acceleration_norm_terms_dependent_variable,
                       get_docstring("PropagationDependentVariables.spherical_harmonic_acceleration_norm_terms_type").c_str())
                .value("body_fixed_relative_cartesian_position_type",
                       tp::PropagationDependentVariables::body_fixed_relative_cartesian_position,
                       get_docstring("PropagationDependentVariables.body_fixed_relative_cartesian_position_type").c_str())
                .value("body_fixed_relative_spherical_position_type",
                       tp::PropagationDependentVariables::body_fixed_relative_spherical_position,
                       get_docstring("PropagationDependentVariables.body_fixed_relative_spherical_position_type").c_str())
                .value("total_gravity_field_variation_acceleration_type",
                       tp::PropagationDependentVariables::total_gravity_field_variation_acceleration,
                       get_docstring("PropagationDependentVariables.total_gravity_field_variation_acceleration_type").c_str())
                .value("single_gravity_field_variation_acceleration_type",
                       tp::PropagationDependentVariables::single_gravity_field_variation_acceleration,
                       get_docstring("PropagationDependentVariables.single_gravity_field_variation_acceleration_type").c_str())
                .value("single_gravity_field_variation_acceleration_terms_type",
                       tp::PropagationDependentVariables::single_gravity_field_variation_acceleration_terms,
                       get_docstring("PropagationDependentVariables.single_gravity_field_variation_acceleration_terms_type").c_str())
                .value("acceleration_partial_wrt_body_translational_state_type",
                       tp::PropagationDependentVariables::acceleration_partial_wrt_body_translational_state,
                       get_docstring("PropagationDependentVariables.acceleration_partial_wrt_body_translational_state_type").c_str())
                .value("local_dynamic_pressure_type",
                       tp::PropagationDependentVariables::local_dynamic_pressure_dependent_variable,
                       get_docstring("PropagationDependentVariables.local_dynamic_pressure_type").c_str())
//                .value("local_aerodynamic_heat_rate_type",
//                       tp::PropagationDependentVariables::local_aerodynamic_heat_rate_dependent_variable,
//                       get_docstring("PropagationDependentVariables.local_aerodynamic_heat_rate_type").c_str())
                .value("euler_angles_to_body_fixed_type",
                       tp::PropagationDependentVariables::euler_angles_to_body_fixed_313,
                       get_docstring("PropagationDependentVariables.euler_angles_to_body_fixed_type").c_str())
                .value("current_body_mass_type",
                       tp::PropagationDependentVariables::current_body_mass_dependent_variable,
                       get_docstring("PropagationDependentVariables.current_body_mass_type").c_str())
                .value("radiation_pressure_coefficient_type",
                       tp::PropagationDependentVariables::radiation_pressure_coefficient_dependent_variable,
                       get_docstring("PropagationDependentVariables.radiation_pressure_coefficient_type").c_str())
                .value("custom_type",
                       tp::PropagationDependentVariables::custom_dependent_variable,
                       get_docstring("PropagationDependentVariables.custom_type").c_str())
                .value("gravity_field_potential_type",
                       tp::PropagationDependentVariables::gravity_field_potential_dependent_variable)
                .value("gravity_field_laplacian_of_potential_type",
                       tp::PropagationDependentVariables::gravity_field_laplacian_of_potential_dependent_variable)
                .export_values();


        //////////////////////////////////////////////////////////////////////////////////////
        /// CLASSES //////////////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////

        py::class_<tp::VariableSettings,
                std::shared_ptr<tp::VariableSettings>>
                (m, "VariableSettings",
                 get_docstring("VariableSettings").c_str());

        py::class_<tp::SingleDependentVariableSaveSettings,
                std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
                tp::VariableSettings>(m, "SingleDependentVariableSaveSettings",
                                      get_docstring("SingleDependentVariableSaveSettings").c_str());
//            .def(py::init<
//                 const tp::PropagationDependentVariables,
//                 const std::string &,
//                 const std::string &,
//                 const int>(),
//                 py::arg("dependent_variable_type"),
//                 py::arg("associated_body"),
//                 py::arg("secondary_body") = "",
//                 py::arg("component_idx") = -1);

        py::class_<
                tp::SingleAccelerationDependentVariableSaveSettings,
                std::shared_ptr<tp::SingleAccelerationDependentVariableSaveSettings>,
                tp::SingleDependentVariableSaveSettings>(m, "SingleAccelerationDependentVariableSaveSettings",
                                                         get_docstring(
                                                                 "SingleAccelerationDependentVariableSaveSettings").c_str());

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


        //////////////////////////////////////////////////////////////////////////////////////
        /// FREE FUNCTIONS ///////////////////////////////////////////////////////////////////
        //////////////////////////////////////////////////////////////////////////////////////

        m.def("create",
              &tp::createDependentVariableSaveSettings,
              py::arg("dependent_variable_list"),
              py::arg("print_variable_indices") = true);

        m.def("mach_number",
              &tp::machNumberDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("mach_number").c_str());

        m.def("altitude",
              &tp::altitudeDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("altitude").c_str());

        m.def("airspeed",
              &tp::airspeedDependentVariable,
              py::arg("body"),
              py::arg("body_with_atmosphere"),
              get_docstring("airspeed").c_str());

        m.def("density",
              &tp::densityDependentVariable,
              py::arg("body"),
              py::arg("body_with_atmosphere"),
              get_docstring("density").c_str());

        m.def("temperature",
              &tp::localTemperatureDependentVariable,
              py::arg("body"),
              get_docstring("temperature").c_str());

        m.def("dynamic_pressure",
              &tp::localDynamicPressureDependentVariable,
              py::arg("body"),
              get_docstring("dynamic_pressure").c_str());

//        m.def("local_aerodynamic_heat_rate",
//              &tp::localAerodynamicHeatRateDependentVariable,
//              py::arg("body"),
//              get_docstring("local_aerodynamic_heat_rate").c_str());

        m.def("local_aerodynamic_g_load",
              &tp::totalAerodynamicGLoadDependentVariable,
              py::arg("body"),
              get_docstring("local_aerodynamic_g_load").c_str());

        m.def("relative_speed",
              &tp::relativeSpeedDependentVariable,
              py::arg("body"),
              py::arg("relative_body"),
              get_docstring("relative_speed").c_str());

        m.def("relative_position",
              &tp::relativePositionDependentVariable,
              py::arg("body"),
              py::arg("relative_body"),
              get_docstring("relative_position").c_str());

        m.def("relative_distance",
              &tp::relativeDistanceDependentVariable,
              py::arg("body"),
              py::arg("relative_body"),
              get_docstring("relative_distance").c_str());

        m.def("relative_velocity",
              &tp::relativeVelocityDependentVariable,
              py::arg("body"),
              py::arg("relative_body"),
              get_docstring("relative_velocity").c_str());

        m.def("keplerian_state",
              &tp::keplerianStateDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("keplerian_state").c_str());

        m.def("modified_equinoctial_state",
              &tp::modifiedEquinoctialStateDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("modified_equinoctial_state").c_str());

        m.def("single_acceleration",
              &tp::singleAccelerationDependentVariable,
              py::arg("acceleration_type"),
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              get_docstring("single_acceleration").c_str());

        m.def("single_acceleration_norm",
              &tp::singleAccelerationNormDependentVariable,
              py::arg("acceleration_type"),
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              get_docstring("single_acceleration_norm").c_str());

        m.def("total_acceleration_norm",
              &tp::totalAccelerationNormDependentVariable,
              py::arg("body"),
              get_docstring("total_acceleration_norm").c_str());

        m.def("total_acceleration",
              &tp::totalAccelerationDependentVariable,
              py::arg("body"),
              get_docstring("total_acceleration").c_str());

        m.def("single_torque_norm",
              &tp::singleTorqueNormVariable,
              py::arg("torque_type"),
              py::arg("body_undergoing_torque"),
              py::arg("body_exerting_torque"),
              get_docstring("single_torque_norm").c_str());

        m.def("single_torque",
              &tp::singleTorqueVariable,
              py::arg("torque_type"),
              py::arg("body_undergoing_torque"),
              py::arg("body_exerting_torque"),
              get_docstring("single_torque").c_str());

        m.def("total_torque_norm",
              &tp::totalTorqueNormDependentVariable,
              py::arg("body"),
              get_docstring("total_torque_norm").c_str());

        m.def("total_torque",
              &tp::totalTorqueDependentVariable,
              py::arg("body"),
              get_docstring("total_torque").c_str());

        m.def("spherical_harmonic_terms_acceleration",
              &tp::sphericalHarmonicAccelerationTermsDependentVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              py::arg("component_indices"),
              get_docstring("spherical_harmonic_terms_acceleration").c_str());

        m.def("spherical_harmonic_terms_acceleration_norm",
              &tp::sphericalHarmonicAccelerationTermsNormDependentVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              py::arg("component_indices"),
              get_docstring("spherical_harmonic_terms_acceleration_norm").c_str());

        m.def("aerodynamic_force_coefficients",
              &tp::aerodynamicForceCoefficientDependentVariable,
              py::arg("body"),
              get_docstring("aerodynamic_force_coefficients").c_str());

        m.def("aerodynamic_moment_coefficients",
              &tp::aerodynamicMomentCoefficientDependentVariable,
              py::arg("body"),
              get_docstring("aerodynamic_moment_coefficients").c_str());

        m.def("latitude",
              &tp::latitudeDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("latitude").c_str());

        m.def("geodetic_latitude",
              &tp::geodeticLatitudeDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("geodetic_latitude").c_str());

        m.def("longitude",
              &tp::longitudeDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("longitude").c_str());

        m.def("heading_angle",
              &tp::headingDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("heading_angle").c_str());

        m.def("flight_path_angle",
              &tp::flightPathAngleDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("flight_path_angle").c_str());

        m.def("angle_of_attack",
              &tp::angleOfAttackDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("angle_of_attack").c_str());

        m.def("sideslip_angle",
              &tp::sideslipAngleDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("sideslip_angle").c_str());

        m.def("bank_angle",
              &tp::bankAngleDependentVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("bank_angle").c_str());

        m.def("radiation_pressure",
              &tp::radiationPressureDependentVariable,
              py::arg("body"),
              py::arg("radiating_body"),
              get_docstring("radiation_pressure").c_str());

        m.def("total_gravity_field_variation_acceleration",
              &tp::totalGravityFieldVariationAccelerationContributionVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              get_docstring("total_gravity_field_variation_acceleration").c_str());

        m.def("single_gravity_field_variation_acceleration",
              &tp::singleGravityFieldVariationAccelerationContributionVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              py::arg("deformation_type"),
              py::arg("identifier") = "",
              get_docstring("single_gravity_field_variation_acceleration").c_str());

        m.def("single_per_term_gravity_field_variation_acceleration",
              &tp::singleGravityFieldVariationSeparateTermsAccelerationContributionVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"),
              py::arg("component_indices"),
              py::arg("deformation_type"),
              py::arg("identifier") = "",
              get_docstring("single_per_term_gravity_field_variation_acceleration").c_str());

        m.def("body_fixed_airspeed_velocity",
              &tp::bodyFixedAirspeedBasedVelocityVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("body_fixed_airspeed_velocity").c_str());

        m.def("body_fixed_groundspeed_velocity",
              &tp::bodyFixedGroundspeedBasedVelocityVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("body_fixed_groundspeed_velocity").c_str());

        m.def("inertial_to_body_fixed_rotation_frame",
              &tp::inertialToBodyFixedRotationMatrixVariable,
              py::arg("body"),
              get_docstring("inertial_to_body_fixed_rotation_frame").c_str());

        m.def("tnw_to_inertial_rotation_matrix",
              &tp::tnwToInertialFrameRotationMatrixVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("tnw_to_inertial_rotation_matrix").c_str());

        m.def("rsw_to_inertial_rotation_matrix",
              &tp::tnwToInertialFrameRotationMatrixVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("rsw_to_inertial_rotation_matrix").c_str());

        m.def("inertial_to_body_fixed_313_euler_angles",
              &tp::eulerAnglesToBodyFixed313Variable,
              py::arg("body"),
              get_docstring("inertial_to_body_fixed_313_euler_angles").c_str());

        m.def("intermediate_aerodynamic_rotation_matrix_variable",
              &tp::intermediateAerodynamicRotationMatrixVariable,
              py::arg("body"),
              py::arg("base_frame"),
              py::arg("target_frame"),
              get_docstring("intermediate_aerodynamic_rotation_matrix_variable").c_str());

        m.def("periapsis_altitude",
              &tp::periapsisAltitudeVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("periapsis_altitude").c_str());

        m.def("control_surface_deflection",
              &tp::controlSurfaceDeflectionDependentVariable,
              py::arg("body"),
              py::arg("control_surface"),
              get_docstring("control_surface_deflection").c_str());

        m.def("central_body_fixed_spherical_position",
              &tp::centralBodyFixedSphericalPositionVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("central_body_fixed_spherical_position").c_str());

        m.def("central_body_fixed_cartesian_position",
              &tp::centralBodyFixedCartesianPositionVariable,
              py::arg("body"),
              py::arg("central_body"),
              get_docstring("central_body_fixed_cartesian_position").c_str());

        m.def("body_mass",
              &tp::bodyMassVariable,
              py::arg("body"),
              get_docstring("body_mass").c_str());

        m.def("radiation_pressure_coefficient",
              &tp::radiationPressureCoefficientVariable,
              py::arg("body"),
              py::arg("emitting_body"),
              get_docstring("radiation_pressure_coefficient").c_str());

//        m.def("stagnation_point_heat_flux",
//              &tp::stagnationPointHeatFluxDependentVariable,
//              py::arg("body"),
//              get_docstring("stagnation_point_heat_flux").c_str());

        m.def("total_mass_rate",
              &tp::totalMassRateDependentVariable,
              py::arg("body"),
              get_docstring("total_mass_rate").c_str());

        m.def("custom",
              &tp::customDependentVariable,
              py::arg("custom_function"),
              py::arg("variable_size"),
              get_docstring("custom").c_str());

        m.def("gravity_field_potential",
              &tp::gravityFieldPotentialDependentVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"));

        m.def("gravity_field_laplacian_of_potential",
              &tp::gravityFieldLaplacianOfPotentialDependentVariable,
              py::arg("body_undergoing_acceleration"),
              py::arg("body_exerting_acceleration"));

    }

}// namespace dependent_variable
}// namespace propagation_setup
}// namespace numerical_simulation
}// namespace tudatpy
