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

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tp = tudat::propagators;

namespace tudatpy {

void expose_propagation_setup(py::module &m) {

  // Type registration
  //  py::class_<tp::SingleStateTypeDerivative<double, double>,
  //             std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>
  //      SingleStateTypeDerivative(m, "SingleStateTypeDerivative");

  // accelerationSettings.h
  py::class_<tss::AccelerationSettings,
             std::shared_ptr<tss::AccelerationSettings>>(m,
                                                         "AccelerationSettings")
      .def(py::init<const tudat::basic_astrodynamics::AvailableAcceleration>(),
           py::arg("acceleration_type"));

  py::class_<tss::SphericalHarmonicAccelerationSettings,
             std::shared_ptr<tss::SphericalHarmonicAccelerationSettings>,
             tss::AccelerationSettings>(m,
                                        "SphericalHarmonicAccelerationSettings")
      .def(py::init<const int, const int>(), py::arg("maximum_degree"),
           py::arg("maximum_order"));

  //            basic_astrodynamics::AccelerationMap
  //            createAccelerationModelsMap(
  //                    const NamedBodyMap& bodyMap,
  //                    const SelectedAccelerationMap&
  //                    selectedAccelerationPerBody, const std::vector<
  //                    std::string >& propagatedBodies, const std::vector<
  //                    std::string >& centralBodies );

  // createAccelerationModels.cpp
  // createAccelerationModelsMap (overload 1)
  m.def("create_acceleration_models_dict",
        py::overload_cast<
            const tss::NamedBodyMap &, const tss::SelectedAccelerationMap &,
            const std::vector<std::string> &, const std::vector<std::string> &>(
            &tss::createAccelerationModelsMap),
        py::arg("body_dict"), py::arg("selected_acceleration_per_body"),
        py::arg("propagated_bodies"), py::arg("central_bodies"));

  //            basic_astrodynamics::AccelerationMap
  //            createAccelerationModelsMap(
  //                    const NamedBodyMap& bodyMap,
  //                    const SelectedAccelerationMap&
  //                    selectedAccelerationPerBody, const std::map<
  //                    std::string, std::string >& centralBodies );

  // createAccelerationModelsMap (overload 2)
  m.def("create_acceleration_models_dict",
        py::overload_cast<const tss::NamedBodyMap &,
                          const tss::SelectedAccelerationMap &,
                          const std::map<std::string, std::string> &>(
            &tss::createAccelerationModelsMap),
        py::arg("body_dict"), py::arg("selected_acceleration_per_body"),
        py::arg("central_bodies"));

  py::class_<
      tp::SingleArcDynamicsSimulator<double, double>,
      std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double>>>(m, "SingleArcDynamicsSimulator")
      .def(py::init<
               const tudat::simulation_setup::NamedBodyMap &,
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
      .def("get_equations_of_motion_numerical_solution_raw",
           &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
      .def("get_dependent_variable_history",
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
      std::shared_ptr<tp::PropagatorSettings<double>>>
      PropagatorSettings_(m, "PropagatorSettings");

  py::class_<
      tp::SingleArcPropagatorSettings<double>,
      std::shared_ptr<tp::SingleArcPropagatorSettings<double>>,
      tp::PropagatorSettings<double>>
      SingleArcPropagatorSettings_(m, "SingleArcPropagatorSettings");

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
          py::arg("accelerations_map"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_body_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
          py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
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
          py::arg("acceleration_settings_map"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_body_states"),
          py::arg("termination_settings"),
          py::arg("propagator") = tp::cowell,
          py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
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
          py::arg("accelerations_map"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_body_states"),
          py::arg("end_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
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
          py::arg("acceleration_settings_map"),
          py::arg("bodies_to_integrate"),
          py::arg("initial_body_states"),
          py::arg("end_time"),
          py::arg("propagator") = tp::cowell,
          py::arg("dependent_variables_to_save") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
          py::arg("print_interval") = TUDAT_NAN);

  //        py::enum_<tp::VariableType>(m, "VariableType")
  //                .value("independent_variable", tp::VariableType::independentVariable)
  //                .value("cpu_time_variable", tp::VariableType::cpuTimeVariable)
  //                .value("state_variable", tp::VariableType::stateVariable)
  //                .value("dependent_variable", tp::VariableType::dependentVariable)
  //                .export_values();

  //        enum PropagationDependentVariables
  //        {
  //            mach_number_dependent_variable = 0,
  //            altitude_dependent_variable = 1,
  //            airspeed_dependent_variable = 2,
  //            local_density_dependent_variable = 3,
  //            relative_speed_dependent_variable = 4,
  //            relative_position_dependent_variable = 5,
  //            relative_distance_dependent_variable = 6,
  //            relative_velocity_dependent_variable = 7,
  //            radiation_pressure_dependent_variable = 8,
  //            total_acceleration_norm_dependent_variable = 9,
  //            single_acceleration_norm_dependent_variable = 10,
  //            total_acceleration_dependent_variable = 11,
  //            single_acceleration_dependent_variable = 12,
  //            aerodynamic_force_coefficients_dependent_variable = 13,
  //            aerodynamic_moment_coefficients_dependent_variable = 14,
  //            rotation_matrix_to_body_fixed_frame_variable = 15,
  //            intermediate_aerodynamic_rotation_matrix_variable = 16,
  //            relative_body_aerodynamic_orientation_angle_variable = 17,
  //            body_fixed_airspeed_based_velocity_variable = 18,
  //            total_aerodynamic_g_load_variable = 19,
  //            stagnation_point_heat_flux_dependent_variable = 20,
  //            local_temperature_dependent_variable = 21,
  //            geodetic_latitude_dependent_variable = 22,
  //            control_surface_deflection_dependent_variable = 23,
  //            total_mass_rate_dependent_variables = 24,
  //            lvlh_to_inertial_frame_rotation_dependent_variable = 25,
  //            periapsis_altitude_dependent_variable = 26,
  //            total_torque_norm_dependent_variable = 27,
  //            single_torque_norm_dependent_variable = 28,
  //            total_torque_dependent_variable = 29,
  //            single_torque_dependent_variable = 30,
  //            body_fixed_groundspeed_based_velocity_variable = 31,
  //            keplerian_state_dependent_variable = 32,
  //            modified_equinocial_state_dependent_variable = 33,
  //            spherical_harmonic_acceleration_terms_dependent_variable = 34,
  //            body_fixed_relative_cartesian_position = 35,
  //            body_fixed_relative_spherical_position = 36,
  //            total_gravity_field_variation_acceleration = 37,
  //            single_gravity_field_variation_acceleration = 38,
  //            single_gravity_field_variation_acceleration_terms = 39,
  //            acceleration_partial_wrt_body_translational_state = 40,
  //            local_dynamic_pressure_dependent_variable = 41,
  //            local_aerodynamic_heat_rate_dependent_variable = 42
  //        };

  py::enum_<tp::PropagationDependentVariables>(m, "PropagationDependentVariables")
      // C++ legacy variable names.
      .value("mach_number_dependent_variable", tp::PropagationDependentVariables::mach_number_dependent_variable)
      .value("altitude_dependent_variable", tp::PropagationDependentVariables::altitude_dependent_variable)
      .value("airspeed_dependent_variable", tp::PropagationDependentVariables::airspeed_dependent_variable)
      .value("local_density_dependent_variable", tp::PropagationDependentVariables::local_density_dependent_variable)
      .value("relative_speed_dependent_variable", tp::PropagationDependentVariables::relative_speed_dependent_variable)
      .value("relative_position_dependent_variable", tp::PropagationDependentVariables::relative_position_dependent_variable)
      .value("relative_distance_dependent_variable", tp::PropagationDependentVariables::relative_distance_dependent_variable)
      .value("relative_velocity_dependent_variable", tp::PropagationDependentVariables::relative_velocity_dependent_variable)
      .value("radiation_pressure_dependent_variable", tp::PropagationDependentVariables::radiation_pressure_dependent_variable)
      .value("total_acceleration_norm_dependent_variable", tp::PropagationDependentVariables::total_acceleration_norm_dependent_variable)
      .value("single_acceleration_norm_dependent_variable", tp::PropagationDependentVariables::single_acceleration_norm_dependent_variable)
      .value("total_acceleration_dependent_variable", tp::PropagationDependentVariables::total_acceleration_dependent_variable)
      .value("single_acceleration_dependent_variable", tp::PropagationDependentVariables::single_acceleration_dependent_variable)
      .value("aerodynamic_force_coefficients_dependent_variable", tp::PropagationDependentVariables::aerodynamic_force_coefficients_dependent_variable)
      .value("aerodynamic_moment_coefficients_dependent_variable", tp::PropagationDependentVariables::aerodynamic_moment_coefficients_dependent_variable)
      .value("rotation_matrix_to_body_fixed_frame_variable", tp::PropagationDependentVariables::rotation_matrix_to_body_fixed_frame_variable)
      .value("intermediate_aerodynamic_rotation_matrix_variable", tp::PropagationDependentVariables::intermediate_aerodynamic_rotation_matrix_variable)
      .value("relative_body_aerodynamic_orientation_angle_variable", tp::PropagationDependentVariables::relative_body_aerodynamic_orientation_angle_variable)
      .value("body_fixed_airspeed_based_velocity_variable", tp::PropagationDependentVariables::body_fixed_airspeed_based_velocity_variable)
      .value("total_aerodynamic_g_load_variable", tp::PropagationDependentVariables::total_aerodynamic_g_load_variable)
      .value("stagnation_point_heat_flux_dependent_variable", tp::PropagationDependentVariables::stagnation_point_heat_flux_dependent_variable)
      .value("local_temperature_dependent_variable", tp::PropagationDependentVariables::local_temperature_dependent_variable)
      .value("geodetic_latitude_dependent_variable", tp::PropagationDependentVariables::geodetic_latitude_dependent_variable)
      .value("control_surface_deflection_dependent_variable", tp::PropagationDependentVariables::control_surface_deflection_dependent_variable)
      .value("total_mass_rate_dependent_variables", tp::PropagationDependentVariables::total_mass_rate_dependent_variables)
      .value("lvlh_to_inertial_frame_rotation_dependent_variable", tp::PropagationDependentVariables::lvlh_to_inertial_frame_rotation_dependent_variable)
      .value("periapsis_altitude_dependent_variable", tp::PropagationDependentVariables::periapsis_altitude_dependent_variable)
      .value("total_torque_norm_dependent_variable", tp::PropagationDependentVariables::total_torque_norm_dependent_variable)
      .value("single_torque_norm_dependent_variable", tp::PropagationDependentVariables::single_torque_norm_dependent_variable)
      .value("total_torque_dependent_variable", tp::PropagationDependentVariables::total_torque_dependent_variable)
      .value("single_torque_dependent_variable", tp::PropagationDependentVariables::single_torque_dependent_variable)
      .value("body_fixed_groundspeed_based_velocity_variable", tp::PropagationDependentVariables::body_fixed_groundspeed_based_velocity_variable)
      .value("keplerian_state_dependent_variable", tp::PropagationDependentVariables::keplerian_state_dependent_variable)
      .value("modified_equinocial_state_dependent_variable", tp::PropagationDependentVariables::modified_equinocial_state_dependent_variable)
      .value("spherical_harmonic_acceleration_terms_dependent_variable", tp::PropagationDependentVariables::spherical_harmonic_acceleration_terms_dependent_variable)
      .value("body_fixed_relative_cartesian_position", tp::PropagationDependentVariables::body_fixed_relative_cartesian_position)
      .value("body_fixed_relative_spherical_position", tp::PropagationDependentVariables::body_fixed_relative_spherical_position)
      .value("total_gravity_field_variation_acceleration", tp::PropagationDependentVariables::total_gravity_field_variation_acceleration)
      .value("single_gravity_field_variation_acceleration", tp::PropagationDependentVariables::single_gravity_field_variation_acceleration)
      .value("single_gravity_field_variation_acceleration_terms", tp::PropagationDependentVariables::single_gravity_field_variation_acceleration_terms)
      .value("acceleration_partial_wrt_body_translational_state", tp::PropagationDependentVariables::acceleration_partial_wrt_body_translational_state)
      .value("local_dynamic_pressure_dependent_variable", tp::PropagationDependentVariables::local_dynamic_pressure_dependent_variable)
      .value("local_aerodynamic_heat_rate_dependent_variable", tp::PropagationDependentVariables::local_aerodynamic_heat_rate_dependent_variable)

      // Proposed changes / aliases.
      .value("mach_n", tp::PropagationDependentVariables::mach_number_dependent_variable)
      .value("altitude", tp::PropagationDependentVariables::altitude_dependent_variable)
      .value("airspeed", tp::PropagationDependentVariables::airspeed_dependent_variable)
      .value("local_density", tp::PropagationDependentVariables::local_density_dependent_variable)
      .value("rel_speed", tp::PropagationDependentVariables::relative_speed_dependent_variable)
      .value("rel_position", tp::PropagationDependentVariables::relative_position_dependent_variable)
      .value("rel_distance", tp::PropagationDependentVariables::relative_distance_dependent_variable)
      .value("rel_velocity", tp::PropagationDependentVariables::relative_velocity_dependent_variable)
      .value("radiation_pressure", tp::PropagationDependentVariables::radiation_pressure_dependent_variable)
      .value("total_accel_norm", tp::PropagationDependentVariables::total_acceleration_norm_dependent_variable)
      .value("single_accel_norm", tp::PropagationDependentVariables::single_acceleration_norm_dependent_variable)
      .value("total_accel", tp::PropagationDependentVariables::total_acceleration_dependent_variable)
      .value("single_accel", tp::PropagationDependentVariables::single_acceleration_dependent_variable)
      .value("aero_f_coeffs", tp::PropagationDependentVariables::aerodynamic_force_coefficients_dependent_variable)
      .value("aero_m_coeffs", tp::PropagationDependentVariables::aerodynamic_moment_coefficients_dependent_variable)
      .value("bf_rotation_matrix", tp::PropagationDependentVariables::rotation_matrix_to_body_fixed_frame_variable)
      .value("intermediate_aero_rotation_matrix", tp::PropagationDependentVariables::intermediate_aerodynamic_rotation_matrix_variable)
      .value("relative_body_aero_orientation_angle", tp::PropagationDependentVariables::relative_body_aerodynamic_orientation_angle_variable)
      .value("body_fixed_airspeed_based_vel", tp::PropagationDependentVariables::body_fixed_airspeed_based_velocity_variable)
      .value("total_aerodynamic_g_load", tp::PropagationDependentVariables::total_aerodynamic_g_load_variable)
      .value("stagnation_point_heat_flux", tp::PropagationDependentVariables::stagnation_point_heat_flux_dependent_variable)
      .value("local_temperature", tp::PropagationDependentVariables::local_temperature_dependent_variable)
      .value("geodetic_latitude", tp::PropagationDependentVariables::geodetic_latitude_dependent_variable)
      .value("control_surface_deflection", tp::PropagationDependentVariables::control_surface_deflection_dependent_variable)
      .value("total_mass_rate_dependent_vars", tp::PropagationDependentVariables::total_mass_rate_dependent_variables)
      .value("lvlh_to_inertial_frame_rotation", tp::PropagationDependentVariables::lvlh_to_inertial_frame_rotation_dependent_variable)
      .value("periapsis_altitude", tp::PropagationDependentVariables::periapsis_altitude_dependent_variable)
      .value("total_torque_norm", tp::PropagationDependentVariables::total_torque_norm_dependent_variable)
      .value("single_torque_norm", tp::PropagationDependentVariables::single_torque_norm_dependent_variable)
      .value("total_torque", tp::PropagationDependentVariables::total_torque_dependent_variable)
      .value("single_torque", tp::PropagationDependentVariables::single_torque_dependent_variable)
      .value("bf_groundspeed_based_velocity", tp::PropagationDependentVariables::body_fixed_groundspeed_based_velocity_variable)
      .value("kep_state", tp::PropagationDependentVariables::keplerian_state_dependent_variable)
      .value("mee_state", tp::PropagationDependentVariables::modified_equinocial_state_dependent_variable)
      .value("sh_accel_terms", tp::PropagationDependentVariables::spherical_harmonic_acceleration_terms_dependent_variable)
      .value("bf_rel_cartesian_pos", tp::PropagationDependentVariables::body_fixed_relative_cartesian_position)
      .value("bf_rel_spherical_pos", tp::PropagationDependentVariables::body_fixed_relative_spherical_position)
      .value("total_gravity_field_variation_accel", tp::PropagationDependentVariables::total_gravity_field_variation_acceleration)
      .value("single_gravity_field_variation_accel", tp::PropagationDependentVariables::single_gravity_field_variation_acceleration)
      .value("single_gravity_field_variation_accel_terms", tp::PropagationDependentVariables::single_gravity_field_variation_acceleration_terms)
      .value("accel_partial_wrt_body_translational_state", tp::PropagationDependentVariables::acceleration_partial_wrt_body_translational_state)
      .value("local_dynamic_pressure", tp::PropagationDependentVariables::local_dynamic_pressure_dependent_variable)
      .value("local_aerodynamic_heat_rate", tp::PropagationDependentVariables::local_aerodynamic_heat_rate_dependent_variable)
      .export_values();

  py::class_<tp::VariableSettings,
             std::shared_ptr<tp::VariableSettings>>
      VariableSettings_(m, "VariableSettings");

  py::class_<tp::SingleDependentVariableSaveSettings,
             std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
             tp::VariableSettings>(m, "SingleDependentVariableSaveSettings")
      .def(py::init<
               const tp::PropagationDependentVariables,
               const std::string &,
               const std::string &,
               const int>(),
           py::arg("dependent_variable_type"),
           py::arg("associated_body"),
           py::arg("secondary_body") = "",
           py::arg("component_idx") = -1);

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
           py::arg("dependent_variable_settings"),
           py::arg("limit_value"),
           py::arg("use_as_lower_limit"),
           py::arg("terminate_exactly_on_final_condition") = false,
           py::arg("termination_root_finder_settings") = nullptr);
}
}// namespace tudatpy