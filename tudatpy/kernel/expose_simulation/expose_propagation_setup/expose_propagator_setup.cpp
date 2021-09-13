/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagator_setup.h"

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
namespace simulation {
namespace propagation_setup {
namespace propagator {

    void expose_propagator_setup(py::module &m) {

        // ENUMS
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

        py::enum_<tp::RotationalPropagatorType>(m, "RotationalPropagatorType")
                .value("undefined_rotational_propagator",
                       tp::RotationalPropagatorType::undefined_rotational_propagator)
                .value("quaternions",
                       tp::RotationalPropagatorType::quaternions)
                .value("modified_rodrigues_parameters",
                       tp::RotationalPropagatorType::modified_rodrigues_parameters)
                .value("exponential_map",
                       tp::RotationalPropagatorType::exponential_map)
                .export_values();

        // TODO: why is this enum defined here and not in Tudat?
          enum PropagationTerminationTypes
          {
            time_stopping_condition = 0,
            cpu_time_stopping_condition = 1,
            dependent_variable_stopping_condition = 2,
            hybrid_stopping_condition = 3,
            custom_stopping_condition = 4
          };
        // TODO: expose enum
//        py::enum_<tss::PropagationTerminationTypes,
//                std::shared_ptr<>>;

        py::enum_<tp::IntegratedStateType>(m, "StateType")
                .value("hybrid_type", tp::IntegratedStateType::hybrid)
                .value("translational_type", tp::IntegratedStateType::translational_state)
                .value("rotational_type", tp::IntegratedStateType::rotational_state)
                .value("mass_type", tp::IntegratedStateType::body_mass_state)
                .value("custom_type", tp::IntegratedStateType::custom_state)
                .export_values();


        //        py::enum_<tp::VariableType>(m, "VariableType")
        //                .value("independent_variable", tp::VariableType::independentVariable)
        //                .value("cpu_time_variable", tp::VariableType::cpuTimeVariable)
        //                .value("state_variable", tp::VariableType::stateVariable)
        //                .value("dependent_variable", tp::VariableType::dependentVariable)
        //                .export_values();


        // CLASSES

        // NOTE: this class does not strictly belong to "propagator settings".
        py::class_<
                tp::SingleArcDynamicsSimulator<double, double>,
                std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double>>>(m,
                                                                                 "SingleArcDynamicsSimulator")
//                .def(py::init<
//                             const tudat::simulation_setup::SystemOfBodies &,
//                             const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double>>,
//                             const std::shared_ptr<tp::PropagatorSettings<double>>,
//                             const bool,
//                             const bool,
//                             const bool,
//                             const bool,
//                             const std::chrono::steady_clock::time_point,
//                             const std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>> &,
//                             const bool>(),
//                     py::arg("body_map"),
//                     py::arg("integrator_settings"),
//                     py::arg("propagator_settings"),
//                     py::arg("are_equations_of_motion_to_be_integrated") = true,
//                     py::arg("clear_numerical_solutions") = false,
//                     py::arg("set_integrated_result") = false,
//                     py::arg("print_number_of_function_evaluations") = false,
//                     py::arg("initial_clock_time") = std::chrono::steady_clock::now(),
//                     py::arg("state_derivative_models") =
//                             std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>(),
//                     py::arg("print_dependent_variable_data") = true)
                .def("integrate_equations_of_motion",
                     &tp::SingleArcDynamicsSimulator<double, double>::integrateEquationsOfMotion,
                     py::arg("initial_states"))
                .def_property_readonly("equations_of_motion_numerical_solution",
                     &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
                .def_property_readonly("state_history",
                                       &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
                .def_property_readonly("equations_of_motion_numerical_solution_raw",
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
                     &tp::SingleArcDynamicsSimulator<double, double>::processNumericalEquationsOfMotionSolution)
                .def("suppress_dependent_variable_terminal_printing",
                     &tp::SingleArcDynamicsSimulator<double, double>::suppressDependentVariableDataPrinting)
                .def("enable_dependent_variable_terminal_printing",
                     &tp::SingleArcDynamicsSimulator<double, double>::enableDependentVariableDataPrinting);

        py::class_<tp::DependentVariableSaveSettings,
                std::shared_ptr<tp::DependentVariableSaveSettings>>(m, "DependentVariableSaveSettings",
                                                                    get_docstring(
                                                                            "DependentVariableSaveSettings").c_str());
//                .def(py::init<
//                             const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings>>,
//                             const bool>(),
//                     py::arg("dependent_variables"),
//                     py::arg("print_dependent_variable_types") = true);

        py::class_<
                tp::PropagatorSettings<double>,
                std::shared_ptr<tp::PropagatorSettings<double>>>(m, "PropagatorSettings",
                                                                 get_docstring("PropagatorSettings").c_str())
                .def("reset_initial_states", &tp::PropagatorSettings<double>::resetInitialStates,
                     py::arg("initial_states"));

        py::class_<
                tp::MultiArcPropagatorSettings<double>,
                std::shared_ptr<tp::MultiArcPropagatorSettings<double>>,
                tp::PropagatorSettings<double>>(m, "MultiArcPropagatorSettings",
                                                get_docstring("MultiArcPropagatorSettings").c_str());

        py::class_<
                tp::SingleArcPropagatorSettings<double>,
                std::shared_ptr<tp::SingleArcPropagatorSettings<double>>,
                tp::PropagatorSettings<double>>(m, "SingleArcPropagatorSettings",
                                                get_docstring("SingleArcPropagatorSettings").c_str())
                .def_property("termination_settings",
                              &tp::SingleArcPropagatorSettings<double>::getTerminationSettings,
                              &tp::SingleArcPropagatorSettings<double>::resetTerminationSettings);

        py::class_<
                tp::TranslationalStatePropagatorSettings<double>,
                std::shared_ptr<tp::TranslationalStatePropagatorSettings<double>>,
                tp::SingleArcPropagatorSettings<double>>(m, "TranslationalStatePropagatorSettings",
                                                         get_docstring("TranslationalStatePropagatorSettings").c_str())
//                .def(// ctor 1
//                        py::init<
//                                const std::vector<std::string> &,
//                                const tba::AccelerationMap &,
//                                const std::vector<std::string> &,
//                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                const std::shared_ptr<tp::PropagationTerminationSettings>,
//                                const tp::TranslationalPropagatorType,
//                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
//                                const double>(),
//                        py::arg("central_bodies"),
//                        py::arg("acceleration_models"),
//                        py::arg("bodies_to_integrate"),
//                        py::arg("initial_states"),
//                        py::arg("termination_settings"),
//                        py::arg("propagator") = tp::TranslationalPropagatorType::cowell,
//                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                        py::arg("print_interval") = TUDAT_NAN)
//                .def(// ctor 2
//                        py::init<const std::vector<std::string> &,
//                                const tss::SelectedAccelerationMap &,
//                                const std::vector<std::string> &,
//                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                const std::shared_ptr<tp::PropagationTerminationSettings>,
//                                const tp::TranslationalPropagatorType,
//                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
//                                const double>(),
//                        py::arg("central_bodies"),
//                        py::arg("acceleration_settings"),
//                        py::arg("bodies_to_integrate"),
//                        py::arg("initial_states"),
//                        py::arg("termination_settings"),
//                        py::arg("propagator") = tp::cowell,
//                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                        py::arg("print_interval") = TUDAT_NAN)
//                .def(// ctor 3
//                        py::init<const std::vector<std::string> &,
//                                const tba::AccelerationMap &,
//                                const std::vector<std::string> &,
//                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                const double,
//                                const tp::TranslationalPropagatorType,
//                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
//                                const double>(),
//                        py::arg("central_bodies"),
//                        py::arg("acceleration_models"),
//                        py::arg("bodies_to_integrate"),
//                        py::arg("initial_states"),
//                        py::arg("termination_time"),
//                        py::arg("propagator") = tp::cowell,
//                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                        py::arg("print_interval") = TUDAT_NAN)
//                .def(// ctor 4
//                        py::init<const std::vector<std::string> &,
//                                const tss::SelectedAccelerationMap &,
//                                const std::vector<std::string> &,
//                                const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
//                                const double,
//                                const tp::TranslationalPropagatorType,
//                                const std::shared_ptr<tp::DependentVariableSaveSettings>,
//                                const double>(),
//                        py::arg("central_bodies"),
//                        py::arg("acceleration_settings"),
//                        py::arg("bodies_to_integrate"),
//                        py::arg("initial_states"),
//                        py::arg("termination_time"),
//                        py::arg("propagator") = tp::cowell,
//                        py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
//                        py::arg("print_interval") = TUDAT_NAN)
                .def("recreate_state_derivative_models",
                     &tp::TranslationalStatePropagatorSettings<double>::resetIntegratedStateModels,
                     py::arg("bodies"))
                .def("get_propagated_state_size",
                     &tp::TranslationalStatePropagatorSettings<double>::getPropagatedStateSize)
                .def("reset_and_recreate_acceleration_models",
                     &tp::TranslationalStatePropagatorSettings<double>::resetAccelerationModelsMap,
                     py::arg("new_acceleration_settings"),
                     py::arg("bodies"))
                .def_property_readonly("acceleration_settings",
                                       &tp::TranslationalStatePropagatorSettings<double>::getAccelerationSettingsMap);


        py::class_<
                tp::MultiTypePropagatorSettings<double>,
                std::shared_ptr<tp::MultiTypePropagatorSettings<double>>,
                tp::SingleArcPropagatorSettings<double>>(m, "MultiTypePropagatorSettings",
                                                         get_docstring("MultiTypePropagatorSettings").c_str())
                .def("reset_initial_states", &tp::MultiTypePropagatorSettings<double>::resetInitialStates,
                     py::arg("initial_states"))
                .def("recreate_state_derivative_models",
                     &tp::MultiTypePropagatorSettings<double>::resetIntegratedStateModels,
                     py::arg("bodies"))
                .def("single_type_settings", &tp::MultiTypePropagatorSettings<double>::getSingleTypePropagatorSettings,
                     py::arg("state_type"))
                .def_property_readonly("propagator_settings_per_type",
                                       &tp::MultiTypePropagatorSettings<double>::getPropagatorSettingsMap);

        py::class_<
                tp::RotationalStatePropagatorSettings<double>,
                std::shared_ptr<tp::RotationalStatePropagatorSettings<double>>,
                tp::SingleArcPropagatorSettings<double>>(m, "RotationalStatePropagatorSettings",
                                                         get_docstring("RotationalStatePropagatorSettings").c_str());

        py::class_<
                tp::MassPropagatorSettings<double>,
                std::shared_ptr<tp::MassPropagatorSettings<double>>,
                tp::SingleArcPropagatorSettings<double>>(m, "MassPropagatorSettings",
                                                         get_docstring("MassPropagatorSettings").c_str());

        // FREE FUNCTIONS

        // NOTE: the following 5 functions do not strictly belong to the actual "propagator settings".

        m.def("get_initial_state_of_bodies",
              py::overload_cast<const std::vector<std::string> &,
                      const std::vector<std::string> &,
                      const tss::SystemOfBodies &,
                      const double>(
                      &tp::getInitialStatesOfBodies<>),
              py::arg("bodies_to_propagate"),
              py::arg("central_bodies"),
              py::arg("body_system"),
              py::arg("initial_time"),
              get_docstring("get_initial_state_of_bodies").c_str());

        m.def("combine_initial_states",
              &tp::createCombinedInitialState<double>,
              py::arg("propagator_settings_per_type"),
              get_docstring("combine_initial_states").c_str());

        // First overload
        m.def("create_acceleration_models",
              py::overload_cast<const tss::SystemOfBodies &,
                      const tss::SelectedAccelerationMap &,
                      const std::map<std::string, std::string> &>(
                      &tss::createAccelerationModelsMap),
              py::arg("body_system"),
              py::arg("selected_acceleration_per_body"),
              py::arg("central_bodies"),
              get_docstring("create_acceleration_models", 0).c_str());

        // Second overload
        m.def("create_acceleration_models",
              py::overload_cast<const tss::SystemOfBodies &,
                      const tss::SelectedAccelerationMap &,
                      const std::vector<std::string> &,
                      const std::vector<std::string> &>(
                      &tss::createAccelerationModelsMap),
              py::arg("body_system"),
              py::arg("selected_acceleration_per_body"),
              py::arg("bodies_to_propagate"),
              py::arg("central_bodies"),
              get_docstring("create_acceleration_models", 1).c_str());

        m.def("create_torque_models",
              &tss::createTorqueModelsMap,
              py::arg("body_system"),
              py::arg("selected_torque_per_body"),
              py::arg("bodies_to_propagate"),
              get_docstring("create_torque_models").c_str());


        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tba::AccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
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
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 0).c_str());


        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tba::AccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const tp::TranslationalPropagatorType,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> > &,
                      const double>(&tp::translationalStatePropagatorSettings<double>),
              py::arg("central_bodies"),
              py::arg("acceleration_models"),
              py::arg("bodies_to_integrate"),
              py::arg("initial_states"),
              py::arg("termination_settings"),
              py::arg("propagator") = tp::cowell,
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 1).c_str());


        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tba::AccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
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
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 2).c_str());


        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tba::AccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const double,
                      const tp::TranslationalPropagatorType,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> > &,
                      const double>(&tp::translationalStatePropagatorSettings<double>),
              py::arg("central_bodies"),
              py::arg("acceleration_models"),
              py::arg("bodies_to_integrate"),
              py::arg("initial_states"),
              py::arg("termination_time"),
              py::arg("propagator") = tp::cowell,
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 3).c_str());

        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tss::SelectedAccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
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
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 4).c_str());

        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tss::SelectedAccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
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
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 5).c_str());

        m.def("translational",
              py::overload_cast<
                      const std::vector<std::string> &,
                      const tss::SelectedAccelerationMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const tp::TranslationalPropagatorType,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> > &,
                      const double>(&tp::translationalStatePropagatorSettings<double>),
              py::arg("central_bodies"),
              py::arg("acceleration_settings"),
              py::arg("bodies_to_integrate"),
              py::arg("initial_states"),
              py::arg("termination_settings"),
              py::arg("propagator") = tp::cowell,
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("translational", 6).c_str());

        m.def("rotational",
              py::overload_cast<
                      const tba::TorqueModelMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const tp::RotationalPropagatorType,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
                      const double>(&tp::rotationalStatePropagatorSettings<double>),
              py::arg("torque_models"),
              py::arg("bodies_to_integrate"),
              py::arg("initial_states"),
              py::arg("termination_settings"),
              py::arg("propagator") = tp::quaternions,
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("rotational", 0).c_str());

        m.def("rotational",
              py::overload_cast<
                      const tss::SelectedTorqueMap &,
                      const std::vector<std::string> &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const tp::RotationalPropagatorType,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
                      const double>(&tp::rotationalStatePropagatorSettings<double>),
              py::arg("torque_settings"),
              py::arg("bodies_to_integrate"),
              py::arg("initial_states"),
              py::arg("termination_settings"),
              py::arg("propagator") = tp::quaternions,
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("rotational", 1).c_str());


        m.def("mass",
              py::overload_cast<
                      const std::vector<std::string>,
                      const std::map<std::string, std::shared_ptr<tba::MassRateModel> > &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::shared_ptr<tp::DependentVariableSaveSettings>,
                      const double>(&tp::massPropagatorSettings<double>),
              py::arg("bodies_with_mass_to_propagate"),
              py::arg("mass_rate_models"),
              py::arg("initial_body_masses"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("mass", 0).c_str());

        m.def("mass",
              py::overload_cast<
                      const std::vector<std::string>,
                      const std::map<std::string, std::vector<std::shared_ptr<tba::MassRateModel> > > &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::shared_ptr<tp::DependentVariableSaveSettings>,
                      const double>(&tp::massPropagatorSettings<double>),
              py::arg("bodies_with_mass_to_propagate"),
              py::arg("mass_rate_models"),
              py::arg("initial_body_masses"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("mass", 1).c_str());


        m.def("mass",
              py::overload_cast<
                      const std::vector<std::string>,
                      const tss::SelectedMassRateModelMap &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::shared_ptr<tp::DependentVariableSaveSettings>,
                      const double>(&tp::massPropagatorSettings<double>),
              py::arg("bodies_with_mass_to_propagate"),
              py::arg("mass_rate_settings"),
              py::arg("initial_body_masses"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("mass", 2).c_str());


        m.def("mass",
              py::overload_cast<
                      const std::vector<std::string>,
                      const std::map<std::string, std::vector<std::shared_ptr<tba::MassRateModel> > > &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
                      const double>(&tp::massPropagatorSettings<double>),
              py::arg("bodies_with_mass_to_propagate"),
              py::arg("mass_rate_models"),
              py::arg("initial_body_masses"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("mass", 3).c_str());


        m.def("mass",
              py::overload_cast<
                      const std::vector<std::string>,
                      const tss::SelectedMassRateModelMap &,
                      const Eigen::Matrix<double, Eigen::Dynamic, 1> &,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
                      const double>(&tp::massPropagatorSettings<double>),
              py::arg("bodies_with_mass_to_propagate"),
              py::arg("mass_rate_settings"),
              py::arg("initial_body_masses"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("mass", 4).c_str());


        m.def("multitype",
              py::overload_cast<
                      const std::vector<std::shared_ptr<tp::SingleArcPropagatorSettings<double> > >,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::shared_ptr<tp::DependentVariableSaveSettings>,
                      const double>(&tp::multiTypePropagatorSettings<double>),
              py::arg("propagator_settings_list"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::shared_ptr<tp::DependentVariableSaveSettings>(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("multitype", 0).c_str());


        m.def("multitype",
              py::overload_cast<
                      const std::vector<std::shared_ptr<tp::SingleArcPropagatorSettings<double> > >,
                      const std::shared_ptr<tp::PropagationTerminationSettings>,
                      const std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >,
                      const double>(&tp::multiTypePropagatorSettings<double>),
              py::arg("propagator_settings_list"),
              py::arg("termination_settings"),
              py::arg("output_variables") = std::vector<std::shared_ptr<tp::SingleDependentVariableSaveSettings> >(),
              py::arg("print_interval") = TUDAT_NAN,
              get_docstring("multitype", 1).c_str());

        m.def("multi_arc",
              &tp::multiArcPropagatorSettings<double>,
              py::arg("single_arc_settings"),
              py::arg("transfer_state_to_next_arc") = false,
              get_docstring("multi_arc").c_str());

        m.def("hybrid_arc",
              &tp::hybridArcPropagatorSettings<double>,
              py::arg("single_arc_settings"),
              py::arg("multi_arc_settings"),
              get_docstring("hybrid_arc").c_str());

        py::class_<tp::PropagationTerminationSettings,
                std::shared_ptr<tp::PropagationTerminationSettings>>
                PropagationTerminationSettings_(m, "PropagationTerminationSettings",
                                                get_docstring("PropagationTerminationSettings").c_str());

        py::class_<
                tp::PropagationDependentVariableTerminationSettings,
                std::shared_ptr<tp::PropagationDependentVariableTerminationSettings>,
                tp::PropagationTerminationSettings>(m, "PropagationDependentVariableTerminationSettings",
                                                    get_docstring(
                                                            "PropagationDependentVariableTerminationSettings").c_str());

//                .def(py::init<
//                             const std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
//                             const double,
//                             const bool,
//                             const bool,
//                             const std::shared_ptr<tudat::root_finders::RootFinderSettings>>(),
//                     py::arg("dependent_variadble_settings"),
//                     py::arg("limit_value"),
//                     py::arg("use_as_lower_limit"),
//                     py::arg("terminate_exactly_on_final_condition") = false,
//                     py::arg("termination_root_finder_settings") = nullptr);

        m.def("time_termination",
              &tp::propagationTimeTerminationSettings,
              py::arg("termination_time"),
              py::arg("terminate_exactly_on_final_condition") = false,
              get_docstring("time_termination").c_str());

        m.def("cpu_time_termination",
              &tp::propagationCPUTimeTerminationSettings,
              py::arg("cpu_termination_time"),
              get_docstring("cpu_time_termination").c_str());

        m.def("dependent_variable_termination",
              &tp::propagationDependentVariableTerminationSettings,
              py::arg("dependent_variable_settings"),
              py::arg("limit_value"),
              py::arg("use_as_lower_limit"),
              py::arg("terminate_exactly_on_final_condition") = false,
              py::arg("termination_root_finder_settings") = nullptr,
              get_docstring("dependent_variable_termination").c_str());

        m.def("custom_termination",
              &tp::popagationCustomTerminationSettings,
              py::arg("custom_condition"),
              get_docstring("custom_termination").c_str());

        m.def("hybrid_termination",
              &tp::propagationHybridTerminationSettings,
              py::arg("termination_settings"),
              py::arg("fulfill_single_condition"),
              get_docstring("hybrid_termination").c_str());

    }

}// namespace propagator
}// namespace propagation_setup
}// namespace simulation
}// namespace tudatpy
