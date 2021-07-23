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

#include "expose_propagation_setup/expose_dependent_variable_setup.h"
#include "expose_propagation_setup/expose_torque_setup.h"
#include "expose_propagation_setup/expose_acceleration_setup.h"
#include "expose_propagation_setup/expose_integrator_setup.h"
#include "expose_propagation_setup/expose_mass_rate_setup.h"
#include "expose_propagation_setup/expose_propagator_setup.h"


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


void expose_propagation_setup(py::module &m) {

    py::enum_<tss::ThrustFrames>(m, "ThrustFrames")
            .value("unspecified_thrust_frame", tss::ThrustFrames::unspecified_thrust_frame)
            .value("inertial_thtust_frame", tss::ThrustFrames::inertial_thtust_frame)
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
    // createTorqueModels.cpp
    //////////////////////////////////////////////////////////////////////////////

    m.def("create_torque_models",// overload [1/2]
              &tss::createTorqueModelsMap,
          py::arg("body_system"),
          py::arg("selected_acceleration_per_body"),
          py::arg("bodies_to_propagate"));


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
                 const std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>&,
                 const bool >(),
                 py::arg("body_map"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("are_equations_of_motion_to_be_integrated") = true,
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_number_of_function_evaluations") = false,
                 py::arg("initial_clock_time") = std::chrono::steady_clock::now(),
                 py::arg("state_derivative_models") =
            std::vector<std::shared_ptr<tp::SingleStateTypeDerivative<double, double>>>(),
                 py::arg("print_dependent_variable_data" )= true )
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
                 &tp::SingleArcDynamicsSimulator<double, double>::processNumericalEquationsOfMotionSolution)
            .def("suppress_dependent_variable_terminal_printing",
                 &tp::SingleArcDynamicsSimulator<double, double>::suppressDependentVariableDataPrinting)
            .def("enable_dependent_variable_terminal_printing",
                 &tp::SingleArcDynamicsSimulator<double, double>::enableDependentVariableDataPrinting);


    //        py::enum_<tp::VariableType>(m, "VariableType")
    //                .value("independent_variable", tp::VariableType::independentVariable)
    //                .value("cpu_time_variable", tp::VariableType::cpuTimeVariable)
    //                .value("state_variable", tp::VariableType::stateVariable)
    //                .value("dependent_variable", tp::VariableType::dependentVariable)
    //                .export_values();



    auto acceleration_setup = m.def_submodule("acceleration");
    expose_acceleration_setup(acceleration_setup);

    auto torque_setup = m.def_submodule("torque");
    expose_torque_setup(torque_setup);

    auto integrator_setup = m.def_submodule("integrator");
    expose_integrator_setup(integrator_setup);

    auto propagator_setup = m.def_submodule("propagator");
    expose_propagator_setup(propagator_setup);

    auto mass_setup = m.def_submodule("mass");
    expose_mass_rate_setup(mass_setup);

    auto dependent_variable_setup = m.def_submodule("dependent_variable");
    expose_dependent_variable_setup(dependent_variable_setup);
}


}// namespace propagation_setup
}// namespace simulation
}// namespace tudatpy
