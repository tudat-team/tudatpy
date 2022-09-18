/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/docstrings.h"

#include "expose_numerical_simulation.h"

#include "expose_numerical_simulation/expose_environment_setup.h"
#include "expose_numerical_simulation/expose_estimation_setup.h"
#include "expose_numerical_simulation/expose_propagation_setup.h"

#include "expose_numerical_simulation/expose_environment.h"
#include "expose_numerical_simulation/expose_estimation.h"
#include "expose_numerical_simulation/expose_propagation.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;

namespace tudatpy {
namespace numerical_simulation {

void expose_numerical_simulation(py::module &m) {


    auto environment_submodule = m.def_submodule("environment");
    environment::expose_environment(environment_submodule);

    auto propagation_submodule = m.def_submodule("propagation");
    propagation::expose_propagation(propagation_submodule);

    auto estimation_submodule = m.def_submodule("estimation");
    estimation::expose_estimation(estimation_submodule);

    auto environment_setup_submodule = m.def_submodule("environment_setup");
    environment_setup::expose_environment_setup(environment_setup_submodule);

    auto propagation_setup_submodule = m.def_submodule("propagation_setup");
    propagation_setup::expose_propagation_setup(propagation_setup_submodule);

    auto estimation_setup_submodule = m.def_submodule("estimation_setup");
    estimation_setup::expose_estimation_setup(estimation_setup_submodule);

    m.def("get_integrated_type_and_body_list",
          &tp::getIntegratedTypeAndBodyList<double>,
          py::arg("propagator_settings") );

    m.def("get_single_integration_size",
          &tp::getSingleIntegrationSize,
          py::arg("state_type") );

    py::class_<
            tp::SingleArcDynamicsSimulator<double, double>,
            std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double>>>(m,
                                                                             "SingleArcSimulator",
                                                                             get_docstring("SingleArcSimulator").c_str())
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies &,
                 const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<double>>,
                 const std::shared_ptr<tp::PropagatorSettings<double>>,
                 const bool,
                 const bool,
                 const bool,
                 const bool,
                 const bool,
                 const bool>(),
                 py::arg("bodies"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("are_equations_of_motion_to_be_integrated") = true,
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_number_of_function_evaluations") = false,
                 py::arg("print_dependent_variable_data") = true,
                 py::arg("print_state_data") = true,
                 get_docstring("SingleArcSimulator.ctor").c_str())
            .def("integrate_equations_of_motion",
                 &tp::SingleArcDynamicsSimulator<double, double>::integrateEquationsOfMotion,
                 py::arg("initial_states"),
                 get_docstring("SingleArcSimulator.integrate_equations_of_motion").c_str())
            .def_property_readonly("state_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution,
                                   get_docstring("SingleArcSimulator.state_history").c_str())
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw,
                                   get_docstring("SingleArcSimulator.unprocessed_state_history").c_str())
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory,
                                   get_docstring("SingleArcSimulator.dependent_variable_history").c_str())
            .def_property_readonly("cumulative_computation_time_history",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistory,
                                   get_docstring("SingleArcSimulator.cumulative_computation_time_history").c_str())
            .def_property_readonly("cumulative_number_of_function_evaluations",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeNumberOfFunctionEvaluations,
                                   get_docstring("SingleArcSimulator.cumulative_number_of_function_evaluations").c_str())
            //          .def("manually_set_and_process_raw_numerical_equations_of_motion_solution",
            //               &tp::SingleArcDynamicsSimulator<double, double>::manuallySetAndProcessRawNumericalEquationsOfMotionSolution,
            //               py::arg("equations_of_motion_numerical_solution"),
            //               py::arg("dependent_variable_history"),
            //               py::arg("process_solution"),
            //               get_docstring("manually_set_and_process_raw_numerical_equations_of_motion_solution").c_str())
            .def_property_readonly("integrator_settings",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getIntegratorSettings,
                                   get_docstring("SingleArcSimulator.integrator_settings").c_str())
            .def_property_readonly("state_derivative_function",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getStateDerivativeFunction,
                                   get_docstring("SingleArcSimulator.state_derivative_function").c_str())
            .def_property_readonly("environment_updater",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getEnvironmentUpdater,
                                   get_docstring("SingleArcSimulator.environment_updater").c_str())
            //          .def_property_readonly("dynamics_state_derivative",
            //                                 &tp::SingleArcDynamicsSimulator<double, double>::getDynamicsStateDerivative,
            //                                 get_docstring("dynamics_state_derivative").c_str())
            //          .def_property_readonly("propagation_termination_condition",
            //                                 &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationCondition,
            //                                 get_docstring("propagation_termination_condition").c_str())
            //          .def_property_readonly("integrated_state_processors",
            //                                 &tp::SingleArcDynamicsSimulator<double, double>::getIntegratedStateProcessors,
            //                                 get_docstring("integrated_state_processors").c_str())
            .def_property_readonly("propagation_termination_details",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getPropagationTerminationReason,
                                   get_docstring("SingleArcSimulator.propagation_termination_reason").c_str())
            .def_property_readonly("integration_completed_successfully",
                                   &tp::SingleArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully,
                                   get_docstring("SingleArcSimulator.integration_completed_successfully").c_str())
            .def_property_readonly("dependent_variable_ids",
                                   &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableIds,
                                   get_docstring("SingleArcSimulator.dependent_variable_ids").c_str());
    //          .def_property_readonly("initial_propagation_time",
    //                                 &tp::SingleArcDynamicsSimulator<double, double>::getInitialPropagationTime,
    //                                 get_docstring("initial_propagation_time").c_str());
    //          .def("reset_initial_propagation_time",
    //               &tp::SingleArcDynamicsSimulator<double, double>::resetInitialPropagationTime,
    //               py::arg("new_initial_propagation_time"),
    //               get_docstring("reset_initial_propagation_time").c_str())
    //          .def_property_readonly("dependent_variables_functions",
    //                                 &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariablesFunctions,
    //                                 get_docstring("dependent_variables_functions").c_str())
    //          .def("reset_propagation_termination_conditions",
    //               &tp::SingleArcDynamicsSimulator<double, double>::resetPropagationTerminationConditions,
    //               get_docstring("reset_propagation_termination_conditions").c_str())
    //          .def("process_numerical_equations_of_motion_solution",
    //               &tp::SingleArcDynamicsSimulator<double, double>::processNumericalEquationsOfMotionSolution,
    //               get_docstring("process_numerical_equations_of_motion_solution").c_str())
    //          .def("suppress_dependent_variable_terminal_printing",
    //               &tp::SingleArcDynamicsSimulator<double, double>::suppressDependentVariableDataPrinting,
    //               get_docstring("suppress_dependent_variable_terminal_printing").c_str())
    //          .def("enable_dependent_variable_terminal_printing",
    //               &tp::SingleArcDynamicsSimulator<double, double>::enableDependentVariableDataPrinting,
    //               get_docstring("enable_dependent_variable_terminal_printing").c_str());



    //TODO: Remove variationalOnlyIntegratorSettings
    py::class_<
            tp::SingleArcVariationalEquationsSolver<double, double>,
            std::shared_ptr<tp::SingleArcVariationalEquationsSolver<double, double>>>(m, "SingleArcVariationalSimulator",
                                                                                      get_docstring("SingleArcVariationalSimulator").c_str() )
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies&,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings<double>>,
                 const std::shared_ptr< tp::PropagatorSettings<double>>,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const bool,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >,
                 const bool,
                 const bool,
                 const bool >(),
                 py::arg("bodies"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("estimated_parameters"),
                 py::arg("integrate_equations_concurrently") = true,
                 py::arg("variational_only_integrator_settings") = std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >( ),
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("integrate_on_creation") = true,
                 py::arg("set_integrated_result") = false,
                 get_docstring("SingleArcVariationalSimulator.ctor").c_str() )
            .def("integrate_equations_of_motion_only",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::integrateDynamicalEquationsOfMotionOnly,
                 py::arg("initial_states"),
                 get_docstring("SingleArcVariationalSimulator.integrate_equations_of_motion_only").c_str() )
            .def("integrate_full_equations",
                 &tp::SingleArcVariationalEquationsSolver<double, double>::integrateVariationalAndDynamicalEquations,
                 py::arg("initial_states"),
                 py::arg("integrate_equations_concurrently") = true,
                 get_docstring("SingleArcVariationalSimulator.integrate_full_equations").c_str() )
            .def_property("parameter_vector",
                          &tp::SingleArcVariationalEquationsSolver<double, double>::getParametersToEstimate,
                          &tp::SingleArcVariationalEquationsSolver<double, double>::resetParameterEstimate,
                          get_docstring("SingleArcVariationalSimulator.parameter_vector").c_str() )
            .def_property_readonly("variational_equations_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getNumericalVariationalEquationsSolution,
                                   get_docstring("SingleArcVariationalSimulator.variational_equations_history").c_str() )
            .def_property_readonly("state_transition_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getStateTransitionMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_transition_matrix_history").c_str() )
            .def_property_readonly("sensitivity_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getSensitivityMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.sensitivity_matrix_history").c_str() )
            .def_property_readonly("state_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getEquationsOfMotionSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_history").c_str() )
            .def_property_readonly("dynamics_simulator",
                                   &tp::SingleArcVariationalEquationsSolver<double, double>::getDynamicsSimulator,
                                   get_docstring("SingleArcVariationalSimulator.dynamics_simulator").c_str() );

    py::class_<
            tss::OrbitDeterminationManager<double, double>,
            std::shared_ptr<tss::OrbitDeterminationManager<double, double>>>(m, "Estimator",
                                                                             get_docstring("Estimator").c_str() )
            .def(py::init<const tss::SystemOfBodies&,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                 const std::shared_ptr< tni::IntegratorSettings< double > >,
                 const std::shared_ptr< tp::PropagatorSettings< double > >,
                 const bool >( ),
                 py::arg("bodies"),
                 py::arg("estimated_parameters"),
                 py::arg("observation_settings"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("integrate_on_creation") = true,
                 get_docstring("Estimator.ctor").c_str() )
            .def_property_readonly("observation_simulators",
                                   &tss::OrbitDeterminationManager<double, double>::getObservationSimulators,
                                   get_docstring("Estimator.observation_simulators").c_str() )
            .def_property_readonly("observation_managers",
                                   &tss::OrbitDeterminationManager<double, double>::getObservationManagers,
                                   get_docstring("Estimator.observation_managers").c_str() )
            .def_property_readonly("state_transition_interface",
                                   &tss::OrbitDeterminationManager<double, double>::getStateTransitionAndSensitivityMatrixInterface,
                                   get_docstring("Estimator.state_transition_interface").c_str() )
            .def("perform_estimation",
                 &tss::OrbitDeterminationManager<double, double>::estimateParameters,
                 py::arg( "estimation_input" ),
                 py::arg( "convergence_checker" ) = std::make_shared< tss::EstimationConvergenceChecker >( ),
                 get_docstring("Estimator.perform_estimation").c_str() )
            .def_property_readonly("variational_solver",
                                   &tss::OrbitDeterminationManager<double, double>::getVariationalEquationsSolver,
                                   get_docstring("Estimator.variational_solver").c_str() );
};

}// namespace numerical_simulation
}// namespace tudatpy
