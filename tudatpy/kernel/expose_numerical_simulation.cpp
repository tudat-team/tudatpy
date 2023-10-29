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
#include "tudatpy/scalarTypes.h"

#include "expose_numerical_simulation.h"

#include "expose_numerical_simulation/expose_environment_setup.h"
#include "expose_numerical_simulation/expose_estimation_setup.h"
#include "expose_numerical_simulation/expose_propagation_setup.h"

#include "expose_numerical_simulation/expose_environment.h"
#include "expose_numerical_simulation/expose_estimation.h"
#include "expose_numerical_simulation/expose_propagation.h"

#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/dateTime.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

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
          &tp::getIntegratedTypeAndBodyList<double,TIME_TYPE>,
          py::arg("propagator_settings") );

    m.def("get_single_integration_size",
          &tp::getSingleIntegrationSize,
          py::arg("state_type") );

    m.def("create_dynamics_simulator",
          &tss::createDynamicsSimulator<double,TIME_TYPE>,
          py::arg("bodies"),
          py::arg("propagator_settings"),
          py::arg("simulate_dynamics_on_creation") = true,
          get_docstring("create_dynamics_simulator").c_str() );

    py::class_<
            tudat::Time >(
                m,"Time", get_docstring("Time").c_str())
            .def(py::init<
                 const int,
                 const long double>(),
                 py::arg("full_periods"),
                 py::arg("seconds_into_full_period") )
            .def(py::self + py::self)
            .def(py::self + double())
            .def(double() + py::self)
            .def(py::self += py::self)
            .def(py::self += double())
            .def(py::self - py::self)
            .def(py::self - double())
            .def(py::self -= py::self)
            .def(py::self -= double())
            .def(double() - py::self)
            .def(py::self * double())
            .def(double() * py::self)
            .def(py::self *= double())
            .def(py::self / double())
            .def(py::self /= double())
            .def(py::self == py::self)
            .def(double() == py::self)
            .def(py::self == double())
            .def(py::self != py::self)
            .def(py::self != double())
            .def(double() != py::self)
            .def(py::self < py::self)
            .def(py::self < double())
            .def(double() < py::self)
            .def(py::self > py::self)
            .def(py::self > double())
            .def(double() > py::self)
            .def(py::self <= py::self)
            .def(py::self <= double())
            .def(double() <= py::self)
            .def(py::self >= py::self)
            .def(double() >= py::self)
            .def(py::self >= double());



    m.def("create_variational_equations_solver",
          &tss::createVariationalEquationsSolver<double,TIME_TYPE>,
          py::arg("bodies"),
          py::arg("propagator_settings"),
          py::arg("parameters_to_estimate"),
          py::arg("simulate_dynamics_on_creation") = true,
          get_docstring("create_dynamics_simulator").c_str() );


    py::class_<
            tp::SingleArcDynamicsSimulator<double, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcDynamicsSimulator<double, TIME_TYPE>>>(
                m,"SingleArcSimulator", get_docstring("SingleArcSimulator").c_str())
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies &,
                 const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<TIME_TYPE>>,
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
                 py::arg("print_state_data") = true)

            .def_property_readonly("state_history",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getEquationsOfMotionNumericalSolution)
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getEquationsOfMotionNumericalSolutionRaw)
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getDependentVariableHistory)
            .def_property_readonly("cumulative_computation_time_history",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getCumulativeComputationTimeHistory)
            .def_property_readonly("cumulative_number_of_function_evaluations",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getCumulativeNumberOfFunctionEvaluations)
            //          .def("manually_set_and_process_raw_numerical_equations_of_motion_solution",
            //               &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::manuallySetAndProcessRawNumericalEquationsOfMotionSolution,
            //               py::arg("equations_of_motion_numerical_solution"),
            //               py::arg("dependent_variable_history"),
            //               py::arg("process_solution"),
            //               get_docstring("manually_set_and_process_raw_numerical_equations_of_motion_solution").c_str())
            .def_property_readonly("integrator_settings",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getIntegratorSettings,
                                   get_docstring("SingleArcSimulator.integrator_settings").c_str())
            .def_property_readonly("state_derivative_function",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getStateDerivativeFunction,
                                   get_docstring("SingleArcSimulator.state_derivative_function").c_str())
            .def_property_readonly("environment_updater",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getEnvironmentUpdater,
                                   get_docstring("SingleArcSimulator.environment_updater").c_str())
            .def_property_readonly("propagation_results",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getPropagationResults,
                                   get_docstring("SingleArcSimulator.propagation_results").c_str())
            //          .def_property_readonly("dynamics_state_derivative",
            //                                 &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getDynamicsStateDerivative,
            //                                 get_docstring("dynamics_state_derivative").c_str())
            //          .def_property_readonly("propagation_termination_condition",
            //                                 &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getPropagationTerminationCondition,
            //                                 get_docstring("propagation_termination_condition").c_str())
            //          .def_property_readonly("integrated_state_processors",
            //                                 &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getIntegratedStateProcessors,
            //                                 get_docstring("integrated_state_processors").c_str())
            .def_property_readonly("propagation_termination_details",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getPropagationTerminationReason)
            .def_property_readonly("integration_completed_successfully",
                                   &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::integrationCompletedSuccessfully);
    //          .def_property_readonly("dependent_variable_ids",
    //                                 &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getDependentVariableIds,
    //                                 get_docstring("SingleArcSimulator.dependent_variable_ids").c_str());
    //          .def_property_readonly("initial_propagation_time",
    //                                 &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getInitialPropagationTime,
    //                                 get_docstring("initial_propagation_time").c_str());
    //          .def("reset_initial_propagation_time",
    //               &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::resetInitialPropagationTime,
    //               py::arg("new_initial_propagation_time"),
    //               get_docstring("reset_initial_propagation_time").c_str())
    //          .def_property_readonly("dependent_variables_functions",
    //                                 &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::getDependentVariablesFunctions,
    //                                 get_docstring("dependent_variables_functions").c_str())
    //          .def("reset_propagation_termination_conditions",
    //               &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::resetPropagationTerminationConditions,
    //               get_docstring("reset_propagation_termination_conditions").c_str())
    //          .def("process_numerical_equations_of_motion_solution",
    //               &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::processNumericalEquationsOfMotionSolution,
    //               get_docstring("process_numerical_equations_of_motion_solution").c_str())
    //          .def("suppress_dependent_variable_terminal_printing",
    //               &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::suppressDependentVariableDataPrinting,
    //               get_docstring("suppress_dependent_variable_terminal_printing").c_str())
    //          .def("enable_dependent_variable_terminal_printing",
    //               &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::enableDependentVariableDataPrinting,
    //               get_docstring("enable_dependent_variable_terminal_printing").c_str());


    //TODO: Remove variationalOnlyIntegratorSettings
    py::class_<
            tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>>>(m, "SingleArcVariationalSimulator",
                                                                                      get_docstring("SingleArcVariationalSimulator").c_str() )
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies&,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings<TIME_TYPE>>,
                 const std::shared_ptr< tp::PropagatorSettings<double>>,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const bool,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >,
                 const bool,
                 const bool,
                 const bool,
                 const bool>(),
                 py::arg("bodies"),
                 py::arg("integrator_settings"),
                 py::arg("propagator_settings"),
                 py::arg("estimated_parameters"),
                 py::arg("integrate_equations_concurrently") = true,
                 py::arg("variational_only_integrator_settings") = std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< double > >( ),
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("integrate_on_creation") = true,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_dependent_variable_data") = true,
                 get_docstring("SingleArcVariationalSimulator.ctor").c_str() )
            .def("integrate_equations_of_motion_only",
                 &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::integrateDynamicalEquationsOfMotionOnly,
                 py::arg("initial_states"),
                 get_docstring("SingleArcVariationalSimulator.integrate_equations_of_motion_only").c_str() )
            .def("integrate_full_equations",
                 &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::integrateVariationalAndDynamicalEquations,
                 py::arg("initial_states"),
                 py::arg("integrate_equations_concurrently") = true,
                 get_docstring("SingleArcVariationalSimulator.integrate_full_equations").c_str() )
            .def_property("parameter_vector",
                          &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::getParametersToEstimate,
                          &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::resetParameterEstimate,
                          get_docstring("SingleArcVariationalSimulator.parameter_vector").c_str() )
            .def_property_readonly("variational_equations_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::getNumericalVariationalEquationsSolution,
                                   get_docstring("SingleArcVariationalSimulator.variational_equations_history").c_str() )
            .def_property_readonly("state_transition_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::getStateTransitionMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_transition_matrix_history").c_str() )
            .def_property_readonly("sensitivity_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::getSensitivityMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.sensitivity_matrix_history").c_str() )
            .def_property_readonly("state_history",
                                   &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::getEquationsOfMotionSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_history").c_str() )
            .def_property_readonly("dynamics_simulator",
                                   &tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>::getDynamicsSimulator,
                                   get_docstring("SingleArcVariationalSimulator.dynamics_simulator").c_str() );

    py::class_<
            tss::OrbitDeterminationManager<double, TIME_TYPE>,
            std::shared_ptr<tss::OrbitDeterminationManager<double, TIME_TYPE>>>(m, "Estimator",
                                                                             get_docstring("Estimator").c_str() )
            .def(py::init<
                 const tss::SystemOfBodies&,
                 const std::shared_ptr< tep::EstimatableParameterSet< double > >,
                 const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                 const std::shared_ptr< tp::PropagatorSettings< double > >,
                 const bool >( ),
                 py::arg("bodies"),
                 py::arg("estimated_parameters"),
                 py::arg("observation_settings"),
                 py::arg("propagator_settings"),
                 py::arg("integrate_on_creation") = true,
                 get_docstring("Estimator.ctor").c_str() )
            .def_property_readonly("observation_simulators",
                                   &tss::OrbitDeterminationManager<double, TIME_TYPE>::getObservationSimulators,
                                   get_docstring("Estimator.observation_simulators").c_str() )
            .def_property_readonly("observation_managers",
                                   &tss::OrbitDeterminationManager<double, TIME_TYPE>::getObservationManagers,
                                   get_docstring("Estimator.observation_managers").c_str() )
            .def_property_readonly("state_transition_interface",
                                   &tss::OrbitDeterminationManager<double, TIME_TYPE>::getStateTransitionAndSensitivityMatrixInterface,
                                   get_docstring("Estimator.state_transition_interface").c_str() )
            .def("perform_estimation",
                 &tss::OrbitDeterminationManager<double, TIME_TYPE>::estimateParameters,
                 py::arg( "estimation_input" ),
                 get_docstring("Estimator.perform_estimation").c_str() )
            .def("compute_covariance",
                 &tss::OrbitDeterminationManager<double, TIME_TYPE>::computeCovariance,
                 py::arg( "covariance_analysis_input" ),
                 get_docstring("Estimator.compute_covariance").c_str() )
            .def_property_readonly("variational_solver",
                                   &tss::OrbitDeterminationManager<double, TIME_TYPE>::getVariationalEquationsSolver,
                                   get_docstring("Estimator.variational_solver").c_str() );
};

}// namespace numerical_simulation
}// namespace tudatpy
