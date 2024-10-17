/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "docstrings.h"
#include "scalarTypes.h"

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
          &tp::getIntegratedTypeAndBodyList<STATE_SCALAR_TYPE,TIME_TYPE>,
          py::arg("propagator_settings") );

    m.def("get_single_integration_size",
          &tp::getSingleIntegrationSize,
          py::arg("state_type") );

    m.def("create_dynamics_simulator",
          &tss::createDynamicsSimulator<STATE_SCALAR_TYPE,TIME_TYPE>,
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
            .def(py::init<
                 const double>(),
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
          &tss::createVariationalEquationsSolver<STATE_SCALAR_TYPE,TIME_TYPE>,
          py::arg("bodies"),
          py::arg("propagator_settings"),
          py::arg("parameters_to_estimate"),
          py::arg("simulate_dynamics_on_creation") = true,
          get_docstring("create_dynamics_simulator").c_str() );



    py::class_<
        tp::DynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
        std::shared_ptr<tp::DynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>>(
        m,"DynamicsSimulator", get_docstring("DynamicsSimulator").c_str());

    py::class_<
            tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>,
            tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE> >(
                m,"SingleArcSimulator", get_docstring("SingleArcSimulator").c_str())
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies &,
                 const std::shared_ptr<tudat::numerical_integrators::IntegratorSettings<TIME_TYPE>>,
                 const std::shared_ptr<tp::PropagatorSettings<STATE_SCALAR_TYPE>>,
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
            .def_property_readonly("bodies",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getSystemOfBodies)
            .def_property_readonly("state_history",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getEquationsOfMotionNumericalSolution)
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getEquationsOfMotionNumericalSolutionRaw)
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getDependentVariableHistory)
            .def_property_readonly("cumulative_computation_time_history",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getCumulativeComputationTimeHistory)
            .def_property_readonly("cumulative_number_of_function_evaluations",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getCumulativeNumberOfFunctionEvaluations)
            .def_property_readonly("integrator_settings",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getIntegratorSettings,
                                   get_docstring("SingleArcSimulator.integrator_settings").c_str())
            .def_property_readonly("state_derivative_function",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getStateDerivativeFunction,
                                   get_docstring("SingleArcSimulator.state_derivative_function").c_str())
            .def_property_readonly("environment_updater",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getEnvironmentUpdater,
                                   get_docstring("SingleArcSimulator.environment_updater").c_str())
            .def_property_readonly("propagation_results",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getPropagationResults,
                                   get_docstring("SingleArcSimulator.propagation_results").c_str())
            .def_property_readonly("propagation_termination_details",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getPropagationTerminationReason)
            .def_property_readonly("integration_completed_successfully",
                                   &tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::integrationCompletedSuccessfully);

    py::class_<
        tp::MultiArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
        std::shared_ptr<tp::MultiArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>,
        tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE>>(
        m,"MultiArcDynamicsSimulator", get_docstring("MultiArcDynamicsSimulator").c_str())
            .def_property_readonly("propagation_results",
                                   &tp::MultiArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getMultiArcPropagationResults,
                                   get_docstring("MultiArcDynamicsSimulator.propagation_results").c_str());

    py::class_<
        tp::HybridArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
        std::shared_ptr<tp::HybridArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>,
        tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE>>(
        m,"HybridArcDynamicsSimulator", get_docstring("HybridArcDynamicsSimulator").c_str())
            .def_property_readonly("propagation_results",
                                   &tp::HybridArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>::getHybridArcPropagationResults,
                                   get_docstring("HybridArcDynamicsSimulator.propagation_results").c_str());


    //TODO: Remove variationalOnlyIntegratorSettings
    py::class_<
            tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "SingleArcVariationalSimulator",
                                                                                      get_docstring("SingleArcVariationalSimulator").c_str() )
            .def(py::init<
                 const tudat::simulation_setup::SystemOfBodies&,
                 const std::shared_ptr< tudat::numerical_integrators::IntegratorSettings<TIME_TYPE>>,
                 const std::shared_ptr< tp::PropagatorSettings<STATE_SCALAR_TYPE>>,
                 const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
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
                 py::arg("variational_only_integrator_settings") = std::shared_ptr< tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >( ),
                 py::arg("clear_numerical_solutions") = false,
                 py::arg("integrate_on_creation") = true,
                 py::arg("set_integrated_result") = false,
                 py::arg("print_dependent_variable_data") = true,
                 get_docstring("SingleArcVariationalSimulator.ctor").c_str() )
            .def("integrate_equations_of_motion_only",
                 &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::integrateDynamicalEquationsOfMotionOnly,
                 py::arg("initial_states"),
                 get_docstring("SingleArcVariationalSimulator.integrate_equations_of_motion_only").c_str() )
            .def("integrate_full_equations",
                 &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::integrateVariationalAndDynamicalEquations,
                 py::arg("initial_states"),
                 py::arg("integrate_equations_concurrently") = true,
                 get_docstring("SingleArcVariationalSimulator.integrate_full_equations").c_str() )
            .def_property("parameter_vector",
                          &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getParametersToEstimate,
                          &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::resetParameterEstimate,
                          get_docstring("SingleArcVariationalSimulator.parameter_vector").c_str() )
            .def_property_readonly("variational_equations_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getNumericalVariationalEquationsSolution,
                                   get_docstring("SingleArcVariationalSimulator.variational_equations_history").c_str() )
            .def_property_readonly("state_transition_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getStateTransitionMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_transition_matrix_history").c_str() )
            .def_property_readonly("sensitivity_matrix_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getSensitivityMatrixSolution,
                                   get_docstring("SingleArcVariationalSimulator.sensitivity_matrix_history").c_str() )
            .def_property_readonly("state_history",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getEquationsOfMotionSolution,
                                   get_docstring("SingleArcVariationalSimulator.state_history").c_str() )
            .def_property_readonly("dynamics_simulator",
                                   &tp::SingleArcVariationalEquationsSolver<STATE_SCALAR_TYPE, TIME_TYPE>::getDynamicsSimulator,
                                   get_docstring("SingleArcVariationalSimulator.dynamics_simulator").c_str() );

    py::class_<
            tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "Estimator",
                                                                             get_docstring("Estimator").c_str() )
            .def(py::init<
                 const tss::SystemOfBodies&,
                 const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                 const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
                 const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                 const bool >( ),
                 py::arg("bodies"),
                 py::arg("estimated_parameters"),
                 py::arg("observation_settings"),
                 py::arg("propagator_settings"),
                 py::arg("integrate_on_creation") = true,
                 get_docstring("Estimator.ctor").c_str() )
            .def_property_readonly("observation_simulators",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationSimulators,
                                   get_docstring("Estimator.observation_simulators").c_str() )
            .def_property_readonly("observation_managers",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationManagers,
                                   get_docstring("Estimator.observation_managers").c_str() )
            .def_property_readonly("state_transition_interface",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getStateTransitionAndSensitivityMatrixInterface,
                                   get_docstring("Estimator.state_transition_interface").c_str() )
            .def("perform_estimation",
                 &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::estimateParameters,
                 py::arg( "estimation_input" ),
                 get_docstring("Estimator.perform_estimation").c_str() )
            .def("compute_covariance",
                 &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::computeCovariance,
                 py::arg( "covariance_analysis_input" ),
                 get_docstring("Estimator.compute_covariance").c_str() )
            .def_property_readonly("variational_solver",
                                   &tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>::getVariationalEquationsSolver,
                                   get_docstring("Estimator.variational_solver").c_str() );
};

}// namespace numerical_simulation
}// namespace tudatpy
