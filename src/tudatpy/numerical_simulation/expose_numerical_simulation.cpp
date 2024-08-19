/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define PYBIND11_DETAILED_ERROR_MESSAGES

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/basic_astro/dateTime.h>
#include <tudat/basics/timeType.h>
#include <tudat/simulation/environment_setup.h>
#include <tudat/simulation/estimation_setup.h>
#include <tudat/simulation/propagation_setup.h>

// #include <tudat/math/integrators/createNumericalIntegrator.h>

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {
    namespace numerical_simulation {

        PYBIND11_MODULE(expose_numerical_simulation, m) {
            // Expose IntegratorSettings locally to avoid circular imports
            py::class_<tni::IntegratorSettings<TIME_TYPE>,
                       std::shared_ptr<tni::IntegratorSettings<TIME_TYPE>>>(
                m, "IntegratorSettings", py::module_local(),
                R"doc(Functional base class to define settings for integrators.

                Class to define settings for numerical integrators, for instance for use in numerical integration of equations of motion/variational equations. This class can be used for simple integrators such as fixed step RK and Euler. Integrators that require more settings to define have their own derived class.)doc")
                .def_readwrite("initial_time",
                               &tni::IntegratorSettings<
                                   TIME_TYPE>::initialTimeDeprecated_);

            // py::module_::import(
            //     "tudatpy.numerical_simulation.propagation_setup.integrator."
            //     "expose_integrator");
            m.def("get_integrated_type_and_body_list",
                  &tp::getIntegratedTypeAndBodyList<double, TIME_TYPE>,
                  py::arg("propagator_settings"));

            m.def("get_single_integration_size", &tp::getSingleIntegrationSize,
                  py::arg("state_type"));

            m.def("create_dynamics_simulator",
                  &tss::createDynamicsSimulator<double, TIME_TYPE>,
                  py::arg("bodies"), py::arg("propagator_settings"),
                  py::arg("simulate_dynamics_on_creation") = true,
                  R"doc(Function to create object that propagates the dynamics.

	Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
	Depending on the specific input type (e.g. which factory function from the :ref:`\`\`propagator\`\`` module was used),
	a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
	function. When using default settings, calling this function will automatically propagate the dynamics.


	:param bodies:
		Object defining the physical environment, with all
		properties of artificial and natural bodies.

	:param propagator_settings:
		Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

	:param simulate_dynamics_on_creation:
		Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
		can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function

	:return:
		Object that propagates the dynamics, and processes the results.
)doc");

            py::class_<tudat::Time>(m, "Time", get_docstring("Time").c_str())
                .def(py::init<const int, const long double>(),
                     py::arg("full_periods"),
                     py::arg("seconds_into_full_period"))
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
                  &tss::createVariationalEquationsSolver<double, TIME_TYPE>,
                  py::arg("bodies"), py::arg("propagator_settings"),
                  py::arg("parameters_to_estimate"),
                  py::arg("simulate_dynamics_on_creation") = true,
                  R"doc(Function to create object that propagates the dynamics.

	Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
	Depending on the specific input type (e.g. which factory function from the :ref:`\`\`propagator\`\`` module was used),
	a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
	function. When using default settings, calling this function will automatically propagate the dynamics.


	:param bodies:
		Object defining the physical environment, with all
		properties of artificial and natural bodies.

	:param propagator_settings:
		Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

	:param simulate_dynamics_on_creation:
		Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
		can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function

	:return:
		Object that propagates the dynamics, and processes the results.
)doc");


            py::class_<
                tp::DynamicsSimulator<double, TIME_TYPE>,
                std::shared_ptr<tp::DynamicsSimulator<double, TIME_TYPE>>>(
                m, "DynamicsSimulator",
                get_docstring("DynamicsSimulator").c_str());

            py::class_<tp::SingleArcDynamicsSimulator<double, TIME_TYPE>,
                       std::shared_ptr<
                           tp::SingleArcDynamicsSimulator<double, TIME_TYPE>>,
                       tp::DynamicsSimulator<double, TIME_TYPE>>(
                m, "SingleArcSimulator",
                R"doc(Class for consolidating single arc dynamics simulation functionality.

	Class for consolidating all functionality required to perform single arc dynamics simulations.

)doc")
                .def(py::init<
                         const tudat::simulation_setup::SystemOfBodies &,
                         const std::shared_ptr<
                             tudat::numerical_integrators::IntegratorSettings<
                                 TIME_TYPE>>,
                         const std::shared_ptr<tp::PropagatorSettings<double>>,
                         const bool, const bool, const bool, const bool,
                         const bool, const bool>(),
                     py::arg("bodies"), py::arg("integrator_settings"),
                     py::arg("propagator_settings"),
                     py::arg("are_equations_of_motion_to_be_integrated") = true,
                     py::arg("clear_numerical_solutions") = false,
                     py::arg("set_integrated_result") = false,
                     py::arg("print_number_of_function_evaluations") = false,
                     py::arg("print_dependent_variable_data") = true,
                     py::arg("print_state_data") = true)
                .def_property_readonly(
                    "state_history",
                    &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::
                        getEquationsOfMotionNumericalSolution)
                .def_property_readonly(
                    "unprocessed_state_history",
                    &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::
                        getEquationsOfMotionNumericalSolutionRaw)
                .def_property_readonly(
                    "dependent_variable_history",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getDependentVariableHistory)
                .def_property_readonly(
                    "cumulative_computation_time_history",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getCumulativeComputationTimeHistory)
                .def_property_readonly(
                    "cumulative_number_of_function_evaluations",
                    &tp::SingleArcDynamicsSimulator<double, TIME_TYPE>::
                        getCumulativeNumberOfFunctionEvaluations)
                .def_property_readonly(
                    "integrator_settings",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getIntegratorSettings,
                    R"doc(Settings to create the numerical integrator that is to be used
for the integration of the equations of motion

	)doc")
                .def_property_readonly(
                    "state_derivative_function",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getStateDerivativeFunction,
                    R"doc(Function that performs a single state derivative function evaluation. This function takes the numerically propagated
state, and current independent variable (time) as input, and returns the derivative of the state that is then used
by the numerical integration routine. Typically, this function is NOT used directly by users.

	)doc")
                .def_property_readonly(
                    "environment_updater",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getEnvironmentUpdater,
                    R"doc(# Object used in the propagation to update the environment, it uses the current time and numerically calculated state
to update the translational state, rotational state, flight conditions, etc. of all bodies in the simulation to be
consistent with this time and state.  Typically, this class is NOT used directly by users.

	)doc")
                .def_property_readonly(
                    "propagation_results",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getPropagationResults,
                    R"doc(This function retrieves all the results of the numerical propagation, stored
in a single wrapper object

	)doc")
                .def_property_readonly(
                    "propagation_termination_details",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::getPropagationTerminationReason)
                .def_property_readonly(
                    "integration_completed_successfully",
                    &tp::SingleArcDynamicsSimulator<
                        double, TIME_TYPE>::integrationCompletedSuccessfully);

            py::class_<tp::MultiArcDynamicsSimulator<double, TIME_TYPE>,
                       std::shared_ptr<
                           tp::MultiArcDynamicsSimulator<double, TIME_TYPE>>,
                       tp::DynamicsSimulator<double, TIME_TYPE>>(
                m, "MultiArcDynamicsSimulator",
                get_docstring("MultiArcDynamicsSimulator").c_str())
                .def_property_readonly(
                    "propagation_results",
                    &tp::MultiArcDynamicsSimulator<
                        double, TIME_TYPE>::getMultiArcPropagationResults,
                    get_docstring(
                        "MultiArcDynamicsSimulator.propagation_results")
                        .c_str());

            py::class_<tp::HybridArcDynamicsSimulator<double, TIME_TYPE>,
                       std::shared_ptr<
                           tp::HybridArcDynamicsSimulator<double, TIME_TYPE>>,
                       tp::DynamicsSimulator<double, TIME_TYPE>>(
                m, "HybridArcDynamicsSimulator",
                get_docstring("HybridArcDynamicsSimulator").c_str())
                .def_property_readonly(
                    "propagation_results",
                    &tp::HybridArcDynamicsSimulator<
                        double, TIME_TYPE>::getHybridArcPropagationResults,
                    get_docstring(
                        "HybridArcDynamicsSimulator.propagation_results")
                        .c_str());


            // TODO: Remove variationalOnlyIntegratorSettings
            py::class_<
                tp::SingleArcVariationalEquationsSolver<double, TIME_TYPE>,
                std::shared_ptr<tp::SingleArcVariationalEquationsSolver<
                    double, TIME_TYPE>>>(
                m, "SingleArcVariationalSimulator",
                R"doc(Class for consolidating single arc variational dynamics functionality.

	Class for consolidating all functionality required to perform single arc variational dynamics simulations.

)doc")
                .def(py::init<
                         const tudat::simulation_setup::SystemOfBodies &,
                         const std::shared_ptr<
                             tudat::numerical_integrators::IntegratorSettings<
                                 TIME_TYPE>>,
                         const std::shared_ptr<tp::PropagatorSettings<double>>,
                         const std::shared_ptr<
                             tep::EstimatableParameterSet<double>>,
                         const bool,
                         const std::shared_ptr<tudat::numerical_integrators::
                                                   IntegratorSettings<double>>,
                         const bool, const bool, const bool, const bool>(),
                     py::arg("bodies"), py::arg("integrator_settings"),
                     py::arg("propagator_settings"),
                     py::arg("estimated_parameters"),
                     py::arg("integrate_equations_concurrently") = true,
                     py::arg("variational_only_integrator_settings") =
                         std::shared_ptr<tni::IntegratorSettings<double>>(),
                     py::arg("clear_numerical_solutions") = false,
                     py::arg("integrate_on_creation") = true,
                     py::arg("set_integrated_result") = false,
                     py::arg("print_dependent_variable_data") = true,
                     R"doc(Class constructor.

	Constructor through which the user can create instances of this class.
	Defines environment, propagation and integrations models, as well as a number of settings related
	to the (estimatable) parameters, w.r.t. which the variational equations are defined.

	.. note:: When using default settings, creating an object of
	          this type automatically triggers the propagation


	:param bodies:
		Object defining the physical environment, with all
		properties of artificial and natural bodies.

	:param integrator_settings:
		Settings to create the numerical integrator that is to be used for the integration of the equations of motion.

	:param propagator_settings:
		Settings to create the propagator that is to be used for the propagation of the dynamics.

	:param estimated_parameters:
		Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
		linked to the environment and acceleration settings of the simulation.

	:param integrate_equations_concurrently:
		Boolean defining whether equations of motion and variational equations are to be propagated concurrently
		(if true) or sequentially (of false).

	:param variational_only_integrator_settings:
		Settings to create the numerical integrator that is to be used for integration the variational equations.
		If none is given (default), the numerical integration settings are taken to be the same as the ones applied
		in the integration of the equations of motions (specified by the `integrator_settings` parameter).

	:param clear_numerical_solutions:
		Boolean to determine whether to clear the raw numerical solution member variables
		and to reset the state transition interface after propagation.

	:param integrate_on_creation:
		Boolean defining whether the propagation should be performed immediately (default), or at a later time
		(when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).

	:param set_integrated_result:
		Boolean to determine whether to automatically use the integrated results to set ephemerides for the
		propagated bodies.

)doc")
                .def(
                    "integrate_equations_of_motion_only",
                    &tp::SingleArcVariationalEquationsSolver<
                        double,
                        TIME_TYPE>::integrateDynamicalEquationsOfMotionOnly,
                    py::arg("initial_states"),
                    R"doc(Function to trigger the integration of the (regular) equations of motion.


	Function to trigger the integration only of the (regular) equations of motion, resulting in a `state_history`.
	This step does not yet use variational dynamics. In order to also solve the variational equations,
	use the :func:`integrate_full_equations` member function.

	:return:
		Creates / modifies the `state_history` property of the :class:`~tudatpy.numerical_simulation.SingleArcVariationalSolver` object.
)doc")
                .def(
                    "integrate_full_equations",
                    &tp::SingleArcVariationalEquationsSolver<
                        double,
                        TIME_TYPE>::integrateVariationalAndDynamicalEquations,
                    py::arg("initial_states"),
                    py::arg("integrate_equations_concurrently") = true,
                    R"doc(Function to trigger the integration of variational and dynamical equations (equations of motion).


	Function to trigger the integration of the (regular) equations of motion as well as the variational equations,
	solving for `state_history` and `variational_equations_history`
	(in its two components `state_transition_matrix_history` & `sensitivity_matrix_history`).


	:param initial_states:
		Initial state to be used for the parameters in the equations of motion.

	:param integrate_equations_concurrently:
		Boolean defining whether equations of motion and variational equations are to be propagated concurrently
		(if true) or sequentially (of false).

	:return:
		Creates / modifies the properties of the VariationalSolver object.
)doc")
                .def_property("parameter_vector",
                              &tp::SingleArcVariationalEquationsSolver<
                                  double, TIME_TYPE>::getParametersToEstimate,
                              &tp::SingleArcVariationalEquationsSolver<
                                  double, TIME_TYPE>::resetParameterEstimate,
                              R"doc(Consolidated set of (estimatable) parameters
w.r.t. the variational dynamics in the Variational Simulator are defined.

	)doc")
                .def_property_readonly(
                    "variational_equations_history",
                    &tp::SingleArcVariationalEquationsSolver<
                        double,
                        TIME_TYPE>::getNumericalVariationalEquationsSolution,
                    R"doc(List containing the solution of the variational equations, i.e. the
state transition matrix history (first entry) and sensitivity matrix history (second vector entry).

	)doc")
                .def_property_readonly(
                    "state_transition_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver<
                        double, TIME_TYPE>::getStateTransitionMatrixSolution,
                    R"doc(State transition matrix history, given as epoch with propagation epochs as keys.
This is (alongside the `sensitivity_matrix_history`) the solution of the variational equations.

	)doc")
                .def_property_readonly(
                    "sensitivity_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver<
                        double, TIME_TYPE>::getSensitivityMatrixSolution,
                    R"doc(Sensitivity matrix history, given as epoch with propagation epochs as keys.
This is (alongside the `state_transition_matrix_history`) the solution of the variational equations.

	)doc")
                .def_property_readonly(
                    "state_history",
                    &tp::SingleArcVariationalEquationsSolver<
                        double, TIME_TYPE>::getEquationsOfMotionSolution,
                    R"doc(State history, given as epoch with propagation epochs as keys.
This is the solution of the (propagated) equations of motion, describing the states along which
the variational dynamics are solved.

	)doc")
                .def_property_readonly(
                    "dynamics_simulator",
                    &tp::SingleArcVariationalEquationsSolver<
                        double, TIME_TYPE>::getDynamicsSimulator,
                    R"doc(Simulator object containing all functionality for solving of the (regular) equations of motion.

	)doc");

            py::class_<tss::OrbitDeterminationManager<double, TIME_TYPE>,
                       std::shared_ptr<
                           tss::OrbitDeterminationManager<double, TIME_TYPE>>>(
                m, "Estimator",
                R"doc(Class for consolidating all estimation functionality.

	Class for consolidating all functionality required to perform an estimation.

)doc")
                .def(py::init<
                         const tss::SystemOfBodies &,
                         const std::shared_ptr<
                             tep::EstimatableParameterSet<double>>,
                         const std::vector<
                             std::shared_ptr<tom::ObservationModelSettings>> &,
                         const std::shared_ptr<tp::PropagatorSettings<double>>,
                         const bool>(),
                     py::arg("bodies"), py::arg("estimated_parameters"),
                     py::arg("observation_settings"),
                     py::arg("propagator_settings"),
                     py::arg("integrate_on_creation") = true,
                     R"doc(Class constructor.

	Constructor through which the user can create instances of this class.
	Defines environment, propagation and integrations models, as well as a number of settings related
	to the estimatable parameters and observation settings.

	.. note:: When using default settings, creating an object of
	          this type automatically triggers the propagation


	:param bodies:
		Object defining the physical environment, with all
		properties of artificial and natural bodies.

	:param estimated_parameters:
		Object defining a consolidated set of estimatable parameters,
		linked to the environment and acceleration settings of the simulation.

	:param observation_settings:
		List of settings objects, each object defining the observation model settings for one
		combination of observable and link geometry that is to be simulated.

	:param integrator_settings:
		Settings to create the numerical integrator that is to be
		used for the integration of the equations of motion

	:param propagator_settings:
		Settings to create the propagator that is to be
		used for the propagation of dynamics

	:param integrate_on_creation:
		Boolean defining whether the propagation should be
		performed immediately (default), or at a later time
		(when calling the :func:`perform_estimation` member function.

)doc")
                .def_property_readonly(
                    "observation_simulators",
                    &tss::OrbitDeterminationManager<
                        double, TIME_TYPE>::getObservationSimulators,
                    R"doc(Observation simulators contained in the Estimator object. A single observation simulator hosts
the functionality for simulating a given observable over the defined link geometry.

	)doc")
                .def_property_readonly(
                    "observation_managers",
                    &tss::OrbitDeterminationManager<
                        double, TIME_TYPE>::getObservationManagers,
                    R"doc(Observation managers contained in the Estimator object. A single observation manager can simulate observations and
calculate observation partials for all link ends involved in the given observable type.

	)doc")
                .def_property_readonly(
                    "state_transition_interface",
                    &tss::OrbitDeterminationManager<double, TIME_TYPE>::
                        getStateTransitionAndSensitivityMatrixInterface,
                    R"doc(State transition and sensitivity matrix interface, setting the variational equations/dynamics in the
Estimator object.

	)doc")
                .def("perform_estimation",
                     &tss::OrbitDeterminationManager<
                         double, TIME_TYPE>::estimateParameters,
                     py::arg("estimation_input"),
                     R"doc(Function to trigger the parameter estimation.


	Function to trigger the parameter estimation. Much of the process and requirements are similar to those described in the
	:func:`~tudatpy.numerical_simulation.Estimator.compute_covariance` function. This function uses an iterative least-squares
	estimate process to fit the data (inside ``estimation_input``) to the model defined by the inputs to the ``Estimator`` constructor.s


	:param estimation_input:
		Object consolidating all relevant settings for the estimation
		This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.

	:return:
		Object containing all outputs from the estimation process.
)doc")
                .def(
                    "compute_covariance",
                    &tss::OrbitDeterminationManager<
                        double, TIME_TYPE>::computeCovariance,
                    py::arg("covariance_analysis_input"),
                    R"doc(Function to perform a covariance analysis for the given observations and parameters


	Function to perform a covariance analysis for the given observations and parameters. The observations are provided through the
	``covariance_analysis_input`` input, as are the weights :math:`\mathbf{W}` and inverse a priori covariance :math:`(\mathbf{P}_{0})^{-1}`.
	Calling this function uses the environment and propagator settings provided to the constructor of this `Estimator` class to simulate
	the dynamics of any relevant bodies for the observations (and associated variational equations). The observations are then
	computed using the observation models created by the settings provided to the constructor of this `Estimator` class, as is the
	associated design matrix :math:`\mathbf{H}`. This function then produces the covariance :math:`\mathbf{P}` (omitting the normalization used
	internally for numerical stability)

	.. math::
	   \mathbf{P}=\left(\mathbf{H}^{T}\mathbf{W}\mathbf{H}+(\mathbf{P}_{0})^{-1}\right)^{-1}

	Note that, although the actual observations are formally not required for a covariance analysis, all additional data (e.g. observation time, type, link ends, etc.)
	are. And, as such, the ``covariance_analysis_input`` does require the full set of observations and associated information, for consistency purposes (e.g., same input as
	``perform_estimation`` function) .


	:param covariance_analysis_input:
		Object consolidating all relevant settings for the covariance analysis
		This includes foremost the simulated observations, as well as a priori information about the estimatable parameters

	:return:
		Object containing all outputs from the estimation process.
)doc")
                .def_property_readonly(
                    "variational_solver",
                    &tss::OrbitDeterminationManager<
                        double, TIME_TYPE>::getVariationalEquationsSolver,
                    R"doc(Variational equations solver, which is used to manage and execute the numerical integration of
equations of motion and variational equations/dynamics in the Estimator object.

	)doc");
        }

    }  // namespace numerical_simulation
}  // namespace tudatpy
