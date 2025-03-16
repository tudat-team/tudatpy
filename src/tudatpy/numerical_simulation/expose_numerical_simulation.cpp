/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/operators.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/timeType.h"
#include "tudat/simulation/environment_setup.h"
#include "tudat/simulation/estimation_setup.h"
#include "tudat/simulation/propagation_setup.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tep = tudat::estimatable_parameters;
namespace tom = tudat::observation_models;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy
{

namespace numerical_simulation
{

PYBIND11_MODULE( expose_numerical_simulation, m )
{
    m.def( "get_integrated_type_and_body_list",
           &tp::getIntegratedTypeAndBodyList< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "propagator_settings" ) );

    m.def( "get_single_integration_size", &tp::getSingleIntegrationSize, py::arg( "state_type" ) );

    m.def( "create_dynamics_simulator",
           &tss::createDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ),
           py::arg( "simulate_dynamics_on_creation" ) = true,
           R"doc(

Function to create object that propagates the dynamics.

Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
Depending on the specific input type (e.g. which function from the :ref:`\`\`propagator\`\`` module was used),
a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
function. When using default settings, calling this function will automatically propagate the dynamics.


Parameters
----------
bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object defining the physical environment, with all
    properties of artificial and natural bodies.

propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
    Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

simulate_dynamics_on_creation : Bool, default=True
    Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
    can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function

Returns
-------
:class:`~tudatpy.numerical_simulation.Simulator`
    Object that propagates the dynamics, and processes the results.






    )doc" );

    py::class_< tudat::Time >( m, "Time", R"doc(No documentation found.)doc" )
            .def( py::init< const int, const long double >( ),
                  py::arg( "full_periods" ),
                  py::arg( "seconds_into_full_period" ) )
            .def( py::init< const double >( ), py::arg( "seconds_into_full_period" ) )
            .def( "to_float",
                  &tudat::Time::getSeconds< double >,
                  R"doc(No documentation found.)doc" )
            .def( py::self + py::self )
            .def( py::self + double( ) )
            .def( double( ) + py::self )
            .def( py::self += py::self )
            .def( py::self += double( ) )
            .def( py::self - py::self )
            .def( py::self - double( ) )
            .def( py::self -= py::self )
            .def( py::self -= double( ) )
            .def( double( ) - py::self )
            .def( py::self * double( ) )
            .def( double( ) * py::self )
            .def( py::self *= double( ) )
            .def( py::self / double( ) )
            .def( py::self /= double( ) )
            .def( py::self == py::self )
            .def( double( ) == py::self )
            .def( py::self == double( ) )
            .def( py::self != py::self )
            .def( py::self != double( ) )
            .def( double( ) != py::self )
            .def( py::self < py::self )
            .def( py::self < double( ) )
            .def( double( ) < py::self )
            .def( py::self > py::self )
            .def( py::self > double( ) )
            .def( double( ) > py::self )
            .def( py::self <= py::self )
            .def( py::self <= double( ) )
            .def( double( ) <= py::self )
            .def( py::self >= py::self )
            .def( double( ) >= py::self )
            .def( py::self >= double( ) );

    m.def( "create_variational_equations_solver",
           &tss::createVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ),
           py::arg( "parameters_to_estimate" ),
           py::arg( "simulate_dynamics_on_creation" ) = true,
           R"doc(

Function to create object that propagates the dynamics.

Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
Depending on the specific input type (e.g. which function from the :ref:`\`\`propagator\`\`` module was used),
a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
function. When using default settings, calling this function will automatically propagate the dynamics.


Parameters
----------
bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object defining the physical environment, with all
    properties of artificial and natural bodies.

propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
    Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

simulate_dynamics_on_creation : Bool, default=True
    Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
    can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function

Returns
-------
:class:`~tudatpy.numerical_simulation.Simulator`
    Object that propagates the dynamics, and processes the results.






    )doc" );

    py::class_< tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m, "DynamicsSimulator", R"doc(No documentation found.)doc" );

    py::class_< tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                         "SingleArcSimulator",
                                                                         R"doc(

        Class for consolidating single arc dynamics simulation functionality.

        Class for consolidating all functionality required to perform single arc dynamics simulations.





     )doc" )
            .def( py::init< const tudat::simulation_setup::SystemOfBodies &,
                            const std::shared_ptr<
                                    tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const bool,
                            const bool,
                            const bool,
                            const bool,
                            const bool,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "integrator_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "are_equations_of_motion_to_be_integrated" ) = true,
                  py::arg( "clear_numerical_solutions" ) = false,
                  py::arg( "set_integrated_result" ) = false,
                  py::arg( "print_number_of_function_evaluations" ) = false,
                  py::arg( "print_dependent_variable_data" ) = true,
                  py::arg( "print_state_data" ) = true )
            .def_property_readonly(
                    "bodies",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getSystemOfBodies )
            .def_property_readonly(
                    "state_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getEquationsOfMotionNumericalSolution )
            .def_property_readonly(
                    "unprocessed_state_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getEquationsOfMotionNumericalSolutionRaw )
            .def_property_readonly(
                    "dependent_variable_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getDependentVariableHistory )
            .def_property_readonly(
                    "cumulative_computation_time_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getCumulativeComputationTimeHistory )
            .def_property_readonly(
                    "cumulative_number_of_function_evaluations",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getCumulativeNumberOfFunctionEvaluations )
            .def_property_readonly(
                    "integrator_settings",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getIntegratorSettings,
                    R"doc(

        Settings to create the numerical integrator that is to be used
        for the integration of the equations of motion


        :type: IntegratorSettings
     )doc" )
            .def_property_readonly(
                    "state_derivative_function",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getStateDerivativeFunction,
                    R"doc(

        **read-only**

        Function that performs a single state derivative function evaluation. This function takes the numerically propagated
        state, and current independent variable (time) as input, and returns the derivative of the state that is then used
        by the numerical integration routine. Typically, this function is NOT used directly by users.


        :type: Callable[[float, numpy.ndarray], numpy.ndarray]
     )doc" )
            .def_property_readonly(
                    "environment_updater",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getEnvironmentUpdater,
                    R"doc(

        **read-only**

        Object used in the propagation to update the environment, it uses the current time and numerically calculated state
        to update the translational state, rotational state, flight conditions, etc. of all bodies in the simulation to be
        consistent with this time and state.  Typically, this class is NOT used directly by users, but can be useful in specific situations.


        :type: EnvironmentUpdater
     )doc" )
            .def_property_readonly(
                    "propagation_results",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getPropagationResults,
                    R"doc(

        **read-only**

        This function retrieves all the results of the numerical propagation, stored
        in a single wrapper object


        :type: SingleArcSimulationResults
     )doc" )
            .def_property_readonly(
                    "propagation_termination_details",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getPropagationTerminationReason )
            .def_property_readonly(
                    "integration_completed_successfully",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            integrationCompletedSuccessfully );

    py::class_< tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "MultiArcDynamicsSimulator", R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "propagation_results",
                    &tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getMultiArcPropagationResults,
                    R"doc(No documentation found.)doc" );

    py::class_< tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "HybridArcDynamicsSimulator", R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "propagation_results",
                    &tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getHybridArcPropagationResults,
                    R"doc(No documentation found.)doc" );

    // TODO: Remove variationalOnlyIntegratorSettings
    py::class_< tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr<
                        tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "SingleArcVariationalSimulator",
            R"doc(

        Class for consolidating single arc variational dynamics functionality.

        Class for consolidating all functionality required to perform single arc variational dynamics simulations.





     )doc" )
            .def( py::init< const tudat::simulation_setup::SystemOfBodies &,
                            const std::shared_ptr<
                                    tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const std::shared_ptr<
                                    tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const bool,
                            const std::shared_ptr<
                                    tudat::numerical_integrators::IntegratorSettings< double > >,
                            const bool,
                            const bool,
                            const bool,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "integrator_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "integrate_equations_concurrently" ) = true,
                  py::arg( "variational_only_integrator_settings" ) = std::shared_ptr<
                          tudat::numerical_integrators::IntegratorSettings< TIME_TYPE > >( ),
                  py::arg( "clear_numerical_solutions" ) = false,
                  py::arg( "integrate_on_creation" ) = true,
                  py::arg( "set_integrated_result" ) = false,
                  py::arg( "print_dependent_variable_data" ) = true,
                  R"doc(

        Class constructor.

        Constructor through which the user can create instances of this class.
        Defines environment, propagation and integrations models, as well as a number of settings related
        to the (estimatable) parameters, w.r.t. which the variational equations are defined.

        .. note:: When using default settings, creating an object of
                  this type automatically triggers the propagation


        Parameters
        ----------
        bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
            Object defining the physical environment, with all
            properties of artificial and natural bodies.

        integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`
            Settings to create the numerical integrator that is to be used for the integration of the equations of motion.

        propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
            Settings to create the propagator that is to be used for the propagation of the dynamics.

        estimated_parameters : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
            Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
            linked to the environment and acceleration settings of the simulation.

        integrate_equations_concurrently : Bool, default = True
            Boolean defining whether equations of motion and variational equations are to be propagated concurrently
            (if true) or sequentially (of false).

        variational_only_integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`, default = []
            Settings to create the numerical integrator that is to be used for integration the variational equations.
            If none is given (default), the numerical integration settings are taken to be the same as the ones applied
            in the integration of the equations of motions (specified by the `integrator_settings` parameter).

        clear_numerical_solutions : Bool, default = False
            Boolean to determine whether to clear the raw numerical solution member variables
            and to reset the state transition interface after propagation.

        integrate_on_creation : Bool, default = True
            Boolean defining whether the propagation should be performed immediately (default), or at a later time
            (when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).

        set_integrated_result : Bool, default = True
            Boolean to determine whether to automatically use the integrated results to set ephemerides for the
            propagated bodies.





    )doc" )
            .def( "integrate_equations_of_motion_only",
                  &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                          integrateDynamicalEquationsOfMotionOnly,
                  py::arg( "initial_states" ),
                  R"doc(

        Function to trigger the integration of the (regular) equations of motion.


        Function to trigger the integration only of the (regular) equations of motion, resulting in a `state_history`.
        This step does not yet use variational dynamics. In order to also solve the variational equations,
        use the :func:`integrate_full_equations` member function.

        Returns
        -------
        None
            Creates / modifies the `state_history` property of the :class:`~tudatpy.numerical_simulation.SingleArcVariationalSolver` object.





    )doc" )
            .def( "integrate_full_equations",
                  &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                          integrateVariationalAndDynamicalEquations,
                  py::arg( "initial_states" ),
                  py::arg( "integrate_equations_concurrently" ) = true,
                  R"doc(

        Function to trigger the integration of variational and dynamical equations (equations of motion).


        Function to trigger the integration of the (regular) equations of motion as well as the variational equations,
        solving for `state_history` and `variational_equations_history`
        (in its two components `state_transition_matrix_history` & `sensitivity_matrix_history`).


        Parameters
        ----------
        initial_states : numpy.ndarray([m, 1])
            Initial state to be used for the parameters in the equations of motion.

        integrate_equations_concurrently : Bool, default = True
            Boolean defining whether equations of motion and variational equations are to be propagated concurrently
            (if true) or sequentially (of false).

        Returns
        -------
        None
            Creates / modifies the properties of the VariationalSolver object.





    )doc" )
            .def_property(
                    "parameter_vector",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE,
                                                              TIME_TYPE >::getParametersToEstimate,
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE,
                                                              TIME_TYPE >::resetParameterEstimate,
                    R"doc(

        Consolidated set of (estimatable) parameters
        w.r.t. the variational dynamics in the Variational Simulator are defined.


        :type: :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
     )doc" )
            .def_property_readonly(
                    "variational_equations_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getNumericalVariationalEquationsSolution,
                    R"doc(

        **read-only**

        List containing the solution of the variational equations, i.e. the
        state transition matrix history (first entry) and sensitivity matrix history (second vector entry).


        :type: list[ dict[float, numpy.ndarray] ]
     )doc" )
            .def_property_readonly(
                    "state_transition_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getStateTransitionMatrixSolution,
                    R"doc(

        **read-only**

        State transition matrix history, given as epoch with propagation epochs as keys.
        This is (alongside the `sensitivity_matrix_history`) the solution of the variational equations.


        :type: dict[float, numpy.ndarray]
     )doc" )
            .def_property_readonly(
                    "sensitivity_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getSensitivityMatrixSolution,
                    R"doc(

        **read-only**

        Sensitivity matrix history, given as epoch with propagation epochs as keys.
        This is (alongside the `state_transition_matrix_history`) the solution of the variational equations.


        :type: dict[float, numpy.ndarray]
     )doc" )
            .def_property_readonly(
                    "state_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getEquationsOfMotionSolution,
                    R"doc(

        **read-only**

        State history, given as epoch with propagation epochs as keys.
        This is the solution of the (propagated) equations of motion, describing the states along which
        the variational dynamics are solved.


        :type: dict[float, numpy.ndarray]
     )doc" )
            .def_property_readonly(
                    "dynamics_simulator",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE,
                                                              TIME_TYPE >::getDynamicsSimulator,
                    R"doc(

        **read-only**

        Simulator object containing all functionality for solving of the (regular) equations of motion.


        :type: :class:`~tudatpy.numerical_simulation.SingleArcSimulator`
     )doc" );

    py::class_< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m, "Estimator", R"doc(

        Class for consolidating all estimation functionality.

        Class for consolidating all functionality required to perform an estimation.





     )doc" )
            .def( py::init< const tss::SystemOfBodies &,
                            const std::shared_ptr<
                                    tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const std::vector< std::shared_ptr< tom::ObservationModelSettings > > &,
                            const std::shared_ptr< tp::PropagatorSettings< STATE_SCALAR_TYPE > >,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "observation_settings" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "integrate_on_creation" ) = true,
                  R"doc(

        Class constructor.

        Constructor through which the user can create instances of this class.
        Defines environment, propagation and integrations models, as well as a number of settings related
        to the estimatable parameters and observation settings.

        .. note:: When using default settings, creating an object of
                  this type automatically triggers the propagation


        Parameters
        ----------
        bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
            Object defining the physical environment, with all
            properties of artificial and natural bodies.

        estimated_parameters : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
            Object defining a consolidated set of estimatable parameters,
            linked to the environment and acceleration settings of the simulation.

        observation_settings : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings`
            List of settings objects, each object defining the observation model settings for one
            combination of observable and link geometry that is to be simulated.

        integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`
            Settings to create the numerical integrator that is to be
            used for the integration of the equations of motion

        propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
            Settings to create the propagator that is to be
            used for the propagation of dynamics

        integrate_on_creation : Bool, default = True
            Boolean defining whether the propagation should be
            performed immediately (default), or at a later time
            (when calling the :func:`perform_estimation` member function.





    )doc" )
            .def_property_readonly(
                    "observation_simulators",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getObservationSimulators,
                    R"doc(

        **read-only**

        Observation simulators contained in the Estimator object. A single observation simulator hosts
        the functionality for simulating a given observable over the defined link geometry.


        :type: list[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
     )doc" )
            .def_property_readonly(
                    "observation_managers",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getObservationManagers,
                    R"doc(

        **read-only**

        Observation managers contained in the Estimator object. A single observation manager can simulate observations and
        calculate observation partials for all link ends involved in the given observable type.


        :type: dict[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`, :class:`~tudatpy.numerical_simulation.estimation.ObservationManager` ]
     )doc" )
            .def_property_readonly(
                    "state_transition_interface",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getStateTransitionAndSensitivityMatrixInterface,
                    R"doc(

        **read-only**

        State transition and sensitivity matrix interface, setting the variational equations/dynamics in the
        Estimator object.


        :type: :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
     )doc" )
            .def( "perform_estimation",
                  &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                   TIME_TYPE >::estimateParameters,
                  py::arg( "estimation_input" ),
                  R"doc(

        Function to trigger the parameter estimation.


        Function to trigger the parameter estimation. Much of the process and requirements are similar to those described in the
        :func:`~tudatpy.numerical_simulation.Estimator.compute_covariance` function. This function uses an iterative least-squares
        estimate process to fit the data (inside ``estimation_input``) to the model defined by the inputs to the ``Estimator`` constructor.s


        Parameters
        ----------
        estimation_input : :class:`~tudatpy.numerical_simulation.estimation.EstimationInput`
            Object consolidating all relevant settings for the estimation
            This includes foremost the simulated observations, as well as a priori information about the estimatable parameters and convergence criteria for the least squares estimation.

        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.EstimationOutput`
            Object containing all outputs from the estimation process.





    )doc" )
            .def( "compute_covariance",
                  &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                   TIME_TYPE >::computeCovariance,
                  py::arg( "covariance_analysis_input" ),
                  R"doc(

        Function to perform a covariance analysis for the given observations and parameters


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


        Parameters
        ----------
        covariance_analysis_input : :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput`
            Object consolidating all relevant settings for the covariance analysis
            This includes foremost the simulated observations, as well as a priori information about the estimatable parameters

        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.vOutput`
            Object containing all outputs from the estimation process.





    )doc" )
            .def_property_readonly(
                    "variational_solver",
                    &tss::OrbitDeterminationManager< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getVariationalEquationsSolver,
                    R"doc(

        **read-only**

        Variational equations solver, which is used to manage and execute the numerical integration of
        equations of motion and variational equations/dynamics in the Estimator object.


        :type: :class:`~tudatpy.numerical_simulation.SingleArcVariationalSolver`
     )doc" );
};

}  // namespace numerical_simulation
}  // namespace tudatpy
