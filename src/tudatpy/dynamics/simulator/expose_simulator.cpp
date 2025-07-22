/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_simulator.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/propagators/stateTransitionMatrixInterface.h"
#include "tudat/simulation/propagation_setup/dynamicsSimulator.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tep = tudat::estimatable_parameters;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void expose_simulator( py::module& m )
{

    /*!
     *************** DYNAMICS SIMULATOR ***************
     */

    py::class_< tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "DynamicsSimulator",
            R"doc(Base class for propagation of dynamics (with derived classes implementing single-, multi- or hybrid-arc.)doc" );

    py::class_< tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >(
            m, "SingleArcSimulator", R"doc(

         Class for propagation of single arc dynamics.

         Class for propagation of single arc dynamics, from propagation settings and environment models, typically
         instantiated using :func:`~create_dynamics_simulator` function. See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagating_dynamics.html>`_ for more details.

      )doc" )
            .def_property_readonly( "bodies",
                                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                                     TIME_TYPE >::getSystemOfBodies,
                                    R"doc(
         **read-only**

         Object storing the set of bodies that comprise the physical environment.

         :type: SystemOfBodies
      )doc" )
            .def_property_readonly(
                    "propagator_settings",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getPropagatorSettings,
                    R"doc(
         **read-only**

         Settings for the propagation used to initialize the simulator

         :type: SingleArcPropagatorProcessingSettings
      )doc" )
            .def( "integrate_equations_of_motion",
                  &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::integrate,
                  py::arg( "initial_state" ),
                  R"doc(

 Function to reintegrate the equations of motion with a new initial state.

 Parameters
 ----------
 initial_state : numpy.ndarray
        New initial state from which the dynamics is to be propagated

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

         Object containing all the results of the numerical propagation, stored
         in a single wrapper object.


         :type: SingleArcSimulationResults
      )doc" )
        .def_property_readonly(
            "state_history_time_object",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
            getEquationsOfMotionNumericalSolution,
            R"doc(
         **read-only**

         Shorthand for propagation_results.state_history_time_object

         :type: dict[float, numpy.ndarray]
      )doc" )
        .def_property_readonly(
            "state_history",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
            getEquationsOfMotionNumericalSolutionDouble,
            R"doc(
         **read-only**

         Shorthand for propagation_results.state_history

         :type: dict[float, numpy.ndarray]
      )doc" )
        .def_property_readonly(
            "unprocessed_state_history_time_object",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
            getEquationsOfMotionNumericalSolutionRaw,
            R"doc(
         **read-only**

         Shorthand for propagation_results.unprocessed_state_history_time_object

         :type: dict[float, numpy.ndarray]
      )doc" )
        .def_property_readonly(
            "unprocessed_state_history",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
            getEquationsOfMotionNumericalSolutionRawDouble,
            R"doc(
         **read-only**

         Shorthand for propagation_results.unprocessed_state_history

         :type: dict[float, numpy.ndarray]
      )doc" )
        .def_property_readonly(
            "dependent_variable_history",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                TIME_TYPE >::getDependentVariableHistoryDouble,
            R"doc(
         **read-only**

         Shorthand for propagation_results.dependent_variable_history

         :type: dict[float, numpy.ndarray]
      )doc" )
        .def_property_readonly(
            "dependent_variable_history_time_object",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                TIME_TYPE >::getDependentVariableHistory,
            R"doc(
         **read-only**

         Shorthand for propagation_results.dependent_variable_history_time_object

         :type: dict[float, numpy.ndarray]
      )doc" )
        .def_property_readonly(
            "cumulative_computation_time_history",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
            getCumulativeComputationTimeHistoryDouble,
            R"doc(
         **read-only**

         Shorthand for propagation_results.cumulative_computation_time_history

         :type: dict[float, float]
      )doc" )
        .def_property_readonly(
            "cumulative_number_of_function_evaluations",
            &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
            getCumulativeNumberOfFunctionEvaluationsDouble,
            R"doc(
         **read-only**

         Shorthand for propagation_results.cumulative_number_of_function_evaluations

         :type: dict[float, int]
      )doc" )
            .def_property_readonly(
                    "propagation_termination_details",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getPropagationTerminationReason,
                    R"doc(
         **read-only**

         Shorthand for propagation_results.termination_details

         :type: PropagationTerminationDetails
      )doc" )
            .def_property_readonly(
                    "integration_completed_successfully",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::integrationCompletedSuccessfully,
                    R"doc(
         **read-only**

         Shorthand for propagation_results.integration_completed_successfully

         :type: bool
      )doc" );

    py::class_< tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                         "MultiArcSimulator",
                                                                         R"doc(

         Class for propagation of multi-arc dynamics.

         Class for propagation of multi-arc dynamics, from propagation settings and environment models, typically
         instantiated using :func:`~create_dynamics_simulator` function. See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagating_dynamics.html>`_ for more details.

      )doc" )
            .def( "integrate_equations_of_motion",
                  &tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::integrate,
                  py::arg( "initial_state" ),
                  R"doc(

 Function to reintegrate the equations of motion with a new initial state.

 Parameters
 ----------
 initial_state : numpy.ndarray
        New initial state from which the dynamics is to be propagated, consisting of concatenated initial states of each arc

      )doc" )
            .def_property_readonly(
                    "propagation_results",
                    &tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getMultiArcPropagationResults,
                    R"doc(

         **read-only**

         Object containing all the results of the numerical propagation, stored
         in a single wrapper object


         :type: MultiArcSimulationResults
         )doc" )
            .def_property_readonly(
                    "single_arc_simulators",
                    &tp::MultiArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                    TIME_TYPE >::getSingleArcDynamicsSimulators,
                    R"doc(

         **read-only**

         List of single-arc simulators, that were used to propagate the dynamics on the constituent single arcs


         :type: list[SingleArcSimulator]
         )doc" );

    py::class_< tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::DynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                         "HybridArcSimulator",
                                                                         R"doc(

         Class for propagation of hybrid-arc dynamics.

         Class for propagation of hybrid-arc dynamics, from propagation settings and environment models, typically
         instantiated using :func:`~create_dynamics_simulator` function. See `user guide <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagating_dynamics.html>`_ for more details.

      )doc" )
            .def( "integrate_equations_of_motion",
                  &tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::integrate,
                  py::arg( "initial_state" ),
                  R"doc(

 Function to reintegrate the equations of motion with a new initial state.

 Parameters
 ----------
 initial_state : numpy.ndarray
        New initial state from which the dynamics is to be propagated, consisting of concatenated initial state single-arc portion of dynamics, followed by concatenated initial states of the constituent arcs of the multi-arc portion

      )doc" )
            .def_property_readonly(
                    "single_arc_simulator",
                    &tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getSingleArcDynamicsSimulator,
                    R"doc(

         **read-only**

         Object used to propagate the single-arc portion of the hybrid-arc dynamics


         :type: SingleArcSimulator
         )doc" )
            .def_property_readonly(
                    "multi_arc_simulator",
                    &tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getMultiArcDynamicsSimulator,
                    R"doc(

         **read-only**

         Object used to propagate the multi-arc portion of the hybrid-arc dynamics

         :type: MultiArcSimulator
         )doc" )
            .def_property_readonly(
                    "propagation_results",
                    &tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getHybridArcPropagationResults,
                    R"doc(

         **read-only**

         Object containing all the results of the numerical propagation, stored
         in a single wrapper object


         :type: HybridArcSimulationResults
         )doc" )
            .def_property_readonly(
                    "propagation_results",
                    &tp::HybridArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getHybridArcPropagationResults,
                    R"doc(

         **read-only**

         Object containing all the results of the numerical propagation, stored
         in a single wrapper object


         :type: HybridArcSimulationResults
         )doc" );

    m.def( "create_dynamics_simulator",
           &tss::createDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ),
           py::arg( "simulate_dynamics_on_creation" ) = true,
           R"doc(

    Function to create object that propagates the dynamics.

    Function to create object that propagates the dynamics, as specified by propagator settings, and the physical environment.
    Depending on the specific input type (e.g. which function from the :ref:`propagator` module was used),
    a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`
    function.

    .. note::

        When using default settings, calling this function will automatically propagate the dynamics.


    Parameters
    ----------
    bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
        Object defining the physical environment, with all
        properties of artificial and natural bodies.

    propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
        Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

    simulate_dynamics_on_creation : bool, default=True
        Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
        can be propagated at a later time by calling the :func:`~tudatpy.dynamics.simulator.DynamicsSimulator.integrate_equations_of_motion` function

    Returns
    -------
    :class:`~tudatpy.dynamics.simulator.DynamicsSimulator`
        Object that propagates the dynamics, and processes the results.
        Depending on the ``propagator_settings``, this object can be a single-, multi- or hybrid-arc simulator.






        )doc" );


    /*!
     *************** VARIATIONAL EQUATIONS SOLVER ***************
     */

    py::class_< tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "VariationalSimulator",
            R"doc(

         Base class for variational equations propagation.

         Base class for variational equations propagation.
         Derived classes :class:`~SingleArcVariationalSimulator`, :class:`~MultiArcVariationalSimulator` and
         :class:`~HybridArcVariationalSimulator` implement single-, multi- and hybrid-arc functionality, respectively."

      )doc" );

    py::class_< tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr<
                        tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > > >(
            m,
            "SingleArcVariationalSimulator",
            R"doc(

         Class for single arc variational equations propagation.

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
         bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
             Object defining the physical environment, with all
             properties of artificial and natural bodies.

         integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`
             Settings to create the numerical integrator that is to be used for the integration of the equations of motion.

         propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
             Settings to create the propagator that is to be used for the propagation of the dynamics.

         estimated_parameters : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
             Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
             linked to the environment and acceleration settings of the simulation.

         integrate_equations_concurrently : bool, default = True
             Boolean defining whether equations of motion and variational equations are to be propagated concurrently
             (if true) or sequentially (of false).

         variational_only_integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`, default = []
             Settings to create the numerical integrator that is to be used for integration the variational equations.
             If none is given (default), the numerical integration settings are taken to be the same as the ones applied
             in the integration of the equations of motions (specified by the `integrator_settings` parameter).

         clear_numerical_solutions : bool, default = False
             Boolean to determine whether to clear the raw numerical solution member variables
             and to reset the state transition interface after propagation.

         integrate_on_creation : bool, default = True
             Boolean defining whether the propagation should be performed immediately (default), or at a later time
             (when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).

         set_integrated_result : bool, default = True
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
             Creates / modifies the `state_history` property of the :class:`~tudatpy.dynamics.simulator.SingleArcVariationalSimulator` object.





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

         integrate_equations_concurrently : bool, default = True
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


         :type: :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
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


         :type: :class:`~tudatpy.dynamics.simulator.SingleArcSimulator`
      )doc" );

    m.def( "create_variational_equations_solver",
           &tss::createVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "propagator_settings" ),
           py::arg( "parameters_to_estimate" ),
           py::arg( "simulate_dynamics_on_creation" ) = true,
           R"doc(

 Function to create object that propagates the dynamics and variational equations.

 Function to create object that propagates the dynamics and variational equations, as specified by propagator settings, the physical environment, and a set of parameters for which to compute the partials.
 Depending on the specific input type (e.g. which function from the :ref:`propagator` module was used to define the propagator settings),
 a single-, multi- or hybrid-arc variational solver is created. The environment is typically created by the :func:`~tudatpy.dynamics.environment_setup.create_system_of_bodies`
 function. When using default settings, calling this function will automatically propagate the dynamics.

 Parameters
 ----------
 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
  Object defining the physical environment, with all
  properties of artificial and natural bodies.

 propagator_settings : :class:`~tudatpy.dynamics.propagation_setup.propagator.PropagatorSettings`
  Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

 parameters_to_estimate : :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
  Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
  linked to the environment and acceleration settings of the simulation.

 simulate_dynamics_on_creation : bool, default=True
  Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics
  can be propagated at a later time by calling the :func:`~tudatpy.dynamics.DynamicsSimulator.integrate_equations_of_motion` function

 Returns
 -------
 :class:`~tudatpy.dynamics.simulator.VariationalSimulator`
  Object that propagates the dynamics, and processes the results.

  )doc" );


    /*!
     *************** STATE TRANSITION INTERFACE ***************
     */
         
    py::class_< tp::CombinedStateTransitionAndSensitivityMatrixInterface,
                std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > >(
            m,
            "CombinedStateTransitionAndSensitivityMatrixInterface",
            R"doc(

         Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.

         Class establishing an interface to the State Transition and Sensitivity Matrices.
         Instances of this class are instantiated automatically upon creation of :class:`~tudatpy.estimation.estimation_analysis.Estimator` objects,
         using the simulation information in the observation, propagation and integration settings that the :class:`~tudatpy.estimation.estimation_analysis.Estimator` instance is linked to.





      )doc" )
            .def( "state_transition_sensitivity_at_epoch",
                  &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                          getCombinedStateTransitionAndSensitivityMatrix,
                  py::arg( "time" ),
                  py::arg( "add_central_body_dependency" ) = true,
                  py::arg( "arc_defining_bodies" ) = std::vector< std::string >( ),
                  R"doc(

         Function to get the concatenated state transition and sensitivity matrix at a given time.

         Function to get the concatenated state transition and sensitivity matrix at a given time.
         Entries corresponding to parameters which are not active at the current arc are omitted.


         Parameters
         ----------
         time : float
             Time at which concatenated state transition and sensitivity matrix are to be retrieved.
         Returns
         -------
         numpy.ndarray[numpy.float64[m, n]]
             Concatenated state transition and sensitivity matrix at a given time.





     )doc" )
            .def( "full_state_transition_sensitivity_at_epoch",
                  &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                          getFullCombinedStateTransitionAndSensitivityMatrix,
                  py::arg( "time" ),
                  py::arg( "add_central_body_dependency" ) = true,
                  py::arg( "arc_defining_bodies" ) = std::vector< std::string >( ),
                  R"doc(


         Parameters
         ----------
         time : float
             Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
         Returns
         -------
         numpy.ndarray[numpy.float64[m, n]]
             Full concatenated state transition and sensitivity matrix at a given time.





     )doc" )
            .def_property_readonly( "state_transition_size",
                                    &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                                            getStateTransitionMatrixSize,
                                    R"doc(

         **read-only**

         Size of the (square) state transition matrix.

         :type: int
      )doc" )
            .def_property_readonly( "sensitivity_size",
                                    &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                                            getSensitivityMatrixSize,
                                    R"doc(

         **read-only**

         Number of columns in the sensitivity matrix.

         :type: int
      )doc" )
            .def_property_readonly( "full_parameter_size",
                                    &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                                            getFullParameterVectorSize,
                                    R"doc(

         **read-only**

         Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.

         :type: int
      )doc" );

}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
