/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_simulator_bindings.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/variationalEquationsSolver.h"
#include "tudat/simulation/estimation_setup/createNumericalSimulator.h"

namespace py = pybind11;
namespace tp = tudat::propagators;
namespace tep = tudat::estimatable_parameters;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace dynamics
{
namespace simulator
{

void expose_simulator_variational_bindings( py::module& m )
{
    py::class_< tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                     "VariationalSimulator",
                                                                                                     R"doc(

         Base class for variational equations propagation.

         Base class for variational equations propagation.
         Derived classes :class:`~SingleArcVariationalSimulator`, :class:`~MultiArcVariationalSimulator` and
         :class:`~HybridArcVariationalSimulator` implement single-, multi- and hybrid-arc functionality, respectively.
      )doc" )
            .def_property_readonly( "state_transition_interface",
                                    &tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getStateTransitionMatrixInterface,
                                    R"doc(

         **read-only**

         State transition interface that includes the state transition and sensitivity matrices, which can be used for instance to propagate a covariance, see :func:`~tudatpy.estimation.estimation_analysis.propagate_covariance`.


         :type: CombinedStateTransitionAndSensitivityMatrixInterface
      )doc" );

    py::class_< tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                  "SingleArcVariationalSimulator",
                                                                                  R"doc(

         Class for single arc variational equations propagation.

      )doc" )
            .def( py::init< const tudat::simulation_setup::SystemOfBodies&,
                            const std::shared_ptr< tp::SingleArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                            const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const bool,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "integrate_equations_concurrently" ) = true,
                  py::arg( "integrate_on_creation" ) = true,
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
                  &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::integrateDynamicalEquationsOfMotionOnly,
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
                  &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::integrateVariationalAndDynamicalEquations,
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
            .def_property_readonly( "parameter_vector",
                                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getParametersToEstimate,
                                    R"doc(

         **read-only**

         Object containing the set of (estimatable) parameters
         w.r.t. the variational dynamics in this object are estimated


         :type: :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet`
      )doc" )
            .def_property_readonly(
                    "variational_equations_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getNumericalVariationalEquationsSolution,
                    R"doc(

         **read-only**

         Deprecated shorthand for the state-transition and sensitivity-matrix histories. Prefer the corresponding
         attributes of :attr:`variational_propagation_results`.

         :type: list[dict[float, numpy.ndarray]]
      )doc" )
            .def_property_readonly(
                    "state_transition_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getStateTransitionMatrixSolution,
                    R"doc(

         **read-only**

         Deprecated shorthand for ``variational_propagation_results.state_transition_matrix_history``.

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly( "sensitivity_matrix_history",
                                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getSensitivityMatrixSolution,
                                    R"doc(

         **read-only**

         Deprecated shorthand for ``variational_propagation_results.sensitivity_matrix_history``.

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly(
                    "state_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getEquationsOfMotionSolutionDouble,
                    R"doc(

         **read-only**

         Deprecated shorthand for ``variational_propagation_results.dynamics_results.state_history``.

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly(
                    "variational_propagation_results",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getVariationalPropagationResults,
                    R"doc(

        **read-only**

        Object containing all the results of the numerical propagation of both the dynamics and the variational equations, stored
        in a single wrapper object.

        Examples
        --------

        Assuming you have an object of type ``SingleArcVariationalSimulator`` named ``variational_simulator``, the numerical solution
        of the state transition matrix and sensitivity matrix can be retrieved as:

        .. code-block:: python

           # Retrieve object containing all propagation results
           results_container = variational_simulator.variational_propagation_results

           # Retrieve dictionaries with numerically integrated time histories of state transition and sensitivity matrices
           state_transition_matrix_history = results_container.state_transition_matrix_history
           sensitivity_matrix_history = results_container.sensitivity_matrix_history

           # Retrieve object containing propagation results of the dynamics (not including variational equations)
           dynamics_results_container = results_container.dynamics_results

           # Retrieve dictionaries with numerically integrated time history of states and associated dependent variables
           state_history = dynamics_results_container.state_history
           dependent_variable_history = dynamics_results_container.dependent_variable_history

        :type: SingleArcVariationalSimulationResults

      )doc" )

            .def_property_readonly( "dynamics_simulator",
                                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getDynamicsSimulator,
                                    R"doc(

         **read-only**

         Simulator object containing all functionality for solving of the (regular) equations of motion.


         :type: :class:`~tudatpy.dynamics.simulator.SingleArcSimulator`
      )doc" );

    py::class_< tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >,
                tp::VariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE > >( m,
                                                                                  "MultiArcVariationalSimulator",
                                                                                  R"doc(
         Class for multi-arc variational equations propagation.
      )doc" )
            .def( py::init< const tudat::simulation_setup::SystemOfBodies&,
                            const std::shared_ptr< tp::MultiArcPropagatorSettings< STATE_SCALAR_TYPE, TIME_TYPE > >,
                            const std::shared_ptr< tep::EstimatableParameterSet< STATE_SCALAR_TYPE > >,
                            const bool >( ),
                  py::arg( "bodies" ),
                  py::arg( "propagator_settings" ),
                  py::arg( "estimated_parameters" ),
                  py::arg( "integrate_on_creation" ) = false,
                  R"doc(
         Class constructor.
         Constructor through which the user can create instances of this class.
         Defines environment, propagation and integration models, as well as a number of settings related
         to the (estimatable) parameters, w.r.t. which the variational equations are defined.
         .. note:: When using default settings, creating an object of
                   this type does **not** automatically trigger the propagation
                   (``integrate_on_creation`` defaults to ``False`` for multi-arc).
                   Call :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` explicitly.
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
         propagation_start_times : list[float]
             List of start times for the separate arcs.
         integrate_equations_concurrently : bool, default = True
             Boolean defining whether equations of motion and variational equations are to be propagated concurrently
             (if true) or sequentially (if false).
         variational_only_integrator_settings : :class:`~tudatpy.dynamics.propagation_setup.integrator.IntegratorSettings`, default = []
             Settings to create the numerical integrator that is to be used for integrating the variational equations only.
             If none is given (default), the numerical integration settings are taken to be the same as the ones applied
             in the integration of the equations of motion (specified by the ``integrator_settings`` parameter).
         clear_numerical_solutions : bool, default = True
             Boolean to determine whether to clear the raw numerical solution member variables
             and to reset the state transition interface after propagation.
         integrate_on_creation : bool, default = False
             Boolean defining whether the propagation should be performed immediately (default), or at a later time
             (when calling the :func:`integrate_full_equations` or :func:`integrate_equations_of_motion_only` member function).
         reset_multi_arc_dynamics_after_propagation : bool, default = True
             Boolean defining whether to reset the multi-arc dynamics simulator state after each propagation.
         set_dependent_variables_interface : bool, default = False
             Boolean defining whether to set up the dependent variables interface after propagation.
     )doc" )
            .def( "integrate_equations_of_motion_only",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE,
                                                                   TIME_TYPE >::integrateDynamicalEquationsOfMotionOnly ),
                  py::arg( "initial_states" ),
                  R"doc(
         Function to trigger the integration of the (regular) equations of motion.
         Function to trigger the integration only of the (regular) equations of motion for all arcs, resulting in
         arc-wise state history entries. This step does not yet use variational dynamics. In order to also solve
         the variational equations, use the :func:`integrate_full_equations` member function.
         Parameters
         ----------
         initial_states : numpy.ndarray([m, 1])
             Concatenated initial states for all arcs, in the same order as defined in the estimated parameters.
         Returns
         -------
         None
             Creates / modifies the arc-wise state history entries of the
             :attr:`~tudatpy.dynamics.simulator.MultiArcVariationalSimulator.variational_propagation_results` object.
     )doc" )
            .def( "integrate_full_equations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&, const bool >(
                          &tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE,
                                                                   TIME_TYPE >::integrateVariationalAndDynamicalEquations ),
                  py::arg( "initial_states" ),
                  py::arg( "integrate_equations_concurrently" ) = true,
                  R"doc(
         Function to trigger the integration of variational and dynamical equations (equations of motion).
         Function to trigger the integration of the (regular) equations of motion as well as the variational equations
         for all arcs, solving for arc-wise ``state_history`` and arc-wise ``variational_equations_history``
         (in its two components ``state_transition_matrix_history`` & ``sensitivity_matrix_history``).
         Parameters
         ----------
         initial_states : numpy.ndarray([m, 1])
             Concatenated initial states for all arcs to be used for the equations of motion, in the same order
             as defined in the estimated parameters.
         integrate_equations_concurrently : bool, default = True
             Boolean defining whether equations of motion and variational equations are to be propagated concurrently
             (if true) or sequentially (if false).
         Returns
         -------
         None
             Creates / modifies the properties of the MultiArcVariationalSimulator object.
     )doc" )
            .def_property_readonly( "arc_wise_parameters_to_estimate",
                                    &tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getArcWiseParametersToEstimate,
                                    R"doc(
         **read-only**
         List of objects, one per arc, each containing the set of (estimatable) parameters
         w.r.t. which the variational dynamics in that arc are estimated.
         :type: list[ :class:`~tudatpy.dynamics.parameters.EstimatableParameterSet` ]
      )doc" )
            .def_property_readonly(
                    "variational_equations_history",
                    &tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getNumericalVariationalEquationsSolution,
                    R"doc(

         **read-only**

         Deprecated shorthand for the variational-equation histories of all arcs. Prefer the individual results in
         :attr:`variational_propagation_results`.

         :type: list[list[dict[float, numpy.ndarray]]]
      )doc" )
            .def_property_readonly(
                    "variational_propagation_results",
                    &tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getMultiArcVariationalPropagationResults,
                    R"doc(
        **read-only**
        Object containing all the results of the numerical propagation of both the dynamics and the variational equations
        for all arcs, stored in a single wrapper object.
        Examples
        --------
        Assuming you have an object of type ``MultiArcVariationalSimulator`` named ``variational_simulator``, the numerical
        solution of the state transition matrix and sensitivity matrix for each arc can be retrieved as:
        .. code-block:: python
           # Retrieve object containing all propagation results (all arcs)
           results_container = variational_simulator.variational_propagation_results
           # Retrieve list of single-arc result objects
           single_arc_results = results_container.single_arc_results
           # For a given arc index i, retrieve state transition and sensitivity matrix histories
           state_transition_matrix_history = single_arc_results[i].state_transition_matrix_history
           sensitivity_matrix_history = single_arc_results[i].sensitivity_matrix_history
           # Retrieve object containing propagation results of the dynamics for arc i (not including variational equations)
           dynamics_results_container = single_arc_results[i].dynamics_results
           # Retrieve dictionaries with numerically integrated time history of states and associated dependent variables
           state_history = dynamics_results_container.state_history
           dependent_variable_history = dynamics_results_container.dependent_variable_history
        :type: MultiArcVariationalSimulationResults
      )doc" )
            .def_property_readonly( "dynamics_simulator",
                                    &tp::MultiArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::getDynamicsSimulatorBase,
                                    R"doc(
         **read-only**
         Simulator object containing all functionality for solving of the (regular) equations of motion across all arcs.
         :type: :class:`~tudatpy.dynamics.simulator.MultiArcSimulator`
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
}

}  // namespace simulator
}  // namespace dynamics
}  // namespace tudatpy
