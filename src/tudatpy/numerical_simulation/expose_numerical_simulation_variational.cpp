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
#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "expose_numerical_simulation.h"
#include "scalarTypes.h"
#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/basics/timeType.h"

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

void expose_numerical_simulation_variational( py::module &m )
{
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

         integrate_equations_concurrently : bool, default = True
             Boolean defining whether equations of motion and variational equations are to be propagated concurrently
             (if true) or sequentially (of false).

         variational_only_integrator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.integrator.IntegratorSettings`, default = []
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
             Creates / modifies the `state_history` property of the :class:`~tudatpy.numerical_simulation.SingleArcVariationalSimulator` object.





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


         :type: :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
      )doc" ) // TIME_TODO
            .def_property_readonly(
                    "variational_equations_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getNumericalVariationalEquationsSolution,
                    R"doc(

         **read-only**

         List containing the solution of the variational equations, i.e. the
         state transition matrix history (first entry) and sensitivity matrix history (second vector entry).


         :type: list[ dict[float, numpy.ndarray] ]
      )doc" ) // TIME_TODO
            .def_property_readonly(
                    "state_transition_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getStateTransitionMatrixSolution,
                    R"doc(

         **read-only**

         State transition matrix history, given as epoch with propagation epochs as keys.
         This is (alongside the `sensitivity_matrix_history`) the solution of the variational equations.


         :type: dict[float, numpy.ndarray]
      )doc" ) // TIME_TODO
            .def_property_readonly(
                    "sensitivity_matrix_history",
                    &tp::SingleArcVariationalEquationsSolver< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getSensitivityMatrixSolution,
                    R"doc(

         **read-only**

         Sensitivity matrix history, given as epoch with propagation epochs as keys.
         This is (alongside the `state_transition_matrix_history`) the solution of the variational equations.


         :type: dict[float, numpy.ndarray]
      )doc" ) // TIME_TODO
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
 a single-, multi- or hybrid-arc variational solver is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
 function. When using default settings, calling this function will automatically propagate the dynamics.

 Parameters
 ----------
 bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
  Object defining the physical environment, with all
  properties of artificial and natural bodies.

 propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
  Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

 parameters_to_estimate : :class:`~tudatpy.numerical_simulation.estimation.EstimatableParameterSet`
  Object defining a consolidated set of (estimatable) parameters (w.r.t. variational equations are defined),
  linked to the environment and acceleration settings of the simulation.

 simulate_dynamics_on_creation : bool, default=True
  Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics
  can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.Simulator.integrate_equations_of_motion` function

 Returns
 -------
 :class:`~tudatpy.numerical_simulation.VariationalSimulator`
  Object that propagates the dynamics, and processes the results.

  )doc" );
}

}  // namespace numerical_simulation
}  // namespace tudatpy
