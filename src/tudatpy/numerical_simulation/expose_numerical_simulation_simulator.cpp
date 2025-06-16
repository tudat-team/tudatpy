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

void expose_numerical_simulation_simulator( py::module &m )
{
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
                    "state_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getEquationsOfMotionNumericalSolution,
                    R"doc(
         **read-only**

         Shorthand for propagation_results.state_history

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly(
                    "unprocessed_state_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getEquationsOfMotionNumericalSolutionRaw,
                    R"doc(
         **read-only**

         Shorthand for propagation_results.unprocessed_state_history

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly(
                    "dependent_variable_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE,
                                                     TIME_TYPE >::getDependentVariableHistory,
                    R"doc(
         **read-only**

         Shorthand for propagation_results.dependent_variable_history

         :type: dict[float, numpy.ndarray]
      )doc" )
            .def_property_readonly(
                    "cumulative_computation_time_history",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getCumulativeComputationTimeHistory,
                    R"doc(
         **read-only**

         Shorthand for propagation_results.cumulative_computation_time_history

         :type: dict[float, float]
      )doc" )
            .def_property_readonly(
                    "cumulative_number_of_function_evaluations",
                    &tp::SingleArcDynamicsSimulator< STATE_SCALAR_TYPE, TIME_TYPE >::
                            getCumulativeNumberOfFunctionEvaluations,
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
    a single-, multi- or hybrid-arc simulator is created. The environment is typically created by the :func:`~tudatpy.numerical_simulation.environment_setup.create_system_of_bodies`
    function.

    .. note::

        When using default settings, calling this function will automatically propagate the dynamics.


    Parameters
    ----------
    bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
        Object defining the physical environment, with all
        properties of artificial and natural bodies.

    propagator_settings : :class:`~tudatpy.numerical_simulation.propagation_setup.propagator.PropagatorSettings`
        Settings to be used for the numerical propagation (dynamics type, termination conditions, integrator, etc.)

    simulate_dynamics_on_creation : bool, default=True
        Boolean defining whether to propagate the dynamics upon creation of the Simulator. If false, the dynamics c
        can be propagated at a later time by calling the :class:`~tudatpy.numerical_simulation.DynamicsSimulator.integrate_equations_of_motion` function

    Returns
    -------
    :class:`~tudatpy.numerical_simulation.DynamicsSimulator`
        Object that propagates the dynamics, and processes the results.
        Depending on the ``propagator_settings``, this object can be a single-, multi- or hybrid-arc simulator.






        )doc" );
};

}  // namespace numerical_simulation
}  // namespace tudatpy
