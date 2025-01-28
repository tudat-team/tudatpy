/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_numerical_simulation.h"

#include <pybind11/operators.h>
#include <pybind11/stl.h>

#include "expose_numerical_simulation/expose_environment.h"
#include "expose_numerical_simulation/expose_environment_setup.h"
#include "expose_numerical_simulation/expose_estimation.h"
#include "expose_numerical_simulation/expose_estimation_setup.h"
#include "expose_numerical_simulation/expose_propagation.h"
#include "expose_numerical_simulation/expose_propagation_setup.h"
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

namespace tudatpy {


    namespace numerical_simulation {

        void expose_numerical_simulation_simulator(py::module &m) {

            m.def(
                "get_integrated_type_and_body_list",
                &tp::getIntegratedTypeAndBodyList<STATE_SCALAR_TYPE, TIME_TYPE>,
                py::arg("propagator_settings"));

            m.def("get_single_integration_size", &tp::getSingleIntegrationSize,
                  py::arg("state_type"));

            m.def("create_dynamics_simulator",
                  &tss::createDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
                  py::arg("bodies"), py::arg("propagator_settings"),
                  py::arg("simulate_dynamics_on_creation") = true,
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






    )doc");


            py::class_<tp::DynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
                       std::shared_ptr<tp::DynamicsSimulator<STATE_SCALAR_TYPE,
                                                             TIME_TYPE>>>(
                m, "DynamicsSimulator", R"doc(No documentation found.)doc");

            py::class_<
                tp::SingleArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
                std::shared_ptr<tp::SingleArcDynamicsSimulator<
                    STATE_SCALAR_TYPE, TIME_TYPE>>,
                tp::DynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>(
                m, "SingleArcSimulator", R"doc(

        Class for consolidating single arc dynamics simulation functionality.

        Class for consolidating all functionality required to perform single arc dynamics simulations.





     )doc")
                .def(py::init<const tudat::simulation_setup::SystemOfBodies &,
                              const std::shared_ptr<
                                  tudat::numerical_integrators::
                                      IntegratorSettings<TIME_TYPE>>,
                              const std::shared_ptr<
                                  tp::PropagatorSettings<STATE_SCALAR_TYPE>>,
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
                    "bodies",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE, TIME_TYPE>::getSystemOfBodies)
                .def_property_readonly(
                    "state_history",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getEquationsOfMotionNumericalSolution)
                .def_property_readonly(
                    "unprocessed_state_history",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getEquationsOfMotionNumericalSolutionRaw)
                .def_property_readonly(
                    "dependent_variable_history",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getDependentVariableHistory)
                .def_property_readonly(
                    "cumulative_computation_time_history",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getCumulativeComputationTimeHistory)
                .def_property_readonly(
                    "cumulative_number_of_function_evaluations",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getCumulativeNumberOfFunctionEvaluations)
                .def_property_readonly(
                    "integrator_settings",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE, TIME_TYPE>::getIntegratorSettings,
                    R"doc(

        Settings to create the numerical integrator that is to be used
        for the integration of the equations of motion


        :type: IntegratorSettings
     )doc")
                .def_property_readonly(
                    "state_derivative_function",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getStateDerivativeFunction,
                    R"doc(

        **read-only**

        Function that performs a single state derivative function evaluation. This function takes the numerically propagated
        state, and current independent variable (time) as input, and returns the derivative of the state that is then used
        by the numerical integration routine. Typically, this function is NOT used directly by users.


        :type: Callable[[float, numpy.ndarray], numpy.ndarray]
     )doc")
                .def_property_readonly(
                    "environment_updater",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE, TIME_TYPE>::getEnvironmentUpdater,
                    R"doc(

        **read-only**

        Object used in the propagation to update the environment, it uses the current time and numerically calculated state
        to update the translational state, rotational state, flight conditions, etc. of all bodies in the simulation to be
        consistent with this time and state.  Typically, this class is NOT used directly by users, but can be useful in specific situations.


        :type: EnvironmentUpdater
     )doc")
                .def_property_readonly(
                    "propagation_results",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE, TIME_TYPE>::getPropagationResults,
                    R"doc(

        **read-only**

        This function retrieves all the results of the numerical propagation, stored
        in a single wrapper object


        :type: SingleArcSimulationResults
     )doc")
                .def_property_readonly(
                    "propagation_termination_details",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getPropagationTerminationReason)
                .def_property_readonly(
                    "integration_completed_successfully",
                    &tp::SingleArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::integrationCompletedSuccessfully);

            py::class_<
                tp::MultiArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
                std::shared_ptr<tp::MultiArcDynamicsSimulator<STATE_SCALAR_TYPE,
                                                              TIME_TYPE>>,
                tp::DynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>(
                m, "MultiArcDynamicsSimulator",
                R"doc(No documentation found.)doc")
                .def_property_readonly(
                    "propagation_results",
                    &tp::MultiArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getMultiArcPropagationResults,
                    R"doc(No documentation found.)doc");

            py::class_<
                tp::HybridArcDynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>,
                std::shared_ptr<tp::HybridArcDynamicsSimulator<
                    STATE_SCALAR_TYPE, TIME_TYPE>>,
                tp::DynamicsSimulator<STATE_SCALAR_TYPE, TIME_TYPE>>(
                m, "HybridArcDynamicsSimulator",
                R"doc(No documentation found.)doc")
                .def_property_readonly(
                    "propagation_results",
                    &tp::HybridArcDynamicsSimulator<
                        STATE_SCALAR_TYPE,
                        TIME_TYPE>::getHybridArcPropagationResults,
                    R"doc(No documentation found.)doc");

        };

    }  // namespace numerical_simulation
}  // namespace tudatpy
