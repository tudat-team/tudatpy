/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_propagation.h"

#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/propagators.h>

#include "scalarTypes.h"


namespace py = pybind11;

namespace ta = tudat::aerodynamics;
namespace tp = tudat::propagators;
namespace tpr = tudat::propulsion;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;

namespace tudat {

    namespace aerodynamics {

        class PyAerodynamicGuidance : public ta::AerodynamicGuidance {
           public:
            /* Inherit the constructors */
            using AerodynamicGuidance::AerodynamicGuidance;

            using AerodynamicGuidance::currentAngleOfAttack_;
            using AerodynamicGuidance::currentAngleOfSideslip_;
            using AerodynamicGuidance::currentBankAngle_;

            void updateGuidance(const double currentTime) override {
                PYBIND11_OVERLOAD_PURE(void, AerodynamicGuidance,
                                       updateGuidance, currentTime);
            }
        };

    }  // namespace aerodynamics

}  // namespace tudat


namespace tudatpy {
    namespace numerical_simulation {
        namespace propagation {


            void expose_propagation(py::module &m) {
                py::class_<ta::AerodynamicGuidance, ta::PyAerodynamicGuidance,
                           std::shared_ptr<ta::AerodynamicGuidance>>(
                    m, "AerodynamicGuidance")
                    .def(py::init<>())
                    .def("updateGuidance",
                         &ta::AerodynamicGuidance::updateGuidance,
                         py::arg("current_time"))
                    .def_readwrite(
                        "angle_of_attack",
                        &ta::PyAerodynamicGuidance::currentAngleOfAttack_)
                    .def_readwrite(
                        "bank_angle",
                        &ta::PyAerodynamicGuidance::currentBankAngle_)
                    .def_readwrite(
                        "sideslip_angle",
                        &ta::PyAerodynamicGuidance::currentAngleOfSideslip_);


                py::class_<tba::TorqueModel, std::shared_ptr<tba::TorqueModel>>(
                    m, "TorqueModel");


                m.def("get_single_integration_size",
                      &tp::getSingleIntegrationSize, py::arg("state_type"));

                m.def("get_single_integration_differential_equation_order",
                      &tp::getSingleIntegrationDifferentialEquationOrder,
                      py::arg("state_type"));

                m.def("get_generalized_acceleration_size",
                      &tp::getGeneralizedAccelerationSize,
                      py::arg("state_type"));

                m.def("get_state_of_bodies",
                      py::overload_cast<const std::vector<std::string> &,
                                        const std::vector<std::string> &,
                                        const tss::SystemOfBodies &,
                                        const TIME_TYPE>(
                          &tp::getInitialStatesOfBodies<TIME_TYPE,
                                                        STATE_SCALAR_TYPE>),
                      py::arg("bodies_to_propagate"), py::arg("central_bodies"),
                      py::arg("body_system"), py::arg("initial_time"),
                      R"doc(

Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time.

Function to get the translational states of a set of bodies, with respect to some set of central bodies, at the requested time. This function
is typically used to extract an initial state for a propagation of a set of bodies, for which the initial state is extracted from the
existing ephemerides of the bodies.


Parameters
----------
bodies_to_propagate : list[str]
    List of names of bodies for which the state is to be extracted
central_bodies : list[str]
    List of central bodies, w.r.t. which the states are to be computed (in the same order as ``bodies_to_propagate``)
bodies_to_propagate : SystemOfBodies
    System of bodies that define the environment
initial_time : float
    Time at which the states are to be extracted from the environment
Returns
-------
numpy.ndarray
    Vector of size :math:`6\times N`, with the translational states of each entry of body from
    ``bodies_to_propagate`` w.r.t. the corresponding central body in ``central_bodies``.







    )doc");


                m.def("get_initial_state_of_bodies",
                      py::overload_cast<const std::vector<std::string> &,
                                        const std::vector<std::string> &,
                                        const tss::SystemOfBodies &,
                                        const TIME_TYPE>(
                          &tp::getInitialStatesOfBodies<TIME_TYPE,
                                                        STATE_SCALAR_TYPE>),
                      py::arg("bodies_to_propagate"), py::arg("central_bodies"),
                      py::arg("body_system"), py::arg("initial_time"));

                m.def(
                    "get_initial_state_of_body",  // overload [2/2]
                    py::overload_cast<const std::string &, const std::string &,
                                      const tss::SystemOfBodies &,
                                      const TIME_TYPE>(
                        &tp::getInitialStateOfBody<TIME_TYPE,
                                                   STATE_SCALAR_TYPE>),
                    py::arg("body_to_propagate"), py::arg("central_body"),
                    py::arg("bodies"), py::arg("initial_time"));

                m.def(
                    "get_initial_rotational_state_of_body",
                    py::overload_cast<const std::string &, const std::string &,
                                      const tss::SystemOfBodies &,
                                      const TIME_TYPE>(
                        &tp::getInitialRotationalStateOfBody<
                            TIME_TYPE, STATE_SCALAR_TYPE>),
                    py::arg("body_to_propagate"), py::arg("base_orientation"),
                    py::arg("bodies"), py::arg("initial_time"));

                py::class_<
                    tp::DampedInitialRotationalStateResults<TIME_TYPE,
                                                            STATE_SCALAR_TYPE>,
                    std::shared_ptr<tp::DampedInitialRotationalStateResults<
                        TIME_TYPE, STATE_SCALAR_TYPE>>>(
                    m, "RotationalProperModeDampingResults",
                    R"doc(

        Object that stores the results of the algorithm to damp the proper mode of rotational dynamics for an initial state,
        as computed by the :func:`~get_damped_proper_mode_initial_rotational_state` function






     )doc")
                    .def_readwrite(
                        "damped_initial_state",
                        &tp::DampedInitialRotationalStateResults<
                            TIME_TYPE, STATE_SCALAR_TYPE>::initialState_,
                        R"doc(

        Initial state produced by the damping algorithm, for which the signature of the proper mode should be
        removed (or at least, substantially reduced). Note that this initial state corresponds to the *full* state vector
        that is provided to the ``get_damped_proper_mode_initial_rotational_state`` function (e.g. is size 7
        for rotational dynamics of a single body, size 13 for coupled orbital-rotational dynamics of a single body, etc.)


        :type: numpy.ndarray
     )doc")
                    .def_readwrite("forward_backward_states",
                                   &tp::DampedInitialRotationalStateResults<
                                       TIME_TYPE, STATE_SCALAR_TYPE>::
                                       forwardBackwardPropagatedStates_,
                                   R"doc(

        Data structure that contains the full state histories used by the damping algorithm. The contents are are as follows:

        * The :math:`i^{th}` entry of the list corresponds to the :math:`i^{th}` iteration of the forward-backward propagation
        * Each tuple in the list contains two dictionaries, the first one corresponding to the forward propagation results, the seconds one to the backward propagation results


        :type: list[tuple[dict[float,numpy.ndarray],dict[float,numpy.ndarray]]]
     )doc")
                    .def_readwrite("forward_backward_dependent_variables",
                                   &tp::DampedInitialRotationalStateResults<
                                       TIME_TYPE, STATE_SCALAR_TYPE>::
                                       forwardBackwardDependentVariables_,
                                   R"doc(

        As ``forward_backward_states``, but for the dependent variables.


        :type: list[tuple[dict[float,numpy.ndarray],dict[float,numpy.ndarray]]]
     )doc");

                m.def("get_damped_proper_mode_initial_rotational_state",
                      py::overload_cast<
                          const tss::SystemOfBodies &,
                          const std::shared_ptr<tp::SingleArcPropagatorSettings<
                              STATE_SCALAR_TYPE, TIME_TYPE>>,
                          const double, const std::vector<double>, const bool>(
                          &tp::getZeroProperModeRotationalStateWithStruct<
                              TIME_TYPE, STATE_SCALAR_TYPE>),
                      py::arg("bodies"), py::arg("propagator_settings"),
                      py::arg("body_mean_rotational_rate"),
                      py::arg("dissipation_times"),
                      py::arg("propagate_undamped") = true,
                      R"doc(

Function to compute an initial rotational state for which the proper mode of rotation is damped.

Function to compute an initial rotational state for which the proper mode of rotation is damped, using the algorithm
used by Rambaux et al. (2010) to compute an initial rotational state for Phobos. This algorithm propagates the
dynamics of the system a number of times, with the settings specified by the user and a specific modification to
damp the proper mode. Since a number of propagations are performed by this function, it may take some time to run.
Specifically, the algorithm works as follows:

* Introduce a damping torque (see below) to damp the proper mode, with damping time :math:`\tau_{d}`
* Propagate the dynamics forward in time for a period of :math:`10\tau_{d}`
* Remove the virtual torque, and propagate the dynamics back to the initial time :math:`t_{0}`
* Repeat the above for the list of damping times provided by the user

The state after the final backwards propagation to :math:`t_{0}` is provided as output by this function, to be
used as damped initial state. The output from this function also provides the user access to the full state history
and dependent variable history of the forward and backward propagations, to allow a user to track and validate
the progress of the algorithm.

The damping torque :math:`\Gamma` is defined as follows:

.. math::
   \boldsymbol{\Gamma}= -\frac{1}{\tau_{d}}\mathbf{I}\begin{pmatrix}\omega_{x}\\ \omega_{y}\\ \omega_{x}-\omega_{p} \end{pmatrix}

where :math:\mathbf{I}` is the body's inertia tensor (in its body-fixed frame), :math:`\tau_{d}` the damping time of the
current propagation, and :math:`\omega_{x}, \omega_{y}, \omega_{z}` the body's current rotation about its
body-fixed, x-, y- and z-axes, respectively. The damping torque is implemented to damp out all rotations along
the body-fixed x- and y-axes, and any deviations from constant rotation with frequency :\omega_{p}: about the body-fixed z-axis.

.. note:: The mean rotation rate of the body :math:`\omega_{p}` is a user-defined input, and must be tuned to the dynamics of the system.


Parameters
----------
bodies : SystemOfBodies
    Set of body objects that defines the environment
propagator_settings : SingleArcPropagatorSettings
    Propagator settings for the dynamics of which the initial rotational state is to be damped. These propagator
    settings must be for rotational dynamics only, or for multi-type rotational dynamics that contains rotational
    dynamics for a single body (e.g. translational-rotational dynamics for a single body)

body_mean_rotational_rate : float
    Mean rotational rate :math:`\omega_{p}` to which the damping algorithm will force the body-fixed rotation about its z-axis.
dissipation_times : list[ float ]
    List of damping times :math:`\tau_{d}` for which the algorithm is to be run. Note that this list should be organized in ascending order for the algorithm to perform properly
propagate_undamped : bool, default = True
    Boolean defining whether the first forward/backward propagation performed by the damping algorithm has damping turned off (damping turned off if True, damping turned on if False).
    Propagating without any damping before starting the damping algorithm is useful for verification purposes, but not required for the algorithm itself.

Returns
-------
DampedInitialRotationalStateResults
    Object that contains the results of the damping algorithm (final damped rotational state, and forward/backward propagation results).






    )doc");

                m.def("combine_initial_states",
                      &tp::createCombinedInitialState<STATE_SCALAR_TYPE,
                                                      TIME_TYPE>,
                      py::arg("propagator_settings_per_type"),
                      R"doc(

Function to retrieve the initial state for a list of propagator settings.

Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
organized in a single container.


Parameters
----------
propagator_settings_per_type : dict
    Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
Returns
-------
numpy.ndarray
    Vector of initial states, sorted in order of IntegratedStateType, and then in the order of the vector of SingleArcPropagatorSettings of given type.






    )doc");

                py::class_<
                    tba::AccelerationModel<Eigen::Vector3d>,
                    std::shared_ptr<tba::AccelerationModel<Eigen::Vector3d>>>(
                    m, "AccelerationModel");

                py::class_<tba::MassRateModel,
                           std::shared_ptr<tba::MassRateModel>>(
                    m, "MassRateModel");


                py::enum_<tp::PropagationTerminationReason>(
                    m, "PropagationTerminationReason",
                    R"doc(

        Enumeration of types of termination of propagation.

     )doc")
                    .value(
                        "propagation_never_run",
                        tp::PropagationTerminationReason::propagation_never_run)
                    .value("unknown_reason",
                           tp::PropagationTerminationReason::
                               unknown_propagation_termination_reason)
                    .value("termination_condition_reached",
                           tp::PropagationTerminationReason::
                               termination_condition_reached)
                    .value("runtime_error_caught_in_propagation",
                           tp::PropagationTerminationReason::
                               runtime_error_caught_in_propagation)
                    .value("nan_or_inf_detected_in_state",
                           tp::PropagationTerminationReason::
                               nan_or_inf_detected_in_state)
                    .export_values();

                py::class_<tp::PropagationTerminationDetails,
                           std::shared_ptr<tp::PropagationTerminationDetails>>(
                    m, "PropagationTerminationDetails",
                    R"doc(

        Class that provides information on the reason for the
        termination of the propagation.






     )doc")
                    .def_property_readonly("termination_reason",
                                           &tp::PropagationTerminationDetails::
                                               getPropagationTerminationReason,
                                           R"doc(

        Enum defining the reason the propagation was terminated


        :type: PropagationTerminationReason
     )doc")
                    .def_property_readonly("terminated_on_exact_condition",
                                           &tp::PropagationTerminationDetails::
                                               getTerminationOnExactCondition,
                                           R"doc(

        Boolean defining whether the propagation was terminated on an *exact* final condition,
        or once the propagation went *past* the determined final condition. The choice of behaviour is
        defined by the termination settings provided as input to the Simulator object. This variable only
        has a meaningful definition if the ``termination_reason`` has value ``termination_condition_reached``


        :type: bool
     )doc");

                py::class_<
                    tp::PropagationTerminationDetailsFromHybridCondition,
                    std::shared_ptr<
                        tp::PropagationTerminationDetailsFromHybridCondition>,
                    tp::PropagationTerminationDetails>(
                    m, "PropagationTerminationDetailsFromHybridCondition",
                    R"doc(

        Class that provides information on the reason for the termination of the propagation, for hybrid termination conditions


        Derived class from :class:`PropagationTerminationDetails` that provides information on the reason for the termination of the propagation,
        for the case of hybrid termination conditions (defined using the :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.hybrid_termination`)
        function





     )doc")
                    .def_property_readonly(
                        "was_condition_met_when_stopping",
                        &tp::PropagationTerminationDetailsFromHybridCondition::
                            getWasConditionMetWhenStopping,
                        R"doc(

        List of booleans defining, per entry in ``termination_settings`` when calling :func:`~tudatpy.numerical_simulation.propagation_setup.propagator.hybrid_termination`,
        whether the corresponding entry of the hybrid termination settings was met or not


        :type: list[bool]
     )doc");

                py::class_<tp::DependentVariablesInterface<TIME_TYPE>,
                           std::shared_ptr<
                               tp::DependentVariablesInterface<TIME_TYPE>>>(
                    m, "DependentVariablesInterface",
                    R"doc(No documentation found.)doc");

                py::class_<tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>,
                           std::shared_ptr<tp::SimulationResults<
                               STATE_SCALAR_TYPE, TIME_TYPE>>>(
                    m, "SimulationResults",
                    R"doc(

        Base class for objects that store all results of a numerical propagation.

        Base class for objects that store all results of a numerical propagation. Derived class are implemented for single-, multi- and hybrid-arc propagation of both dynamics and variational equations





     )doc")
                    .def_property_readonly(
                        "dependent_variable_interface",
                        &tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>::
                            getDependentVariablesInterface,
                        R"doc(No documentation found.)doc");

                py::class_<tp::SingleArcSimulationResults<STATE_SCALAR_TYPE,
                                                          TIME_TYPE>,
                           std::shared_ptr<tp::SingleArcSimulationResults<
                               STATE_SCALAR_TYPE, TIME_TYPE>>,
                           tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "SingleArcSimulationResults",
                    R"doc(

        Class that stores all the results (including logging data) of a single-arc propagation






     )doc")
                    .def_property_readonly(
                        "state_history",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getEquationsOfMotionNumericalSolution,
                        R"doc(

        **read-only**

        Numerical solution of the equations of motion as key-value pairs. The key denotes the epoch. The value contains the
        numerically calculated state at this epoch. For this function, the states are always converted to so-called
        'processed' formulations (e.g. Cartesian states for translational dynamics), see `here <https://docs.tudat.space/en/latest/_src_user_guide/state_propagation/propagation_setup/processed_propagated_elements.html>`_
        for details. For the history of the states that were actually propagated, use the ``unprocessed_state_history``.

        .. note:: The propagated state at each epoch contains the state types in the following order: Translational ( **T** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
                  When propagating two bodies, an example of what the output state would look like is for instance:
                  [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ] The specifics can be retrieved using the :attr:`state_ids` attribute of this class

        .. note:: For propagation of translational dynamics using cowell
                  propagator, the conventional and propagated
                  coordinates are identical.


        :type: dict[float, numpy.ndarray]
     )doc")
                    .def_property_readonly(
                        "state_history_float",
                        &tp::SingleArcSimulationResults<STATE_SCALAR_TYPE,
                                                        TIME_TYPE>::
                            getEquationsOfMotionNumericalSolutionDouble)
                    .def_property_readonly(
                        "state_history_float_split",
                        &tp::SingleArcSimulationResults<STATE_SCALAR_TYPE,
                                                        TIME_TYPE>::
                            getEquationsOfMotionNumericalSolutionDoubleSplit)
                    .def_property_readonly(
                        "unprocessed_state_history",
                        &tp::SingleArcSimulationResults<STATE_SCALAR_TYPE,
                                                        TIME_TYPE>::
                            getEquationsOfMotionNumericalSolutionRaw,
                        R"doc(

        **read-only**

        Numerical solution of the equations of motion as key-value pairs, without any processing applied. The key denotes the epoch. The value contains the
        numerically calculated state at this epoch. This attribute contains the states of the propagated bodies expressed in the
        "raw" form in which the propagation took place. For instance, when using a Gauss-Kepler propagation scheme, this
        attribute will contain the numerically propagated Keplerian elements at each time epoch


        :type: dict[float, numpy.ndarray]
     )doc")
                    .def_property_readonly(
                        "dependent_variable_history",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getDependentVariableHistory,
                        R"doc(

        **read-only**

        Dependent variables computed during the propagation as key-value pairs.
        The vector of all dependent variables concatenated into a single vector as value, with the epoch as key.
        They order of the concatenated dependent variables in a single value is provided by the ``dependent_variable_ids`` attribute of this object.


        :type: dict[float, numpy.ndarray]
     )doc")
                    .def_property_readonly(
                        "cumulative_computation_time_history",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getCumulativeComputationTimeHistory,
                        R"doc(

        **read-only**

        History of cumulative computation time in seconds needed during the propagation as key-value
        pairs. At each epoch (key) the computation time (value) in seconds is the total computation time
        used up to and including that time step. This includes the total time up to and including the current time step,
        since the beginning of the (single-arc) propagation.


        :type: dict[float, float]
     )doc")
                    .def_property_readonly(
                        "cumulative_computation_time_history",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getCumulativeComputationTimeHistory,
                        R"doc(

        **read-only**

        History of cumulative computation time in seconds needed during the propagation as key-value
        pairs. At each epoch (key) the computation time (value) in seconds is the total computation time
        used up to and including that time step. This includes the total time up to and including the current time step,
        since the beginning of the (single-arc) propagation.


        :type: dict[float, float]
     )doc")
                    .def_property_readonly(
                        "cumulative_number_of_function_evaluations_history",
                        &tp::SingleArcSimulationResults<STATE_SCALAR_TYPE,
                                                        TIME_TYPE>::
                            getCumulativeNumberOfFunctionEvaluations,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "total_computation_time",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getTotalComputationRuntime,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "total_number_of_function_evaluations",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getTotalNumberOfFunctionEvaluations,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "termination_details",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getPropagationTerminationReason,
                        R"doc(

        **read-only**

        Object describing the details of the event that triggered the termination of the last propagation.


        :type: PropagationTerminationDetails
     )doc")
                    .def_property_readonly(
                        "integration_completed_successfully",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::integrationCompletedSuccessfully,
                        R"doc(

        **read-only**

        Boolean defining whether the last propagation was finished
        successfully, as defined by the termination conditions, or if
        it was terminated prematurely (for instance due to an
        exception, or an Inf/NaN state entry being detected).


        :type: bool
     )doc")
                    .def_property_readonly(
                        "dependent_variable_ids",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getDependentVariableId,
                        R"doc(

        **read-only**

        Key-value container with the starting entry of the dependent variables saved (key), along with associated ID (value).


        :type: dict[[int,int], str]
     )doc")
                    .def_property_readonly(
                        "ordered_dependent_variable_settings",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getOrderedDependentVariableSettings,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "unordered_dependent_variable_settings",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getOriginalDependentVariableSettings,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "processed_state_ids",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getProcessedStateIds,
                        R"doc(

        **read-only**

        Key-value container with the starting entry of the states (key), along with associated ID (value).


        :type: dict[[int,int] str]
     )doc")
                    .def_property_readonly(
                        "propagated_state_ids",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getPropagatedStateIds,
                        R"doc(

        **read-only**

        Key-value container with the starting entry of the states (key), along with associated ID (value).


        :type: dict[[int,int] str]
     )doc")
                    .def_property_readonly(
                        "initial_and_final_times",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getArcInitialAndFinalTime,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "propagated_state_vector_length",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getPropagatedStateSize,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "propagation_is_performed",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getPropagationIsPerformed,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "solution_is_cleared",
                        &tp::SingleArcSimulationResults<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getSolutionIsCleared,
                        R"doc(No documentation found.)doc");

                py::class_<
                    tp::SingleArcVariationalSimulationResults<STATE_SCALAR_TYPE,
                                                              TIME_TYPE>,
                    std::shared_ptr<tp::SingleArcVariationalSimulationResults<
                        STATE_SCALAR_TYPE, TIME_TYPE>>,
                    tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "SingleArcVariationalSimulationResults",
                    R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "state_transition_matrix_history",
                        &tp::SingleArcVariationalSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getStateTransitionSolution,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "sensitivity_matrix_history",
                        &tp::SingleArcVariationalSimulationResults<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSensitivitySolution,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "dynamics_results",
                        &tp::SingleArcVariationalSimulationResults<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getDynamicsResults,
                        R"doc(No documentation found.)doc");

                py::class_<tp::MultiArcSimulationResults<
                               tp::SingleArcSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>,
                           std::shared_ptr<tp::MultiArcSimulationResults<
                               tp::SingleArcSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>>,
                           tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "MultiArcSimulationResults",
                    R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSingleArcResults,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "arc_start_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getArcStartTimes,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "arc_end_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getArcEndTimes,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "propagation_is_performed",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getPropagationIsPerformed,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "solution_is_cleared",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSolutionIsCleared,
                        R"doc(No documentation found.)doc");

                py::class_<tp::MultiArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>,
                           std::shared_ptr<tp::MultiArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>>,
                           tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "MultiArcVariationalSimulationResults",
                    R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE, TIME_TYPE>::getSingleArcResults,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "arc_start_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE, TIME_TYPE>::getArcStartTimes,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "arc_end_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE, TIME_TYPE>::getArcEndTimes,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "propagation_is_performed",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getPropagationIsPerformed,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "solution_is_cleared",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE, TIME_TYPE>::getSolutionIsCleared,
                        R"doc(No documentation found.)doc");

                py::class_<tp::HybridArcSimulationResults<
                               tp::SingleArcSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>,
                           std::shared_ptr<tp::HybridArcSimulationResults<
                               tp::SingleArcSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>>,
                           tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "HybridArcSimulationResults",
                    R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSingleArcResults,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "multi_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcSimulationResults, STATE_SCALAR_TYPE,
                            TIME_TYPE>::getMultiArcResults,
                        R"doc(No documentation found.)doc");

                py::class_<tp::HybridArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>,
                           std::shared_ptr<tp::HybridArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               STATE_SCALAR_TYPE, TIME_TYPE>>,
                           tp::SimulationResults<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "HybridArcVariationalSimulationResults",
                    R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE, TIME_TYPE>::getSingleArcResults,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "multi_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults,
                            STATE_SCALAR_TYPE, TIME_TYPE>::getMultiArcResults,
                        R"doc(No documentation found.)doc");

                py::class_<tpr::ThrustMagnitudeWrapper,
                           std::shared_ptr<tpr::ThrustMagnitudeWrapper>>(
                    m, "ThrustMagnitudeWrapper");

                py::class_<tpr::ConstantThrustMagnitudeWrapper,
                           std::shared_ptr<tpr::ConstantThrustMagnitudeWrapper>,
                           tpr::ThrustMagnitudeWrapper>(
                    m, "ConstantThrustMagnitudeWrapper")
                    .def_property("constant_thrust_magnitude",
                                  &tpr::ConstantThrustMagnitudeWrapper::
                                      getConstantThrustForceMagnitude,
                                  &tpr::ConstantThrustMagnitudeWrapper::
                                      resetConstantThrustForceMagnitude);


                py::class_<tpr::CustomThrustMagnitudeWrapper,
                           std::shared_ptr<tpr::CustomThrustMagnitudeWrapper>,
                           tpr::ThrustMagnitudeWrapper>(
                    m, "CustomThrustMagnitudeWrapper")
                    .def_property("custom_thrust_magnitude", nullptr,
                                  &tpr::CustomThrustMagnitudeWrapper::
                                      resetThrustMagnitudeFunction);
            }
        }  // namespace propagation
    }      // namespace numerical_simulation
}  // namespace tudatpy
