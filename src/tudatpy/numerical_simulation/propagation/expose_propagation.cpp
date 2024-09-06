/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/propagators.h>

#include "tudatpy/scalarTypes.h"

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


            PYBIND11_MODULE(expose_propagation, m) {
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
                          &tp::getInitialStatesOfBodies<TIME_TYPE, double>),
                      py::arg("bodies_to_propagate"), py::arg("central_bodies"),
                      py::arg("body_system"), py::arg("initial_time"),
R"doc(Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.

	Function to get the states of a set of bodies, with respect to some set of central bodies, at the requested time.


	:param bodies_to_propagate:
		List of bodies to be propagated.
	:param central_bodies:
		List of central bodies, each referred to a body being propagated (in the same order).
	:param bodies_to_propagate:
		System of bodies used in the propagation.
	:param initial_time:
		Initial time of the propagation.
	:return:
		Time at which the states should be retrieved.
)doc");


                m.def("get_initial_state_of_bodies",
                      py::overload_cast<const std::vector<std::string> &,
                                        const std::vector<std::string> &,
                                        const tss::SystemOfBodies &,
                                        const TIME_TYPE>(
                          &tp::getInitialStatesOfBodies<TIME_TYPE, double>),
                      py::arg("bodies_to_propagate"), py::arg("central_bodies"),
                      py::arg("body_system"), py::arg("initial_time"));

                m.def(
                    "get_initial_state_of_body",  // overload [2/2]
                    py::overload_cast<const std::string &, const std::string &,
                                      const tss::SystemOfBodies &,
                                      const TIME_TYPE>(
                        &tp::getInitialStateOfBody<TIME_TYPE, double>),
                    py::arg("body_to_propagate"), py::arg("central_body"),
                    py::arg("bodies"), py::arg("initial_time"));

                m.def(
                    "get_initial_rotational_state_of_body",
                    py::overload_cast<const std::string &, const std::string &,
                                      const tss::SystemOfBodies &,
                                      const TIME_TYPE>(
                        &tp::getInitialRotationalStateOfBody<TIME_TYPE,
                                                             double>),
                    py::arg("body_to_propagate"), py::arg("base_orientation"),
                    py::arg("bodies"), py::arg("initial_time"));

                py::class_<
                    tp::DampedInitialRotationalStateResults<TIME_TYPE, double>,
                    std::shared_ptr<tp::DampedInitialRotationalStateResults<
                        TIME_TYPE, double>>>(
                    m, "RotationalProperModeDampingResults",
"")
                    .def_readwrite("damped_initial_state",
                                   &tp::DampedInitialRotationalStateResults<
                                       TIME_TYPE, double>::initialState_,
"")
                    .def_readwrite(
                        "forward_backward_states",
                        &tp::DampedInitialRotationalStateResults<
                            TIME_TYPE,
                            double>::forwardBackwardPropagatedStates_,
"")
                    .def_readwrite(
                        "forward_backward_dependent_variables",
                        &tp::DampedInitialRotationalStateResults<
                            TIME_TYPE,
                            double>::forwardBackwardDependentVariables_,
"");

                m.def(
                    "get_damped_proper_mode_initial_rotational_state",
                    py::overload_cast<
                        const tss::SystemOfBodies &,
                        const std::shared_ptr<
                            tp::SingleArcPropagatorSettings<double, TIME_TYPE>>,
                        const double, const std::vector<double>, const bool>(
                        &tp::getZeroProperModeRotationalStateWithStruct<
                            TIME_TYPE, double>),
                    py::arg("bodies"), py::arg("propagator_settings"),
                    py::arg("body_mean_rotational_rate"),
                    py::arg("dissipation_times"),
                    py::arg("propagate_undamped") = true,
"");

                m.def("combine_initial_states",
                      &tp::createCombinedInitialState<double, TIME_TYPE>,
                      py::arg("propagator_settings_per_type"),
R"doc(Function to retrieve the initial state for a list of propagator settings.

	Function to retrieve the initial state for a list of propagator settings. This way, the initial state for
	different quantities to be propagated (e.g., translational state, rotational state, mass) are retrieved and
	organized in a single container.


	:param propagator_settings_per_type:
		Propagator settings where the type of propagation is reported as key and the respective list of propagator settings as value.
	:return:
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
R"doc(Enumeration of types of termination of propagation.


	:member propagation_never_run:
	:member unknown_reason:
	:member termination_condition_reached:
	:member runtime_error_caught_in_propagation:
	:member nan_or_inf_detected_in_state:
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
R"doc(Object that provides information on the reason for the
termination of the propagation.


)doc")
                    .def_property_readonly(
                        "termination_reason",
                        &tp::PropagationTerminationDetails::
                            getPropagationTerminationReason,
R"doc(Enum defining the reason the propagation was terminated

	)doc")
                    .def_property_readonly(
                        "terminated_on_exact_condition",
                        &tp::PropagationTerminationDetails::
                            getTerminationOnExactCondition,
R"doc(Boolean defining whether the propagation was terminated on an *exact* final condition,
or once the propagation went *past* the determined final condition. The choice of behaviour is
defined by the termination settings provided as input to the Simulator object. This variable only
has a meaningful definition if the ``termination_reason`` has value ``termination_condition_reached``

	)doc");

                py::class_<
                    tp::PropagationTerminationDetailsFromHybridCondition,
                    std::shared_ptr<
                        tp::PropagationTerminationDetailsFromHybridCondition>,
                    tp::PropagationTerminationDetails>(
                    m, "PropagationTerminationDetailsFromHybridCondition",
"")
                    .def_property_readonly(
                        "was_condition_met_when_stopping",
                        &tp::PropagationTerminationDetailsFromHybridCondition::
                            getWasConditionMetWhenStopping,
"");

                py::class_<tp::DependentVariablesInterface<TIME_TYPE>,
                           std::shared_ptr<
                               tp::DependentVariablesInterface<TIME_TYPE>>>(
                    m, "DependentVariablesInterface",
"");

                py::class_<
                    tp::SimulationResults<double, TIME_TYPE>,
                    std::shared_ptr<tp::SimulationResults<double, TIME_TYPE>>>(
                    m, "SimulationResults",
R"doc(	Base class for objects that store all results of a numerical propagation. Derived class are implemented for single-, multi- and hybrid-arc propagation of botj dynamics and variational equations

)doc")
                    .def_property_readonly(
                        "dependent_variable_interface",
                        &tp::SimulationResults<
                            double, TIME_TYPE>::getDependentVariablesInterface,
"");

                py::class_<tp::SingleArcSimulationResults<double, TIME_TYPE>,
                           std::shared_ptr<tp::SingleArcSimulationResults<
                               double, TIME_TYPE>>,
                           tp::SimulationResults<double, TIME_TYPE>>(
                    m, "SingleArcSimulationResults",
R"doc(Class that stores all the results (including logging data) of a single-arc propagation


)doc")
                    .def_property_readonly(
                        "state_history",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getEquationsOfMotionNumericalSolution,
R"doc(Numerical solution of the equations of motion as key-value pairs. The key denotes the epoch. The value contains the
numerically calculated state at this epoch. For this function, the states are always converted to so-called
'conventional' formulations (e.g. Cartesian states for translational dynamics), see `here <https://tudat-space.readthedocs.io/en/latest/_src_api/propagation_setup/settings/conventional_vs_propagated_coordinates.html>`_
for details. For the history of the states that were actually propagated, use the ``unprocessed_state_history``.

.. note:: The propagated state at each epoch contains the state types in the following order: Translational ( **C** ), Rotational ( **R** ), Mass ( **M** ), and Custom ( **C** ).
          When propagating two bodies, an example of what the output state would look like is for instance:
          [ **T** Body 1, **T** Body 2, **R** Body 1, **R** Body 2, **M** Body 1, **M** Body 2 ] The specifics can be retrieved using the :attr:`state_ids` attribute of this class

.. note:: For propagation of translational dynamics using cowell
          propagator, the conventional and propagated
          coordinates are identical.

	)doc")
                    .def_property_readonly(
                        "unprocessed_state_history",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getEquationsOfMotionNumericalSolutionRaw,
R"doc(Numerical solution of the equations of motion as key-value pairs, without any processing applied. The key denotes the epoch. The value contains the
numerically calculated state at this epoch. This attribute contains the states of the propagated bodies expressed in the
"raw" form in which the propagation took place. For instance, when using a Gauss-Kepler propagation scheme, this
attribute will contain the numerically propagated Keplerian elements at each time epoch

	)doc")
                    .def_property_readonly(
                        "dependent_variable_history",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getDependentVariableHistory,
R"doc(Dependent variables computed during the propagation as key-value pairs.
They are returned in the order with the same order of the DependentVariableSaveSettings object as values,
as value, with the epoch as key.

	)doc")
                    .def_property_readonly(
                        "cumulative_computation_time_history",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getCumulativeComputationTimeHistory,
R"doc(History of cumulative computation time in seconds needed during the propagation as key-value
pairs. At each epoch (key) the computation time (value) in seconds is the total computation time
used up to and including that time step. This includes the total time up to and including the current time step,
since the beginning of the (single-arc) propagation.

	)doc")
                    .def_property_readonly(
                        "cumulative_computation_time_history",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getCumulativeComputationTimeHistory,
R"doc(History of cumulative computation time in seconds needed during the propagation as key-value
pairs. At each epoch (key) the computation time (value) in seconds is the total computation time
used up to and including that time step. This includes the total time up to and including the current time step,
since the beginning of the (single-arc) propagation.

	)doc")
                    .def_property_readonly(
                        "cumulative_number_of_function_evaluations_history",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getCumulativeNumberOfFunctionEvaluations,
"")
                    .def_property_readonly(
                        "total_computation_time",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getTotalComputationRuntime,
"")
                    .def_property_readonly(
                        "total_number_of_function_evaluations",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getTotalNumberOfFunctionEvaluations,
"")
                    .def_property_readonly(
                        "termination_details",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getPropagationTerminationReason,
R"doc(Object describing the details of the event that triggered the termination of the last propagation.

	)doc")
                    .def_property_readonly(
                        "integration_completed_successfully",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            integrationCompletedSuccessfully,
R"doc(Boolean defining whether the last propagation was finished
successfully, as defined by the termination conditions, or if
it was terminated prematurely (for instance due to an
exception, or an Inf/NaN state entry being detected).

	)doc")
                    .def_property_readonly(
                        "dependent_variable_ids",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getDependentVariableId,
R"doc(Key-value container with the starting entry of the dependent variables saved (key), along with associated ID (value).

	)doc")
                    .def_property_readonly(
                        "ordered_dependent_variable_settings",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getOrderedDependentVariableSettings,
"")
                    .def_property_readonly(
                        "unordered_dependent_variable_settings",
                        &tp::SingleArcSimulationResults<double, TIME_TYPE>::
                            getOriginalDependentVariableSettings,
"")
                    .def_property_readonly(
                        "processed_state_ids",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getProcessedStateIds,
R"doc(Key-value container with the starting entry of the states (key), along with associated ID (value).

	)doc")
                    .def_property_readonly(
                        "propagated_state_ids",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getPropagatedStateIds,
R"doc(Key-value container with the starting entry of the states (key), along with associated ID (value).

	)doc")
                    .def_property_readonly(
                        "initial_and_final_times",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getArcInitialAndFinalTime,
"")
                    .def_property_readonly(
                        "propagated_state_vector_length",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getPropagatedStateSize,
"")
                    .def_property_readonly(
                        "propagation_is_performed",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getPropagationIsPerformed,
"")
                    .def_property_readonly(
                        "solution_is_cleared",
                        &tp::SingleArcSimulationResults<
                            double, TIME_TYPE>::getSolutionIsCleared,
"");

                py::class_<
                    tp::SingleArcVariationalSimulationResults<double,
                                                              TIME_TYPE>,
                    std::shared_ptr<tp::SingleArcVariationalSimulationResults<
                        double, TIME_TYPE>>,
                    tp::SimulationResults<double, TIME_TYPE>>(
                    m, "SingleArcVariationalSimulationResults",
"")
                    .def_property_readonly(
                        "state_transition_matrix_history",
                        &tp::SingleArcVariationalSimulationResults<
                            double, TIME_TYPE>::getStateTransitionSolution,
"")
                    .def_property_readonly(
                        "sensitivity_matrix_history",
                        &tp::SingleArcVariationalSimulationResults<
                            double, TIME_TYPE>::getSensitivitySolution,
"")
                    .def_property_readonly(
                        "dynamics_results",
                        &tp::SingleArcVariationalSimulationResults<
                            double, TIME_TYPE>::getDynamicsResults,
"");

                py::class_<
                    tp::MultiArcSimulationResults<
                        tp::SingleArcSimulationResults, double, TIME_TYPE>,
                    std::shared_ptr<tp::MultiArcSimulationResults<
                        tp::SingleArcSimulationResults, double, TIME_TYPE>>,
                    tp::SimulationResults<double, TIME_TYPE>>(
                    m, "MultiArcSimulationResults",
"")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getSingleArcResults,
"")
                    .def_property_readonly(
                        "arc_start_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getArcStartTimes,
"")
                    .def_property_readonly(
                        "arc_end_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getArcEndTimes,
"")
                    .def_property_readonly(
                        "propagation_is_performed",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getPropagationIsPerformed,
"")
                    .def_property_readonly(
                        "solution_is_cleared",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getSolutionIsCleared,
"");

                py::class_<tp::MultiArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               double, TIME_TYPE>,
                           std::shared_ptr<tp::MultiArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               double, TIME_TYPE>>,
                           tp::SimulationResults<double, TIME_TYPE>>(
                    m, "MultiArcVariationalSimulationResults",
"")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getSingleArcResults,
"")
                    .def_property_readonly(
                        "arc_start_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getArcStartTimes,
"")
                    .def_property_readonly(
                        "arc_end_times",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getArcEndTimes,
"")
                    .def_property_readonly(
                        "propagation_is_performed",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getPropagationIsPerformed,
"")
                    .def_property_readonly(
                        "solution_is_cleared",
                        &tp::MultiArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getSolutionIsCleared,
"");

                py::class_<
                    tp::HybridArcSimulationResults<
                        tp::SingleArcSimulationResults, double, TIME_TYPE>,
                    std::shared_ptr<tp::HybridArcSimulationResults<
                        tp::SingleArcSimulationResults, double, TIME_TYPE>>,
                    tp::SimulationResults<double, TIME_TYPE>>(
                    m, "HybridArcSimulationResults",
"")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getSingleArcResults,
"")
                    .def_property_readonly(
                        "multi_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcSimulationResults, double,
                            TIME_TYPE>::getMultiArcResults,
"");

                py::class_<tp::HybridArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               double, TIME_TYPE>,
                           std::shared_ptr<tp::HybridArcSimulationResults<
                               tp::SingleArcVariationalSimulationResults,
                               double, TIME_TYPE>>,
                           tp::SimulationResults<double, TIME_TYPE>>(
                    m, "HybridArcVariationalSimulationResults",
"")
                    .def_property_readonly(
                        "single_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getSingleArcResults,
"")
                    .def_property_readonly(
                        "multi_arc_results",
                        &tp::HybridArcSimulationResults<
                            tp::SingleArcVariationalSimulationResults, double,
                            TIME_TYPE>::getMultiArcResults,
"");

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
    }  // namespace numerical_simulation
}  // namespace tudatpy
