/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/propagators.h>


#include "expose_propagation.h"


namespace py = pybind11;

namespace ta = tudat::aerodynamics;
namespace tp = tudat::propagators;
namespace tpr = tudat::propulsion;
namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;

namespace tudat
{

namespace aerodynamics
{

class PyAerodynamicGuidance : public ta::AerodynamicGuidance {
public:
    /* Inherit the constructors */
    using AerodynamicGuidance::AerodynamicGuidance;

    using AerodynamicGuidance::currentAngleOfAttack_;
    using AerodynamicGuidance::currentAngleOfSideslip_;
    using AerodynamicGuidance::currentBankAngle_;

    void updateGuidance( const double currentTime ) override {
        PYBIND11_OVERLOAD_PURE(void, AerodynamicGuidance, updateGuidance, currentTime ); }
};

}

}


namespace tudatpy {
namespace numerical_simulation {
namespace propagation {


void expose_propagation(py::module &m)
{



    py::class_<ta::AerodynamicGuidance, ta::PyAerodynamicGuidance,
            std::shared_ptr< ta::AerodynamicGuidance > >(m, "AerodynamicGuidance")
            .def(py::init<>())
            .def("updateGuidance", &ta::AerodynamicGuidance::updateGuidance, py::arg("current_time") )
            .def_readwrite("angle_of_attack", &ta::PyAerodynamicGuidance::currentAngleOfAttack_)
            .def_readwrite("bank_angle", &ta::PyAerodynamicGuidance::currentBankAngle_)
            .def_readwrite("sideslip_angle", &ta::PyAerodynamicGuidance::currentAngleOfSideslip_);



    py::class_<
            tba::TorqueModel,
            std::shared_ptr<tba::TorqueModel> >
            (m, "TorqueModel");


    m.def("get_single_integration_size",
          &tp::getSingleIntegrationSize,
          py::arg("state_type"));

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
              &tp::getInitialStatesOfBodies<TIME_TYPE,double>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"),
          get_docstring("get_state_of_bodies").c_str());


    m.def("get_initial_state_of_bodies",
          py::overload_cast<const std::vector<std::string> &,
          const std::vector<std::string> &,
          const tss::SystemOfBodies &,
          const TIME_TYPE>(
              &tp::getInitialStatesOfBodies<TIME_TYPE,double>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"));

    m.def("get_initial_state_of_body",// overload [2/2]
          py::overload_cast<const std::string&,
          const std::string&,
          const tss::SystemOfBodies&,
          const TIME_TYPE>(
              &tp::getInitialStateOfBody<TIME_TYPE,double>),
          py::arg("body_to_propagate"),
          py::arg("central_body"),
          py::arg("bodies"),
          py::arg("initial_time"));

    m.def("get_initial_rotational_state_of_body",
          py::overload_cast<const std::string&,
          const std::string&,
          const tss::SystemOfBodies&,
          const TIME_TYPE>(
              &tp::getInitialRotationalStateOfBody<TIME_TYPE,double>),
          py::arg("body_to_propagate"),
          py::arg("base_orientation"),
          py::arg("bodies"),
          py::arg("initial_time"));
    
    py::class_<
            tp::DampedInitialRotationalStateResults<TIME_TYPE, double>,
            std::shared_ptr<tp::DampedInitialRotationalStateResults<TIME_TYPE, double>>>(m, "RotationalProperModeDampingResults",
                                                                 get_docstring("RotationalProperModeDampingResults").c_str())
            .def_readwrite("damped_initial_state", &tp::DampedInitialRotationalStateResults<TIME_TYPE,double>::initialState_,
                           get_docstring("RotationalProperModeDampingResults.damped_initial_state").c_str())
            .def_readwrite("forward_backward_states", &tp::DampedInitialRotationalStateResults<TIME_TYPE,double>::forwardBackwardPropagatedStates_,
                           get_docstring("RotationalProperModeDampingResults.forward_backward_states").c_str())
            .def_readwrite("forward_backward_dependent_variables", &tp::DampedInitialRotationalStateResults<TIME_TYPE,double>::forwardBackwardDependentVariables_ ,
                           get_docstring("RotationalProperModeDampingResults.forward_backward_dependent_variables").c_str());

    m.def("get_damped_proper_mode_initial_rotational_state",
          py::overload_cast<
          const tss::SystemOfBodies&,
          const std::shared_ptr< tp::SingleArcPropagatorSettings< double, TIME_TYPE > >,
          const double,
          const std::vector< double >,
          const bool >( &tp::getZeroProperModeRotationalStateWithStruct< TIME_TYPE, double > ),
          py::arg("bodies"),
          py::arg("propagator_settings"),
          py::arg("body_mean_rotational_rate"),
          py::arg("dissipation_times"),
          py::arg("propagate_undamped") = true ,
          get_docstring("get_damped_proper_mode_initial_rotational_state").c_str());

    m.def("combine_initial_states",
          &tp::createCombinedInitialState<double,TIME_TYPE>,
          py::arg("propagator_settings_per_type"),
          get_docstring("combine_initial_states").c_str());

    py::class_<
            tba::AccelerationModel<Eigen::Vector3d>,
            std::shared_ptr<tba::AccelerationModel<Eigen::Vector3d>>>(m, "AccelerationModel");

    py::class_<
            tba::MassRateModel,
            std::shared_ptr<tba::MassRateModel>>(m, "MassRateModel");



    py::enum_<tp::PropagationTerminationReason>(m, "PropagationTerminationReason",
                                                get_docstring("PropagationTerminationReason").c_str())
            .value("propagation_never_run",
                   tp::PropagationTerminationReason::propagation_never_run)
            .value("unknown_reason",
                   tp::PropagationTerminationReason::unknown_propagation_termination_reason)
            .value("termination_condition_reached",
                   tp::PropagationTerminationReason::termination_condition_reached)
            .value("runtime_error_caught_in_propagation",
                   tp::PropagationTerminationReason::runtime_error_caught_in_propagation)
            .value("nan_or_inf_detected_in_state",
                   tp::PropagationTerminationReason::nan_or_inf_detected_in_state)
            .export_values();

    py::class_<
            tp::PropagationTerminationDetails,
            std::shared_ptr<tp::PropagationTerminationDetails>>(
                    m,
                    "PropagationTerminationDetails",
                    get_docstring("PropagationTerminationDetails").c_str())
            .def_property_readonly(
                    "termination_reason",
                    &tp::PropagationTerminationDetails::getPropagationTerminationReason,
                    get_docstring("PropagationTerminationDetails.termination_reason").c_str())
            .def_property_readonly("terminated_on_exact_condition",
                    &tp::PropagationTerminationDetails::getTerminationOnExactCondition,
                    get_docstring("PropagationTerminationDetails.terminated_on_exact_condition").c_str());

    py::class_<
            tp::PropagationTerminationDetailsFromHybridCondition,
            std::shared_ptr<tp::PropagationTerminationDetailsFromHybridCondition>,
            tp::PropagationTerminationDetails>(
                    m,
                    "PropagationTerminationDetailsFromHybridCondition",
                    get_docstring("PropagationTerminationDetailsFromHybridCondition").c_str())
            .def_property_readonly(
                    "was_condition_met_when_stopping",
                    &tp::PropagationTerminationDetailsFromHybridCondition::getWasConditionMetWhenStopping,
                    get_docstring("PropagationTerminationDetailsFromHybridCondition.was_condition_met_when_stopping").c_str());

    py::class_<
            tp::DependentVariablesInterface<TIME_TYPE>,
            std::shared_ptr<tp::DependentVariablesInterface<TIME_TYPE>>>(m, "DependentVariablesInterface",
                                                                         get_docstring("DependentVariablesInterface").c_str());

    py::class_<
            tp::SimulationResults<double, TIME_TYPE>,
            std::shared_ptr<tp::SimulationResults<double, TIME_TYPE>>>(m, "SimulationResults",
                                                                       get_docstring("SimulationResults").c_str())
            .def_property_readonly("dependent_variable_interface",
                                   &tp::SimulationResults<double, TIME_TYPE>::getDependentVariablesInterface,
                                   get_docstring("SimulationResults.dependent_variable_interface").c_str() );

    py::class_<
            tp::SingleArcSimulationResults<double, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcSimulationResults<double, TIME_TYPE>>,
            tp::SimulationResults<double, TIME_TYPE> >(m, "SingleArcSimulationResults"
                                                          "cd",
                                                       get_docstring("SingleArcSimulationResults").c_str())
            .def_property_readonly("state_history",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getEquationsOfMotionNumericalSolution,
                                   get_docstring("SingleArcSimulationResults.state_history").c_str() )
            .def_property_readonly("state_history_float",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getEquationsOfMotionNumericalSolutionDouble)
            .def_property_readonly("state_history_float_split",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getEquationsOfMotionNumericalSolutionDoubleSplit)
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getEquationsOfMotionNumericalSolutionRaw,
                                   get_docstring("SingleArcSimulationResults.unprocessed_state_history").c_str() )
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getDependentVariableHistory,
                                   get_docstring("SingleArcSimulationResults.dependent_variable_history").c_str() )
            .def_property_readonly("cumulative_computation_time_history",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getCumulativeComputationTimeHistory,
                                   get_docstring("SingleArcSimulationResults.cumulative_computation_time_history").c_str() )
            .def_property_readonly("cumulative_computation_time_history",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getCumulativeComputationTimeHistory,
                                   get_docstring("SingleArcSimulationResults.cumulative_computation_time_history").c_str() )
            .def_property_readonly("cumulative_number_of_function_evaluations_history",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getCumulativeNumberOfFunctionEvaluations,
                                   get_docstring("SingleArcSimulationResults.cumulative_number_of_function_evaluations_history").c_str() )
            .def_property_readonly("total_computation_time",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getTotalComputationRuntime,
                                   get_docstring("SingleArcSimulationResults.total_computation_time").c_str() )
            .def_property_readonly("total_number_of_function_evaluations",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getTotalNumberOfFunctionEvaluations,
                                   get_docstring("SingleArcSimulationResults.total_number_of_function_evaluations").c_str() )
            .def_property_readonly("termination_details",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getPropagationTerminationReason,
                                   get_docstring("SingleArcSimulationResults.termination_details").c_str() )
            .def_property_readonly("integration_completed_successfully",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::integrationCompletedSuccessfully,
                                   get_docstring("SingleArcSimulationResults.integration_completed_successfully").c_str() )
            .def_property_readonly("dependent_variable_ids",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getDependentVariableId,
                                   get_docstring("SingleArcSimulationResults.dependent_variable_ids").c_str() )
            .def_property_readonly("ordered_dependent_variable_settings",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getOrderedDependentVariableSettings,
                                   get_docstring("SingleArcSimulationResults.ordered_dependent_variable_settings").c_str() )
            .def_property_readonly("unordered_dependent_variable_settings",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getOriginalDependentVariableSettings,
                                   get_docstring("SingleArcSimulationResults.unordered_dependent_variable_settings").c_str() )
            .def_property_readonly("processed_state_ids",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getProcessedStateIds,
                                   get_docstring("SingleArcSimulationResults.state_ids").c_str() )
            .def_property_readonly("propagated_state_ids",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getPropagatedStateIds,
                                   get_docstring("SingleArcSimulationResults.state_ids").c_str() )
            .def_property_readonly("initial_and_final_times",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getArcInitialAndFinalTime,
                                   get_docstring("SingleArcSimulationResults.initial_and_final_times").c_str() )
            .def_property_readonly("propagated_state_vector_length",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getPropagatedStateSize,
                                   get_docstring("SingleArcSimulationResults.propagated_state_vector_length").c_str() )
            .def_property_readonly("propagation_is_performed",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getPropagationIsPerformed,
                                   get_docstring("SingleArcSimulationResults.propagation_is_performed").c_str() )
            .def_property_readonly("solution_is_cleared",
                                   &tp::SingleArcSimulationResults<double, TIME_TYPE>::getSolutionIsCleared,
                                   get_docstring("SingleArcSimulationResults.solution_is_cleared").c_str() );

    py::class_<
            tp::SingleArcVariationalSimulationResults<double, TIME_TYPE>,
            std::shared_ptr<tp::SingleArcVariationalSimulationResults<double, TIME_TYPE>>,
            tp::SimulationResults<double, TIME_TYPE> >(m, "SingleArcVariationalSimulationResults",
                                                       get_docstring("SingleArcVariationalSimulationResults").c_str())
            .def_property_readonly("state_transition_matrix_history",
                                   &tp::SingleArcVariationalSimulationResults<double, TIME_TYPE>::getStateTransitionSolution,
                                   get_docstring("SingleArcVariationalSimulationResults.state_transition_matrix_history").c_str() )
            .def_property_readonly("sensitivity_matrix_history",
                                   &tp::SingleArcVariationalSimulationResults<double, TIME_TYPE>::getSensitivitySolution,
                                   get_docstring("SingleArcVariationalSimulationResults.sensitivity_matrix_history").c_str() )
            .def_property_readonly("dynamics_results",
                                   &tp::SingleArcVariationalSimulationResults<double, TIME_TYPE>::getDynamicsResults,
                                   get_docstring("SingleArcVariationalSimulationResults.dynamics_results").c_str() );

    py::class_<
            tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>,
            std::shared_ptr<tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>>,
            tp::SimulationResults<double, TIME_TYPE> >(m, "MultiArcSimulationResults",
                                                       get_docstring("MultiArcSimulationResults").c_str())
            .def_property_readonly("single_arc_results",
                                   &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getSingleArcResults,
                                   get_docstring("MultiArcSimulationResults.single_arc_results").c_str() )
            .def_property_readonly("arc_start_times",
                                   &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getArcStartTimes,
                                   get_docstring("MultiArcSimulationResults.arc_start_times").c_str() )
            .def_property_readonly("arc_end_times",
                                   &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getArcEndTimes,
                                   get_docstring("MultiArcSimulationResults.arc_end_times").c_str() )
            .def_property_readonly("propagation_is_performed",
                                   &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getPropagationIsPerformed,
                                   get_docstring("MultiArcSimulationResults.propagation_is_performed").c_str() )
            .def_property_readonly("solution_is_cleared",
                                   &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getSolutionIsCleared,
                                   get_docstring("MultiArcSimulationResults.solution_is_cleared").c_str() );

    py::class_<
            tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>,
            std::shared_ptr<tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>>,
            tp::SimulationResults<double, TIME_TYPE> >(m, "MultiArcVariationalSimulationResults",
                                                       get_docstring("MultiArcVariationalSimulationResults").c_str())
            .def_property_readonly("single_arc_results",
                                   &tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getSingleArcResults,
                                   get_docstring("MultiArcVariationalSimulationResults.single_arc_results").c_str() )
            .def_property_readonly("arc_start_times",
                                   &tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getArcStartTimes,
                                   get_docstring("MultiArcVariationalSimulationResults.arc_start_times").c_str() )
            .def_property_readonly("arc_end_times",
                                   &tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getArcEndTimes,
                                   get_docstring("MultiArcVariationalSimulationResults.arc_end_times").c_str() )
            .def_property_readonly("propagation_is_performed",
                                   &tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getPropagationIsPerformed,
                                   get_docstring("MultiArcVariationalSimulationResults.propagation_is_performed").c_str() )
            .def_property_readonly("solution_is_cleared",
                                   &tp::MultiArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getSolutionIsCleared,
                                   get_docstring("MultiArcVariationalSimulationResults.solution_is_cleared").c_str() );

    py::class_<
            tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>,
            std::shared_ptr<tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>>,
            tp::SimulationResults<double, TIME_TYPE> >(m, "HybridArcSimulationResults",
                                                       get_docstring("HybridArcSimulationResults").c_str())
            .def_property_readonly("single_arc_results",
                                   &tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getSingleArcResults,
                                   get_docstring("HybridArcSimulationResults.single_arc_results").c_str() )
            .def_property_readonly("multi_arc_results",
                                   &tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, TIME_TYPE>::getMultiArcResults,
                                   get_docstring("HybridArcSimulationResults.arc_start_times").c_str() );

        py::class_<
            tp::HybridArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>,
            std::shared_ptr<tp::HybridArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>>,
            tp::SimulationResults<double, TIME_TYPE> >(m, "HybridArcVariationalSimulationResults",
                                                       get_docstring("HybridArcVariationalSimulationResults").c_str())
            .def_property_readonly("single_arc_results",
                                   &tp::HybridArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getSingleArcResults,
                                   get_docstring("HybridArcVariationalSimulationResults.single_arc_results").c_str() )
            .def_property_readonly("multi_arc_results",
                                   &tp::HybridArcSimulationResults<tp::SingleArcVariationalSimulationResults, double, TIME_TYPE>::getMultiArcResults,
                                   get_docstring("HybridArcVariationalSimulationResults.arc_start_times").c_str() );

    py::class_<
            tpr::ThrustMagnitudeWrapper,
            std::shared_ptr< tpr::ThrustMagnitudeWrapper > >(m, "ThrustMagnitudeWrapper" );

    py::class_<
            tpr::ConstantThrustMagnitudeWrapper,
            std::shared_ptr< tpr::ConstantThrustMagnitudeWrapper >,
            tpr::ThrustMagnitudeWrapper >(m, "ConstantThrustMagnitudeWrapper" )
            .def_property("constant_thrust_magnitude", &tpr::ConstantThrustMagnitudeWrapper::getConstantThrustForceMagnitude,
                          &tpr::ConstantThrustMagnitudeWrapper::resetConstantThrustForceMagnitude );


    py::class_<
        tpr::CustomThrustMagnitudeWrapper,
        std::shared_ptr< tpr::CustomThrustMagnitudeWrapper >,
        tpr::ThrustMagnitudeWrapper >(m, "CustomThrustMagnitudeWrapper" )
        .def_property("custom_thrust_magnitude", nullptr,
                      &tpr::CustomThrustMagnitudeWrapper::resetThrustMagnitudeFunction );
}
}// namespace propagation
}// namespace numerical_simulation
}// namespace tudatpy
