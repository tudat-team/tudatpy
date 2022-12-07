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

#include <tudat/astro/aerodynamics/aerodynamicGuidance.h>
#include <tudat/astro/basic_astro.h>
#include <tudat/astro/propagators.h>


#include "expose_propagation.h"


namespace py = pybind11;

namespace ta = tudat::aerodynamics;
namespace tp = tudat::propagators;
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


void expose_propagation(py::module &m) {



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
          const double>(
              &tp::getInitialStatesOfBodies<>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"),
          get_docstring("get_state_of_bodies").c_str());


    m.def("get_initial_state_of_bodies",
          py::overload_cast<const std::vector<std::string> &,
          const std::vector<std::string> &,
          const tss::SystemOfBodies &,
          const double>(
              &tp::getInitialStatesOfBodies<>),
          py::arg("bodies_to_propagate"),
          py::arg("central_bodies"),
          py::arg("body_system"),
          py::arg("initial_time"));

    m.def("get_initial_state_of_body",// overload [2/2]
          py::overload_cast<const std::string&,
          const std::string&,
          const tss::SystemOfBodies&,
          const double>(
              &tp::getInitialStateOfBody<>),
          py::arg("body_to_propagate"),
          py::arg("central_body"),
          py::arg("bodies"),
          py::arg("initial_time"));

    m.def("get_initial_rotational_state_of_body",
          py::overload_cast<const std::string&,
          const std::string&,
          const tss::SystemOfBodies&,
          const double>(
              &tp::getInitialRotationalStateOfBody<>),
          py::arg("body_to_propagate"),
          py::arg("base_orientation"),
          py::arg("bodies"),
          py::arg("initial_time"));

    m.def("get_zero_proper_mode_rotational_state",
          py::overload_cast<
          const tss::SystemOfBodies&,
          const std::shared_ptr< tni::IntegratorSettings< double > >,
          const std::shared_ptr< tp::SingleArcPropagatorSettings< double > >,
          const double,
          const std::vector< double >,
          const bool >( &tp::getZeroProperModeRotationalState< > ),
          py::arg("bodies"),
          py::arg("integrator_settings"),
          py::arg("propagator_settings"),
          py::arg("body_mean_rotational_rate"),
          py::arg("dissipation_times"),
          py::arg("propagate_undamped") = true );

    m.def("combine_initial_states",
          &tp::createCombinedInitialState<double>,
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
                    get_docstring("PropagationTerminationDetailsFromHybridCondition.termination_reason").c_str());

    py::class_<
            tp::SingleArcSimulationResults<double, double>,
            std::shared_ptr<tp::SingleArcSimulationResults<double, double>>>(m, "SingleArcSimulationResults",
                                                                             get_docstring("SingleArcSimulationResults").c_str())
            .def_property_readonly("state_history",
                                   &tp::SingleArcSimulationResults<double, double>::getEquationsOfMotionNumericalSolution,
                                   get_docstring("SingleArcSimulationResults.state_history").c_str() )
            .def_property_readonly("unprocessed_state_history",
                                   &tp::SingleArcSimulationResults<double, double>::getEquationsOfMotionNumericalSolutionRaw,
                                   get_docstring("SingleArcSimulationResults.unprocessed_state_history").c_str() )
            .def_property_readonly("dependent_variable_history",
                                   &tp::SingleArcSimulationResults<double, double>::getDependentVariableHistory,
                                   get_docstring("SingleArcSimulationResults.dependent_variable_history").c_str() )
            .def_property_readonly("cumulative_computation_time_history",
                                   &tp::SingleArcSimulationResults<double, double>::getCumulativeComputationTimeHistory,
                                   get_docstring("SingleArcSimulationResults.cumulative_computation_time_history").c_str() )
            .def_property_readonly("cumulative_number_of_function_evaluations",
                                   &tp::SingleArcSimulationResults<double, double>::getCumulativeNumberOfFunctionEvaluations,
                                   get_docstring("SingleArcSimulationResults.cumulative_number_of_function_evaluations").c_str() )
            .def_property_readonly("termination_details",
                                   &tp::SingleArcSimulationResults<double, double>::getPropagationTerminationReason,
                                   get_docstring("SingleArcSimulationResults.termination_details").c_str() )
            .def_property_readonly("integration_completed_successfully",
                                   &tp::SingleArcSimulationResults<double, double>::integrationCompletedSuccessfully,
                                   get_docstring("SingleArcSimulationResults.integration_completed_successfully").c_str() )
            .def_property_readonly("dependent_variable_ids",
                                   &tp::SingleArcSimulationResults<double, double>::getDependentVariableId,
                                   get_docstring("SingleArcSimulationResults.dependent_variable_ids").c_str() )
            .def_property_readonly("state_ids",
                                   &tp::SingleArcSimulationResults<double, double>::getStateIds,
                                   get_docstring("SingleArcSimulationResults.state_ids").c_str() );


}
}// namespace propagation
}// namespace numerical_simulation
}// namespace tudatpy
