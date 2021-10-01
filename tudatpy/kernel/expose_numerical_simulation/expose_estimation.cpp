/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation.h"

#include "tudat/astro/propagators/propagateCovariance.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;


namespace tudatpy {
namespace numerical_simulation {
namespace estimation {


void expose_estimation(py::module &m) {

   /*!
    *************** OBSERVATIONS ***************
    */

    py::class_< tom::ObservationCollection<>,
            std::shared_ptr<tom::ObservationCollection<>>>(m, "ObservationCollection");

    py::class_<tom::ObservationViabilityCalculator,
               std::shared_ptr<tom::ObservationViabilityCalculator>>(m, "ObservationViabilityCalculator")
            .def("is_observation_viable", &tom::ObservationViabilityCalculator::isObservationViable,
                 py::arg( "link_end_states" ),
                 py::arg( "link_end_times" ) );

    py::class_<tom::ObservationSimulatorBase<double,double>,
               std::shared_ptr<tom::ObservationSimulatorBase<double,double>>>(m, "ObservationSimulator");

    py::class_<tom::ObservationSimulator<1,double,double>,
               std::shared_ptr<tom::ObservationSimulator<1,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_1");

    py::class_<tom::ObservationSimulator<2,double,double>,
               std::shared_ptr<tom::ObservationSimulator<2,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_2");

    py::class_<tom::ObservationSimulator<3,double,double>,
               std::shared_ptr<tom::ObservationSimulator<3,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_3");

    py::class_<tom::ObservationSimulator<6,double,double>,
               std::shared_ptr<tom::ObservationSimulator<6,double,double>>,
               tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_6");

    /*!
     *************** STATE TRANSITION INTERFACE ***************
     */

    py::class_<
            tp::CombinedStateTransitionAndSensitivityMatrixInterface,
            std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                m, "CombinedStateTransitionAndSensitivityMatrixInterface")
            .def("state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"))
            .def("full_state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getFullCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"))
            .def_property_readonly(
                "state_transition_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getStateTransitionMatrixSize)
            .def_property_readonly(
                "sensitivity_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getSensitivityMatrixSize)
            .def_property_readonly(
                "full_parameter_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getFullParameterVectorSize);;

    /*!
     *************** COVARIANCE ***************
     */
    m.def("propagate_covariance",
          py::overload_cast<
          const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >(
              &tp::propagateCovariance ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times") );

    m.def("propagate_formal_errors",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateFormalErrors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times") );

    m.def("propagate_covariance",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateCovariance ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times") );

    /*!
     *************** PARAMETERS ***************
     */
    py::class_<tep::EstimatableParameterSet<double>,
            std::shared_ptr<tep::EstimatableParameterSet<double>>>(m, "EstimatableParameterSet")
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize )
            .def_property_readonly( "initial_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalSingleArcStateParameterSize )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalMultiArcStateParameterSize )
            .def_property_readonly( "constraints_size",
                                    &tep::EstimatableParameterSet<double>::getConstraintSize )
            .def_property_readonly( "values",
                                    &tep::EstimatableParameterSet<double>::getFullParameterValues< double > )
            .def_property( "parameter_vector",
                           &tep::EstimatableParameterSet<double>::getFullParameterValues< double >,
                           &tep::EstimatableParameterSet<double>::resetParameterValues< double > )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet<double>::getIndicesForParameterType,
                  py::arg("parameter_type") );

    m.def("create_parameters_to_estimate",
          &tss::createParametersToEstimate< double >,
          py::arg("parameter_settings"),
          py::arg("bodies"),
          py::arg("propagator_settings") =
            std::shared_ptr< tp::PropagatorSettings< double > >( ) );

    /*!
     *************** ESTIMATION ***************
     */
    py::class_<
            tss::EstimationConvergenceChecker,
            std::shared_ptr<tss::EstimationConvergenceChecker>>(m, "EstimationConvergenceChecker")
            .def(py::init< const unsigned int,
                 const double,
                 const double,
                 const int >( ),
                 py::arg( "maximum_iterations" ) = 5,
                 py::arg( "minimum_residual_change" ) = 0.0,
                 py::arg( "minimum_residual" ) = 0.0,
                 py::arg( "number_of_iterations_without_improvement" ) = 2 );

    py::class_<
            tss::PodInput<double, double>,
            std::shared_ptr<tss::PodInput<double, double>>>(m, "PodInput")
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection< > >&,
                 const int,
                 const Eigen::MatrixXd,
                 const Eigen::VectorXd >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "parameter_size" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "apriori_parameter_correction" ) = Eigen::VectorXd( 0 ) )
            .def( "set_constant_weight",
                  &tss::PodInput<double, double>::setConstantWeightsMatrix,
                  py::arg( "weight" ) )
            .def( "set_constant_weight_per_observable",
                  &tss::PodInput<double, double>::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ) )
//            .def( "set_constant_weight_per_observable_and_link_end",
//                  &tss::PodInput<double, double>::setConstantPerObservableAndLinkEndsWeights,
//                  py::arg( "weight_per_observable_and_link" ) )
            .def( "define_estimation_settings",
                  &tss::PodInput<double, double>::defineEstimationSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_residuals_and_parameters_per_iteration" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false );

    py::class_<
            tss::PodOutput<double, double>,
            std::shared_ptr<tss::PodOutput<double, double>>>(m, "PodOutput")
            .def_property_readonly("inverse_covariance",
                                   &tss::PodOutput<double, double>::getUnnormalizedInverseCovarianceMatrix)
            .def_property_readonly("covariance",
                                   &tss::PodOutput<double, double>::getUnnormalizedCovarianceMatrix)
            .def_property_readonly("formal_errors",
                                   &tss::PodOutput<double, double>::getFormalErrorVector)
            .def_property_readonly("correlations",
                                   &tss::PodOutput<double, double>::getCorrelationMatrix)
            .def_property_readonly("residual_history",
                                   &tss::PodOutput<double, double>::getResidualHistoryMatrix)
            .def_property_readonly("parameter_history",
                                   &tss::PodOutput<double, double>::getParameterHistoryMatrix)
            .def_property_readonly("design_matrix",
                                   &tss::PodOutput<double, double>::getUnnormalizedInformationMatrix)
            .def_property_readonly("normalized_design_matrix",
                                   &tss::PodOutput<double, double>::getNormalizedInformationMatrix)
            .def_property_readonly("weighted_design_matrix",
                                   &tss::PodOutput<double, double>::getUnnormalizedWeightedInformationMatrix)
            .def_property_readonly("weighted_normalized_design_matrix",
                                   &tss::PodOutput<double, double>::getNormalizedWeightedInformationMatrix);



}

}
}
}// namespace tudatpy
