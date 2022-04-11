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

#include "tudatpy/docstrings.h"

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
     *************** PARAMETERS ***************
     */


    py::class_<tep::EstimatableParameterSet<double>,
            std::shared_ptr<tep::EstimatableParameterSet<double>>>(m, "EstimatableParameterSet",
                                                                   get_docstring("EstimatableParameterSet").c_str() )
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize,
                                    get_docstring("EstimatableParameterSet.parameter_set_size").c_str() )
            .def_property_readonly( "initial_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_states_size").c_str() )
            .def_property_readonly( "initial_single_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalSingleArcStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_single_arc_states_size").c_str() )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalMultiArcStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_multi_arc_states_size").c_str() )
            .def_property_readonly( "constraints_size",
                                    &tep::EstimatableParameterSet<double>::getConstraintSize,
                                    get_docstring("EstimatableParameterSet.constraints_size").c_str() )
            .def_property( "parameter_vector",
                           &tep::EstimatableParameterSet<double>::getFullParameterValues< double >,
                           &tep::EstimatableParameterSet<double>::resetParameterValues< double >,
                           get_docstring("EstimatableParameterSet.parameter_vector").c_str() )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet<double>::getIndicesForParameterType,
                  py::arg("parameter_type"),
                  get_docstring("EstimatableParameterSet.indices_for_parameter_type").c_str() );

    /*!
     *************** OBSERVATIONS ***************
     */

    py::class_<tom::ObservationViabilityCalculator,
            std::shared_ptr<tom::ObservationViabilityCalculator>>(m, "ObservationViabilityCalculator",
                                                                  get_docstring("ObservationViabilityCalculator").c_str() )
            .def("is_observation_viable", &tom::ObservationViabilityCalculator::isObservationViable,
                 py::arg( "link_end_states" ),
                 py::arg( "link_end_times" ),
                 get_docstring("ObservationViabilityCalculator.is_observation_viable").c_str() );

    py::class_<tom::ObservationSimulatorBase<double,double>,
            std::shared_ptr<tom::ObservationSimulatorBase<double,double>>>(m, "ObservationSimulator",
                                                                           get_docstring("ObservationSimulator").c_str() );

    py::class_<tom::ObservationSimulator<1,double,double>,
            std::shared_ptr<tom::ObservationSimulator<1,double,double>>,
            tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_1",
                                                          get_docstring("ObservationSimulator_1").c_str() );

    py::class_<tom::ObservationSimulator<2,double,double>,
            std::shared_ptr<tom::ObservationSimulator<2,double,double>>,
            tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_2",
                                                          get_docstring("ObservationSimulator_2").c_str() );

    py::class_<tom::ObservationSimulator<3,double,double>,
            std::shared_ptr<tom::ObservationSimulator<3,double,double>>,
            tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_3",
                                                          get_docstring("ObservationSimulator_3").c_str() );

    py::class_<tom::ObservationSimulator<6,double,double>,
            std::shared_ptr<tom::ObservationSimulator<6,double,double>>,
            tom::ObservationSimulatorBase<double,double>>(m, "ObservationSimulator_6",
                                                          get_docstring("ObservationSimulator_6").c_str() );

    m.def("simulate_observations",
          &tss::simulateObservations< >,
          py::arg("simulation_settings"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("simulate_observations").c_str() );

    m.def("compute_target_angles_and_range",
          &tss::getTargetAnglesAndRange,
          py::arg("bodies"),
          py::arg("station_id" ),
          py::arg("target_body" ),
          py::arg("observation_times"),
          py::arg("is_station_transmitting"),
          get_docstring("compute_target_angles_and_range").c_str() );


    py::class_< tom::ObservationCollection<>,
            std::shared_ptr<tom::ObservationCollection<>>>(m, "ObservationCollection",
                                                           get_docstring("ObservationCollection").c_str() )
            .def_property_readonly("concatenated_times", &tom::ObservationCollection<>::getConcatenatedTimeVector,
                                   get_docstring("ObservationCollection.concatenated_times").c_str() )
            .def_property_readonly("concatenated_observations", &tom::ObservationCollection<>::getObservationVector,
                                   get_docstring("ObservationCollection.concatenated_observations").c_str() );

    py::class_< tom::SingleObservationSet<>,
            std::shared_ptr<tom::SingleObservationSet<>>>(m, "SingleObservationSet",
                                                          get_docstring("SingleObservationSet").c_str() );


    /*!
     *************** STATE TRANSITION INTERFACE ***************
     */

    py::class_<
            tp::CombinedStateTransitionAndSensitivityMatrixInterface,
            std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                m, "CombinedStateTransitionAndSensitivityMatrixInterface",
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface").c_str() )
            .def("state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"),
                 get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_sensitivity_at_epoch").c_str() )
            .def("full_state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getFullCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"),
                 get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.full_state_transition_sensitivity_at_epoch").c_str() )
            .def_property_readonly(
                "state_transition_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getStateTransitionMatrixSize,
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_size").c_str() )
            .def_property_readonly(
                "sensitivity_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getSensitivityMatrixSize,
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.sensitivity_size").c_str() )
            .def_property_readonly(
                "full_parameter_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getFullParameterVectorSize,
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.full_parameter_size").c_str() );

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
          py::arg("output_times"),
          get_docstring("propagate_covariance").c_str() );


    m.def("propagate_formal_errors",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateFormalErrors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_formal_errors").c_str() );



    /*!
     *************** ESTIMATION ***************
     */

    py::class_<
            tss::EstimationConvergenceChecker,
            std::shared_ptr<tss::EstimationConvergenceChecker>>(m, "EstimationConvergenceChecker",
                                                                get_docstring("EstimationConvergenceChecker").c_str() );
    //       .def(py::init< const unsigned int,
    //            const double,
    //            const double,
    //            const int >( ),
    //            py::arg( "maximum_iterations" ) = 5,
    //            py::arg( "minimum_residual_change" ) = 0.0,
    //            py::arg( "minimum_residual" ) = 0.0,
    //            py::arg( "number_of_iterations_without_improvement" ) = 2 );


    m.def("estimation_convergence_checker",
          &tss::estimationConvergenceChecker,
          py::arg("maximum_iterations") = 5,
          py::arg("minimum_residual_change") = 0.0,
          py::arg("minimum_residual") = 0.0,
          py::arg("number_of_iterations_without_improvement") = 2,
          get_docstring("estimation_convergence_checker").c_str() );


    py::class_<
            tss::CovarianceAnalysisInput<double, double>,
            std::shared_ptr<tss::CovarianceAnalysisInput<double, double>>>(m, "CovarianceAnalysisInput",
                                                                           get_docstring("CovarianceAnalysisInput").c_str() )
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection< > >&,
                 const int,
                 const Eigen::MatrixXd >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "parameter_size" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 get_docstring("CovarianceAnalysisInput.ctor").c_str() )
            .def( "set_constant_weight",
                  &tss::CovarianceAnalysisInput<double, double>::setConstantWeightsMatrix,
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_weight").c_str() )
            .def( "set_constant_weight_per_observable",
                  &tss::CovarianceAnalysisInput<double, double>::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_weight_per_observable").c_str() )
            .def( "define_covariance_settings",
                  &tss::CovarianceAnalysisInput<double, double>::defineCovarianceSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false,
                  get_docstring("CovarianceAnalysisInput.define_covariance_settings").c_str() )
            .def_property_readonly("weight_matrix_diagonal",
                                   &tss::CovarianceAnalysisInput<double, double>::getWeightsMatrixDiagonals);

    py::class_<
            tss::EstimationInput<double, double>,
            std::shared_ptr<tss::EstimationInput<double, double>>,
            tss::CovarianceAnalysisInput<double, double>>(m, "EstimationInput",
                                                           get_docstring("EstimationInput").c_str() )
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection< > >&,
                 const int,
                 const Eigen::MatrixXd,
                 const Eigen::VectorXd >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "parameter_size" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "apriori_parameter_deviation" ) = Eigen::VectorXd::Zero( 0 ),
                 get_docstring("EstimationInput.ctor").c_str() )
            .def( "define_estimation_settings",
                  &tss::EstimationInput<double, double>::defineEstimationSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_residuals_and_parameters_per_iteration" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false,
                  get_docstring("EstimationInput.define_estimation_settings").c_str() );

    m.attr("PodInput") = m.attr("EstimationInput");


    py::class_<
            tss::CovarianceAnalysisOutput<double, double>,
            std::shared_ptr<tss::CovarianceAnalysisOutput<double, double>>>(m, "CovarianceAnalysisOutput",
                                                                            get_docstring("CovarianceAnalysisOutput").c_str() )
            .def_property_readonly("inverse_covariance",
                                   &tss::CovarianceAnalysisOutput<double, double>::getUnnormalizedInverseCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.inverse_covariance").c_str() )
            .def_property_readonly("covariance",
                                   &tss::CovarianceAnalysisOutput<double, double>::getUnnormalizedCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.covariance").c_str() )
            .def_property_readonly("formal_errors",
                                   &tss::CovarianceAnalysisOutput<double, double>::getFormalErrorVector,
                                   get_docstring("CovarianceAnalysisOutput.formal_errors").c_str() )
            .def_property_readonly("correlations",
                                   &tss::CovarianceAnalysisOutput<double, double>::getCorrelationMatrix,
                                   get_docstring("CovarianceAnalysisOutput.correlations").c_str() )
            .def_property_readonly("design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, double>::getUnnormalizedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.design_matrix").c_str() )
            .def_property_readonly("normalized_design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, double>::getNormalizedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.normalized_design_matrix").c_str() )
            .def_property_readonly("weighted_design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, double>::getUnnormalizedWeightedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.weighted_design_matrix").c_str() )
            .def_property_readonly("weighted_normalized_design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, double>::getNormalizedWeightedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.weighted_normalized_design_matrix").c_str() )
            .def_readonly("normalization_terms",
                          &tss::CovarianceAnalysisOutput<double, double>::designMatrixTransformationDiagonal_,
                          get_docstring("CovarianceAnalysisOutput.normalization_terms").c_str() );


    py::class_<
            tss::EstimationOutput<double, double>,
            std::shared_ptr<tss::EstimationOutput<double, double>>,
            tss::CovarianceAnalysisOutput<double, double>>(m, "EstimationOutput",
                                            get_docstring("EstimationOutput").c_str() )
            .def_property_readonly("residual_history",
                                   &tss::EstimationOutput<double, double>::getResidualHistoryMatrix,
                                   get_docstring("EstimationOutput.residual_history").c_str() )
            .def_property_readonly("parameter_history",
                                   &tss::EstimationOutput<double, double>::getParameterHistoryMatrix,
                                   get_docstring("EstimationOutput.parameter_history").c_str() )
            .def_readonly("final_residuals",
                          &tss::EstimationOutput<double, double>::residuals_,
                          get_docstring("EstimationOutput.final_residuals").c_str() );

    m.attr("PodOutput") = m.attr("EstimationOutput");



}

}
}
}// namespace tudatpy
