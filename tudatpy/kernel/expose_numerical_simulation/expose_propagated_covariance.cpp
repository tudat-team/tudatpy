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

#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"

#include "docstrings.h"
#include "scalarTypes.h"

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
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;


namespace tudat
{

namespace propagators
{

std::map< double, Eigen::MatrixXd > propagateCovarianceRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput< STATE_SCALAR_TYPE, TIME_TYPE > > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance(
                propagatedCovariance, estimationOutput->getUnnormalizedCovarianceMatrix( ),
                orbitDeterminationManager->getStateTransitionAndSensitivityMatrixInterface( ), evaluationTimes );

    tss::SystemOfBodies bodies = orbitDeterminationManager->getBodies( );

    std::shared_ptr< tep::EstimatableParameterSet<STATE_SCALAR_TYPE> > parameterSet = orbitDeterminationManager->getParametersToEstimate( );

    std::map< int, std::shared_ptr< tep::EstimatableParameter< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > > > > initialStates =
        parameterSet->getInitialStateParameters( );
    std::map< std::pair< std::string, std::string >, std::vector< int > > transformationList;
    for( auto it : initialStates )
    {
        if( std::dynamic_pointer_cast< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second ) )
        {
            std::shared_ptr< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > > currentInitialState =
                    std::dynamic_pointer_cast< tep::InitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second );
            transformationList[ std::make_pair( currentInitialState->getParameterName( ).second.first, currentInitialState->getCentralBody( ) ) ].push_back(
                        it.first );

        }
        else if( std::dynamic_pointer_cast< tep::ArcWiseInitialTranslationalStateParameter< STATE_SCALAR_TYPE > >( it.second ) )
        {
            throw std::runtime_error( "Error, multi-arc not yet supported in automatic covariance conversion" );
        }
    }

    Eigen::Matrix3d currentInertialToRswPosition;
    Eigen::Matrix6d currentInertialToRswState;
    Eigen::MatrixXd currentFullInertialToRswState = Eigen::MatrixXd::Zero( 6, 6 );

    std::map< double, Eigen::MatrixXd > propagatedRswCovariance;
    for( auto it : propagatedCovariance )
    {
        double currentTime = static_cast< double >( it.first );
        Eigen::MatrixXd currentCovariance = it.second;
        currentFullInertialToRswState.setZero( );

        for( auto it_body : transformationList )
        {
            Eigen::Vector6d relativeState =
                    bodies.getBody( it_body.first.first )->getStateInBaseFrameFromEphemeris( currentTime ) -
                    bodies.getBody( it_body.first.second )->getStateInBaseFrameFromEphemeris( currentTime );
            currentInertialToRswPosition = trf::getInertialToRswSatelliteCenteredFrameRotationMatrix( relativeState );
            currentInertialToRswState.block( 0, 0, 3, 3 ) = currentInertialToRswPosition;
            currentInertialToRswState.block( 3, 3, 3, 3 ) = currentInertialToRswPosition;
            for( unsigned int j = 0; j < it_body.second.size( ); j++ )
            {
                int currentStartIndex = it_body.second.at( j );
                currentFullInertialToRswState.block( currentStartIndex, currentStartIndex, 6, 6 ) =  currentInertialToRswState;
            }
        }
        propagatedRswCovariance[ currentTime ] = currentFullInertialToRswState * currentCovariance * currentFullInertialToRswState.transpose( );

    }
    return propagatedRswCovariance;
}


std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE> > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedRswCovariance = propagateCovarianceRsw(
            estimationOutput, orbitDeterminationManager, evaluationTimes );

    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedRswCovariance ),
                           utilities::createVectorFromMapValues(
                               propagatedRswCovariance ) );
}

std::map< double, Eigen::VectorXd > propagateFormalErrorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE> > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;

    propagatedCovariance = propagateCovarianceRsw(
                estimationOutput, orbitDeterminationManager, evaluationTimes );
    tp::convertCovarianceHistoryToFormalErrorHistory( propagatedFormalErrors, propagatedCovariance );

    return propagatedFormalErrors;
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE> > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors =
            propagateFormalErrorsRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );
    tp::propagateFormalErrorsRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedFormalErrors ),
                           utilities::createVectorFromMapValues(
                               propagatedFormalErrors ) );
}



std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance(
                propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedCovariance ),
                           utilities::createVectorFromMapValues(
                               propagatedCovariance ) );
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;
    tp::propagateFormalErrors( propagatedFormalErrors, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedFormalErrors ),
                           utilities::createVectorFromMapValues(
                               propagatedFormalErrors ) );
}

}

namespace simulation_setup
{

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > getTargetAnglesAndRangeVector(
        const simulation_setup::SystemOfBodies& bodies,
        const std::pair< std::string, std::string > groundStationId,
        const std::string& targetBody,
        const std::vector< double > times,
        const bool transmittingToTarget )
{
    std::map< double, Eigen::VectorXd > targetAnglesAndRange = getTargetAnglesAndRange(
                bodies, groundStationId, targetBody, times, transmittingToTarget );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               targetAnglesAndRange ),
                           utilities::createVectorFromMapValues(
                               targetAnglesAndRange ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSetWithoutDependentVariables(
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType > observationTimes,
            const tom::LinkEndType referenceLinkEnd,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr )
{
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >(
            observableType, linkEnds, observations, observationTimes, referenceLinkEnd,
            std::vector< Eigen::VectorXd >( ), nullptr, ancilliarySettings );
}


}

}

namespace tudatpy {
namespace numerical_simulation {
namespace estimation {



void expose_propagated_covariance(py::module &m) {


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
                 py::arg("add_central_body_dependency") = true,
                 py::arg("arc_defining_bodies" ) = std::vector< std::string >( ),
                 get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_sensitivity_at_epoch").c_str() )
            .def("full_state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getFullCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"),
                 py::arg("add_central_body_dependency") = true,
                 py::arg("arc_defining_bodies" ) = std::vector< std::string >( ),
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

    m.def("propagate_covariance_rsw_split_output",
          &tp::propagateCovarianceVectorsRsw,
          py::arg("covariance_output"),
          py::arg("estimator"),
          py::arg("output_times"),
          get_docstring("propagate_covariance_rsw_split_output").c_str( ) );


    m.def("propagate_formal_errors_rsw_split_output",
          &tp::propagateFormalErrorVectorsRsw,
          py::arg("covariance_output"),
          py::arg("estimator"),
          py::arg("output_times"),
          get_docstring("propagate_formal_errors_rsw_split_output").c_str( ) );

    m.def("propagate_covariance_split_output",
          py::overload_cast<
          const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >(
              &tp::propagateCovarianceVectors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_covariance_split_output").c_str() );

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

    m.def("propagate_formal_errors_split_output",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateFormalErrorVectors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_formal_errors").c_str() );


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

    m.def("estimation_convergence_checker",
          &tss::estimationConvergenceChecker,
          py::arg("maximum_iterations") = 5,
          py::arg("minimum_residual_change") = 0.0,
          py::arg("minimum_residual") = 0.0,
          py::arg("number_of_iterations_without_improvement") = 2,
          get_docstring("estimation_convergence_checker").c_str() );


    py::class_<
            tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "CovarianceAnalysisInput",
                                                                           get_docstring("CovarianceAnalysisInput").c_str() )
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE> >&,
                 const Eigen::MatrixXd,
                 const Eigen::MatrixXd  >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "consider_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 get_docstring("CovarianceAnalysisInput.ctor").c_str() )
            .def( "set_constant_weight",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantWeightsMatrix,
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_weight").c_str() )
            .def( "set_weights_from_observation_collection",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setWeightsFromObservationCollection,
                  get_docstring("CovarianceAnalysisInput.set_weights_from_observation_collection").c_str() )
            .def( "set_constant_single_observable_weight",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantSingleObservableWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_weight").c_str() )
            .def( "set_constant_single_observable_vector_weight",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantSingleObservableVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_vector_weight").c_str() )
            .def( "set_constant_single_observable_and_link_end_weight",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_and_link_end_weight").c_str() )
            .def( "set_constant_single_observable_and_link_end_vector_weight",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantSingleObservableAndLinkEndsVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_and_link_end_vector_weight").c_str() )
            .def( "set_total_single_observable_and_link_end_vector_weight",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setTabulatedSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight_vector" ),
                  get_docstring("CovarianceAnalysisInput.set_total_single_observable_and_link_end_vector_weight").c_str() )
            .def( "set_constant_weight_per_observable",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_weight_per_observable").c_str() )
            .def( "set_constant_vector_weight_per_observable",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setConstantPerObservableVectorWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_vector_weight_per_observable").c_str() )
            .def( "define_covariance_settings",
                  &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::defineCovarianceSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "limit_condition_number_for_warning" ) = 1.0E8,
                  get_docstring("CovarianceAnalysisInput.define_covariance_settings").c_str() )
            .def_property("weight_matrix_diagonal",
                          &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::getWeightsMatrixDiagonals,
                          &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>::setWeightsMatrixDiagonals,
                          get_docstring("CovarianceAnalysisInput.weight_matrix_diagonal").c_str() );

    py::class_<
            tss::EstimationInput<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tss::EstimationInput<STATE_SCALAR_TYPE, TIME_TYPE>>,
            tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "EstimationInput",
                                                          get_docstring("EstimationInput").c_str() )
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE> >&,
                 const Eigen::MatrixXd,
                 std::shared_ptr< tss::EstimationConvergenceChecker >,
                 const Eigen::MatrixXd,
                 const Eigen::VectorXd,
                 const bool >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "convergence_checker" ) = std::make_shared< tss::EstimationConvergenceChecker >( ),
                 py::arg( "consider_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "consider_parameters_deviations" ) = Eigen::VectorXd::Zero( 0 ),
                 py::arg( "apply_final_parameter_correction" ) = true,
                 get_docstring("EstimationInput.ctor").c_str() )
            .def( "define_estimation_settings",
                  &tss::EstimationInput<STATE_SCALAR_TYPE, TIME_TYPE>::defineEstimationSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_residuals_and_parameters_per_iteration" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false,
                  py::arg( "limit_condition_number_for_warning" ) = 1.0E8,
                  py::arg( "condition_number_warning_each_iteration" ) = true,
                  get_docstring("EstimationInput.define_estimation_settings").c_str() );

    m.attr("PodInput") = m.attr("EstimationInput");


    py::class_<
            tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "CovarianceAnalysisOutput",
                                                                            get_docstring("CovarianceAnalysisOutput").c_str() )
            .def_property_readonly("inverse_covariance",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getUnnormalizedInverseCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.inverse_covariance").c_str() )
            .def_property_readonly("covariance",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getUnnormalizedCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.covariance").c_str() )
            .def_property_readonly("inverse_normalized_covariance",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getNormalizedInverseCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.inverse_normalized_covariance").c_str() )
            .def_property_readonly("normalized_covariance",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getNormalizedCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.normalized_covariance").c_str() )
            .def_property_readonly("formal_errors",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getFormalErrorVector,
                                   get_docstring("CovarianceAnalysisOutput.formal_errors").c_str() )
            .def_property_readonly("correlations",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getCorrelationMatrix,
                                   get_docstring("CovarianceAnalysisOutput.correlations").c_str() )
            .def_property_readonly("design_matrix",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getUnnormalizedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.design_matrix").c_str() )
            .def_property_readonly("normalized_design_matrix",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getNormalizedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.normalized_design_matrix").c_str() )
            .def_property_readonly("weighted_design_matrix",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getUnnormalizedWeightedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.weighted_design_matrix").c_str() )
            .def_property_readonly("weighted_normalized_design_matrix",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getNormalizedWeightedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.weighted_normalized_design_matrix").c_str() )
            .def_property_readonly("consider_covariance_contribution",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getConsiderCovarianceContribution,
                                   get_docstring("CovarianceAnalysisOutput.consider_covariance_contribution").c_str() )
            .def_property_readonly("normalized_covariance_with_consider_parameters",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getNormalizedCovarianceWithConsiderParameters,
                                   get_docstring("CovarianceAnalysisOutput.normalized_covariance_with_consider_parameters").c_str() )
            .def_property_readonly("unnormalized_covariance_with_consider_parameters",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getUnnormalizedCovarianceWithConsiderParameters,
                                   get_docstring("CovarianceAnalysisOutput.unnormalized_covariance_with_consider_parameters").c_str() )
            .def_property_readonly("normalized_design_matrix_consider_parameters",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getNormalizedDesignMatrixConsiderParameters,
                                   get_docstring("CovarianceAnalysisOutput.normalized_design_matrix_consider_parameters").c_str() )
            .def_property_readonly("consider_normalization_factors",
                                   &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getConsiderNormalizationFactors,
                                   get_docstring("CovarianceAnalysisOutput.consider_normalization_factors").c_str() )
            .def_readonly("normalization_terms",
                          &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>::designMatrixTransformationDiagonal_,
                          get_docstring("CovarianceAnalysisOutput.normalization_terms").c_str() );


    py::class_<
            tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>>,
            tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "EstimationOutput",
                                                           get_docstring("EstimationOutput").c_str() )
            .def_property_readonly("residual_history",
                                   &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getResidualHistoryMatrix,
                                   get_docstring("EstimationOutput.residual_history").c_str() )
            .def_property_readonly("parameter_history",
                                   &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getParameterHistoryMatrix,
                                   get_docstring("EstimationOutput.parameter_history").c_str() )
            .def_property_readonly("simulation_results_per_iteration",
                                   &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::getSimulationResults,
                                   get_docstring("EstimationOutput.simulation_results_per_iteration").c_str() )
            .def_readonly("final_residuals",
                          &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::residuals_,
                          get_docstring("EstimationOutput.final_residuals").c_str() )
            .def_readonly("final_parameters",
                          &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::parameterEstimate_,
                          get_docstring("EstimationOutput.final_parameters").c_str() )
            .def_readonly("best_iteration",
                          &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::bestIteration_,
                          get_docstring("EstimationOutput.best_iteration").c_str() );

    m.attr("PodOutput") = m.attr("EstimationOutput");



}

}
}
}// namespace tudatpy
