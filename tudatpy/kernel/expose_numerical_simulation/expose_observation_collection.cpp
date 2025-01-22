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



void expose_observation_collection(py::module &m) {



    py::class_< tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "ObservationCollection",
                                                           get_docstring("ObservationCollection").c_str() )
            .def(py::init< std::vector< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > > >(),
                 py::arg("observation_sets") )
            .def_property_readonly("concatenated_times", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getConcatenatedTimeVector,
                                   get_docstring("ObservationCollection.concatenated_times").c_str() )
            .def_property_readonly("concatenated_float_times", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getConcatenatedDoubleTimeVector,
                                   get_docstring("ObservationCollection.concatenated_times").c_str() )
            .def_property_readonly("concatenated_weights", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getUnparsedConcatenatedWeights,
                                   get_docstring("ObservationCollection.concatenated_weights").c_str() )
            .def_property_readonly("concatenated_observations", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationVector,
                                   get_docstring("ObservationCollection.concatenated_observations").c_str() )
            .def_property_readonly("concatenated_link_definition_ids",
                                   py::overload_cast< >( &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getConcatenatedLinkEndIds ),
                                   get_docstring("ObservationCollection.concatenated_link_definition_ids").c_str() )
            .def_property_readonly("link_definition_ids", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getInverseLinkEndIdentifierMap,
                                   get_docstring("ObservationCollection.link_definition_ids").c_str() )
            .def_property_readonly("observable_type_start_index_and_size", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationTypeStartAndSize,
                                   get_docstring("ObservationCollection.observable_type_start_index_and_size").c_str() )
            .def_property_readonly("observation_set_start_index_and_size", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getObservationSetStartAndSizePerLinkEndIndex,
                                   get_docstring("ObservationCollection.observation_set_start_index_and_size").c_str() )
            .def_property_readonly("observation_vector_size", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getTotalObservableSize,
                                   get_docstring("ObservationCollection.observation_vector_size").c_str() )
            .def_property_readonly("sorted_observation_sets", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getSortedObservationSets,
                                   get_docstring("ObservationCollection.sorted_observation_sets").c_str() )
            .def_property_readonly("link_ends_per_observable_type", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getLinkEndsPerObservableType,
                                   get_docstring("ObservationCollection.link_ends_per_observable_type").c_str() )
            .def_property_readonly("link_definitions_per_observable", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getLinkDefinitionsPerObservable,
                                   get_docstring("ObservationCollection.link_definitions_per_observable").c_str() )
            .def_property_readonly("time_bounds", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getTimeBounds,
                                   get_docstring("ObservationCollection.time_bounds").c_str() )
            .def_property_readonly("sorted_per_set_time_bounds", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getSortedObservationSetsTimeBounds,
                                   get_docstring("ObservationCollection.time_bounds").c_str() )
            .def("set_observations", py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                 &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::setObservations ),
                 py::arg("observations"), get_docstring("set_observations").c_str( ) )
            .def("set_observations", py::overload_cast<
                const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                const std::shared_ptr< tom::ObservationCollectionParser > >(
                &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::setObservations ),
                py::arg("observations"), py::arg("observation_parser"), get_docstring("set_observations").c_str( ) )
            .def("set_observations", py::overload_cast<
                    const std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                    &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::setObservations ),
                 py::arg("observations_per_parser"), get_docstring("set_observations").c_str( ) )
            .def("set_residuals", py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                 &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals ),
                 py::arg("residuals"), get_docstring("set_residuals").c_str( ) )
            .def("set_residuals", py::overload_cast<
                 const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                 const std::shared_ptr< tom::ObservationCollectionParser > >(
                 &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals ),
                 py::arg("residuals"), py::arg("observation_parser"), get_docstring("set_residuals").c_str( ) )
            .def("set_residuals", py::overload_cast<
                    const std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                    &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals ),
                 py::arg("residuals_per_parser"), get_docstring("set_residuals").c_str( ) )
            .def("get_link_definitions_for_observables", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getLinkDefinitionsForSingleObservable,
                 py::arg( "observable_type" ),
                 get_docstring("ObservationCollection.get_link_definitions_for_observables").c_str() )
            .def("get_single_link_and_type_observations", &tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>::getSingleLinkAndTypeObservationSets,
                                   py::arg( "observable_type" ),
                                   py::arg( "link_definition" ),
                                   get_docstring("ObservationCollection.get_single_link_and_type_observations").c_str() )
            .def( "get_observable_types", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableTypes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observable_types" ).c_str() )
            .def( "get_bodies_in_link_ends", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getBodiesInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_bodies_in_link_ends" ).c_str() )
            .def( "get_reference_points_in_link_ends", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getReferencePointsInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_reference_points_in_link_ends" ).c_str() )
            .def( "get_time_bounds_list", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsList,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_time_bounds_list" ).c_str() )
            .def( "get_time_bounds_per_set", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSet,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_time_bounds_per_set" ).c_str() )
            .def( "get_observations", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observations" ).c_str() )
            .def( "get_concatenated_observations", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_observations" ).c_str() )
            .def( "get_observation_times", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observation_times" ).c_str() )
            .def( "get_concatenated_observation_times", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_observation_times" ).c_str() )
            .def( "get_concatenated_float_observation_times", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_float_observation_times" ).c_str() )
            .def( "get_observations_and_times", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observations_and_times" ).c_str() )
            .def( "get_concatenated_observations_and_times", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_observations_and_times" ).c_str() )
            .def( "get_concatenated_link_definition_ids",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                  py::arg( "observation_parser" ),
                  get_docstring( "get_concatenated_link_definition_ids" ).c_str() )
            .def( "get_weights", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_weights" ).c_str() )
            .def( "get_concatenated_weights", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_weights" ).c_str() )
            .def( "get_residuals", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_residuals" ).c_str() )
            .def( "get_concatenated_residuals", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_residuals" ).c_str() )
            .def( "get_rms_residuals", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_rms_residuals" ).c_str() )
            .def( "get_mean_residuals", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_mean_residuals" ).c_str() )
            .def( "get_computed_observations", &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_computed_observations" ).c_str() )
            .def( "get_concatenated_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_computed_observations" ).c_str() )
            .def( "set_constant_weight", py::overload_cast<
                    const double,
                    const std::shared_ptr< tom::ObservationCollectionParser > >(
                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_constant_weight" ).c_str() )
            .def( "set_constant_weight", py::overload_cast<
                    const Eigen::VectorXd,
                    const std::shared_ptr< tom::ObservationCollectionParser > >(
                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_constant_weight" ).c_str() )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, double > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  get_docstring( "set_constant_weight_per_observation_parser" ).c_str() )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  get_docstring( "set_constant_weight_per_observation_parser" ).c_str() )
            .def( "set_tabulated_weights",
                  py::overload_cast<
                  const Eigen::VectorXd,
                  const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_tabulated_weights" ).c_str() )
            .def( "set_tabulated_weights",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  get_docstring( "set_tabulated_weights" ).c_str() )
             .def( "append",
                   &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::appendObservationCollection,
                   py::arg( "observation_collection_to_append" ) )
        .def( "filter_observations",
                  py::overload_cast<
                  const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                  std::shared_ptr< tom::ObservationFilterBase > >&,
                  const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg("save_filtered_observations") = true,
                  get_docstring( "filter_observations" ).c_str( ) )
            .def( "filter_observations",
                  py::overload_cast<
                  std::shared_ptr< tom::ObservationFilterBase >,
                  std::shared_ptr< tom::ObservationCollectionParser >,
                  const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  py::arg("save_filtered_observations") = true,
                  get_docstring( "filter_observations" ).c_str(  ) )
            .def( "split_observation_sets",
                  py::overload_cast<
                  std::shared_ptr< tom::ObservationSetSplitterBase >,
                  std::shared_ptr< tom::ObservationCollectionParser > >(
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::splitObservationSets ),
                  py::arg( "observation_set_splitter" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "split_observation_sets" ).c_str( ) )
            .def( "get_single_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservationSets,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_single_observation_sets" ).c_str( ) )
            .def( "print_observation_sets_start_and_size",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::printObservationSetsStartAndSize,
                  get_docstring( "print_observation_sets_start_index_and_size" ).c_str( ) )
            .def( "remove_single_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationSets ),
                  py::arg("observation_parser"),
                  get_docstring( "remove_single_observation_sets" ).c_str( ) )
            .def("set_reference_point",
                 py::overload_cast< tss::SystemOfBodies&, const Eigen::Vector3d&, const std::string&,
                 const std::string&, const tom::LinkEndType, const std::shared_ptr< tom::ObservationCollectionParser > >(
                         &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE>::setReferencePoint ),
                 py::arg( "bodies" ),
                 py::arg( "antenna_position" ),
                 py::arg( "antenna_name" ),
                 py::arg( "spacecraft_name" ),
                 py::arg( "link_end_type" ),
                 py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                 get_docstring("set_reference_point").c_str( ) )
            .def("set_reference_points",
                  py::overload_cast< tss::SystemOfBodies&, const std::map< double, Eigen::Vector3d >&, const std::string&,
                  const tom::LinkEndType, const std::shared_ptr< tom::ObservationCollectionParser > >(
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoints ),
                  py::arg("bodies"),
                  py::arg("antenna_switch_history"),
                  py::arg("spacecraft_name"),
                  py::arg("link_end_type"),
                  py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_reference_points" ).c_str( ) )
            .def("set_reference_point",
                 py::overload_cast< tss::SystemOfBodies&, const std::shared_ptr< te::Ephemeris >, const std::string&,
                 const std::string&, const tom::LinkEndType, const std::shared_ptr< tom::ObservationCollectionParser > >(
                 &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                 py::arg("bodies"),
                 py::arg("antenna_body_fixed_ephemeris"),
                 py::arg("antenna_name"),
                 py::arg("spacecraft_name"),
                 py::arg("link_end_type"),
                 py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
                 get_docstring( "set_reference_points" ).c_str( ) )
            .def("set_transponder_delay",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTransponderDelay,
                  py::arg( "spacecraft_name" ),
                  py::arg( "transponder_delay" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_transponder_delay" ).c_str( ) )
            .def( "remove_empty_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeEmptySingleObservationSets,
                  get_docstring( "remove_empty_observation_sets" ).c_str( ) )
            .def( "add_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariable,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "bodies" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "add_dependent_variable" ).c_str( ) )
            .def( "dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "dependent_variable" ).c_str( ) )
            .def( "concatenated_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "concatenated_dependent_variable" ).c_str( ) )
            .def( "compatible_dependent_variable_settings",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "compatible_dependent_variable_settings" ).c_str( ) )
            .def( "compatible_dependent_variables_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "compatible_dependent_variables_list" ).c_str( ) )
            .def( "dependent_variable_history_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSet,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "dependent_variable_history_per_set" ).c_str( ) )
            .def( "dependent_variable_history",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "dependent_variable_history" ).c_str( ) );

    m.def("merge_observation_collections",
          &tss::mergeObservationCollections< STATE_SCALAR_TYPE, TIME_TYPE >,
          py::arg("observation_collection_list") );

    m.def( "create_single_observation_set",
           py::overload_cast<
                const tom::ObservableType, const tom::LinkEnds&,
                const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                const std::vector< TIME_TYPE >, const tom::LinkEndType,
                const std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
                &tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings" ),
           get_docstring( "create_single_observation_set" ).c_str( ) );

    m.def("filter_observations",
          py::overload_cast<
                  const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationFilterBase >, const bool >(
                  &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
            py::arg( "original_observation_set" ),
            py::arg( "observation_filter" ),
            py::arg( "save_filtered_observations" ) = false,
            get_docstring( "filter_observations" ).c_str( ) );

    m.def("split_observation_set",
          py::overload_cast<
                  const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationSetSplitterBase >,
                  const bool >(
                  &tom::splitObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
          py::arg( "original_observation_set" ),
          py::arg( "observation_splitter" ),
          py::arg( "print_warning" ) = true,
          get_docstring( "split_observation_set" ).c_str( ) );




}

}
}
}// namespace tudatpy
