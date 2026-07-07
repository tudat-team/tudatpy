/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGERUTILITIESIMPLEMENTATION_H
#define TUDAT_ORBITDETERMINATIONMANAGERUTILITIESIMPLEMENTATION_H

#include <cmath>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/math/basic/leastSquaresTraits.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::getObservationSimulators( ) const
{
    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > > observationSimulators;

    for( typename std::map< observation_models::ObservableType,
                            std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >::
                 const_iterator managerIterator = observationManagers_.begin( );
         managerIterator != observationManagers_.end( );
         managerIterator++ )
    {
        observationSimulators.push_back( managerIterator->second->getObservationSimulator( ) );
    }

    return observationSimulators;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::normalizeAprioriCovariance(
        const Eigen::MatrixXd& inverseAPrioriCovariance,
        const Eigen::VectorXd& normalizationValues )
{
    int numberOfEstimatedParameters = inverseAPrioriCovariance.rows( );
    Eigen::MatrixXd normalizedInverseAprioriCovarianceMatrix =
            Eigen::MatrixXd::Zero( numberOfEstimatedParameters, numberOfEstimatedParameters );

    for( int j = 0; j < numberOfEstimatedParameters; j++ )
    {
        for( int k = 0; k < numberOfEstimatedParameters; k++ )
        {
            normalizedInverseAprioriCovarianceMatrix( j, k ) =
                    inverseAPrioriCovariance( j, k ) / ( normalizationValues( j ) * normalizationValues( k ) );
        }
    }
    return normalizedInverseAprioriCovarianceMatrix;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::normalizeCovariance(
        const Eigen::MatrixXd& covariance,
        const Eigen::VectorXd& normalizationFactors )
{
    int numberParameters = covariance.rows( );
    Eigen::MatrixXd normalizedCovariance = Eigen::MatrixXd::Zero( numberParameters, numberParameters );
    for( int j = 0; j < numberParameters; j++ )
    {
        for( int k = 0; k < numberParameters; k++ )
        {
            normalizedCovariance( j, k ) = covariance( j, k ) * ( normalizationFactors( j ) * normalizationFactors( k ) );
        }
    }
    return normalizedCovariance;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::normalizeDesignMatrix(
        Eigen::MatrixXd& observationMatrix )
{
    Eigen::VectorXd normalizationTerms = Eigen::VectorXd( observationMatrix.cols( ) );

    for( int i = 0; i < observationMatrix.cols( ); i++ )
    {
        Eigen::VectorXd currentVector = observationMatrix.block( 0, i, observationMatrix.rows( ), 1 );
        double minimum = currentVector.minCoeff( );
        double maximum = currentVector.maxCoeff( );
        if( std::fabs( minimum ) > maximum )
        {
            normalizationTerms( i ) = minimum;
        }
        else
        {
            normalizationTerms( i ) = maximum;
        }
        if( normalizationTerms( i ) == 0.0 )
        {
            normalizationTerms( i ) = 1.0;
        }
        currentVector = currentVector / normalizationTerms( i );

        observationMatrix.block( 0, i, observationMatrix.rows( ), 1 ) = currentVector;
    }

    return normalizationTerms;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
template< typename EstimatedDesignMatrixType >
void OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::computeDesignMatrices(
        std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput,
        EstimatedDesignMatrixType& designMatrixEstimatedParameters,
        Eigen::MatrixXd& designMatrixConsiderParameters,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >* residuals )
{
    typedef typename linear_algebra::from_eigen< EstimatedDesignMatrixType >::traits EstimatedDesignMatrixTraits;
    typedef typename linear_algebra::from_eigen< EstimatedDesignMatrixType >::value_type EstimatedDesignMatrixValueType;

    const int totalNumberOfObservations = estimationInput->getObservationCollection( )->getTotalObservableSize( );

    if( residuals != nullptr )
    {
        *residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalNumberOfObservations, 1 );
    }

    std::vector< Eigen::Triplet< EstimatedDesignMatrixValueType > > estimatedParameterTriplets;
    estimatedParameterTriplets.reserve( static_cast< std::size_t >( totalNumberOfObservations ) * 6 );
    if( considerParametersIncluded_ )
    {
        designMatrixConsiderParameters = Eigen::MatrixXd::Zero( totalNumberOfObservations, numberConsiderParameters_ );
    }
    else
    {
        designMatrixConsiderParameters = Eigen::MatrixXd::Zero( 0, 0 );
    }

    typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets sortedObservations =
            estimationInput->getObservationCollection( )->getObservationsSets( );

    for( auto observableIt : sortedObservations )
    {
        observation_models::ObservableType currentObservableType = observableIt.first;

        for( auto linkEndIt : observableIt.second )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndIt.first;
            for( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
            {
                std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > currentObservations =
                        linkEndIt.second.at( i );
                std::pair< int, int > observationIndices =
                        estimationInput->getObservationCollection( )->getObservationSetStartAndSize( ).at( currentObservableType ).at(
                                currentLinkEnds ).at( i );

                if( observationIndices.second > 0 )
                {
                    const std::vector< TimeType >& observationTimes = currentObservations->getObservationTimesReference( );
                    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentObservationVector =
                            residuals == nullptr ? Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >( ) :
                                                   currentObservations->getObservationsVector( );
                    const int rowsPerObservation = observationIndices.second / static_cast< int >( observationTimes.size( ) );
                    if( rowsPerObservation * static_cast< int >( observationTimes.size( ) ) != observationIndices.second )
                    {
                        throw std::runtime_error( "Error when computing sparse covariance design matrix: observation size is inconsistent." );
                    }

                    const int maximumNumberOfObservationsPerBatch = 100;
                    for( int firstObservation = 0; firstObservation < static_cast< int >( observationTimes.size( ) );
                         firstObservation += maximumNumberOfObservationsPerBatch )
                    {
                        const int currentBatchSize =
                                std::min( maximumNumberOfObservationsPerBatch,
                                          static_cast< int >( observationTimes.size( ) ) - firstObservation );
                        std::vector< TimeType > currentObservationTimes(
                                observationTimes.begin( ) + firstObservation,
                                observationTimes.begin( ) + firstObservation + currentBatchSize );

                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
                        Eigen::MatrixXd partialsMatrix;
                        observationManagers_.at( currentObservableType )
                                ->computeObservationsWithPartials( currentObservationTimes,
                                                                   currentLinkEnds,
                                                                   currentObservations->getReferenceLinkEnd( ),
                                                                   currentObservations->getAncillarySettings( ),
                                                                   observationsVector,
                                                                   partialsMatrix,
                                                                   residuals != nullptr,
                                                                   true );

                        const int rowOffset = observationIndices.first + firstObservation * rowsPerObservation;
                        if( residuals != nullptr )
                        {
                            residuals->block( rowOffset, 0, partialsMatrix.rows( ), 1 ) =
                                    currentObservationVector.segment( firstObservation * rowsPerObservation, partialsMatrix.rows( ) ) -
                                    observationsVector;
                        }

                        for( int row = 0; row < partialsMatrix.rows( ); row++ )
                        {
                            for( unsigned int column = 0; column < numberEstimatedParameters_; column++ )
                            {
                                const double value = partialsMatrix( row, column );
                                if( value != 0.0 )
                                {
                                    estimatedParameterTriplets.emplace_back( rowOffset + row, column, value );
                                }
                            }
                        }

                        if( considerParametersIncluded_ )
                        {
                            designMatrixConsiderParameters.block(
                                    rowOffset, 0, partialsMatrix.rows( ), numberConsiderParameters_ ) =
                                    partialsMatrix.block( 0, numberEstimatedParameters_, partialsMatrix.rows( ), numberConsiderParameters_ );
                        }
                    }
                }
            }
        }
    }

    designMatrixEstimatedParameters = EstimatedDesignMatrixTraits::from_sparse_range(
            static_cast< unsigned int >( totalNumberOfObservations ),
            static_cast< unsigned int >( numberEstimatedParameters_ ),
            estimatedParameterTriplets );
    EstimatedDesignMatrixTraits::make_compressed( designMatrixEstimatedParameters );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::getNormalizedConsiderCovariance(
        const std::shared_ptr< CovarianceAnalysisInput< ObservationScalarType, TimeType > > estimationInput,
        const Eigen::VectorXd& considerNormalizationTerms,
        Eigen::MatrixXd& normalizedConsiderCovariance )
{
    Eigen::MatrixXd unnormalizedConsiderCovariance = estimationInput->getConsiderCovariance( );
    if( unnormalizedConsiderCovariance.rows( ) == 0 && unnormalizedConsiderCovariance.cols( ) == 0 )
    {
        unnormalizedConsiderCovariance = Eigen::MatrixXd::Zero( numberConsiderParameters_, numberConsiderParameters_ );
    }
    else if( unnormalizedConsiderCovariance.rows( ) != numberConsiderParameters_ &&
             unnormalizedConsiderCovariance.cols( ) == numberConsiderParameters_ )
    {
        throw std::runtime_error( "Error, consider covariance size: [" + std::to_string( unnormalizedConsiderCovariance.rows( ) ) + ", " +
                                  std::to_string( unnormalizedConsiderCovariance.cols( ) ) +
                                  "] does not match number of consider parameters: " + std::to_string( numberConsiderParameters_ ) );
    }
    normalizedConsiderCovariance = normalizeCovariance( unnormalizedConsiderCovariance, considerNormalizationTerms );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::getObservationManager(
        const observation_models::ObservableType observableType ) const
{
    // Check if manager exists for requested observable type.
    if( observationManagers_.count( observableType ) == 0 )
    {
        throw std::runtime_error( "Error when retrieving observation manager of type " + std::to_string( observableType ) +
                                  ", manager not found" );
    }

    return observationManagers_.at( observableType );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::pair< Eigen::MatrixXd, Eigen::MatrixXd >
OrbitDeterminationManager< ObservationScalarType, TimeType, Dummy >::separateEstimatedAndConsiderDesignMatrices(
        const Eigen::MatrixXd& designMatrix,
        const int numberObservations )
{
    Eigen::MatrixXd designMatrixEstimatedParameters = designMatrix.block( 0, 0, numberObservations, numberEstimatedParameters_ );

    Eigen::MatrixXd designMatrixConsiderParameters =
            designMatrix.block( 0, numberEstimatedParameters_, numberObservations, numberConsiderParameters_ );

    return std::make_pair( designMatrixEstimatedParameters, designMatrixConsiderParameters );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERUTILITIESIMPLEMENTATION_H
