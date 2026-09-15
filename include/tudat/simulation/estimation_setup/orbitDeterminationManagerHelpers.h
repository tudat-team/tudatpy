/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ORBITDETERMINATIONMANAGERHELPERS_H
#define TUDAT_ORBITDETERMINATIONMANAGERHELPERS_H

#include <cmath>
#include <limits>
#include <iostream>
#include <stdexcept>
#include <string>

#include "tudat/basics/timeType.h"
#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/astro/propagators/propagateCovariance.h"

namespace tudat
{

namespace simulation_setup
{

//! Express physical-space linear constraints in normalized correction coordinates.
/*!
 * If design-matrix column j is divided by normalizationValues(j), the normalized least-squares variable is the physical
 * correction multiplied by normalizationValues(j). Each constraint column is therefore divided by the same value. Each
 * complete constraint equation is subsequently scaled to have a maximum absolute multiplier of one; the right-hand side is
 * scaled with it, so the physical constraint is unchanged.
 * \param constraintMatrix Multiplier in the physical-space linear constraint equation, modified in place.
 * \param constraintRightHandSide Right-hand side of the physical-space constraint equation, modified in place.
 * \param normalizationValues Design-matrix column normalization values.
 */
inline void normalizeLinearConstraints( Eigen::MatrixXd& constraintMatrix,
                                        Eigen::VectorXd& constraintRightHandSide,
                                        const Eigen::VectorXd& normalizationValues )
{
    if( constraintMatrix.cols( ) != normalizationValues.size( ) )
    {
        throw std::runtime_error( "Error when normalizing estimation constraints, constraint matrix has " +
                                  std::to_string( constraintMatrix.cols( ) ) + " columns but " +
                                  std::to_string( normalizationValues.size( ) ) + " normalization values were provided." );
    }
    if( constraintMatrix.rows( ) != constraintRightHandSide.size( ) )
    {
        throw std::runtime_error( "Error when normalizing estimation constraints, constraint matrix has " +
                                  std::to_string( constraintMatrix.rows( ) ) + " rows but the right-hand side has size " +
                                  std::to_string( constraintRightHandSide.size( ) ) + "." );
    }

    for( int i = 0; i < constraintMatrix.cols( ); i++ )
    {
        if( !std::isfinite( normalizationValues( i ) ) || normalizationValues( i ) == 0.0 )
        {
            throw std::runtime_error( "Error when normalizing estimation constraints, normalization value at index " + std::to_string( i ) +
                                      " is zero or non-finite." );
        }
        constraintMatrix.col( i ) /= normalizationValues( i );
    }

    for( int i = 0; i < constraintMatrix.rows( ); i++ )
    {
        const double constraintScale = constraintMatrix.row( i ).cwiseAbs( ).maxCoeff( );
        if( constraintScale > 0.0 )
        {
            constraintMatrix.row( i ) /= constraintScale;
            constraintRightHandSide( i ) /= constraintScale;
        }
    }
}

template< typename ObservationScalarType >
void wrapObservationResiduals(
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
        const std::pair< int, int > observableResidualStartAndSize,
        const observation_models::ObservableType observableType,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observedObservationBlock =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >( ),
        const observation_models::ResidualWrappingSettings& residualWrappingSettings = observation_models::ResidualWrappingSettings( ) )
{
    // Get the block of residuals for this observable type
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentResidualBlock =
            residuals.block( observableResidualStartAndSize.first, 0, observableResidualStartAndSize.second, 1 );

    const int residualBlockSize = observableResidualStartAndSize.second;
    if( residualBlockSize == 0 )
    {
        return;
    }

    const int singleObservationSize = observation_models::getObservableSize( observableType );
    const int numberOfObservations = residualBlockSize / singleObservationSize;

    if( observableType == observation_models::angular_position && residualWrappingSettings.normalizeRightAscension &&
        observedObservationBlock.rows( ) != residualBlockSize )
    {
        throw std::runtime_error( "Error when wrapping normalized angular position residuals: observed observation block has size " +
                                  std::to_string( residualBlockSize ) + "." );
                                  std::to_string( observedObservationBlock.rows( ) ) + ", expected " +
                                  std::to_string( residualBlockSize ) + "." );
                                  std::to_string( residualBlockSize ) + "." );
    }

    // Determine which components are periodic and their wrapping ranges
    // based on the observable type
    std::vector< observation_models::ResidualWrappingRange > wrappingRanges;
    if( observation_models::isResidualWrappingRequired( observableType ) )
    {
        wrappingRanges = observation_models::getResidualWrappingRanges( observableType );
    }

    // Wrap each periodic component per observation
    for( int currentObservationIndex = 0; currentObservationIndex < numberOfObservations; currentObservationIndex++ )
    {
        for( int currentObservableComponentIndex = 0; currentObservableComponentIndex < singleObservationSize;
             currentObservableComponentIndex++ )
        {
            if( currentObservableComponentIndex < static_cast< int >( wrappingRanges.size( ) ) &&
                wrappingRanges[ currentObservableComponentIndex ].period( ) > 0.0 )
            {
                const int currentResidualComponentIndex = currentObservationIndex * singleObservationSize + currentObservableComponentIndex;
                double residualWrappingPeriod = wrappingRanges[ currentObservableComponentIndex ].period( );
                double residualWrappingCenter = wrappingRanges[ currentObservableComponentIndex ].center( );

                if( observableType == observation_models::angular_position && residualWrappingSettings.normalizeRightAscension &&
                    currentObservableComponentIndex == 0 )
                {
                    const double rightAscensionNormalizationFactor =
                            std::abs( std::cos( static_cast< double >( observedObservationBlock( currentResidualComponentIndex + 1 ) ) ) );
                    if( rightAscensionNormalizationFactor <= 10.0 * std::numeric_limits< double >::epsilon( ) )
                    {
                        continue;
                    }
                    residualWrappingPeriod *= rightAscensionNormalizationFactor;
                    residualWrappingCenter *= rightAscensionNormalizationFactor;
                }

                currentResidualBlock( currentResidualComponentIndex, 0 ) = currentResidualBlock( currentResidualComponentIndex, 0 ) -
                        residualWrappingPeriod *
                                std::round( ( currentResidualBlock( currentResidualComponentIndex, 0 ) - residualWrappingCenter ) /
                                            residualWrappingPeriod );
            }
        }
    }

    residuals.block( observableResidualStartAndSize.first, 0, observableResidualStartAndSize.second, 1 ) = currentResidualBlock;
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulator,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( observationsCollection->getTotalObservableSize( ), 1 );

    typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets sortedObservations =
            observationsCollection->getObservationsSets( );

    // Iterate over all observable types in observationsAndTimes
    for( auto observablesIterator : sortedObservations )
    {
        observation_models::ObservableType currentObservableType = observablesIterator.first;

        // Iterate over all link ends for current observable type in observationsAndTimes
        for( auto dataIterator : observablesIterator.second )
        {
            observation_models::LinkEnds currentLinkEnds = dataIterator.first;
            const observation_models::ResidualWrappingSettings residualWrappingSettings =
                    observationSimulator.at( currentObservableType )->getResidualWrappingSettings( currentLinkEnds );
            for( unsigned int i = 0; i < dataIterator.second.size( ); i++ )
            {
                std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > currentObservations =
                        dataIterator.second.at( i );
                std::pair< int, int > observationIndices =
                        observationsCollection->getObservationSetStartAndSize( ).at( currentObservableType ).at( currentLinkEnds ).at( i );

                // Compute estimated ranges and range partials from current parameter estimate.
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
                observationSimulator.at( currentObservableType )
                        ->computeObservations( currentObservations->getObservationTimes( ),
                                               currentLinkEnds,
                                               currentObservations->getReferenceLinkEnd( ),
                                               currentObservations->getAncillarySettings( ),
                                               observationsVector );

                residuals.block( observationIndices.first, 0, observationIndices.second, 1 ) =
                        ( currentObservations->getObservationsVector( ) - observationsVector );

                wrapObservationResiduals< ObservationScalarType >( residuals,
                                                                   observationIndices,
                                                                   currentObservableType,
                                                                   currentObservations->getObservationsVector( ),
                                                                   residualWrappingSettings );
            }
        }
    }
}

//! Function to calculate the observation partials matrix and residuals
/*!
 *  This function calculates the observation partials matrix and residuals, based on the state transition matrix,
 *  sensitivity matrix and body states resulting from the previous numerical integration iteration.
 *  Partials and observations are calculated by the observationManagers_.
 *  \param observationsAndTimes Observable values and associated time tags, per observable type and set of link ends.
 *  \param parameterVectorSize Length of the vector of estimated parameters
 *  \param totalObservationSize Total number of observations in observationsAndTimes map.
 *  \param residualsAndPartials Pair of residuals of computed w.r.t. input observable values and partials of
 *  observables w.r.t. parameter vector (return by reference).
 */
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateDesignMatrixAndResiduals(
        std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalNumberParameters,
        const int totalObservationSize,
        Eigen::MatrixXd& designMatrix,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
        const bool calculateResiduals = true,
        const bool calculatePartials = true )
{
    if( calculatePartials && totalNumberParameters <= 0 )
    {
        throw std::runtime_error( "Error when computing observation partials; number of parameters is 0 or smaller: " +
                                  std::to_string( totalNumberParameters ) );
    }

    // Initialize return data.
    if( calculatePartials )
    {
        designMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalNumberParameters );
    }

    if( calculateResiduals )
    {
        residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservationSize, 1 );
    }

    typename observation_models::ObservationCollection< ObservationScalarType, TimeType >::SortedObservationSets sortedObservations =
            observationsCollection->getObservationsSets( );

    // Iterate over all observable types in observationsAndTimes
    for( auto observableIt : sortedObservations )
    {
        observation_models::ObservableType currentObservableType = observableIt.first;

        // Iterate over all link ends for current observable type in observationsAndTimes
        for( auto linkEndIt : observableIt.second )
        {
            observation_models::LinkEnds currentLinkEnds = linkEndIt.first;
            observation_models::ResidualWrappingSettings residualWrappingSettings;
            if( calculateResiduals )
            {
                residualWrappingSettings = observationManagers.at( currentObservableType )
                                                   ->getObservationSimulator( )
                                                   ->getResidualWrappingSettings( currentLinkEnds );
            }
            for( unsigned int i = 0; i < linkEndIt.second.size( ); i++ )
            {
                std::shared_ptr< observation_models::SingleObservationSet< ObservationScalarType, TimeType > > currentObservations =
                        linkEndIt.second.at( i );
                std::pair< int, int > observationIndices =
                        observationsCollection->getObservationSetStartAndSize( ).at( currentObservableType ).at( currentLinkEnds ).at( i );

                if( observationIndices.second > 0 )
                {
                    // Compute estimated ranges and range partials from current parameter estimate.
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
                    Eigen::MatrixXd partialsMatrix;
                    observationManagers.at( currentObservableType )
                            ->computeObservationsWithPartials( currentObservations->getObservationTimes( ),
                                                               currentLinkEnds,
                                                               currentObservations->getReferenceLinkEnd( ),
                                                               currentObservations->getAncillarySettings( ),
                                                               observationsVector,
                                                               partialsMatrix,
                                                               calculateResiduals,
                                                               calculatePartials );

                    if( calculatePartials )
                    {
                        // Set current observation partials in matrix of all partials
                        designMatrix.block( observationIndices.first, 0, observationIndices.second, totalNumberParameters ) =
                                partialsMatrix;
                    }

                    // Compute residuals for current link ends and observable type.
                    if( calculateResiduals )
                    {
                        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsVector =
                                currentObservations->getObservationsVector( ) - observationsVector;
                        residuals.block( observationIndices.first, 0, observationIndices.second, 1 ) = ( residualsVector );

                        wrapObservationResiduals< ObservationScalarType >( residuals,
                                                                           observationIndices,
                                                                           currentObservableType,
                                                                           currentObservations->getObservationsVector( ),
                                                                           residualWrappingSettings );
                    }
                }
            }
        }
    }
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateDesignMatrix(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalNumberParameters,
        const int totalObservationSize,
        Eigen::MatrixXd& designMatrix )
{
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > dummyVector;
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( observationsCollection,
                                                                          observationManagers,
                                                                          totalNumberParameters,
                                                                          totalObservationSize,
                                                                          designMatrix,
                                                                          dummyVector,
                                                                          false,
                                                                          true );
}

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
        const std::shared_ptr< observation_models::ObservationCollection< ObservationScalarType, TimeType > > observationsCollection,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalObservationSize,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    Eigen::VectorXd dummyMatrix;
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >(
            observationsCollection, observationManagers, 0, totalObservationSize, dummyMatrix, residuals, true, false );
}

//! Function to propagate full covariance at the initial time to state formal errors at later times
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::map< double, Eigen::MatrixXd > propagateCovarianceFromObjects(
        const std::shared_ptr< CovarianceAnalysisOutput< ObservationScalarType, TimeType > > estimationOutput,
        const std::shared_ptr< propagators::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    Eigen::MatrixXd initialCovariance;
    if( !estimationOutput->considerParametersIncluded_ )
    {
        initialCovariance = estimationOutput->getUnnormalizedCovarianceMatrix( );
    }
    else
    {
        Eigen::MatrixXd parameterCovariance =
                estimationOutput->getUnnormalizedCovarianceMatrix( ) + estimationOutput->getConsiderCovarianceContribution( );
        Eigen::MatrixXd considerCovariance = estimationOutput->considerCovariance_;
        initialCovariance = Eigen::MatrixXd::Zero( parameterCovariance.rows( ) + considerCovariance.rows( ),
                                                   parameterCovariance.cols( ) + considerCovariance.cols( ) );
        initialCovariance.block( 0, 0, parameterCovariance.rows( ), parameterCovariance.cols( ) ) = parameterCovariance;
        initialCovariance.block(
                parameterCovariance.rows( ), parameterCovariance.cols( ), considerCovariance.rows( ), considerCovariance.cols( ) ) =
                considerCovariance;
    }
    return propagators::propagateCovariance( initialCovariance, stateTransitionInterface, evaluationTimes );
}

}  // namespace simulation_setup

}  // namespace tudat

#endif  // TUDAT_ORBITDETERMINATIONMANAGERHELPERS_H
