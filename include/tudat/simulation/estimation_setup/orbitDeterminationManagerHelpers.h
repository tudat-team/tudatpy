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

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
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
    const int residualBlockSize = observableResidualStartAndSize.second;
    if( residualBlockSize == 0 )
    {
        return;
    }

    const int singleObservationSize = observation_models::getObservableSize( observableType );
    const int numberOfObservations = residualBlockSize / singleObservationSize;

    const std::vector< int >& wrappedComponentIndices = observation_models::getResidualWrappingComponentIndices( observableType );
    if( wrappedComponentIndices.empty( ) )
    {
        return;
    }

    const std::vector< observation_models::ResidualWrappingRange > wrappingRanges =
            observation_models::getResidualWrappingRanges( observableType );

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > currentResidualBlock =
            residuals.block( observableResidualStartAndSize.first, 0, residualBlockSize, 1 );

    // Wrap each periodic component per observation
    for( int currentObservationIndex = 0; currentObservationIndex < numberOfObservations; currentObservationIndex++ )
    {
        for( const int currentObservableComponentIndex : wrappedComponentIndices )
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

    residuals.block( observableResidualStartAndSize.first, 0, residualBlockSize, 1 ) = currentResidualBlock;
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
        const bool wrapResiduals = observation_models::isResidualWrappingRequired( currentObservableType );

        // Iterate over all link ends for current observable type in observationsAndTimes
        for( auto dataIterator : observablesIterator.second )
        {
            observation_models::LinkEnds currentLinkEnds = dataIterator.first;
            observation_models::ResidualWrappingSettings residualWrappingSettings;
            if( wrapResiduals )
            {
                residualWrappingSettings = observationSimulator.at( currentObservableType )->getResidualWrappingSettings( currentLinkEnds );
            }
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

                if( wrapResiduals )
                {
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

//! Function to retrieve the weights of a set of observations, checking their size against the number of observations
/*!
 * \param flattenedObservationData Flattened data of the observations for which the weights are to be retrieved
 * \param weightsMatrixDiagonals Diagonal of the weights matrix (returned by reference)
 * \param weightsMatrix Full weights matrix; only set when off-diagonal weights are present (returned by reference)
 * \param hasOffDiagonalWeights Boolean denoting whether off-diagonal weights are present (returned by reference)
 */
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void retrieveObservationWeights(
        const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& flattenedObservationData,
        Eigen::VectorXd& weightsMatrixDiagonals,
        Eigen::SparseMatrix< double >& weightsMatrix,
        bool& hasOffDiagonalWeights )
{
    const int totalNumberOfObservations = static_cast< int >( flattenedObservationData.getObservationVector( ).size( ) );

    weightsMatrixDiagonals = flattenedObservationData.getWeightVector( );
    hasOffDiagonalWeights = flattenedObservationData.hasOffDiagonalWeights( );
    if( hasOffDiagonalWeights )
    {
        weightsMatrix = flattenedObservationData.getSparseWeightMatrix( );
        if( weightsMatrix.rows( ) != totalNumberOfObservations || weightsMatrix.cols( ) != totalNumberOfObservations )
        {
            throw std::runtime_error( "Error when estimating parameters, size of weights matrix (" +
                                      std::to_string( weightsMatrix.rows( ) ) + ", " + std::to_string( weightsMatrix.cols( ) ) +
                                      ") is not compatible with number of observations (" + std::to_string( totalNumberOfObservations ) +
                                      ")" );
        }
    }
    else if( weightsMatrixDiagonals.rows( ) != totalNumberOfObservations )
    {
        throw std::runtime_error( "Error when estimating parameters, size of weights diagonal (" +
                                  std::to_string( weightsMatrixDiagonals.rows( ) ) + ") is not compatible with number of observations (" +
                                  std::to_string( totalNumberOfObservations ) + ")" );
    }
}

//! Function to extract the rows of a matrix or vector that belong to the observations used in the estimation
/*!
 * When outlier rejection is used, residuals and partials are computed for all observations, including the ones that
 * are currently rejected, while the least-squares calculation may only use the observations that are active. The rows
 * of the active observations are a subset of the rows of the complete data, in the same order, so this function
 * copies, per active observation, the rows of that observation from the complete data to the estimation data.
 *
 * \param computationData Flattened data of all observations, defining the row order of computedRows
 * \param estimationData Flattened data of the observations used in the estimation, defining the row order of the output
 * \param computedRows Matrix (or vector) with one row per scalar entry of computationData
 * \return Matrix (or vector) with one row per scalar entry of estimationData
 */
template< typename ObservationScalarType, typename TimeType, typename EigenType >
EigenType extractEstimationObservationRows(
        const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& computationData,
        const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& estimationData,
        const EigenType& computedRows )
{
    const int numberOfEstimationRows = static_cast< int >( estimationData.getObservationVector( ).size( ) );
    EigenType estimationRows = EigenType::Zero( numberOfEstimationRows, computedRows.cols( ) );

    // Each observation occupies one row per component of its observable, so the loop advances by the number of rows
    // of the observation that was just copied.
    const std::vector< unsigned int >& observationIds = estimationData.getObservationIds( );
    for( int estimationRow = 0; estimationRow < numberOfEstimationRows; )
    {
        const unsigned int observationId = observationIds.at( estimationRow );
        const unsigned int numberOfRows = estimationData.getScalarSizeForObservation( observationId );
        const int computationRow = computationData.getFirstVectorRowForObservation( observationId );

        if( computationRow < 0 || computationData.getScalarSizeForObservation( observationId ) != numberOfRows )
        {
            throw std::runtime_error( "Error when extracting observations used in the estimation, observation " +
                                      std::to_string( observationId ) + " is not present in the computed data." );
        }

        estimationRows.block( estimationRow, 0, numberOfRows, computedRows.cols( ) ) =
                computedRows.block( computationRow, 0, numberOfRows, computedRows.cols( ) );

        estimationRow += static_cast< int >( numberOfRows );
    }

    return estimationRows;
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
        std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& observationVectorData,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalNumberParameters,
        Eigen::MatrixXd& designMatrix,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
        const bool calculateResiduals = true,
        const bool calculatePartials = true )
{
    observationDataset->validateObservationVectorData( observationVectorData );
    if( calculatePartials && totalNumberParameters <= 0 )
    {
        throw std::runtime_error( "Error when computing observation partials; number of parameters is 0 or smaller: " +
                                  std::to_string( totalNumberParameters ) );
    }

    const int totalObservationSize = static_cast< int >( observationVectorData.getObservationVector( ).size( ) );
    // Initialize return data.
    if( calculatePartials )
    {
        designMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalNumberParameters );
    }

    if( calculateResiduals )
    {
        residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservationSize, 1 );
    }

    // Iterate over all observation sets; the entries of each set are contiguous in the observation vector data
    for( const unsigned int setId : observationVectorData.getSetIdsInRowOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = observationVectorData.getUniqueObservationIdsForSetInRowOrder( setId );
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                observationVectorData.getSetMetadata( setId );
        const observation_models::ObservableType currentObservableType = metadata.observableType_;
        const observation_models::LinkEnds currentLinkEnds = observationVectorData.getLinkDefinitionForSet( setId ).linkEnds_;
        const std::pair< int, int > observationIndices =
                std::make_pair( observationVectorData.getVectorRow( setObservationIds.front( ), 0 ),
                                static_cast< int >( setObservationIds.size( ) * metadata.observableSize_ ) );
        if( observationVectorData.getVectorRow( setObservationIds.back( ), 0 ) + static_cast< int >( metadata.observableSize_ ) !=
            observationIndices.first + observationIndices.second )
        {
            throw std::runtime_error( "Error when computing design matrix and residuals, observation set is not contiguous." );
        }

        std::vector< TimeType > times;
        for( const unsigned int observationId : setObservationIds )
        {
            times.push_back( observationVectorData.getTimes( ).at( observationVectorData.getVectorRow( observationId, 0 ) ) );
        }

        // The observation manager returns the observations and partials ordered by time, while the observations of a set in an
        // ObservationDataset need not be sorted by time: compute them for the sorted times, and put them back in the order of the set
        std::vector< std::size_t > timeOrder( times.size( ) );
        std::iota( timeOrder.begin( ), timeOrder.end( ), 0 );
        std::stable_sort( timeOrder.begin( ), timeOrder.end( ), [ &times ]( const std::size_t first, const std::size_t second ) {
            return times.at( first ) < times.at( second );
        } );
        std::vector< TimeType > sortedTimes;
        sortedTimes.reserve( times.size( ) );
        for( const std::size_t index : timeOrder )
        {
            sortedTimes.push_back( times.at( index ) );
        }

        // Compute estimated ranges and range partials from current parameter estimate.
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > sortedObservationsVector;
        Eigen::MatrixXd sortedPartialsMatrix;
        observationManagers.at( currentObservableType )
                ->computeObservationsWithPartials( sortedTimes,
                                                   currentLinkEnds,
                                                   metadata.referenceLinkEnd_,
                                                   observationVectorData.getAncillarySettingsForSet( setId ),
                                                   sortedObservationsVector,
                                                   sortedPartialsMatrix,
                                                   calculateResiduals,
                                                   calculatePartials );
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector = sortedObservationsVector;
        Eigen::MatrixXd partialsMatrix = sortedPartialsMatrix;
        const int observableSize = static_cast< int >( metadata.observableSize_ );
        if( ( calculatePartials && sortedPartialsMatrix.rows( ) != static_cast< Eigen::Index >( times.size( ) * observableSize ) ) ||
            ( calculateResiduals && sortedObservationsVector.rows( ) != static_cast< Eigen::Index >( times.size( ) * observableSize ) ) )
        {
            throw std::runtime_error( "Error when computing observations and partials of observation set " + std::to_string( setId ) +
                                      ", the number of computed values is inconsistent with the number of observations." );
        }
        for( std::size_t i = 0; i < timeOrder.size( ); ++i )
        {
            if( calculateResiduals )
            {
                observationsVector.segment( timeOrder.at( i ) * observableSize, observableSize ) =
                        sortedObservationsVector.segment( i * observableSize, observableSize );
            }
            if( calculatePartials )
            {
                partialsMatrix.block( timeOrder.at( i ) * observableSize, 0, observableSize, partialsMatrix.cols( ) ) =
                        sortedPartialsMatrix.block( i * observableSize, 0, observableSize, partialsMatrix.cols( ) );
            }
        }

        if( calculatePartials )
        {
            // Set current observation partials in matrix of all partials
            designMatrix.block( observationIndices.first, 0, observationIndices.second, totalNumberParameters ) = partialsMatrix;
        }

        // Compute residuals for current link ends and observable type.
        if( calculateResiduals )
        {
            const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observedObservations =
                    observationVectorData.getObservationVector( ).segment( observationIndices.first, observationIndices.second );
            residuals.block( observationIndices.first, 0, observationIndices.second, 1 ) = ( observedObservations - observationsVector );

            if( observation_models::isResidualWrappingRequired( currentObservableType ) )
            {
                wrapObservationResiduals< ObservationScalarType >( residuals,
                                                                   observationIndices,
                                                                   currentObservableType,
                                                                   observedObservations,
                                                                   observationManagers.at( currentObservableType )
                                                                           ->getObservationSimulator( )
                                                                           ->getResidualWrappingSettings( currentLinkEnds ) );
            }
        }
    }
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
        // Covariance from a constrained estimation has an additional row and column per constraint
        const int numberOfParameters = estimationOutput->getNormalizationTerms( ).rows( );
        initialCovariance = estimationOutput->getUnnormalizedCovarianceMatrix( ).topLeftCorner( numberOfParameters, numberOfParameters );
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
