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

//! Calculate residuals using an explicit dataset vector mapping.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const observation_models::ObservationVectorData< ObservationScalarType, TimeType >& observationVectorData,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulator,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    observationDataset->validateObservationVectorData( observationVectorData );
    residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( observationVectorData.getObservationVector( ).size( ), 1 );

    for( const unsigned int setId : observationVectorData.getSetIdsInRowOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = observationVectorData.getUniqueObservationIdsForSetInRowOrder( setId );
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                observationVectorData.getSetMetadata( setId );
        const observation_models::ObservableType currentObservableType = metadata.observableType_;
        const observation_models::LinkEnds currentLinkEnds = observationVectorData.getLinkDefinitionForSet( setId ).linkEnds_;
        const unsigned int observableSize = metadata.observableSize_;
        const int currentObservationSize = static_cast< int >( setObservationIds.size( ) * observableSize );

        if( currentObservationSize > 0 )
        {
            std::vector< TimeType > times;
            times.reserve( setObservationIds.size( ) );
            for( const unsigned int observationId : setObservationIds )
            {
                times.push_back( observationVectorData.getTimes( ).at( observationVectorData.getVectorRow( observationId, 0 ) ) );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
            observationSimulator.at( currentObservableType )
                    ->computeObservations( times,
                                           currentLinkEnds,
                                           metadata.referenceLinkEnd_,
                                           observationVectorData.getAncillarySettingsForSet( setId ),
                                           observationsVector );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualBlock =
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( currentObservationSize );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observedObservationBlock =
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( currentObservationSize );
            for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
            {
                const unsigned int observationId = setObservationIds.at( observationIndex );
                for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                {
                    const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                    const int targetRow = observationVectorData.getVectorRow( observationId, componentIndex );
                    observedObservationBlock( sourceRow ) = observationVectorData.getObservationVector( )( targetRow );
                    residualBlock( sourceRow ) = observedObservationBlock( sourceRow ) - observationsVector( sourceRow );
                }
            }
            if( observation_models::isResidualWrappingRequired( currentObservableType ) )
            {
                wrapObservationResiduals< ObservationScalarType >(
                        residualBlock,
                        std::make_pair( 0, currentObservationSize ),
                        currentObservableType,
                        observedObservationBlock,
                        observationSimulator.at( currentObservableType )->getResidualWrappingSettings( currentLinkEnds ) );
            }
            for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
            {
                const unsigned int observationId = setObservationIds.at( observationIndex );
                for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                {
                    const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                    const int targetRow = observationVectorData.getVectorRow( observationId, componentIndex );
                    residuals( targetRow ) = residualBlock( sourceRow );
                }
            }
        }
    }
}

//! Calculate residuals for the dataset's computation vector data.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationSimulatorBase< ObservationScalarType, TimeType > > >&
                observationSimulator,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    calculateResiduals< ObservationScalarType, TimeType >(
            observationDataset, observationDataset->createComputationObservationVectorData( true ), observationSimulator, residuals );
}

//! Calculate residuals for a legacy observation collection.
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
    const auto dataset = observationsCollection->getObservationDataset( );
    calculateResiduals< ObservationScalarType, TimeType >(
            dataset, dataset->createOrderedObservationVectorData( true ), observationSimulator, residuals );
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
//! Calculate design-matrix rows and residuals using a dataset-derived vector mapping.
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
    if( calculatePartials )
    {
        designMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalNumberParameters );
    }

    if( calculateResiduals )
    {
        residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservationSize, 1 );
    }

    for( const unsigned int setId : observationVectorData.getSetIdsInRowOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = observationVectorData.getUniqueObservationIdsForSetInRowOrder( setId );
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                observationVectorData.getSetMetadata( setId );
        const observation_models::ObservableType currentObservableType = metadata.observableType_;
        const observation_models::LinkEnds currentLinkEnds = observationVectorData.getLinkDefinitionForSet( setId ).linkEnds_;
        const unsigned int observableSize = metadata.observableSize_;
        const int currentObservationSize = static_cast< int >( setObservationIds.size( ) * observableSize );

        if( currentObservationSize > 0 )
        {
            std::vector< TimeType > times;
            times.reserve( setObservationIds.size( ) );
            for( const unsigned int observationId : setObservationIds )
            {
                times.push_back( observationVectorData.getTimes( ).at( observationVectorData.getVectorRow( observationId, 0 ) ) );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
            Eigen::MatrixXd partialsMatrix;
            observationManagers.at( currentObservableType )
                    ->computeObservationsWithPartials( times,
                                                       currentLinkEnds,
                                                       metadata.referenceLinkEnd_,
                                                       observationVectorData.getAncillarySettingsForSet( setId ),
                                                       observationsVector,
                                                       partialsMatrix,
                                                       calculateResiduals,
                                                       calculatePartials );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualBlock;
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observedObservationBlock;
            if( calculateResiduals )
            {
                residualBlock = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( currentObservationSize );
                observedObservationBlock = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( currentObservationSize );
            }

            for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
            {
                const unsigned int observationId = setObservationIds.at( observationIndex );
                for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                {
                    const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                    const int targetRow = observationVectorData.getVectorRow( observationId, componentIndex );
                    if( calculatePartials )
                    {
                        designMatrix.row( targetRow ) = partialsMatrix.row( sourceRow );
                    }
                    if( calculateResiduals )
                    {
                        observedObservationBlock( sourceRow ) = observationVectorData.getObservationVector( )( targetRow );
                        residualBlock( sourceRow ) = observedObservationBlock( sourceRow ) - observationsVector( sourceRow );
                    }
                }
            }

            if( calculateResiduals )
            {
                if( observation_models::isResidualWrappingRequired( currentObservableType ) )
                {
                    wrapObservationResiduals< ObservationScalarType >( residualBlock,
                                                                       std::make_pair( 0, currentObservationSize ),
                                                                       currentObservableType,
                                                                       observedObservationBlock,
                                                                       observationManagers.at( currentObservableType )
                                                                               ->getObservationSimulator( )
                                                                               ->getResidualWrappingSettings( currentLinkEnds ) );
                }
                for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
                {
                    const unsigned int observationId = setObservationIds.at( observationIndex );
                    for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                    {
                        const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                        const int targetRow = observationVectorData.getVectorRow( observationId, componentIndex );
                        residuals( targetRow ) = residualBlock( sourceRow );
                    }
                }
            }
        }
    }
}

//! Calculate design-matrix rows and residuals for an observation dataset.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateDesignMatrixAndResiduals(
        std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
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
    const observation_models::ObservationVectorData< ObservationScalarType, TimeType > observationVectorData =
            observationDataset->createObservationVectorData( !calculatePartials );
    if( observationVectorData.getObservationVector( ).size( ) != totalObservationSize )
    {
        throw std::runtime_error(
                "Error when computing observation partials, requested size is inconsistent with observation vector data size." );
    }
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( observationDataset,
                                                                          observationVectorData,
                                                                          observationManagers,
                                                                          totalNumberParameters,
                                                                          designMatrix,
                                                                          residuals,
                                                                          calculateResiduals,
                                                                          calculatePartials );
}

//! Calculate design-matrix rows and residuals for a legacy observation collection.
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
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( observationsCollection->getObservationDataset( ),
                                                                          observationManagers,
                                                                          totalNumberParameters,
                                                                          totalObservationSize,
                                                                          designMatrix,
                                                                          residuals,
                                                                          calculateResiduals,
                                                                          calculatePartials );
}

//! Calculate a design matrix for an observation dataset.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateDesignMatrix(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalNumberParameters,
        const int totalObservationSize,
        Eigen::MatrixXd& designMatrix )
{
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > dummyVector;
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >(
            observationDataset, observationManagers, totalNumberParameters, totalObservationSize, designMatrix, dummyVector, false, true );
}

//! Calculate a design matrix for a legacy observation collection.
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
    calculateDesignMatrix< ObservationScalarType, TimeType >( observationsCollection->getObservationDataset( ),
                                                              observationManagers,
                                                              totalNumberParameters,
                                                              totalObservationSize,
                                                              designMatrix );
}

//! Calculate residuals with observation managers for an observation dataset.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void calculateResiduals(
        const std::shared_ptr< observation_models::ObservationDataset< ObservationScalarType, TimeType > > observationDataset,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalObservationSize,
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals )
{
    Eigen::MatrixXd dummyMatrix;
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >(
            observationDataset, observationManagers, 0, totalObservationSize, dummyMatrix, residuals, true, false );
}

//! Calculate residuals with observation managers for a legacy collection.
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
    calculateResiduals< ObservationScalarType, TimeType >(
            observationsCollection->getObservationDataset( ), observationManagers, totalObservationSize, residuals );
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
