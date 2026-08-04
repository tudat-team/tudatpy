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
#include <iostream>

#include "tudat/basics/timeType.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/astro/observation_models/observationManager.h"
#include "tudat/astro/orbit_determination/podInputOutputTypes.h"
#include "tudat/astro/propagators/propagateCovariance.h"

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType >
void checkObservationResidualDiscontinuities( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residuals,
                                              const std::pair< int, int > observableStartAndSize,
                                              const observation_models::ObservableType observableType )
{
    if( observableType == observation_models::angular_position || observableType == observation_models::euler_angle_313_observable ||
        observableType == observation_models::relative_angular_position )
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualsBlock =
                residuals.block( observableStartAndSize.first, 0, observableStartAndSize.second, 1 );
        for( int i = 1; i < residualsBlock.rows( ); i++ )
        {
            if( std::fabs( residualsBlock( i, 0 ) - residualsBlock( i - 1, 0 ) ) > 6.0 )
            {
                if( residualsBlock( i, 0 ) > 0 )
                {
                    residualsBlock( i, 0 ) = residualsBlock( i, 0 ) - 2.0 * mathematical_constants::PI;
                }
                else
                {
                    residualsBlock( i, 0 ) = residualsBlock( i, 0 ) + 2.0 * mathematical_constants::PI;
                }
            }
            else if( std::fabs( residualsBlock( i, 0 ) - residualsBlock( i - 1, 0 ) ) > 3.0 )
            {
                std::cerr << "Warning, detected jump in observation residual of size "
                          << std::fabs( residualsBlock( i, 0 ) - residualsBlock( i - 1, 0 ) ) << " for observable type " << observableType
                          << std::endl;
            }
        }
        residuals.block( observableStartAndSize.first, 0, observableStartAndSize.second, 1 ) = residualsBlock;
    }
}

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
    const observation_models::FlattenedObservationData< ObservationScalarType, TimeType > flattenedObservationData =
            observationDataset->createComputationFlattenedObservationData( true );
    residuals =
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( flattenedObservationData.getObservationVector( ).size( ), 1 );

    for( const unsigned int setId : flattenedObservationData.getSetIdsInRowOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = flattenedObservationData.getUniqueObservationIdsForSetInRowOrder( setId );
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                observationDataset->getObservationSetMetadata( setId );
        const observation_models::ObservableType currentObservableType = metadata.observableType_;
        const observation_models::LinkEnds currentLinkEnds = observationDataset->getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const unsigned int observableSize = metadata.observableSize_;
        const int currentObservationSize = static_cast< int >( setObservationIds.size( ) * observableSize );

        if( currentObservationSize > 0 )
        {
            std::vector< TimeType > times;
            times.reserve( setObservationIds.size( ) );
            for( const unsigned int observationId : setObservationIds )
            {
                times.push_back( observationDataset->getObservationTime( observationId ) );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
            observationSimulator.at( currentObservableType )
                    ->computeObservations( times,
                                           currentLinkEnds,
                                           metadata.referenceLinkEnd_,
                                           observationDataset->getAncillarySettings( metadata.ancillarySettingsId_ ),
                                           observationsVector );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualBlock =
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( currentObservationSize );
            for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
            {
                const unsigned int observationId = setObservationIds.at( observationIndex );
                for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                {
                    const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                    const int targetRow = flattenedObservationData.getFlattenedRow( observationId, componentIndex );
                    residualBlock( sourceRow ) =
                            flattenedObservationData.getObservationVector( )( targetRow ) - observationsVector( sourceRow );
                }
            }
            checkObservationResidualDiscontinuities< ObservationScalarType >(
                    residualBlock, std::make_pair( 0, currentObservationSize ), currentObservableType );
            for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
            {
                const unsigned int observationId = setObservationIds.at( observationIndex );
                for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                {
                    const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                    const int targetRow = flattenedObservationData.getFlattenedRow( observationId, componentIndex );
                    residuals( targetRow ) = residualBlock( sourceRow );
                }
            }
        }
    }
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
    calculateResiduals< ObservationScalarType, TimeType >(
            observationsCollection->getObservationDataset( ), observationSimulator, residuals );
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
void retrieveObservationWeights( const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData,
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
        const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& computationData,
        const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& estimationData,
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
        const int computationRow = computationData.getFirstFlattenedRowForObservation( observationId );

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
        const observation_models::FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData,
        const std::map< observation_models::ObservableType,
                        std::shared_ptr< observation_models::ObservationManagerBase< ObservationScalarType, TimeType > > >&
                observationManagers,
        const int totalNumberParameters,
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

    const int totalObservationSize = static_cast< int >( flattenedObservationData.getObservationVector( ).size( ) );
    if( calculatePartials )
    {
        designMatrix = Eigen::MatrixXd::Zero( totalObservationSize, totalNumberParameters );
    }

    if( calculateResiduals )
    {
        residuals = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( totalObservationSize, 1 );
    }

    for( const unsigned int setId : flattenedObservationData.getSetIdsInRowOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = flattenedObservationData.getUniqueObservationIdsForSetInRowOrder( setId );
        const observation_models::ObservationSetMetadata< ObservationScalarType, TimeType >& metadata =
                observationDataset->getObservationSetMetadata( setId );
        const observation_models::ObservableType currentObservableType = metadata.observableType_;
        const observation_models::LinkEnds currentLinkEnds = observationDataset->getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const unsigned int observableSize = metadata.observableSize_;
        const int currentObservationSize = static_cast< int >( setObservationIds.size( ) * observableSize );

        if( currentObservationSize > 0 )
        {
            std::vector< TimeType > times;
            times.reserve( setObservationIds.size( ) );
            for( const unsigned int observationId : setObservationIds )
            {
                times.push_back( observationDataset->getObservationTime( observationId ) );
            }

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observationsVector;
            Eigen::MatrixXd partialsMatrix;
            observationManagers.at( currentObservableType )
                    ->computeObservationsWithPartials( times,
                                                       currentLinkEnds,
                                                       metadata.referenceLinkEnd_,
                                                       observationDataset->getAncillarySettings( metadata.ancillarySettingsId_ ),
                                                       observationsVector,
                                                       partialsMatrix,
                                                       calculateResiduals,
                                                       calculatePartials );

            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residualBlock;
            if( calculateResiduals )
            {
                residualBlock = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( currentObservationSize );
            }

            for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
            {
                const unsigned int observationId = setObservationIds.at( observationIndex );
                for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                {
                    const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                    const int targetRow = flattenedObservationData.getFlattenedRow( observationId, componentIndex );
                    if( calculatePartials )
                    {
                        designMatrix.row( targetRow ) = partialsMatrix.row( sourceRow );
                    }
                    if( calculateResiduals )
                    {
                        residualBlock( sourceRow ) =
                                flattenedObservationData.getObservationVector( )( targetRow ) - observationsVector( sourceRow );
                    }
                }
            }

            if( calculateResiduals )
            {
                checkObservationResidualDiscontinuities< ObservationScalarType >(
                        residualBlock, std::make_pair( 0, currentObservationSize ), currentObservableType );
                for( std::size_t observationIndex = 0; observationIndex < setObservationIds.size( ); ++observationIndex )
                {
                    const unsigned int observationId = setObservationIds.at( observationIndex );
                    for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
                    {
                        const int sourceRow = static_cast< int >( observationIndex * observableSize + componentIndex );
                        const int targetRow = flattenedObservationData.getFlattenedRow( observationId, componentIndex );
                        residuals( targetRow ) = residualBlock( sourceRow );
                    }
                }
            }
        }
    }
}

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
    const observation_models::FlattenedObservationData< ObservationScalarType, TimeType > flattenedObservationData =
            observationDataset->createOrderedFlattenedObservationData( true );
    if( flattenedObservationData.getObservationVector( ).size( ) != totalObservationSize )
    {
        throw std::runtime_error( "Error when computing observation partials, requested size is inconsistent with flattened data size." );
    }
    calculateDesignMatrixAndResiduals< ObservationScalarType, TimeType >( observationDataset,
                                                                          flattenedObservationData,
                                                                          observationManagers,
                                                                          totalNumberParameters,
                                                                          designMatrix,
                                                                          residuals,
                                                                          calculateResiduals,
                                                                          calculatePartials );
}

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
