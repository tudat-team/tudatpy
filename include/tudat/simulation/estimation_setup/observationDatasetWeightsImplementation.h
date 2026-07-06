/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

// Apply one compact scalar weight to every observation row in one set.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantSingleObservationScalarWeightForSet( const unsigned int setId,
                                                                                                                   const double weight )
{
    // The compact scalar representation is valid for any observable size.
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        observationWeights_.setScalarWeight( observationId, weight );
    }
}

// Apply one diagonal component-weight vector to every observation row in one set.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantSingleObservationDiagonalWeightForSet(
        const unsigned int setId,
        const Eigen::VectorXd& weight )
{
    // Validate once against the observable size before mutating row weights.
    if( weight.size( ) != static_cast< int >( getObservationSetMetadata( setId ).observableSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset weights, weight size is inconsistent." );
    }

    // Store the same diagonal vector on each observation row in the set.
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        setWeightValue( observationId, weight );
    }
}

// Store one full set-level weight matrix for a single observation set.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightMatrixForSet( const unsigned int setId,
                                                                                          const Eigen::MatrixXd& weightMatrix )
{
    // A set-level block spans all scalar components in the selected set.
    if( weightMatrix.rows( ) != static_cast< int >( getTotalScalarSizeForSet( setId ) ) ||
        weightMatrix.cols( ) != static_cast< int >( getTotalScalarSizeForSet( setId ) ) )
    {
        throw std::runtime_error( "Error when setting dataset set weight matrix, matrix size is inconsistent." );
    }
    observationWeights_.setSetWeightBlock( setId, weightMatrix );
}

// Apply one dense observable-size block to every observation row in one set.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantSingleObservationMatrixWeightForSet(
        const unsigned int setId,
        const Eigen::MatrixXd& weight )
{
    // Validate against the observable dimension once before applying the block repeatedly.
    if( weight.rows( ) != static_cast< int >( getObservationSetMetadata( setId ).observableSize_ ) ||
        weight.cols( ) != static_cast< int >( getObservationSetMetadata( setId ).observableSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset observation weight matrix, matrix size is inconsistent." );
    }

    // Store the same dense per-observation block on every row in the set.
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        setWeightMatrixForObservation( observationId, weight );
    }
}

// Report whether a set has an explicit full set-level weight block.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
bool ObservationDataset< ObservationScalarType, TimeType, Dummy >::hasWeightMatrixForSet( const unsigned int setId ) const
{
    return observationWeights_.hasSetWeightBlock( setId );
}

// Store a dense observable-size block for one observation row.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightMatrixForObservation( const unsigned int observationId,
                                                                                                  const Eigen::MatrixXd& weightMatrix )
{
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );

    // The row scalar size is the observable dimension for this single event.
    if( weightMatrix.rows( ) != static_cast< int >( row.scalarSize_ ) || weightMatrix.cols( ) != static_cast< int >( row.scalarSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset observation weight matrix, matrix size is inconsistent." );
    }
    observationWeights_.setWeightBlock( observationId, weightMatrix );
}

// Report whether one observation row has an explicit dense weight block.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
bool ObservationDataset< ObservationScalarType, TimeType, Dummy >::hasWeightMatrixForObservation( const unsigned int observationId ) const
{
    return observationWeights_.hasObservationWeightBlock( observationId );
}

// Add an already scalar-component-indexed off-diagonal weight block.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::addExtraWeightBlock( const ObservationWeightBlock& weightBlock )
{
    observationWeights_.addExtraWeightBlock( weightBlock );
}

// Store a sparse/dense block between selected scalar components of selected observation rows.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightBlock( const std::vector< unsigned int >& rowObservationIds,
                                                                                   const std::vector< unsigned int >& columnObservationIds,
                                                                                   const Eigen::MatrixXd& weightBlock,
                                                                                   const std::vector< unsigned int >& rowComponents,
                                                                                   const std::vector< unsigned int >& columnComponents )
{
    // Resolve observation ids plus optional component selections to scalar-component ids.
    const std::vector< unsigned int > rowScalarComponentIds =
            getScalarComponentIdsForObservationSelection( rowObservationIds, rowComponents );
    const std::vector< unsigned int > columnScalarComponentIds =
            getScalarComponentIdsForObservationSelection( columnObservationIds, columnComponents );

    // Validate the full request before adding any block, keeping this call exception-safe.
    if( weightBlock.rows( ) != static_cast< int >( rowScalarComponentIds.size( ) ) ||
        weightBlock.cols( ) != static_cast< int >( columnScalarComponentIds.size( ) ) )
    {
        throw std::runtime_error( "Error when setting dataset weight block, matrix size is inconsistent with selected observations." );
    }

    if( rowScalarComponentIds == columnScalarComponentIds )
    {
        if( weightBlock.rows( ) != weightBlock.cols( ) || !weightBlock.isApprox( weightBlock.transpose( ) ) )
        {
            throw std::runtime_error(
                    "Error when setting symmetric dataset weight block, block with identical row and column selection is not symmetric." );
        }
    }

    // Store the requested block and its transpose for distinct row/column selections.
    ObservationWeightBlock datasetWeightBlock;
    datasetWeightBlock.rowScalarComponentIds_ = rowScalarComponentIds;
    datasetWeightBlock.columnScalarComponentIds_ = columnScalarComponentIds;
    datasetWeightBlock.weightBlock_ = weightBlock;
    addExtraWeightBlock( datasetWeightBlock );

    if( rowScalarComponentIds != columnScalarComponentIds )
    {
        ObservationWeightBlock transposedWeightBlock;
        transposedWeightBlock.rowScalarComponentIds_ = columnScalarComponentIds;
        transposedWeightBlock.columnScalarComponentIds_ = rowScalarComponentIds;
        transposedWeightBlock.weightBlock_ = weightBlock.transpose( );
        addExtraWeightBlock( transposedWeightBlock );
    }
}

// Apply one compact scalar weight to all observation rows matching a condition.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantSingleObservationScalarWeight(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const double weight )
{
    // Scalars are dimension-independent, so matching ids can be mutated directly.
    for( const unsigned int observationId : getObservationIdsMatchingCondition( condition ) )
    {
        observationWeights_.setScalarWeight( observationId, weight );
    }
}

// Apply one diagonal component-weight vector to all observation rows matching a condition.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantSingleObservationDiagonalWeight(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const Eigen::VectorXd& weight )
{
    const std::vector< unsigned int > observationIds = getObservationIdsMatchingCondition( condition );

    // Prevalidate every selected row so a mixed-size selection cannot be partly mutated.
    for( const unsigned int observationId : observationIds )
    {
        if( weight.size( ) != static_cast< int >( observationRows_.at( observationId ).scalarSize_ ) )
        {
            throw std::runtime_error( "Error when setting dataset weights by condition, weight size is inconsistent." );
        }
    }

    // Store the same diagonal vector on every validated row.
    for( const unsigned int observationId : observationIds )
    {
        setWeightValue( observationId, weight );
    }
}

// Apply one dense observable-size block to all observation rows matching a condition.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setConstantSingleObservationMatrixWeight(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const Eigen::MatrixXd& weight )
{
    const std::vector< unsigned int > observationIds = getObservationIdsMatchingCondition( condition );

    // Prevalidate dimensions before writing any row-specific block.
    for( const unsigned int observationId : observationIds )
    {
        if( weight.rows( ) != static_cast< int >( observationRows_.at( observationId ).scalarSize_ ) ||
            weight.cols( ) != static_cast< int >( observationRows_.at( observationId ).scalarSize_ ) )
        {
            throw std::runtime_error( "Error when setting dataset weight matrices by condition, matrix size is inconsistent." );
        }
    }

    // Store the same dense block on every validated row.
    for( const unsigned int observationId : observationIds )
    {
        setWeightMatrixForObservation( observationId, weight );
    }
}

// Return all explicitly stored scalar-component-indexed off-diagonal blocks.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::vector< ObservationWeightBlock >& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getExtraWeightBlocks( ) const
{
    return observationWeights_.getExtraWeightBlocks( );
}

// Report whether the dataset contains any scalar-component-indexed off-diagonal blocks.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
bool ObservationDataset< ObservationScalarType, TimeType, Dummy >::hasExtraWeightBlocks( ) const
{
    return observationWeights_.hasExtraWeightBlocks( );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H
