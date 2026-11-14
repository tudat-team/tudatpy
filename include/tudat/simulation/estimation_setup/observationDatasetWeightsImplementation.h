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
#include <numeric>

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::promoteWeightStructureForEntry( const unsigned int firstScalar,
                                                                                                   const unsigned int secondScalar )
{
    const auto& first = getObservationRow( getScalarComponentRow( firstScalar ).observationId_ );
    const auto& second = getObservationRow( getScalarComponentRow( secondScalar ).observationId_ );
    const auto structure = first.setId_ != second.setId_
            ? ObservationWeightStructure::inter_set_weights
            : ( first.observationId_ != second.observationId_ ? ObservationWeightStructure::per_set
                                                              : ObservationWeightStructure::per_observation );
    for( const unsigned int setId : { first.setId_, second.setId_ } )
    {
        auto& current = setMetadata_.at( setId ).weightStructure_;
        current = std::max( current, structure );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::refreshWeightStructures( )
{
    // A former cross-set block may now be wholly within one set. Keep its within-set size conservative.
    for( auto& metadata : setMetadata_ )
    {
        if( metadata.weightStructure_ == ObservationWeightStructure::inter_set_weights )
        {
            metadata.weightStructure_ = ObservationWeightStructure::per_set;
        }
    }
    for( const auto& entry : observationWeights_.getOffDiagonalEntries( ) )
    {
        promoteWeightStructureForEntry( entry.first.first, entry.first.second );
    }
}

// Store one full set-level weight matrix for a single observation set.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightMatrixForSet( const unsigned int setId,
                                                                                          const Eigen::MatrixXd& weightMatrix )
{
    const auto indices = getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} );
    completeWeightData_.reset( );
    observationWeights_.setBlock( indices, indices, weightMatrix );
    auto& structure = setMetadata_.at( setId ).weightStructure_;
    structure = std::max( structure, ObservationWeightStructure::per_set );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightVectorForSet( const unsigned int setId,
                                                                                          const Eigen::VectorXd& weightVector )
{
    completeWeightData_.reset( );
    observationWeights_.setDiagonal( getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ), weightVector );
    auto& structure = setMetadata_.at( setId ).weightStructure_;
    if( structure != ObservationWeightStructure::inter_set_weights )
    {
        structure = ObservationWeightStructure::diagonal;
    }
}

// Store a dense observable-size block for one observation row.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightMatrixForObservation( const unsigned int observationId,
                                                                                                  const Eigen::MatrixXd& weightMatrix )
{
    const auto indices = getScalarComponentIdsForObservationSelection( { observationId }, {} );
    completeWeightData_.reset( );
    observationWeights_.setBlock( indices, indices, weightMatrix );
    auto& structure = setMetadata_.at( getObservationRow( observationId ).setId_ ).weightStructure_;
    structure = std::max( structure, ObservationWeightStructure::per_observation );
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
    const auto rows = getScalarComponentIdsForObservationSelection( rowObservationIds, rowComponents );
    const auto columns = getScalarComponentIdsForObservationSelection( columnObservationIds, columnComponents );
    completeWeightData_.reset( );
    observationWeights_.setBlock( rows, columns, weightBlock );
    // Only update metadata after validation. An existing cross-set connection may have been removed.
    if( std::any_of( rowObservationIds.begin( ), rowObservationIds.end( ), [ this ]( unsigned int id ) {
            return setMetadata_.at( getObservationRow( id ).setId_ ).weightStructure_ == ObservationWeightStructure::inter_set_weights;
        } ) )
    {
        refreshWeightStructures( );
        return;
    }
    const auto& entries = observationWeights_.getOffDiagonalEntries( );
    for( const unsigned int row : rows )
    {
        for( const unsigned int column : columns )
        {
            if( entries.count( std::minmax( row, column ) ) )
            {
                promoteWeightStructureForEntry( row, column );
            }
        }
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
    ObservationWeights::validateDiagonal( Eigen::VectorXd::Constant( 1, weight ) );
    // Scalars are dimension-independent, so matching ids can be mutated directly.
    for( const unsigned int observationId : getObservationIdsMatchingCondition( condition ) )
    {
        setWeightValue( observationId, Eigen::VectorXd::Constant( getObservationRow( observationId ).scalarSize_, weight ) );
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
    ObservationWeights::validateDiagonal( weight );
    const std::vector< unsigned int > observationIds = getObservationIdsMatchingCondition( condition );

    // Prevalidate every selected row so a mixed-size selection cannot be partly mutated.
    for( const unsigned int observationId : observationIds )
    {
        if( weight.size( ) != static_cast< int >( getObservationRow( observationId ).scalarSize_ ) )
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
    ObservationWeights validation;
    validation.appendDiagonal( Eigen::VectorXd::Ones( weight.rows( ) ) );
    std::vector< unsigned int > indices( weight.rows( ) );
    std::iota( indices.begin( ), indices.end( ), 0 );
    validation.setBlock( indices, indices, weight );
    const std::vector< unsigned int > observationIds = getObservationIdsMatchingCondition( condition );

    // Prevalidate dimensions before writing any row-specific block.
    for( const unsigned int observationId : observationIds )
    {
        if( weight.rows( ) != static_cast< int >( getObservationRow( observationId ).scalarSize_ ) ||
            weight.cols( ) != static_cast< int >( getObservationRow( observationId ).scalarSize_ ) ||
            !weight.isApprox( weight.transpose( ) ) )
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

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETWEIGHTSIMPLEMENTATION_H
