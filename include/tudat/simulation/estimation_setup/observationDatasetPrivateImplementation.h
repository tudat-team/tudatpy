/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETPRIVATEIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETPRIVATEIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::ObservationDataset( const ObservationDataset& other ):
    std::enable_shared_from_this< ObservationDataset >( ),
    observationRows_( other.observationRows_ ), rowPositionById_( other.rowPositionById_ ),
    nextObservationId_( other.nextObservationId_ ), setMetadata_( other.setMetadata_ ),
    observationIdsBySet_( other.observationIdsBySet_ ), linkDefinitionRegistry_( other.linkDefinitionRegistry_ ),
    observedValues_( other.observedValues_ ), residualValues_( other.residualValues_ ), observationWeights_( other.observationWeights_ )
{
    for( const auto& settings : other.ancillarySettingsRegistry_ )
    {
        ancillarySettingsRegistry_.push_back( settings ? std::make_shared< ObservationAncillarySimulationSettings >( *settings ) : nullptr );
    }
    for( const auto& bookkeeping : other.dependentVariableLayoutRegistry_ )
    {
        dependentVariableLayoutRegistry_.push_back( bookkeeping ? bookkeeping->clone( ) : nullptr );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::rebuildRowIndex( )
{
    std::unordered_map< unsigned int, std::size_t > index;
    index.reserve( observationRows_.size( ) );
    unsigned int nextScalar = 0;
    for( std::size_t i = 0; i < observationRows_.size( ); ++i )
    {
        const auto& row = observationRows_.at( i );
        if( !index.emplace( row.observationId_, i ).second || row.observationId_ >= nextObservationId_ ||
            row.firstScalarComponent_ != nextScalar || row.scalarSize_ == 0 ||
            row.scalarSize_ != setMetadata_.at( row.setId_ ).observableSize_ ||
            observationIdsBySet_.at( row.setId_ ).at( row.indexInSet_ ) != row.observationId_ )
        {
            throw std::runtime_error( "Observation row identity or component mapping is inconsistent." );
        }
        nextScalar += row.scalarSize_;
    }
    std::size_t groupedRows = 0;
    for( const auto& ids : observationIdsBySet_ ) { groupedRows += ids.size( ); }
    if( nextScalar != observedValues_.size( ) || residualValues_.size( ) != observedValues_.size( ) ||
        observationWeights_.size( ) != observedValues_.size( ) || groupedRows != observationRows_.size( ) )
    {
        throw std::runtime_error( "Observation scalar storage or set membership is inconsistent." );
    }
    rowPositionById_ = std::move( index );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::sortObservationIdsForSet( const unsigned int setId )
{
    auto& ids = observationIdsBySet_.at( setId );
    std::stable_sort( ids.begin( ), ids.end( ), [ this ]( unsigned int lhs, unsigned int rhs ) {
        return getObservationRow( lhs ).time_ < getObservationRow( rhs ).time_;
    } );
    for( std::size_t i = 0; i < ids.size( ); ++i ) { mutableObservationRow( ids.at( i ) ).indexInSet_ = i; }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::retainObservationRows(
        const std::vector< unsigned int >& retainedIds )
{
    const std::unordered_set< unsigned int > selected( retainedIds.begin( ), retainedIds.end( ) );
    if( selected.size( ) != retainedIds.size( ) )
    {
        throw std::runtime_error( "Observation selection contains duplicate identities." );
    }
    const auto scalars = getScalarComponentIdsForObservationSelection( retainedIds, {} );
    auto weights = observationWeights_.restricted( scalars );
    std::vector< ObservationDatasetRow< TimeType > > rows;
    std::vector< ObservationScalarType > observations, residuals;
    std::unordered_map< unsigned int, std::size_t > index;
    auto groups = observationIdsBySet_;
    rows.reserve( retainedIds.size( ) );
    index.reserve( retainedIds.size( ) );
    observations.reserve( scalars.size( ) );
    residuals.reserve( scalars.size( ) );
    for( const unsigned int id : retainedIds )
    {
        auto row = getObservationRow( id );
        row.firstScalarComponent_ = observations.size( );
        const auto& original = getObservationRow( id );
        const auto first = original.firstScalarComponent_;
        observations.insert( observations.end( ), observedValues_.begin( ) + first, observedValues_.begin( ) + first + row.scalarSize_ );
        residuals.insert( residuals.end( ), residualValues_.begin( ) + first, residualValues_.begin( ) + first + row.scalarSize_ );
        index.emplace( id, rows.size( ) );
        rows.push_back( std::move( row ) );
    }
    for( auto& ids : groups )
    {
        ids.erase( std::remove_if( ids.begin( ), ids.end( ), [ &selected ]( unsigned int id ) { return !selected.count( id ); } ), ids.end( ) );
        for( std::size_t i = 0; i < ids.size( ); ++i ) { rows.at( index.at( ids.at( i ) ) ).indexInSet_ = i; }
    }
    // Publish all mutually dependent storage only after successful validation and allocation.
    observationRows_ = std::move( rows );
    rowPositionById_ = std::move( index );
    observationIdsBySet_ = std::move( groups );
    observedValues_ = std::move( observations );
    residualValues_ = std::move( residuals );
    observationWeights_ = std::move( weights );
    ++structuralVersion_;
}



template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setResidualVector(
        const FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector )
{
    validateProjection( flattenedObservationData );
    if( residualVector.size( ) != flattenedObservationData.getObservationVector( ).size( ) )
    {
        throw std::runtime_error(
                "Error when setting dataset residual vector from flattened observation data, input size is inconsistent with flattened "
                "data size." );
    }

    for( int i = 0; i < residualVector.size( ); ++i )
    {
        residualValues_.at( flattenedObservationData.getScalarComponentIds( ).at( i ) ) = residualVector( i );
    }
}













template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::validateObservationSetData(
        const unsigned int setId,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals ) const
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    if( observations.size( ) != times.size( ) )
    {
        throw std::runtime_error( "Error when updating observation dataset, observation and time sizes are inconsistent." );
    }
    if( !dependentVariables.empty( ) && dependentVariables.size( ) != observations.size( ) )
    {
        throw std::runtime_error( "Error when updating observation dataset, dependent variable size is inconsistent." );
    }
    if( !weights.empty( ) && weights.size( ) != observations.size( ) )
    {
        throw std::runtime_error( "Error when updating observation dataset, weight size is inconsistent." );
    }
    if( !residuals.empty( ) && residuals.size( ) != observations.size( ) )
    {
        throw std::runtime_error( "Error when updating observation dataset, residual size is inconsistent." );
    }

    for( const auto& weight : weights )
    {
        ObservationWeights::validateDiagonal( weight );
    }

    int dependentVariableSize = -1;
    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& bookkeeping =
            getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
    if( bookkeeping != nullptr )
    {
        dependentVariableSize = bookkeeping->getTotalDependentVariableSize( );
    }
    else if( !dependentVariables.empty( ) )
    {
        dependentVariableSize = dependentVariables.front( ).size( );
    }

    for( std::size_t i = 0; i < observations.size( ); ++i )
    {
        if( observations.at( i ).size( ) != static_cast< int >( observableSize ) ||
            ( !weights.empty( ) && weights.at( i ).size( ) != static_cast< int >( observableSize ) ) ||
            ( !residuals.empty( ) && residuals.at( i ).size( ) != static_cast< int >( observableSize ) ) ||
            ( !dependentVariables.empty( ) && dependentVariables.at( i ).size( ) != dependentVariableSize ) )
        {
            throw std::runtime_error( "Error when updating observation dataset, scalar component size is inconsistent." );
        }
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::size_t > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTimeSortingPermutation(
        const std::vector< TimeType >& observationTimes )
{
    std::vector< std::size_t > permutation( observationTimes.size( ) );
    for( std::size_t i = 0; i < observationTimes.size( ); ++i )
    {
        permutation.at( i ) = i;
    }

    std::stable_sort( permutation.begin( ), permutation.end( ), [ &observationTimes ]( const std::size_t i, const std::size_t j ) {
        return observationTimes.at( i ) < observationTimes.at( j );
    } );

    return permutation;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
template< typename VectorType >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::reorderVector( std::vector< VectorType >& data,
                                                                                  const std::vector< std::size_t >& permutation )
{
    std::vector< VectorType > reorderedData( data.size( ) );
    for( std::size_t i = 0; i < data.size( ); ++i )
    {
        reorderedData.at( i ) = data.at( permutation.at( i ) );
    }
    data.swap( reorderedData );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
template< typename VectorType >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::removeEntries( std::vector< VectorType >& data,
                                                                                  const std::vector< unsigned int >& indicesToRemove )
{
    for( std::vector< unsigned int >::const_reverse_iterator indexIterator = indicesToRemove.rbegin( );
         indexIterator != indicesToRemove.rend( );
         ++indexIterator )
    {
        data.erase( data.begin( ) + *indexIterator );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationDataset< ObservationScalarType, TimeType, Dummy >::createSetVector(
        const unsigned int setId,
        const std::vector< ObservationScalarType >& scalarValues ) const
{
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > vector =
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( getTotalScalarSizeForSet( setId ) );
    std::size_t currentIndex = 0;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        const ObservationDatasetRow< TimeType >& row = getObservationRow( observationId );
        for( unsigned int i = 0; i < row.scalarSize_; ++i )
        {
            vector( currentIndex++ ) = scalarValues.at( row.firstScalarComponent_ + i );
        }
    }
    return vector;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setObservationValue(
        const unsigned int observationId,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation )
{
    const ObservationDatasetRow< TimeType >& row = getObservationRow( observationId );
    if( observation.size( ) != static_cast< int >( row.scalarSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset observation, scalar size is inconsistent." );
    }
    for( unsigned int i = 0; i < row.scalarSize_; ++i )
    {
        observedValues_.at( row.firstScalarComponent_ + i ) = observation( i );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setResidualValue(
        const unsigned int observationId,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residual )
{
    const ObservationDatasetRow< TimeType >& row = getObservationRow( observationId );
    if( residual.size( ) != static_cast< int >( row.scalarSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset residual, scalar size is inconsistent." );
    }
    for( unsigned int i = 0; i < row.scalarSize_; ++i )
    {
        residualValues_.at( row.firstScalarComponent_ + i ) = residual( i );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setWeightValue( const unsigned int observationId,
                                                                                   const Eigen::VectorXd& weight )
{
    const auto& row = getObservationRow( observationId );
    if( weight.size( ) != static_cast< Eigen::Index >( row.scalarSize_ ) )
    {
        throw std::runtime_error( "Observation weight component size is inconsistent." );
    }
    observationWeights_.setObservationDiagonal( row.firstScalarComponent_, weight );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponentIdsForObservationSelection(
        const std::vector< unsigned int >& observationIds,
        const std::vector< unsigned int >& components ) const
{
    std::vector< unsigned int > scalarComponentIds;
    for( const unsigned int observationId : observationIds )
    {
        const ObservationDatasetRow< TimeType >& row = getObservationRow( observationId );
        if( components.empty( ) )
        {
            for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
            {
                scalarComponentIds.push_back( row.firstScalarComponent_ + componentIndex );
            }
        }
        else
        {
            for( const unsigned int componentIndex : components )
            {
                if( componentIndex >= row.scalarSize_ )
                {
                    throw std::runtime_error(
                            "Error when setting dataset weight block, selected component index is inconsistent with observation "
                            "size." );
                }
                scalarComponentIds.push_back( row.firstScalarComponent_ + componentIndex );
            }
        }
    }
    return scalarComponentIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
int ObservationDataset< ObservationScalarType, TimeType, Dummy >::registerLinkDefinition( const LinkDefinition& linkDefinition )
{
    for( std::size_t i = 0; i < linkDefinitionRegistry_.size( ); ++i )
    {
        if( linkDefinitionRegistry_.at( i ) == linkDefinition )
        {
            return static_cast< int >( i );
        }
    }

    linkDefinitionRegistry_.push_back( linkDefinition );
    return static_cast< int >( linkDefinitionRegistry_.size( ) - 1 );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
int ObservationDataset< ObservationScalarType, TimeType, Dummy >::registerAncillarySettings(
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
{
    for( std::size_t i = 0; i < ancillarySettingsRegistry_.size( ); ++i )
    {
        if( ancillarySettingsRegistry_.at( i ) == ancillarySettings )
        {
            return static_cast< int >( i );
        }
    }

    ancillarySettingsRegistry_.push_back( ancillarySettings );
    return static_cast< int >( ancillarySettingsRegistry_.size( ) - 1 );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
int ObservationDataset< ObservationScalarType, TimeType, Dummy >::registerDependentVariableLayout(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& bookkeeping )
{
    for( std::size_t i = 0; i < dependentVariableLayoutRegistry_.size( ); ++i )
    {
        if( dependentVariableLayoutRegistry_.at( i ) == bookkeeping )
        {
            return static_cast< int >( i );
        }
    }

    dependentVariableLayoutRegistry_.push_back( bookkeeping );
    return static_cast< int >( dependentVariableLayoutRegistry_.size( ) - 1 );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETPRIVATEIMPLEMENTATION_H
