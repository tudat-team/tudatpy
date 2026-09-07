/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETFLATTENINGIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETFLATTENINGIMPLEMENTATION_H

#include <iostream>
#include <string>

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
FlattenedObservationData< ObservationScalarType, TimeType >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createEstimationFlattenedObservationData( const bool includeRejected ) const
{
    return createEstimationProjection( includeRejected );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
FlattenedObservationData< ObservationScalarType, TimeType >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createComputationFlattenedObservationData( const bool includeRejected ) const
{
    return createFlattenedObservationDataFromObservationIds( getAllObservationIds( ), includeRejected );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
FlattenedObservationData< ObservationScalarType, TimeType >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createOrderedFlattenedObservationData( const bool includeInactive ) const
{
    return createFlattenedObservationDataFromObservationIds( getObservationIdsInOrderedFlattenedDataOrder( ), includeInactive );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSetIdsInOrderedFlattenedDataOrder( ) const
{
    std::map< ObservableType, std::map< LinkEnds, std::vector< unsigned int > > > setIdsByObservableAndLinkEnds;
    for( unsigned int setId = 0; setId < setMetadata_.size( ); ++setId )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = setMetadata_.at( setId );
        setIdsByObservableAndLinkEnds[ metadata.observableType_ ][ linkDefinitionRegistry_.at( metadata.linkDefinitionId_ ).linkEnds_ ]
                .push_back( setId );
    }

    std::vector< unsigned int > setIds;
    setIds.reserve( setMetadata_.size( ) );
    for( const auto& observableIterator : setIdsByObservableAndLinkEnds )
    {
        for( const auto& linkEndsIterator : observableIterator.second )
        {
            setIds.insert( setIds.end( ), linkEndsIterator.second.begin( ), linkEndsIterator.second.end( ) );
        }
    }
    return setIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getStructuralVersion( ) const
{
    return structuralVersion_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::resetLinkDefinitionForSet( const unsigned int setId,
                                                                                              const LinkDefinition& linkDefinition )
{
    setMetadata_.at( setId ).linkDefinitionId_ = registerLinkDefinition( linkDefinition );
    ++structuralVersion_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setLinkEndReferencePoint(
        const std::string& bodyName,
        const std::string& referencePointName,
        const LinkEndType linkEndType,
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    for( unsigned int setId = 0; setId < getNumberOfObservationSets( ); ++setId )
    {
        const std::vector< unsigned int >& observationIds = getObservationIdsForSet( setId );
        bool setMatchesCondition = false;
        for( const unsigned int observationId : observationIds )
        {
            if( condition( *this, observationId ) )
            {
                setMatchesCondition = true;
                break;
            }
        }

        if( !setMatchesCondition )
        {
            continue;
        }

        std::map< LinkEndType, LinkEndId > linkEnds = getLinkDefinition( setMetadata_.at( setId ).linkDefinitionId_ ).linkEnds_;
        typename std::map< LinkEndType, LinkEndId >::iterator linkEndIterator = linkEnds.find( linkEndType );
        if( linkEndIterator == linkEnds.end( ) || linkEndIterator->second.bodyName_ != bodyName ||
            linkEndIterator->second.getReferencePointName( ) == referencePointName )
        {
            continue;
        }

        linkEndIterator->second = LinkEndId( linkEndIterator->second.bodyName_, referencePointName );
        resetLinkDefinitionForSet( setId, LinkDefinition( linkEnds ) );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::resetDependentVariableBookkeepingForSet(
        const unsigned int setId,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping )
{
    const auto& metadata = getObservationSetMetadata( setId );
    if( dependentVariableBookkeeping )
    {
        if( dependentVariableBookkeeping->getObservableType( ) != metadata.observableType_ ||
            !( dependentVariableBookkeeping->getLinkEnds( ) == getLinkDefinition( metadata.linkDefinitionId_ ) ) )
        {
            throw std::runtime_error( "Dependent-variable bookkeeping is incompatible with the observation set." );
        }
        for( const auto id : getObservationIdsForSet( setId ) )
        {
            const auto size = getDependentVariables( id ).size( );
            if( size != 0 && size != dependentVariableBookkeeping->getTotalDependentVariableSize( ) )
            {
                throw std::runtime_error( "Dependent-variable layout does not match the stored values." );
            }
        }
    }
    setMetadata_.at( setId ).dependentVariableLayoutId_ = registerDependentVariableLayout( dependentVariableBookkeeping );
    ++structuralVersion_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAllObservationIds( ) const
{
    std::vector< unsigned int > result;
    result.reserve( observationRows_.size( ) );
    for( const auto& row : observationRows_ )
    {
        result.push_back( row.observationId_ );
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationIdsInOrderedFlattenedDataOrder( )
        const
{
    std::vector< unsigned int > observationIds;
    observationIds.reserve( observationRows_.size( ) );
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = observationIdsBySet_.at( setId );
        observationIds.insert( observationIds.end( ), setObservationIds.begin( ), setObservationIds.end( ) );
    }
    return observationIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
FlattenedObservationData< ObservationScalarType, TimeType >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createFlattenedObservationDataFromObservationIds(
        const std::vector< unsigned int >& selectedObservationIds,
        const bool includeInactive ) const
{
    FlattenedObservationData< ObservationScalarType, TimeType > result;
    result.source_ = getLifetimeToken( );
    result.structuralVersion_ = structuralVersion_;
    result.projectionVersion_ = projectionVersion_;
    result.uniqueObservationIdsBySet_.resize( setMetadata_.size( ) );
    std::vector< unsigned int > selected;
    for( const unsigned int id : selectedObservationIds )
    {
        const auto& row = getObservationRow( id );
        if( includeInactive || row.isActive_ )
        {
            selected.push_back( id );
        }
    }
    result.scalarComponentIds_ = getScalarComponentIdsForObservationSelection( selected, {} );
    const auto weights = observationWeights_.restricted( result.scalarComponentIds_ );
    result.weights_ = weights.diagonalVector( );
    result.isDiagonalWeightOnly_ = !weights.hasOffDiagonalWeights( );
    if( weights.hasOffDiagonalWeights( ) )
    {
        result.weightMatrix_ = weights.sparseMatrix( );
    }
    const auto size = result.scalarComponentIds_.size( );
    result.observations_.resize( size );
    result.residuals_.resize( size );
    result.times_.reserve( size );
    result.observationIds_.reserve( size );
    result.setIds_.reserve( size );
    result.rowMapping_.reserve( selected.size( ) );
    unsigned int first = 0;
    for( const unsigned int id : selected )
    {
        const auto& row = getObservationRow( id );
        result.rowMapping_.emplace( id, std::make_pair( first, row.scalarSize_ ) );
        auto& group = result.uniqueObservationIdsBySet_.at( row.setId_ );
        if( group.empty( ) )
        {
            result.setIdsInRowOrder_.push_back( row.setId_ );
            const auto& metadata = getObservationSetMetadata( row.setId_ );
            result.metadataBySet_.emplace( row.setId_, metadata );
            result.linksBySet_.emplace( row.setId_, getLinkDefinition( metadata.linkDefinitionId_ ) );
            const auto ancillary = getAncillarySettings( metadata.ancillarySettingsId_ );
            result.ancillaryBySet_.emplace(
                    row.setId_, ancillary ? std::make_shared< ObservationAncillarySimulationSettings >( *ancillary ) : nullptr );
        }
        group.push_back( id );
        result.dependentVariables_.emplace( id, row.dependentVariableValues_ );
        for( unsigned int component = 0; component < row.scalarSize_; ++component, ++first )
        {
            const unsigned int scalar = row.firstScalarComponent_ + component;
            result.observations_( first ) = observedValues_.at( scalar );
            result.residuals_( first ) = residualValues_.at( scalar );
            result.times_.push_back( row.time_ );
            result.observationIds_.push_back( id );
            result.setIds_.push_back( row.setId_ );
        }
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::pair< int, int > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetStartAndSizeInDatasetOrder( ) const
{
    std::vector< std::pair< int, int > > startAndSize;
    startAndSize.reserve( getNumberOfObservationSets( ) );

    int currentIndex = 0;
    for( unsigned int setId = 0; setId < getNumberOfObservationSets( ); ++setId )
    {
        const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
        startAndSize.push_back( std::make_pair( currentIndex, currentSize ) );
        currentIndex += currentSize;
    }
    return startAndSize;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetIdsForObservableType(
        const ObservableType observableType ) const
{
    std::vector< unsigned int > setIds;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        if( getObservationSetMetadata( setId ).observableType_ == observableType )
        {
            setIds.push_back( setId );
        }
    }
    return setIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTotalScalarSizeForObservableType(
        const ObservableType observableType ) const
{
    std::size_t totalSize = 0;
    for( const unsigned int setId : getObservationSetIdsForObservableType( observableType ) )
    {
        totalSize += getTotalScalarSizeForSet( setId );
    }
    return totalSize;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSingleLinkObservations( const ObservableType observableType,
                                                                                         const LinkDefinition& linkDefinition ) const
{
    return getSingleLinkObservationsAndTimes( observableType, linkDefinition ).first;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< TimeType > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSingleLinkTimes(
        const ObservableType observableType,
        const LinkDefinition& linkDefinition ) const
{
    return getSingleLinkObservationsAndTimes( observableType, linkDefinition ).second;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::pair< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >, std::vector< TimeType > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSingleLinkObservationsAndTimes(
        const ObservableType observableType,
        const LinkDefinition& linkDefinition ) const
{
    const ObservationSelectionCondition< ObservationScalarType, TimeType > condition =
            ObservationSelectionCondition< ObservationScalarType, TimeType >::observableType( observableType ) &&
            ObservationSelectionCondition< ObservationScalarType, TimeType >::linkDefinition( linkDefinition );
    const FlattenedObservationData< ObservationScalarType, TimeType > flattenedData =
            createFlattenedObservationDataFromObservationIds( getObservationIdsMatchingCondition( condition ), true );
    if( flattenedData.getObservationVector( ).size( ) == 0 )
    {
        throw std::runtime_error( "Error when getting single-link observations from dataset, no matching observations found." );
    }
    return std::make_pair( flattenedData.getObservationVector( ), flattenedData.getTimes( ) );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETFLATTENINGIMPLEMENTATION_H
