/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETVECTORDATAIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETVECTORDATAIMPLEMENTATION_H

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
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSetIdsInObservationVectorOrder( ) const
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
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setLinkEndReferencePoint(
        const std::string& bodyName,
        const std::string& referencePointName,
        const LinkEndType linkEndType,
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    const std::vector< unsigned int > matchingSetIds = getSetIds( condition );
    for( const unsigned int setId : std::set< unsigned int >( matchingSetIds.begin( ), matchingSetIds.end( ) ) )
    {
        std::map< LinkEndType, LinkEndId > linkEnds = getLinkDefinition( setMetadata_.at( setId ).linkDefinitionId_ ).linkEnds_;
        typename std::map< LinkEndType, LinkEndId >::iterator linkEndIterator = linkEnds.find( linkEndType );
        if( linkEndIterator == linkEnds.end( ) || linkEndIterator->second.bodyName_ != bodyName ||
            linkEndIterator->second.getReferencePointName( ) == referencePointName )
        {
            continue;
        }
        const long numberOfSelectedObservations = std::count( matchingSetIds.begin( ), matchingSetIds.end( ), setId );
        if( numberOfSelectedObservations != static_cast< long >( getNumberOfObservationsForSet( setId ) ) )
        {
            throw std::runtime_error(
                    "Error when setting link end reference point, the condition selects " + std::to_string( numberOfSelectedObservations ) +
                    " of the " + std::to_string( getNumberOfObservationsForSet( setId ) ) + " observations of set " +
                    std::to_string( setId ) + ", but the reference point can only be changed for complete observation sets." );
        }

        const LinkEndId oldLinkEndId = linkEndIterator->second;
        linkEndIterator->second = LinkEndId( linkEndIterator->second.bodyName_, referencePointName );
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > oldBookkeeping =
                getDependentVariableBookkeeping( setMetadata_.at( setId ).dependentVariableLayoutId_ );
        resetLinkDefinitionForSet( setId, LinkDefinition( linkEnds ) );

        // The dependent variables of the set are defined for the new link ends (their stored values, like the residuals, are not
        // recomputed here)
        if( oldBookkeeping != nullptr )
        {
            resetDependentVariableBookkeepingForSet(
                    setId,
                    oldBookkeeping->cloneWithChangedLinkEnd(
                            LinkDefinition( linkEnds ), linkEndType, oldLinkEndId, linkEndIterator->second ) );
        }
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
    setMetadata_.at( setId ).dependentVariableLayoutId_ = registerEntry( dependentVariableLayoutRegistry_, dependentVariableBookkeeping );
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
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationIdsInObservationVectorOrder( ) const
{
    std::vector< unsigned int > observationIds;
    observationIds.reserve( observationRows_.size( ) );
    for( const unsigned int setId : getSetIdsInObservationVectorOrder( ) )
    {
        const std::vector< unsigned int >& setObservationIds = observationIdsBySet_.at( setId );
        observationIds.insert( observationIds.end( ), setObservationIds.begin( ), setObservationIds.end( ) );
    }
    return observationIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
ObservationVectorData< ObservationScalarType, TimeType >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createObservationVectorData( const bool includeRejected ) const
{
    ObservationVectorData< ObservationScalarType, TimeType > result;
    result.source_ = lifetimeToken_;
    result.structuralVersion_ = structuralVersion_;
    result.vectorDataVersion_ = vectorDataVersion_;
    result.uniqueObservationIdsBySet_.resize( setMetadata_.size( ) );
    // Rejection changes the selection, not the complete weights or their covariance cache.
    if( !completeWeightData_ || completeWeightData_->structuralVersion_ != structuralVersion_ )
    {
        auto data = std::make_shared< typename ObservationVectorData< ObservationScalarType, TimeType >::CompleteWeightData >( );
        data->weights_ = observationWeights_;
        data->setMetadata_ = setMetadata_;
        data->observationIdsBySet_ = observationIdsBySet_;
        data->structuralVersion_ = structuralVersion_;
        data->observationMapping_.reserve( observationRows_.size( ) );
        for( const auto& row : observationRows_ )
        {
            data->observationMapping_.emplace( row.observationId_, std::make_pair( row.firstScalarComponent_, row.setId_ ) );
        }
        completeWeightData_ = std::move( data );
    }
    result.completeWeightData_ = completeWeightData_;
    const std::vector< unsigned int > selected = resolveObservationIds(
            ObservationSelectionCondition< ObservationScalarType, TimeType >::all( ), ObservationOrdering::estimation, includeRejected );
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

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETVECTORDATAIMPLEMENTATION_H
