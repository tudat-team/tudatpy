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
    return createFlattenedObservationDataFromObservationIds( getAllObservationIds( ), includeRejected );
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
    bool hasUpdatedSet = false;
    for( unsigned int setId = 0; setId < getNumberOfObservationSets( ); ++setId )
    {
        const std::vector< unsigned int >& observationIds = getObservationIdsForSet( setId );
        bool setMatchesCondition = observationIds.empty( );
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
        hasUpdatedSet = true;
    }

    if( hasUpdatedSet )
    {
        ++structuralVersion_;
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::resetDependentVariableBookkeepingForSet(
        const unsigned int setId,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping )
{
    setMetadata_.at( setId ).dependentVariableLayoutId_ = registerDependentVariableLayout( dependentVariableBookkeeping );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAllObservationIds( ) const
{
    std::vector< unsigned int > observationIds;
    observationIds.reserve( observationRows_.size( ) );
    for( unsigned int observationId = 0; observationId < observationRows_.size( ); ++observationId )
    {
        observationIds.push_back( observationId );
    }
    return observationIds;
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
    FlattenedObservationData< ObservationScalarType, TimeType > flattenedObservationData;

    std::size_t flattenedDataSize = 0;
    bool materializeWeightMatrix = false;
    std::map< unsigned int, bool > selectedScalarComponentIds;
    std::map< unsigned int, bool > selectedSetIds;
    for( const unsigned int observationId : selectedObservationIds )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        if( includeInactive || row.isActive_ )
        {
            flattenedDataSize += row.scalarSize_;
            selectedSetIds[ row.setId_ ] = true;
            // A row-level block requires a matrix only when it can contribute to the final precedence-resolved projection.
            if( ( !observationWeights_.hasSetWeightBlock( row.setId_ ) ||
                  observationWeights_.hasExplicitObservationWeight( observationId ) ) &&
                !observationWeights_.isObservationWeightDiagonalOnly( observationId, row.scalarSize_ ) )
            {
                materializeWeightMatrix = true;
            }
            for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
            {
                selectedScalarComponentIds[ row.firstScalarComponent_ + componentIndex ] = true;
            }
        }
    }

    for( const auto& selectedSet : selectedSetIds )
    {
        // Set-level blocks provide the lowest-priority matrix layer for all selected components in the set.
        if( observationWeights_.hasSetWeightBlock( selectedSet.first ) &&
            !observationWeights_.isSetWeightBlockDiagonalOnly( selectedSet.first ) )
        {
            materializeWeightMatrix = true;
        }
    }

    for( const ObservationWeightBlock& extraWeightBlock : observationWeights_.getExtraWeightBlocks( ) )
    {
        // Extra scalar-component blocks may connect arbitrary observations; only selected components matter here.
        for( std::size_t i = 0; i < extraWeightBlock.rowScalarComponentIds_.size( ); ++i )
        {
            const unsigned int rowScalarComponentId = extraWeightBlock.rowScalarComponentIds_.at( i );
            if( selectedScalarComponentIds.count( rowScalarComponentId ) == 0 )
            {
                continue;
            }
            for( std::size_t j = 0; j < extraWeightBlock.columnScalarComponentIds_.size( ); ++j )
            {
                const unsigned int columnScalarComponentId = extraWeightBlock.columnScalarComponentIds_.at( j );
                if( selectedScalarComponentIds.count( columnScalarComponentId ) > 0 && rowScalarComponentId != columnScalarComponentId &&
                    extraWeightBlock.weightBlock_( i, j ) != 0.0 )
                {
                    materializeWeightMatrix = true;
                }
            }
        }
    }

    flattenedObservationData.observations_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( flattenedDataSize );
    flattenedObservationData.residuals_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( flattenedDataSize );
    flattenedObservationData.weights_ = Eigen::VectorXd::Zero( flattenedDataSize );
    if( materializeWeightMatrix )
    {
        flattenedObservationData.weightMatrix_.resize( flattenedDataSize, flattenedDataSize );
    }
    flattenedObservationData.times_.reserve( flattenedDataSize );
    flattenedObservationData.observationIds_.reserve( flattenedDataSize );
    flattenedObservationData.setIds_.reserve( flattenedDataSize );
    flattenedObservationData.scalarComponentIds_.reserve( flattenedDataSize );

    std::size_t currentIndex = 0;
    std::map< unsigned int, std::size_t > flattenedDataIndexByScalarComponent;
    std::map< int, bool > observationAlreadyRegistered;
    struct FlattenedWeightEntry {
        double weight_;
        std::string source_;
        unsigned int rowScalarComponentId_;
        unsigned int columnScalarComponentId_;
    };

    std::map< std::size_t, FlattenedWeightEntry > diagonalWeightEntries;
    std::map< std::pair< std::size_t, std::size_t >, FlattenedWeightEntry > sparseWeightEntries;
    const std::string precedenceDescription =
            "set-level weights, then explicit per-observation weights, then extra scalar-component blocks";

    auto warnWeightConflict = [ &precedenceDescription ]( const std::size_t rowIndex,
                                                          const std::size_t columnIndex,
                                                          const unsigned int rowScalarComponentId,
                                                          const unsigned int columnScalarComponentId,
                                                          const FlattenedWeightEntry& previousEntry,
                                                          const std::string& source,
                                                          const double weight ) {
        if( previousEntry.source_ != source && previousEntry.weight_ != weight )
        {
            // Keep warning text stable for tests; large override sets can emit one line per conflicting entry.
            std::cerr << "[WARNING] Conflicting observation weight entry at flattened matrix row " << rowIndex << ", column " << columnIndex
                      << " (scalar component ids " << rowScalarComponentId << ", " << columnScalarComponentId
                      << "): " << previousEntry.source_ << " value " << previousEntry.weight_ << " overwritten by " << source << " value "
                      << weight << ". Precedence is " << precedenceDescription << "." << std::endl;
        }
    };

    auto setFlattenedDiagonalWeightEntry =
            [ &diagonalWeightEntries, &flattenedObservationData, &warnWeightConflict ](
                    const std::size_t rowIndex, const unsigned int scalarComponentId, const double weight, const std::string& source ) {
                const auto existingEntry = diagonalWeightEntries.find( rowIndex );
                if( existingEntry != diagonalWeightEntries.end( ) )
                {
                    warnWeightConflict( rowIndex, rowIndex, scalarComponentId, scalarComponentId, existingEntry->second, source, weight );
                }
                diagonalWeightEntries[ rowIndex ] = { weight, source, scalarComponentId, scalarComponentId };
                flattenedObservationData.weights_( rowIndex ) = weight;
            };

    // Store matrix entries in a map first so later, higher-priority blocks can overwrite previous entries cleanly.
    auto setFlattenedWeightEntry = [ &sparseWeightEntries, &warnWeightConflict ]( const std::size_t rowIndex,
                                                                                  const std::size_t columnIndex,
                                                                                  const unsigned int rowScalarComponentId,
                                                                                  const unsigned int columnScalarComponentId,
                                                                                  const double weight,
                                                                                  const std::string& source,
                                                                                  const bool warnOnConflict = true ) {
        const std::pair< std::size_t, std::size_t > indexPair = std::make_pair( rowIndex, columnIndex );
        const auto existingEntry = sparseWeightEntries.find( indexPair );
        if( warnOnConflict && existingEntry != sparseWeightEntries.end( ) )
        {
            warnWeightConflict(
                    rowIndex, columnIndex, rowScalarComponentId, columnScalarComponentId, existingEntry->second, source, weight );
        }
        sparseWeightEntries[ indexPair ] = { weight, source, rowScalarComponentId, columnScalarComponentId };
    };
    for( const unsigned int observationId : selectedObservationIds )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        if( includeInactive || row.isActive_ )
        {
            if( observationAlreadyRegistered.count( observationId ) == 0 )
            {
                if( flattenedObservationData.uniqueObservationIdsBySet_.size( ) <= row.setId_ )
                {
                    flattenedObservationData.uniqueObservationIdsBySet_.resize( row.setId_ + 1 );
                }
                if( flattenedObservationData.uniqueObservationIdsBySet_.at( row.setId_ ).empty( ) )
                {
                    flattenedObservationData.setIdsInRowOrder_.push_back( row.setId_ );
                }
                flattenedObservationData.uniqueObservationIdsBySet_.at( row.setId_ ).push_back( observationId );
                observationAlreadyRegistered[ observationId ] = true;
            }
            if( flattenedObservationData.firstFlattenedRowByObservation_.size( ) <= observationId )
            {
                flattenedObservationData.firstFlattenedRowByObservation_.resize( observationId + 1, -1 );
                flattenedObservationData.scalarSizeByObservation_.resize( observationId + 1, 0 );
            }
            flattenedObservationData.firstFlattenedRowByObservation_.at( observationId ) = static_cast< int >( currentIndex );
            flattenedObservationData.scalarSizeByObservation_.at( observationId ) = row.scalarSize_;

            for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
            {
                const unsigned int scalarComponentId = row.firstScalarComponent_ + componentIndex;
                flattenedObservationData.observations_( currentIndex ) = observedValues_.at( scalarComponentId );
                flattenedObservationData.residuals_( currentIndex ) = residualValues_.at( scalarComponentId );
                flattenedObservationData.times_.push_back( row.time_ );
                flattenedObservationData.observationIds_.push_back( observationId );
                flattenedObservationData.setIds_.push_back( row.setId_ );
                flattenedObservationData.scalarComponentIds_.push_back( scalarComponentId );
                flattenedDataIndexByScalarComponent[ scalarComponentId ] = currentIndex;
                ++currentIndex;
            }
        }
    }

    for( unsigned int setId = 0; setId < setMetadata_.size( ); ++setId )
    {
        if( observationWeights_.hasSetWeightBlock( setId ) )
        {
            // The set-level block is indexed in local set order and projected to the selected global rows.
            const Eigen::MatrixXd& setWeightBlock = observationWeights_.getSetWeightBlock( setId );
            for( const unsigned int rowObservationId : observationIdsBySet_.at( setId ) )
            {
                const ObservationDatasetRow< TimeType >& row = observationRows_.at( rowObservationId );
                if( includeInactive || row.isActive_ )
                {
                    for( unsigned int rowComponentIndex = 0; rowComponentIndex < row.scalarSize_; ++rowComponentIndex )
                    {
                        const unsigned int rowScalarComponentId = row.firstScalarComponent_ + rowComponentIndex;
                        if( flattenedDataIndexByScalarComponent.count( rowScalarComponentId ) == 0 )
                        {
                            continue;
                        }
                        const std::size_t flattenedRow = flattenedDataIndexByScalarComponent.at( rowScalarComponentId );
                        const std::size_t setLocalRow = row.indexInSet_ * setMetadata_.at( setId ).observableSize_ + rowComponentIndex;
                        setFlattenedDiagonalWeightEntry(
                                flattenedRow, rowScalarComponentId, setWeightBlock( setLocalRow, setLocalRow ), "set-level block" );

                        if( materializeWeightMatrix )
                        {
                            for( const unsigned int columnObservationId : observationIdsBySet_.at( setId ) )
                            {
                                const ObservationDatasetRow< TimeType >& columnRow = observationRows_.at( columnObservationId );
                                if( includeInactive || columnRow.isActive_ )
                                {
                                    for( unsigned int columnComponentIndex = 0; columnComponentIndex < columnRow.scalarSize_;
                                         ++columnComponentIndex )
                                    {
                                        const unsigned int columnScalarComponentId = columnRow.firstScalarComponent_ + columnComponentIndex;
                                        if( flattenedDataIndexByScalarComponent.count( columnScalarComponentId ) == 0 )
                                        {
                                            continue;
                                        }
                                        const std::size_t flattenedColumn =
                                                flattenedDataIndexByScalarComponent.at( columnScalarComponentId );
                                        const std::size_t setLocalColumn =
                                                columnRow.indexInSet_ * setMetadata_.at( setId ).observableSize_ + columnComponentIndex;
                                        setFlattenedWeightEntry( flattenedRow,
                                                                 flattenedColumn,
                                                                 rowScalarComponentId,
                                                                 columnScalarComponentId,
                                                                 setWeightBlock( setLocalRow, setLocalColumn ),
                                                                 "set-level block",
                                                                 false );
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }

    for( const unsigned int observationId : selectedObservationIds )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        if( includeInactive || row.isActive_ )
        {
            const bool insertObservationWeight = !observationWeights_.hasSetWeightBlock( row.setId_ ) ||
                    observationWeights_.hasExplicitObservationWeight( observationId );
            if( !insertObservationWeight )
            {
                continue;
            }

            const Eigen::VectorXd observationWeightVector =
                    observationWeights_.getObservationWeightVector( observationId, row.scalarSize_ );
            const std::string source = "per-observation weight for observation id " + std::to_string( observationId );
            for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
            {
                const unsigned int scalarComponentId = row.firstScalarComponent_ + componentIndex;
                const std::size_t flattenedRow = flattenedDataIndexByScalarComponent.at( scalarComponentId );
                setFlattenedDiagonalWeightEntry( flattenedRow, scalarComponentId, observationWeightVector( componentIndex ), source );
            }

            if( materializeWeightMatrix )
            {
                const Eigen::MatrixXd observationWeightMatrix =
                        observationWeights_.getObservationWeightMatrix( observationId, row.scalarSize_ );
                for( unsigned int rowComponentIndex = 0; rowComponentIndex < row.scalarSize_; ++rowComponentIndex )
                {
                    const unsigned int rowScalarComponentId = row.firstScalarComponent_ + rowComponentIndex;
                    const std::size_t flattenedRow = flattenedDataIndexByScalarComponent.at( rowScalarComponentId );
                    for( unsigned int columnComponentIndex = 0; columnComponentIndex < row.scalarSize_; ++columnComponentIndex )
                    {
                        const unsigned int columnScalarComponentId = row.firstScalarComponent_ + columnComponentIndex;
                        const std::size_t flattenedColumn = flattenedDataIndexByScalarComponent.at( columnScalarComponentId );
                        setFlattenedWeightEntry( flattenedRow,
                                                 flattenedColumn,
                                                 rowScalarComponentId,
                                                 columnScalarComponentId,
                                                 observationWeightMatrix( rowComponentIndex, columnComponentIndex ),
                                                 source,
                                                 rowComponentIndex != columnComponentIndex );
                    }
                }
            }
        }
    }

    for( const ObservationWeightBlock& extraWeightBlock : observationWeights_.getExtraWeightBlocks( ) )
    {
        // Extra blocks are inserted last and can overwrite individual scalar-component entries.
        for( std::size_t i = 0; i < extraWeightBlock.rowScalarComponentIds_.size( ); ++i )
        {
            const unsigned int rowScalarComponentId = extraWeightBlock.rowScalarComponentIds_.at( i );
            if( flattenedDataIndexByScalarComponent.count( rowScalarComponentId ) == 0 )
            {
                continue;
            }
            for( std::size_t j = 0; j < extraWeightBlock.columnScalarComponentIds_.size( ); ++j )
            {
                const unsigned int columnScalarComponentId = extraWeightBlock.columnScalarComponentIds_.at( j );
                if( flattenedDataIndexByScalarComponent.count( columnScalarComponentId ) == 0 )
                {
                    continue;
                }
                if( rowScalarComponentId == columnScalarComponentId )
                {
                    setFlattenedDiagonalWeightEntry( flattenedDataIndexByScalarComponent.at( rowScalarComponentId ),
                                                     rowScalarComponentId,
                                                     extraWeightBlock.weightBlock_( i, j ),
                                                     "extra scalar-component block" );
                }
                if( materializeWeightMatrix )
                {
                    setFlattenedWeightEntry( flattenedDataIndexByScalarComponent.at( rowScalarComponentId ),
                                             flattenedDataIndexByScalarComponent.at( columnScalarComponentId ),
                                             rowScalarComponentId,
                                             columnScalarComponentId,
                                             extraWeightBlock.weightBlock_( i, j ),
                                             "extra scalar-component block",
                                             rowScalarComponentId != columnScalarComponentId );
                }
            }
        }
    }

    if( materializeWeightMatrix )
    {
        // Create the compressed sparse matrix only when an off-diagonal path was detected.
        std::vector< Eigen::Triplet< double > > sparseWeightTriplets;
        sparseWeightTriplets.reserve( sparseWeightEntries.size( ) );
        for( const auto& weightEntry : sparseWeightEntries )
        {
            if( weightEntry.second.weight_ != 0.0 )
            {
                sparseWeightTriplets.emplace_back( weightEntry.first.first, weightEntry.first.second, weightEntry.second.weight_ );
            }
            if( weightEntry.first.first == weightEntry.first.second )
            {
                flattenedObservationData.weights_( weightEntry.first.first ) = weightEntry.second.weight_;
            }
        }
        flattenedObservationData.weightMatrix_.setFromTriplets( sparseWeightTriplets.begin( ), sparseWeightTriplets.end( ) );
        flattenedObservationData.weightMatrix_.makeCompressed( );
        flattenedObservationData.isDiagonalWeightOnly_ = false;
    }
    else
    {
        flattenedObservationData.isDiagonalWeightOnly_ = true;
    }

    return flattenedObservationData;
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
