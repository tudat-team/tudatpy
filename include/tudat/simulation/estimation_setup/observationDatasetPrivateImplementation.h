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
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setResidualVector(
        const FlattenedObservationData< ObservationScalarType, TimeType >& flattenedObservationData,
        const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector )
{
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
unsigned int ObservationDataset< ObservationScalarType, TimeType, Dummy >::invalidObservationId( )
{
    return std::numeric_limits< unsigned int >::max( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::copyObservationStateAndWeightFrom(
        const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
        const unsigned int sourceObservationId,
        const unsigned int targetObservationId )
{
    observationRows_.at( targetObservationId ).isActive_ = sourceDataset.observationRows_.at( sourceObservationId ).isActive_;
    observationRows_.at( targetObservationId ).rejectionReason_ = sourceDataset.observationRows_.at( sourceObservationId ).rejectionReason_;
    observationWeights_.setObservationWeight( targetObservationId,
                                              sourceDataset.observationWeights_.getObservationWeight( sourceObservationId ),
                                              sourceDataset.observationWeights_.hasExplicitObservationWeight( sourceObservationId ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::copySetWeightBlockSubsetFrom(
        const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
        const unsigned int sourceSetId,
        const std::vector< unsigned int >& sourceObservationIds,
        const unsigned int targetSetId )
{
    std::vector< std::size_t > selectedSetScalarIndices;
    for( const unsigned int observationId : sourceObservationIds )
    {
        const ObservationDatasetRow< TimeType >& row = sourceDataset.observationRows_.at( observationId );
        for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
        {
            selectedSetScalarIndices.push_back( row.indexInSet_ * row.scalarSize_ + componentIndex );
        }
    }

    const Eigen::MatrixXd& sourceSetWeightMatrix = sourceDataset.observationWeights_.getSetWeightBlock( sourceSetId );
    Eigen::MatrixXd targetSetWeightMatrix = Eigen::MatrixXd::Zero( selectedSetScalarIndices.size( ), selectedSetScalarIndices.size( ) );
    for( std::size_t rowIndex = 0; rowIndex < selectedSetScalarIndices.size( ); ++rowIndex )
    {
        for( std::size_t columnIndex = 0; columnIndex < selectedSetScalarIndices.size( ); ++columnIndex )
        {
            targetSetWeightMatrix( rowIndex, columnIndex ) =
                    sourceSetWeightMatrix( selectedSetScalarIndices.at( rowIndex ), selectedSetScalarIndices.at( columnIndex ) );
        }
    }
    setWeightMatrixForSet( targetSetId, targetSetWeightMatrix );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::copyRemappedExtraWeightBlocksFrom(
        const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
        const std::map< unsigned int, unsigned int >& scalarComponentIdMap )
{
    for( const ObservationWeightBlock& sourceWeightBlock : sourceDataset.observationWeights_.getExtraWeightBlocks( ) )
    {
        std::vector< std::size_t > selectedRowIndices;
        std::vector< unsigned int > targetRowScalarComponentIds;
        for( std::size_t i = 0; i < sourceWeightBlock.rowScalarComponentIds_.size( ); ++i )
        {
            const unsigned int sourceScalarComponentId = sourceWeightBlock.rowScalarComponentIds_.at( i );
            if( scalarComponentIdMap.count( sourceScalarComponentId ) > 0 )
            {
                selectedRowIndices.push_back( i );
                targetRowScalarComponentIds.push_back( scalarComponentIdMap.at( sourceScalarComponentId ) );
            }
        }

        std::vector< std::size_t > selectedColumnIndices;
        std::vector< unsigned int > targetColumnScalarComponentIds;
        for( std::size_t i = 0; i < sourceWeightBlock.columnScalarComponentIds_.size( ); ++i )
        {
            const unsigned int sourceScalarComponentId = sourceWeightBlock.columnScalarComponentIds_.at( i );
            if( scalarComponentIdMap.count( sourceScalarComponentId ) > 0 )
            {
                selectedColumnIndices.push_back( i );
                targetColumnScalarComponentIds.push_back( scalarComponentIdMap.at( sourceScalarComponentId ) );
            }
        }

        if( selectedRowIndices.empty( ) || selectedColumnIndices.empty( ) )
        {
            continue;
        }

        Eigen::MatrixXd targetWeightBlock = Eigen::MatrixXd::Zero( selectedRowIndices.size( ), selectedColumnIndices.size( ) );
        for( std::size_t rowIndex = 0; rowIndex < selectedRowIndices.size( ); ++rowIndex )
        {
            for( std::size_t columnIndex = 0; columnIndex < selectedColumnIndices.size( ); ++columnIndex )
            {
                targetWeightBlock( rowIndex, columnIndex ) =
                        sourceWeightBlock.weightBlock_( selectedRowIndices.at( rowIndex ), selectedColumnIndices.at( columnIndex ) );
            }
        }

        ObservationWeightBlock targetObservationWeightBlock;
        targetObservationWeightBlock.rowScalarComponentIds_ = targetRowScalarComponentIds;
        targetObservationWeightBlock.columnScalarComponentIds_ = targetColumnScalarComponentIds;
        targetObservationWeightBlock.weightBlock_ = targetWeightBlock;
        observationWeights_.addExtraWeightBlock( targetObservationWeightBlock );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::replaceObservationSetDataWithSourceRows(
        const unsigned int setId,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
        const std::vector< unsigned int >& sourceObservationIdsForReplacement,
        const std::vector< bool >& explicitWeightsForReplacement )
{
    validateObservationSetData( setId, observations, times, dependentVariables, weights, residuals );
    if( !sourceObservationIdsForReplacement.empty( ) && sourceObservationIdsForReplacement.size( ) != observations.size( ) )
    {
        throw std::runtime_error( "Error when rebuilding observation dataset, source-row map size is inconsistent." );
    }
    if( !explicitWeightsForReplacement.empty( ) && explicitWeightsForReplacement.size( ) != observations.size( ) )
    {
        throw std::runtime_error( "Error when rebuilding observation dataset, explicit-weight map size is inconsistent." );
    }

    const ObservationDataset< ObservationScalarType, TimeType > sourceDataset = *this;
    ObservationDataset< ObservationScalarType, TimeType > rebuiltDataset;
    std::map< unsigned int, unsigned int > scalarComponentIdMap;

    for( unsigned int currentSetId = 0; currentSetId < sourceDataset.setMetadata_.size( ); ++currentSetId )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = sourceDataset.setMetadata_.at( currentSetId );
        const unsigned int newSetId = rebuiltDataset.addObservationSet(
                metadata.observableType_,
                sourceDataset.linkDefinitionRegistry_.at( metadata.linkDefinitionId_ ),
                currentSetId == setId ? observations : sourceDataset.getObservationsForSet( currentSetId ),
                currentSetId == setId ? times : sourceDataset.getObservationTimesForSet( currentSetId ),
                metadata.referenceLinkEnd_,
                currentSetId == setId ? dependentVariables : sourceDataset.getDependentVariablesForSet( currentSetId ),
                sourceDataset.dependentVariableLayoutRegistry_.at( metadata.dependentVariableLayoutId_ ),
                sourceDataset.ancillarySettingsRegistry_.at( metadata.ancillarySettingsId_ ),
                currentSetId == setId ? weights : sourceDataset.getWeightsForSet( currentSetId ),
                currentSetId == setId ? residuals : sourceDataset.getResidualsForSet( currentSetId ) );

        std::vector< unsigned int > sourceObservationIdsToPreserve;
        if( currentSetId == setId )
        {
            sourceObservationIdsToPreserve = sourceObservationIdsForReplacement;
        }
        else
        {
            sourceObservationIdsToPreserve = sourceDataset.observationIdsBySet_.at( currentSetId );
        }

        const std::vector< unsigned int >& targetObservationIds = rebuiltDataset.observationIdsBySet_.at( newSetId );
        std::vector< unsigned int > validSourceObservationIdsForSetBlock;
        for( std::size_t i = 0; i < sourceObservationIdsToPreserve.size( ); ++i )
        {
            const unsigned int sourceObservationId = sourceObservationIdsToPreserve.at( i );
            if( sourceObservationId == invalidObservationId( ) )
            {
                if( !explicitWeightsForReplacement.empty( ) && !explicitWeightsForReplacement.at( i ) )
                {
                    rebuiltDataset.observationWeights_.setObservationWeight(
                            targetObservationIds.at( i ),
                            rebuiltDataset.observationWeights_.getObservationWeight( targetObservationIds.at( i ) ),
                            false );
                }
                continue;
            }

            rebuiltDataset.copyObservationStateAndWeightFrom( sourceDataset, sourceObservationId, targetObservationIds.at( i ) );
            validSourceObservationIdsForSetBlock.push_back( sourceObservationId );

            const ObservationDatasetRow< TimeType >& sourceRow = sourceDataset.observationRows_.at( sourceObservationId );
            const ObservationDatasetRow< TimeType >& targetRow = rebuiltDataset.observationRows_.at( targetObservationIds.at( i ) );
            for( unsigned int componentIndex = 0; componentIndex < sourceRow.scalarSize_; ++componentIndex )
            {
                scalarComponentIdMap[ sourceRow.firstScalarComponent_ + componentIndex ] = targetRow.firstScalarComponent_ + componentIndex;
            }
        }

        if( sourceDataset.hasWeightMatrixForSet( currentSetId ) && !validSourceObservationIdsForSetBlock.empty( ) )
        {
            rebuiltDataset.copySetWeightBlockSubsetFrom( sourceDataset, currentSetId, validSourceObservationIdsForSetBlock, newSetId );
        }
    }

    rebuiltDataset.copyRemappedExtraWeightBlocksFrom( sourceDataset, scalarComponentIdMap );

    const std::size_t newStructuralVersion = structuralVersion_ + 1;
    *this = rebuiltDataset;
    structuralVersion_ = newStructuralVersion;
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

    for( std::size_t i = 0; i < observations.size( ); ++i )
    {
        if( observations.at( i ).size( ) != static_cast< int >( observableSize ) ||
            ( !weights.empty( ) && weights.at( i ).size( ) != static_cast< int >( observableSize ) ) ||
            ( !residuals.empty( ) && residuals.at( i ).size( ) != static_cast< int >( observableSize ) ) )
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

    std::sort( permutation.begin( ), permutation.end( ), [ &observationTimes ]( const std::size_t i, const std::size_t j ) {
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
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
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
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
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
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
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
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
    if( weight.size( ) != static_cast< int >( row.scalarSize_ ) )
    {
        throw std::runtime_error( "Error when setting dataset weight, scalar size is inconsistent." );
    }
    observationWeights_.setDiagonalWeightVector( observationId, weight );
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
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
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
