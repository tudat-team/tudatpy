/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETACCESSORSIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETACCESSORSIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSingleDependentVariableForSet(
        const unsigned int setId,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
        const bool returnFirstCompatibleSettings ) const
{
    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
            getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
    if( dependentVariableBookkeeping == nullptr )
    {
        throw std::runtime_error(
                "Error when getting dependent variable from observation dataset, no dependent variable bookkeeping is available." );
    }

    std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > settingsIndicesAndSizes =
            dependentVariableBookkeeping->getSettingsIndicesAndSizes( );

    std::vector< std::pair< int, int > > indicesAndSizes;
    for( const auto& settingsIterator : settingsIndicesAndSizes )
    {
        if( dependentVariableSettings->areSettingsCompatible( settingsIterator.second ) )
        {
            indicesAndSizes.push_back( settingsIterator.first );
        }
    }

    if( indicesAndSizes.empty( ) )
    {
        throw std::runtime_error( "Error when getting dependent variable, no dependent variable values found for given settings." );
    }
    if( indicesAndSizes.size( ) > 1 && !returnFirstCompatibleSettings )
    {
        throw std::runtime_error( "Error when getting dependent variable, multiple dependent variables found for given settings." );
    }

    return getSingleDependentVariableForSet( setId, indicesAndSizes.at( 0 ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getCompatibleDependentVariableSettingsForSet(
        const unsigned int setId,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
{
    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
            getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
    if( dependentVariableBookkeeping == nullptr )
    {
        return std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >( );
    }

    std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > compatibleSettings;
    for( const auto& settings : dependentVariableBookkeeping->getDependentVariableSettings( ) )
    {
        if( dependentVariableSettings->areSettingsCompatible( settings ) )
        {
            compatibleSettings.push_back( settings );
        }
    }
    return compatibleSettings;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::MatrixXd > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAllCompatibleDependentVariablesForSet(
        const unsigned int setId,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
{
    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
            getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
    if( dependentVariableBookkeeping == nullptr )
    {
        return std::vector< Eigen::MatrixXd >( );
    }

    std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > settingsIndicesAndSizes =
            dependentVariableBookkeeping->getSettingsIndicesAndSizes( );

    std::vector< Eigen::MatrixXd > dependentVariablesList;
    for( const auto& settingsIterator : settingsIndicesAndSizes )
    {
        if( dependentVariableSettings->areSettingsCompatible( settingsIterator.second ) )
        {
            dependentVariablesList.push_back( getSingleDependentVariableForSet( setId, settingsIterator.first ) );
        }
    }
    return dependentVariablesList;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setDependentVariablesForSet(
        const unsigned int setId,
        const std::vector< Eigen::VectorXd >& dependentVariables )
{
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );
    if( dependentVariables.size( ) != observationIds.size( ) )
    {
        throw std::runtime_error( "Error when setting dataset dependent variables, size is inconsistent." );
    }
    for( std::size_t i = 0; i < observationIds.size( ); ++i )
    {
        observationRows_.at( observationIds.at( i ) ).dependentVariableValues_ = dependentVariables.at( i );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::clearDependentVariablesForSet( const unsigned int setId )
{
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        observationRows_.at( observationId ).dependentVariableValues_ = Eigen::VectorXd( );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getNumberOfObservationsForSet( const unsigned int setId ) const
{
    return observationIdsBySet_.at( setId ).size( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTotalScalarSizeForSet( const unsigned int setId ) const
{
    return observationIdsBySet_.at( setId ).size( ) * getObservationSetMetadata( setId ).observableSize_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const LinkDefinition& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getLinkDefinition(
        const unsigned int linkDefinitionId ) const
{
    return linkDefinitionRegistry_.at( linkDefinitionId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getNumberOfLinkDefinitions( ) const
{
    return linkDefinitionRegistry_.size( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::shared_ptr< ObservationAncillarySimulationSettings >&
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAncillarySettings( const unsigned int ancillarySettingsId ) const
{
    return ancillarySettingsRegistry_.at( ancillarySettingsId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::shared_ptr< ObservationAncillarySimulationSettings >&
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAncillarySettingsForSet( const unsigned int setId ) const
{
    return getAncillarySettings( getObservationSetMetadata( setId ).ancillarySettingsId_ );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >&
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getDependentVariableBookkeeping(
        const unsigned int dependentVariableLayoutId ) const
{
    return dependentVariableLayoutRegistry_.at( dependentVariableLayoutId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetIdsForDependentVariableSettings(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
{
    std::vector< unsigned int > matchingSetIds;
    for( const unsigned int setId : getSetIdsInOrderedFlattenedDataOrder( ) )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const LinkEnds& linkEnds = getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        if( !getCompatibleDependentVariableSettingsForSet( setId, dependentVariableSettings ).empty( ) &&
            !simulation_setup::createAllCompatibleDependentVariableSettings( metadata.observableType_, linkEnds, dependentVariableSettings )
                     .empty( ) )
        {
            matchingSetIds.push_back( setId );
        }
    }
    return matchingSetIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetIdsWithDependentVariableValues(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
{
    std::vector< unsigned int > matchingSetIds;
    for( const unsigned int setId : getObservationSetIdsForDependentVariableSettings( dependentVariableSettings ) )
    {
        if( !getAllCompatibleDependentVariablesForSet( setId, dependentVariableSettings ).empty( ) )
        {
            matchingSetIds.push_back( setId );
        }
    }
    return matchingSetIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getCompatibleDependentVariableSettingsPerSet(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
{
    std::vector< std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > > compatibleSettingsPerSet;
    for( const unsigned int setId : getObservationSetIdsWithDependentVariableValues( dependentVariableSettings ) )
    {
        const std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > compatibleSettings =
                getCompatibleDependentVariableSettingsForSet( setId, dependentVariableSettings );
        if( !compatibleSettings.empty( ) )
        {
            compatibleSettingsPerSet.push_back( compatibleSettings );
        }
    }
    return compatibleSettingsPerSet;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< std::vector< Eigen::MatrixXd > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getAllCompatibleDependentVariablesPerSet(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
{
    std::vector< std::vector< Eigen::MatrixXd > > dependentVariablesPerSet;
    for( const unsigned int setId : getObservationSetIdsWithDependentVariableValues( dependentVariableSettings ) )
    {
        const std::vector< Eigen::MatrixXd > currentVariables =
                getAllCompatibleDependentVariablesForSet( setId, dependentVariableSettings );
        if( !currentVariables.empty( ) )
        {
            dependentVariablesPerSet.push_back( currentVariables );
        }
    }
    return dependentVariablesPerSet;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::map< TimeType, Eigen::VectorXd > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getDependentVariableHistory(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
        const bool returnFirstCompatibleSettings ) const
{
    std::map< TimeType, Eigen::VectorXd > dependentVariableHistory;
    for( const unsigned int setId : getObservationSetIdsWithDependentVariableValues( dependentVariableSettings ) )
    {
        const Eigen::MatrixXd dependentVariables =
                getSingleDependentVariableForSet( setId, dependentVariableSettings, returnFirstCompatibleSettings );
        const std::vector< TimeType > observationTimes = getObservationTimesForSet( setId );
        for( unsigned int i = 0; i < observationTimes.size( ); ++i )
        {
            dependentVariableHistory[ observationTimes.at( i ) ] =
                    dependentVariables.block( i, 0, 1, dependentVariables.cols( ) ).transpose( );
        }
    }
    return dependentVariableHistory;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::addDependentVariableToSets(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
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

        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const LinkEnds& linkEnds = getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > allSettingsToCreate =
                simulation_setup::createAllCompatibleDependentVariableSettings(
                        metadata.observableType_, linkEnds, dependentVariableSettings );
        if( allSettingsToCreate.empty( ) )
        {
            continue;
        }

        std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > bookkeeping =
                getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
        if( bookkeeping == nullptr )
        {
            bookkeeping =
                    std::make_shared< simulation_setup::ObservationDependentVariableBookkeeping >( metadata.observableType_, linkEnds );
            setMetadata_.at( setId ).dependentVariableLayoutId_ = registerDependentVariableLayout( bookkeeping );
        }
        else if( !getDependentVariablesForSet( setId ).empty( ) )
        {
            throw std::runtime_error(
                    "Error when adding dependent variable settings to observation dataset, dependent-variable values already exist." );
        }
        bookkeeping->addDependentVariables( allSettingsToCreate );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationIdsMatchingCondition(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
{
    std::vector< unsigned int > observationIds;
    for( unsigned int observationId = 0; observationId < observationRows_.size( ); ++observationId )
    {
        if( condition( *this, observationId ) )
        {
            observationIds.push_back( observationId );
        }
    }
    return observationIds;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createNewAndKeep(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
{
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > reducedDataset =
            std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );
    std::map< unsigned int, unsigned int > scalarComponentIdMap;

    for( unsigned int setId = 0; setId < setMetadata_.size( ); ++setId )
    {
        std::vector< unsigned int > selectedObservationIds;
        for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
        {
            if( condition( *this, observationId ) )
            {
                selectedObservationIds.push_back( observationId );
            }
        }

        if( !selectedObservationIds.empty( ) )
        {
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = setMetadata_.at( setId );
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
            std::vector< TimeType > times;
            std::vector< Eigen::VectorXd > dependentVariables;
            std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;
            const bool hasDependentVariables = !getDependentVariablesForSet( setId ).empty( );
            const bool hasSetWeightBlock = hasWeightMatrixForSet( setId );

            for( const unsigned int observationId : selectedObservationIds )
            {
                observations.push_back( getObservationValue( observationId ) );
                times.push_back( getObservationTime( observationId ) );
                if( !hasSetWeightBlock )
                {
                    weights.push_back( getWeightValue( observationId ) );
                }
                residuals.push_back( getResidualValue( observationId ) );
                if( hasDependentVariables )
                {
                    dependentVariables.push_back( getDependentVariables( observationId ) );
                }
            }

            const unsigned int newSetId =
                    reducedDataset->addObservationSet( metadata.observableType_,
                                                       linkDefinitionRegistry_.at( metadata.linkDefinitionId_ ),
                                                       observations,
                                                       times,
                                                       metadata.referenceLinkEnd_,
                                                       dependentVariables,
                                                       dependentVariableLayoutRegistry_.at( metadata.dependentVariableLayoutId_ ),
                                                       ancillarySettingsRegistry_.at( metadata.ancillarySettingsId_ ),
                                                       weights,
                                                       residuals );

            if( hasSetWeightBlock )
            {
                std::vector< std::size_t > selectedSetScalarIndices;
                for( const unsigned int observationId : selectedObservationIds )
                {
                    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
                    for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
                    {
                        selectedSetScalarIndices.push_back( row.indexInSet_ * row.scalarSize_ + componentIndex );
                    }
                }

                const Eigen::MatrixXd sourceSetWeightMatrix = getWeightMatrixForSet( setId );
                Eigen::MatrixXd reducedSetWeightMatrix =
                        Eigen::MatrixXd::Zero( selectedSetScalarIndices.size( ), selectedSetScalarIndices.size( ) );
                for( std::size_t rowIndex = 0; rowIndex < selectedSetScalarIndices.size( ); ++rowIndex )
                {
                    for( std::size_t columnIndex = 0; columnIndex < selectedSetScalarIndices.size( ); ++columnIndex )
                    {
                        reducedSetWeightMatrix( rowIndex, columnIndex ) = sourceSetWeightMatrix(
                                selectedSetScalarIndices.at( rowIndex ), selectedSetScalarIndices.at( columnIndex ) );
                    }
                }
                reducedDataset->setWeightMatrixForSet( newSetId, reducedSetWeightMatrix );
            }

            const std::vector< unsigned int >& newObservationIds = reducedDataset->getObservationIdsForSet( newSetId );
            for( std::size_t i = 0; i < selectedObservationIds.size( ); ++i )
            {
                reducedDataset->observationWeights_.setObservationWeight(
                        newObservationIds.at( i ),
                        observationWeights_.getObservationWeight( selectedObservationIds.at( i ) ),
                        observationWeights_.hasExplicitObservationWeight( selectedObservationIds.at( i ) ) );
            }

            for( std::size_t i = 0; i < selectedObservationIds.size( ); ++i )
            {
                const ObservationDatasetRow< TimeType >& sourceRow = observationRows_.at( selectedObservationIds.at( i ) );
                ObservationDatasetRow< TimeType >& targetRow =
                        reducedDataset->observationRows_.at( reducedDataset->observationIdsBySet_.at( newSetId ).at( i ) );
                targetRow.isActive_ = sourceRow.isActive_;
                targetRow.rejectionReason_ = sourceRow.rejectionReason_;
                for( unsigned int componentIndex = 0; componentIndex < sourceRow.scalarSize_; ++componentIndex )
                {
                    scalarComponentIdMap[ sourceRow.firstScalarComponent_ + componentIndex ] =
                            targetRow.firstScalarComponent_ + componentIndex;
                }
            }
        }
    }

    reducedDataset->copyRemappedExtraWeightBlocksFrom( *this, scalarComponentIdMap );

    return reducedDataset;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETACCESSORSIMPLEMENTATION_H
