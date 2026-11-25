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
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setDependentVariablesForSet(
        const unsigned int setId,
        const std::vector< Eigen::VectorXd >& dependentVariables )
{
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );
    if( dependentVariables.empty( ) )
    {
        clearDependentVariablesForSet( setId );
        return;
    }
    if( dependentVariables.size( ) != observationIds.size( ) )
    {
        throw std::runtime_error( "Error when setting dataset dependent variables, size is inconsistent." );
    }
    int dependentVariableSize = dependentVariables.front( ).size( );
    const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& bookkeeping =
            getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
    if( bookkeeping != nullptr )
    {
        dependentVariableSize = bookkeeping->getTotalDependentVariableSize( );
    }
    for( const Eigen::VectorXd& dependentVariable : dependentVariables )
    {
        if( dependentVariable.size( ) != dependentVariableSize )
        {
            throw std::runtime_error(
                    "Error when setting dataset dependent variables, dependent-variable component size is inconsistent." );
        }
    }
    for( std::size_t i = 0; i < observationIds.size( ); ++i )
    {
        mutableObservationRow( observationIds.at( i ) ).dependentVariableValues_ = dependentVariables.at( i );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::clearDependentVariablesForSet( const unsigned int setId )
{
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        mutableObservationRow( observationId ).dependentVariableValues_ = Eigen::VectorXd( );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::addDependentVariableToSets(
        const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    if( dependentVariableSettings == nullptr )
    {
        throw std::runtime_error( "Error when adding dependent variable settings to observation dataset, settings are null." );
    }
    // Create the new bookkeeping of all sets before changing any of them, so that the dataset is unchanged if an error occurs
    const std::vector< unsigned int > matchingSetIds = getSetIds( condition );
    std::map< unsigned int, std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > > newBookkeepingPerSet;
    for( const unsigned int setId : std::set< unsigned int >( matchingSetIds.begin( ), matchingSetIds.end( ) ) )
    {
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
        }
        else if( !getDependentVariablesForSet( setId ).empty( ) )
        {
            throw std::runtime_error(
                    "Error when adding dependent variable settings to observation dataset, dependent-variable values "
                    "already exist in observation set " +
                    std::to_string( setId ) + "." );
        }
        else
        {
            bookkeeping = bookkeeping->clone( );
        }
        bookkeeping->addDependentVariables( allSettingsToCreate );
        newBookkeepingPerSet[ setId ] = bookkeeping;
    }
    for( const auto& setBookkeeping : newBookkeepingPerSet )
    {
        resetDependentVariableBookkeepingForSet( setBookkeeping.first, setBookkeeping.second );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< unsigned int > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationIdsMatchingCondition(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
{
    std::vector< unsigned int > observationIds;
    for( const auto& row : observationRows_ )
    {
        const unsigned int observationId = row.observationId_;
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
    const auto selected = getObservationIdsMatchingCondition( condition );
    auto result = std::make_shared< ObservationDataset >( *this );
    result->retainObservationRows( selected );
    return result;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETACCESSORSIMPLEMENTATION_H
