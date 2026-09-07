/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETMUTATIONIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETMUTATIONIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

namespace observation_dataset_detail
{

template< typename ObservationScalarType >
Eigen::VectorXd accumulateResidualStatistic( const unsigned int observableSize,
                                             const std::size_t numberOfObservations,
                                             const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
                                             const bool rootMeanSquare )
{
    Eigen::VectorXd statistic = Eigen::VectorXd::Zero( observableSize );
    for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
    {
        for( std::size_t observationIndex = 0; observationIndex < numberOfObservations; ++observationIndex )
        {
            const double residual = static_cast< double >( residuals.at( observationIndex )( componentIndex, 0 ) );
            statistic[ componentIndex ] += rootMeanSquare ? residual * residual : residual;
        }
        statistic[ componentIndex ] /= static_cast< double >( numberOfObservations );
        if( rootMeanSquare )
        {
            statistic[ componentIndex ] = std::sqrt( statistic[ componentIndex ] );
        }
    }
    return statistic;
}

}  // namespace observation_dataset_detail

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
int ObservationDataset< ObservationScalarType, TimeType, Dummy >::addObservationSet(
        const ObservableType observableType,
        const LinkDefinition& linkDefinition,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const LinkEndType referenceLinkEnd,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings,
        const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
        const bool sortObservations,
        const bool eraseDuplicateObservations )
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > preparedObservations = observations;
    std::vector< TimeType > preparedTimes = times;
    std::vector< Eigen::VectorXd > preparedDependentVariables = dependentVariables;
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > preparedWeights = weights;
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > preparedResiduals = residuals;

    const int declaredObservableSize = getObservableSize( observableType );
    if( declaredObservableSize <= 0 )
    {
        throw std::runtime_error( "Error when adding observation set to dataset, observable size is invalid." );
    }
    const unsigned int observableSize = static_cast< unsigned int >( declaredObservableSize );

    if( preparedObservations.size( ) != preparedTimes.size( ) )
    {
        throw std::runtime_error( "Error when adding observation set to dataset, observation and time sizes are inconsistent." );
    }
    if( !preparedWeights.empty( ) && preparedWeights.size( ) != preparedObservations.size( ) )
    {
        throw std::runtime_error( "Error when adding observation set to dataset, weight size is inconsistent." );
    }
    if( !preparedResiduals.empty( ) && preparedResiduals.size( ) != preparedObservations.size( ) )
    {
        throw std::runtime_error( "Error when adding observation set to dataset, residual size is inconsistent." );
    }
    if( !preparedDependentVariables.empty( ) && preparedDependentVariables.size( ) != preparedObservations.size( ) )
    {
        throw std::runtime_error( "Error when adding observation set to dataset, dependent variable size is inconsistent." );
    }

    int dependentVariableSize = -1;
    if( dependentVariableBookkeeping != nullptr )
    {
        if( dependentVariableBookkeeping->getObservableType( ) != observableType ||
            !( dependentVariableBookkeeping->getLinkEnds( ) == linkDefinition ) )
        {
            throw std::runtime_error(
                    "Error when adding observation set to dataset, dependent-variable bookkeeping is incompatible with the set." );
        }
        dependentVariableSize = dependentVariableBookkeeping->getTotalDependentVariableSize( );
    }
    else if( !preparedDependentVariables.empty( ) )
    {
        dependentVariableSize = preparedDependentVariables.front( ).size( );
    }

    for( std::size_t i = 0; i < preparedObservations.size( ); ++i )
    {
        if( preparedObservations.at( i ).size( ) != static_cast< int >( observableSize ) ||
            ( !preparedResiduals.empty( ) && preparedResiduals.at( i ).size( ) != static_cast< int >( observableSize ) ) ||
            ( !preparedWeights.empty( ) && preparedWeights.at( i ).size( ) != static_cast< int >( observableSize ) ) ||
            ( !preparedDependentVariables.empty( ) && preparedDependentVariables.at( i ).size( ) != dependentVariableSize ) )
        {
            throw std::runtime_error( "Error when adding observation set to dataset, scalar component size is inconsistent." );
        }
    }

    if( sortObservations && preparedTimes.size( ) > 1 )
    {
        const std::vector< std::size_t > permutation = getTimeSortingPermutation( preparedTimes );
        reorderVector( preparedObservations, permutation );
        reorderVector( preparedTimes, permutation );
        if( !preparedDependentVariables.empty( ) )
        {
            reorderVector( preparedDependentVariables, permutation );
        }
        if( !preparedWeights.empty( ) )
        {
            reorderVector( preparedWeights, permutation );
        }
        if( !preparedResiduals.empty( ) )
        {
            reorderVector( preparedResiduals, permutation );
        }
    }

    if( eraseDuplicateObservations && preparedTimes.size( ) > 1 )
    {
        std::set< TimeType > retainedTimes;
        std::vector< unsigned int > indicesToRemove;
        for( unsigned int i = 0; i < preparedTimes.size( ); ++i )
        {
            if( !retainedTimes.insert( preparedTimes.at( i ) ).second )
            {
                indicesToRemove.push_back( i );
            }
        }

        if( !indicesToRemove.empty( ) )
        {
            const std::size_t beforeCount = preparedObservations.size( );
            removeEntries( preparedObservations, indicesToRemove );
            removeEntries( preparedTimes, indicesToRemove );
            if( !preparedDependentVariables.empty( ) )
            {
                removeEntries( preparedDependentVariables, indicesToRemove );
            }
            if( !preparedWeights.empty( ) )
            {
                removeEntries( preparedWeights, indicesToRemove );
            }
            if( !preparedResiduals.empty( ) )
            {
                removeEntries( preparedResiduals, indicesToRemove );
            }
            std::cerr << "[WARNING] Detected and removed " << beforeCount - preparedObservations.size( )
                      << "duplicate observations when creating observation dataset" << std::endl;
        }
    }

    for( const auto& weight : weights )
    {
        ObservationWeights::validateDiagonal( weight );
    }

    const unsigned int linkDefinitionId = registerLinkDefinition( linkDefinition );
    const unsigned int ancillarySettingsId = registerAncillarySettings( ancillarySettings );
    const unsigned int dependentVariableLayoutId = registerDependentVariableLayout( dependentVariableBookkeeping );

    const unsigned int setId = setMetadata_.size( );
    setMetadata_.push_back(
            { observableType, linkDefinitionId, referenceLinkEnd, observableSize, ancillarySettingsId, dependentVariableLayoutId } );

    observationIdsBySet_.push_back( std::vector< unsigned int >( ) );
    observationIdsBySet_.back( ).reserve( preparedObservations.size( ) );

    addObservationsToSet( setId, preparedObservations, preparedTimes, preparedDependentVariables,
                          preparedWeights, preparedResiduals, false );
    ++structuralVersion_;
    return setId;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
int ObservationDataset< ObservationScalarType, TimeType, Dummy >::addObservationSetFromDataset(
        const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
        const unsigned int sourceSetId )
{
    if( &sourceDataset == this )
    {
        const ObservationDataset< ObservationScalarType, TimeType > sourceSnapshot( sourceDataset );
        return addObservationSetFromDataset( sourceSnapshot, sourceSetId );
    }

    const ObservationSetMetadata< ObservationScalarType, TimeType >& sourceMetadata =
            sourceDataset.getObservationSetMetadata( sourceSetId );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& sourceBookkeeping =
            sourceDataset.getDependentVariableBookkeeping( sourceMetadata.dependentVariableLayoutId_ );
    const std::shared_ptr< ObservationAncillarySimulationSettings >& sourceAncillarySettings =
            sourceDataset.getAncillarySettings( sourceMetadata.ancillarySettingsId_ );
    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > copiedBookkeeping = sourceBookkeeping == nullptr
            ? nullptr
            : sourceBookkeeping->clone( );
    const std::shared_ptr< ObservationAncillarySimulationSettings > copiedAncillarySettings = sourceAncillarySettings == nullptr
            ? nullptr
            : std::make_shared< ObservationAncillarySimulationSettings >( *sourceAncillarySettings );
    const unsigned int newSetId = addObservationSet( sourceMetadata.observableType_,
                                                     sourceDataset.getLinkDefinition( sourceMetadata.linkDefinitionId_ ),
                                                     sourceDataset.getObservationsForSet( sourceSetId ),
                                                     sourceDataset.getObservationTimesForSet( sourceSetId ),
                                                     sourceMetadata.referenceLinkEnd_,
                                                     sourceDataset.getDependentVariablesForSet( sourceSetId ),
                                                     copiedBookkeeping,
                                                     copiedAncillarySettings,
                                                     sourceDataset.getWeightsForSet( sourceSetId ),
                                                     sourceDataset.getResidualsForSet( sourceSetId ) );

    const auto& sourceIds = sourceDataset.getObservationIdsForSet( sourceSetId );
    const auto& targetIds = getObservationIdsForSet( newSetId );
    for( std::size_t i = 0; i < sourceIds.size( ); ++i )
    {
        const auto& source = sourceDataset.getObservationRow( sourceIds.at( i ) );
        auto& target = mutableObservationRow( targetIds.at( i ) );
        target.isActive_ = source.isActive_;
        target.rejectionReason_ = source.rejectionReason_;
    }
    observationWeights_.copyBlock( sourceDataset.observationWeights_.restricted(
        sourceDataset.getScalarComponentIdsForObservationSelection( sourceIds, {} ) ),
        getScalarComponentIdsForObservationSelection( targetIds, {} ) );
    return newSetId;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
int ObservationDataset< ObservationScalarType, TimeType, Dummy >::addObservationSetWithWeights(
        const ObservableType observableType,
        const LinkDefinition& linkDefinition,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const LinkEndType referenceLinkEnd,
        const ObservationWeightSettings& weightSettings,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
{
    const auto weights = ObservationWeights::forSet( observations.size( ), getObservableSize( observableType ), weightSettings );
    const unsigned int setId = addObservationSet( observableType, linkDefinition, observations, times, referenceLinkEnd,
        dependentVariables, dependentVariableBookkeeping, ancillarySettings, {}, residuals );
    observationWeights_.copyBlock( weights, getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ) );
    return setId;
}

// Replace every vector-valued observation in one set without changing row metadata.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setObservationsForSet(
        const unsigned int setId,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations )
{
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );

    // The replacement is row-for-row; structural changes go through add/remove/rebuild helpers.
    if( observations.size( ) != observationIds.size( ) )
    {
        throw std::runtime_error( "Error when setting dataset observations, number of observations is inconsistent." );
    }

    for( const auto& value : observations )
    {
        if( value.size( ) != static_cast< Eigen::Index >( getObservationSetMetadata( setId ).observableSize_ ) )
        {
            throw std::runtime_error( "Observation component size is inconsistent." );
        }
    }

    // Delegate each row assignment to the scalar-component aware setter.
    for( std::size_t i = 0; i < observations.size( ); ++i )
    {
        setObservationValue( observationIds.at( i ), observations.at( i ) );
    }
}

// Replace every vector-valued residual in one set without changing row metadata.
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::setResidualsForSet(
        const unsigned int setId,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
{
    const std::vector< unsigned int >& observationIds = observationIdsBySet_.at( setId );

    // The replacement is row-for-row; structural changes go through add/remove/rebuild helpers.
    if( residuals.size( ) != observationIds.size( ) )
    {
        throw std::runtime_error( "Error when setting dataset residuals, number of observations is inconsistent." );
    }

    for( const auto& value : residuals )
    {
        if( value.size( ) != static_cast< Eigen::Index >( getObservationSetMetadata( setId ).observableSize_ ) )
        {
            throw std::runtime_error( "Observation component size is inconsistent." );
        }
    }

    // Delegate each row assignment to the scalar-component aware setter.
    for( std::size_t i = 0; i < residuals.size( ); ++i )
    {
        setResidualValue( observationIds.at( i ), residuals.at( i ) );
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::addObservationsToSet(
        const unsigned int setId,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
        const bool sortObservations )
{
    validateObservationSetData( setId, observations, times, dependentVariables, weights, residuals );
    const auto& metadata = getObservationSetMetadata( setId );
    const auto& ids = observationIdsBySet_.at( setId );
    if( !observations.empty( ) && !ids.empty( ) &&
        ( getObservationRow( ids.front( ) ).dependentVariableValues_.size( ) != 0 ) != !dependentVariables.empty( ) )
    {
        throw std::runtime_error( "Dependent variables must be present for all observations in a set or none." );
    }
    if( observations.size( ) > std::numeric_limits< unsigned int >::max( ) - nextObservationId_ ||
        observations.size( ) > ( std::numeric_limits< unsigned int >::max( ) - observedValues_.size( ) ) / metadata.observableSize_ )
    {
        throw std::overflow_error( "Observation identity or scalar storage capacity exceeded." );
    }
    for( std::size_t i = 0; i < observations.size( ); ++i )
    {
        const unsigned int id = nextObservationId_++;
        const unsigned int first = observedValues_.size( );
        rowPositionById_.emplace( id, observationRows_.size( ) );
        observationRows_.push_back( { id, times.at( i ), setId, first, metadata.observableSize_,
            static_cast< unsigned int >( observationIdsBySet_.at( setId ).size( ) ),
            dependentVariables.empty( ) ? Eigen::VectorXd( ) : dependentVariables.at( i ), true, "" } );
        observationIdsBySet_.at( setId ).push_back( id );
        for( unsigned int component = 0; component < metadata.observableSize_; ++component )
        {
            observedValues_.push_back( observations.at( i )( component ) );
            residualValues_.push_back( residuals.empty( ) ? ObservationScalarType( 0 ) : residuals.at( i )( component ) );
        }
        observationWeights_.appendDiagonal( weights.empty( ) ? Eigen::VectorXd::Ones( metadata.observableSize_ ) : weights.at( i ) );
    }
    if( sortObservations )
    {
        sortObservationIdsForSet( setId );
    }
    if( !observations.empty( ) )
    {
        ++structuralVersion_;
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::removeObservationsFromSet( const unsigned int setId,
                                                                                              std::vector< unsigned int > indicesToRemove )
{
    const auto& ids = observationIdsBySet_.at( setId );
    std::unordered_set< unsigned int > removed;
    for( const unsigned int index : indicesToRemove )
    {
        if( index >= ids.size( ) )
        {
            throw std::runtime_error( "Observation removal index is out of bounds." );
        }
        removed.insert( ids.at( index ) );
    }
    std::vector< unsigned int > retained;
    for( const auto& row : observationRows_ )
    {
        if( !removed.count( row.observationId_ ) )
        {
            retained.push_back( row.observationId_ );
        }
    }
    retainObservationRows( retained );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::pair< TimeType, TimeType > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTimeBoundsForSet(
        const unsigned int setId ) const
{
    const std::vector< TimeType > observationTimes = getObservationTimesForSet( setId );
    if( observationTimes.empty( ) )
    {
        return std::make_pair( TUDAT_NAN, TUDAT_NAN );
    }
    return std::make_pair( *std::min_element( observationTimes.begin( ), observationTimes.end( ) ),
                           *std::max_element( observationTimes.begin( ), observationTimes.end( ) ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::pair< TimeType, TimeType > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTimeBounds( ) const
{
    if( observationRows_.empty( ) )
    {
        return std::make_pair( TUDAT_NAN, TUDAT_NAN );
    }

    TimeType startTime = observationRows_.front( ).time_;
    TimeType endTime = observationRows_.front( ).time_;
    for( const ObservationDatasetRow< TimeType >& observationRow : observationRows_ )
    {
        startTime = std::min( startTime, observationRow.time_ );
        endTime = std::max( endTime, observationRow.time_ );
    }
    return std::make_pair( startTime, endTime );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getComputedObservationsForSet( const unsigned int setId ) const
{
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations = getObservationsForSet( setId );
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = getResidualsForSet( setId );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > computedObservations;
    computedObservations.reserve( observations.size( ) );
    for( std::size_t i = 0; i < observations.size( ); ++i )
    {
        computedObservations.push_back( observations.at( i ) - residuals.at( i ) );
    }
    return computedObservations;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getRmsResidualsForSet( const unsigned int setId ) const
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const std::size_t numberOfObservations = getNumberOfObservationsForSet( setId );
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = getResidualsForSet( setId );
    return observation_dataset_detail::accumulateResidualStatistic( observableSize, numberOfObservations, residuals, true );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getMeanResidualsForSet( const unsigned int setId ) const
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const std::size_t numberOfObservations = getNumberOfObservationsForSet( setId );
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = getResidualsForSet( setId );
    return observation_dataset_detail::accumulateResidualStatistic( observableSize, numberOfObservations, residuals, false );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getNumberOfObservationSets( ) const
{
    return setMetadata_.size( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getNumberOfObservations( ) const
{
    return observationRows_.size( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::size_t ObservationDataset< ObservationScalarType, TimeType, Dummy >::getTotalScalarSize( ) const
{
    return observedValues_.size( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > >&
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetMetadata( ) const
{
    return setMetadata_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const ObservationSetMetadata< ObservationScalarType, TimeType >&
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationSetMetadata( const unsigned int setId ) const
{
    return setMetadata_.at( setId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::vector< ObservationDatasetRow< TimeType > >& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationRows( )
        const
{
    return observationRows_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const ObservationDatasetRow< TimeType >& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationRow(
        const unsigned int observationId ) const
{
    return observationRows_.at( rowPositionById_.at( observationId ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< ObservationScalarComponentRow > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponentRows( )
        const
{
    std::vector< ObservationScalarComponentRow > result;
    result.reserve( getTotalScalarSize( ) );
    for( const auto& row : observationRows_ )
    {
        for( unsigned int component = 0; component < row.scalarSize_; ++component )
        {
            result.push_back( { row.observationId_, component } );
        }
    }
    return result;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
ObservationScalarComponentRow ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponentRow(
        const unsigned int scalarComponentId ) const
{
    if( scalarComponentId >= getTotalScalarSize( ) )
    {
        throw std::out_of_range( "Observation scalar component is out of bounds." );
    }
    const auto next = std::upper_bound( observationRows_.begin( ), observationRows_.end( ), scalarComponentId,
        []( unsigned int scalar, const ObservationDatasetRow< TimeType >& row ) { return scalar < row.firstScalarComponent_; } );
    const auto& row = *std::prev( next );
    return { row.observationId_, scalarComponentId - row.firstScalarComponent_ };
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::vector< unsigned int >& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationIdsForSet(
        const unsigned int setId ) const
{
    return observationIdsBySet_.at( setId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationsForSet( const unsigned int setId ) const
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        observations.push_back( getObservationValue( observationId ) );
    }
    return observations;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationValue(
        const unsigned int observationId ) const
{
    const ObservationDatasetRow< TimeType >& row = getObservationRow( observationId );
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > value =
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( row.scalarSize_ );
    for( unsigned int i = 0; i < row.scalarSize_; ++i )
    {
        value( i ) = observedValues_.at( row.firstScalarComponent_ + i );
    }
    return value;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< TimeType > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationTimesForSet(
        const unsigned int setId ) const
{
    std::vector< TimeType > times;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        times.push_back( getObservationRow( observationId ).time_ );
    }
    return times;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
TimeType ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationTime( const unsigned int observationId ) const
{
    return getObservationRow( observationId ).time_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightsForSet(
        const unsigned int setId ) const
{
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        weights.push_back( getWeightValue( observationId ) );
    }
    return weights;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightValue(
        const unsigned int observationId ) const
{
    return observationWeights_.getDiagonal( getScalarComponentIdsForObservationSelection( { observationId }, {} ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightMatrixForObservation(
        const unsigned int observationId ) const
{
    return observationWeights_.restricted( getScalarComponentIdsForObservationSelection( { observationId }, {} ) ).sparseMatrix( ).toDense( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightMatrixForSet( const unsigned int setId ) const
{
    return observationWeights_.restricted( getScalarComponentIdsForObservationSelection( observationIdsBySet_.at( setId ), {} ) ).sparseMatrix( ).toDense( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::getResidualsForSet( const unsigned int setId ) const
{
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        residuals.push_back( getResidualValue( observationId ) );
    }
    return residuals;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getResidualValue(
        const unsigned int observationId ) const
{
    const ObservationDatasetRow< TimeType >& row = getObservationRow( observationId );
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > value =
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( row.scalarSize_ );
    for( unsigned int i = 0; i < row.scalarSize_; ++i )
    {
        value( i ) = residualValues_.at( row.firstScalarComponent_ + i );
    }
    return value;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::vector< Eigen::VectorXd > ObservationDataset< ObservationScalarType, TimeType, Dummy >::getDependentVariablesForSet(
        const unsigned int setId ) const
{
    std::vector< Eigen::VectorXd > dependentVariables;
    bool hasNonEmptyDependentVariables = false;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        const Eigen::VectorXd dependentVariable = getObservationRow( observationId ).dependentVariableValues_;
        if( dependentVariable.size( ) > 0 )
        {
            hasNonEmptyDependentVariables = true;
        }
        dependentVariables.push_back( dependentVariable );
    }
    return hasNonEmptyDependentVariables ? dependentVariables : std::vector< Eigen::VectorXd >( );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getDependentVariables(
        const unsigned int observationId ) const
{
    return getObservationRow( observationId ).dependentVariableValues_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSingleDependentVariableForSet(
        const unsigned int setId,
        const std::pair< int, int >& dependentVariableIndexAndSize ) const
{
    const std::vector< Eigen::VectorXd > observationsDependentVariables = getDependentVariablesForSet( setId );
    if( observationsDependentVariables.empty( ) )
    {
        throw std::runtime_error(
                "Error when retrieving single observation dependent variable, the set has no dependent-variable values." );
    }
    Eigen::MatrixXd singleDependentVariable =
            Eigen::MatrixXd::Zero( getNumberOfObservationsForSet( setId ), dependentVariableIndexAndSize.second );
    for( unsigned int i = 0; i < observationsDependentVariables.size( ); ++i )
    {
        if( dependentVariableIndexAndSize.first + dependentVariableIndexAndSize.second > observationsDependentVariables.at( i ).size( ) )
        {
            throw std::runtime_error(
                    "Error when retrieving single observation dependent variable, required index and size incompatible with dependent "
                    "variables size." );
        }
        Eigen::VectorXd singleDependentVariableVector =
                observationsDependentVariables.at( i ).segment( dependentVariableIndexAndSize.first, dependentVariableIndexAndSize.second );
        singleDependentVariable.block( i, 0, 1, dependentVariableIndexAndSize.second ) = singleDependentVariableVector.transpose( );
    }
    return singleDependentVariable;
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETMUTATIONIMPLEMENTATION_H
