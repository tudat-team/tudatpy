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

template< typename ObservationScalarType, typename TimeType >
int getTotalScalarSizeForNewObservationSet( const ObservableType observableType,
                                            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations )
{
    if( observations.empty( ) )
    {
        return 0;
    }

    const int observableSize = observations.front( ).size( );
    if( observableSize != getObservableSize( observableType ) )
    {
        throw std::runtime_error( "Error when adding observation set with weight block, observable size is inconsistent." );
    }

    for( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observation : observations )
    {
        if( observation.size( ) != observableSize )
        {
            throw std::runtime_error( "Error when adding observation set with weight block, scalar component size is inconsistent." );
        }
    }
    return static_cast< int >( observations.size( ) ) * observableSize;
}

template< typename ObservationScalarType, typename TimeType >
void validateObservationWeightBlock( const ObservableType observableType,
                                     const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                     const Eigen::MatrixXd& weightBlock )
{
    const int observableSize = observations.empty( ) ? getObservableSize( observableType ) : observations.front( ).size( );
    getTotalScalarSizeForNewObservationSet< ObservationScalarType, TimeType >( observableType, observations );
    if( weightBlock.rows( ) != observableSize || weightBlock.cols( ) != observableSize )
    {
        throw std::runtime_error( "Error when adding observation set with weight block, matrix size is inconsistent." );
    }
}

template< typename ObservationScalarType, typename TimeType >
void validateObservationWeightBlocks( const ObservableType observableType,
                                      const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                      const std::vector< Eigen::MatrixXd >& weightBlocks )
{
    if( weightBlocks.size( ) != observations.size( ) )
    {
        throw std::runtime_error( "Error when adding observation set with weight blocks, weight count is inconsistent." );
    }
    getTotalScalarSizeForNewObservationSet< ObservationScalarType, TimeType >( observableType, observations );
    for( std::size_t i = 0; i < observations.size( ); ++i )
    {
        if( weightBlocks.at( i ).rows( ) != observations.at( i ).size( ) || weightBlocks.at( i ).cols( ) != observations.at( i ).size( ) )
        {
            throw std::runtime_error( "Error when adding observation set with weight blocks, matrix size is inconsistent." );
        }
    }
}

template< typename ObservationScalarType, typename TimeType >
void validateSetWeightBlock( const ObservableType observableType,
                             const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                             const Eigen::MatrixXd& setWeightBlock )
{
    const int totalScalarSize = getTotalScalarSizeForNewObservationSet< ObservationScalarType, TimeType >( observableType, observations );
    if( setWeightBlock.rows( ) != totalScalarSize || setWeightBlock.cols( ) != totalScalarSize )
    {
        throw std::runtime_error( "Error when adding observation set with set weight block, matrix size is inconsistent." );
    }
}

template< typename ObservationScalarType, typename TimeType >
void validateObservationWeightSettings( const ObservableType observableType,
                                        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                        const ObservationWeightSettings& weightSettings )
{
    switch( weightSettings.type_ )
    {
        case ObservationWeightSettings::Type::default_weights:
        case ObservationWeightSettings::Type::constant_scalar:
            return;
        case ObservationWeightSettings::Type::scalar_per_observation:
            if( weightSettings.scalarWeights_.size( ) != observations.size( ) )
            {
                throw std::runtime_error( "Error when adding observation set with scalar weights, weight count is inconsistent." );
            }
            return;
        case ObservationWeightSettings::Type::constant_block:
            validateObservationWeightBlock< ObservationScalarType, TimeType >( observableType, observations, weightSettings.weightBlock_ );
            return;
        case ObservationWeightSettings::Type::block_per_observation:
            validateObservationWeightBlocks< ObservationScalarType, TimeType >(
                    observableType, observations, weightSettings.weightBlocks_ );
            return;
        case ObservationWeightSettings::Type::set_block:
            validateSetWeightBlock< ObservationScalarType, TimeType >( observableType, observations, weightSettings.weightBlock_ );
            return;
        default:
            throw std::runtime_error( "Error when adding observation set with weights, weight settings type is invalid." );
    }
}

template< typename DatasetType, typename ObservationScalarType, typename TimeType, typename ApplyWeightsFunction >
int addObservationSetAndApplyWeights(
        DatasetType& dataset,
        const ObservableType observableType,
        const LinkDefinition& linkDefinition,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType >& times,
        const LinkEndType referenceLinkEnd,
        const std::vector< Eigen::VectorXd >& dependentVariables,
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping,
        const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals,
        ApplyWeightsFunction applyWeights )
{
    const unsigned int setId = dataset.addObservationSet( observableType,
                                                          linkDefinition,
                                                          observations,
                                                          times,
                                                          referenceLinkEnd,
                                                          dependentVariables,
                                                          dependentVariableBookkeeping,
                                                          ancillarySettings,
                                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                          residuals );
    applyWeights( setId );
    return setId;
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

    const unsigned int observableSize = preparedObservations.empty( ) ? static_cast< unsigned int >( getObservableSize( observableType ) )
                                                                      : preparedObservations.at( 0 ).size( );
    if( observableSize == 0 )
    {
        throw std::runtime_error( "Error when adding observation set to dataset, observable size is zero." );
    }

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

    for( std::size_t i = 0; i < preparedObservations.size( ); ++i )
    {
        if( preparedObservations.at( i ).size( ) != static_cast< int >( observableSize ) ||
            ( !preparedResiduals.empty( ) && preparedResiduals.at( i ).size( ) != static_cast< int >( observableSize ) ) ||
            ( !preparedWeights.empty( ) && preparedWeights.at( i ).size( ) != static_cast< int >( observableSize ) ) )
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
        std::vector< unsigned int > indicesToRemove;
        for( unsigned int i = 1; i < preparedTimes.size( ); ++i )
        {
            if( preparedTimes.at( i ) == preparedTimes.at( i - 1 ) )
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

    const unsigned int linkDefinitionId = registerLinkDefinition( linkDefinition );
    const unsigned int ancillarySettingsId = registerAncillarySettings( ancillarySettings );
    const unsigned int dependentVariableLayoutId = registerDependentVariableLayout( dependentVariableBookkeeping );

    const unsigned int setId = setMetadata_.size( );
    setMetadata_.push_back(
            { observableType, linkDefinitionId, referenceLinkEnd, observableSize, ancillarySettingsId, dependentVariableLayoutId } );

    observationIdsBySet_.push_back( std::vector< unsigned int >( ) );
    observationIdsBySet_.back( ).reserve( preparedObservations.size( ) );

    for( std::size_t i = 0; i < preparedObservations.size( ); ++i )
    {
        const unsigned int observationId = observationRows_.size( );
        const unsigned int firstScalarComponent = scalarComponentRows_.size( );
        const Eigen::VectorXd dependentVariablesForObservation =
                preparedDependentVariables.empty( ) ? Eigen::VectorXd( ) : preparedDependentVariables.at( i );

        observationRows_.push_back( { preparedTimes.at( i ),
                                      setId,
                                      firstScalarComponent,
                                      observableSize,
                                      static_cast< unsigned int >( i ),
                                      dependentVariablesForObservation,
                                      true,
                                      "" } );
        observationIdsBySet_.back( ).push_back( observationId );

        for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
        {
            scalarComponentRows_.push_back( { observationId, componentIndex } );
            observedValues_.push_back( preparedObservations.at( i )( componentIndex ) );
            residualValues_.push_back( preparedResiduals.empty( ) ? ObservationScalarType( 0 )
                                                                  : preparedResiduals.at( i )( componentIndex ) );
        }

        if( preparedWeights.empty( ) )
        {
            observationWeights_.appendScalarWeight( 1.0 );
        }
        else
        {
            observationWeights_.appendDiagonalWeightVector( preparedWeights.at( i ) );
        }
    }

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
    const ObservationSetMetadata< ObservationScalarType, TimeType >& sourceMetadata =
            sourceDataset.getObservationSetMetadata( sourceSetId );
    const unsigned int newSetId =
            addObservationSet( sourceMetadata.observableType_,
                               sourceDataset.getLinkDefinition( sourceMetadata.linkDefinitionId_ ),
                               sourceDataset.getObservationsForSet( sourceSetId ),
                               sourceDataset.getObservationTimesForSet( sourceSetId ),
                               sourceMetadata.referenceLinkEnd_,
                               sourceDataset.getDependentVariablesForSet( sourceSetId ),
                               sourceDataset.getDependentVariableBookkeeping( sourceMetadata.dependentVariableLayoutId_ ),
                               sourceDataset.getAncillarySettings( sourceMetadata.ancillarySettingsId_ ),
                               sourceDataset.getWeightsForSet( sourceSetId ),
                               sourceDataset.getResidualsForSet( sourceSetId ) );
    if( sourceDataset.hasWeightMatrixForSet( sourceSetId ) )
    {
        setWeightMatrixForSet( newSetId, sourceDataset.getWeightMatrixForSet( sourceSetId ) );
    }
    else
    {
        const std::vector< unsigned int >& sourceObservationIds = sourceDataset.getObservationIdsForSet( sourceSetId );
        const std::vector< unsigned int >& targetObservationIds = getObservationIdsForSet( newSetId );
        for( std::size_t i = 0; i < sourceObservationIds.size( ); ++i )
        {
            observationWeights_.setObservationWeight(
                    targetObservationIds.at( i ), sourceDataset.observationWeights_.getObservationWeight( sourceObservationIds.at( i ) ) );
        }
    }
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
    // Validate the weight policy before adding the set, so invalid weights cannot leave a partial set behind.
    observation_dataset_detail::validateObservationWeightSettings< ObservationScalarType, TimeType >(
            observableType, observations, weightSettings );

    // Add observations through the normal path, then apply the selected weight policy to the new set.
    return observation_dataset_detail::addObservationSetAndApplyWeights(
            *this,
            observableType,
            linkDefinition,
            observations,
            times,
            referenceLinkEnd,
            dependentVariables,
            dependentVariableBookkeeping,
            ancillarySettings,
            residuals,
            [ & ]( const unsigned int setId ) {
                const std::vector< unsigned int >& observationIds = getObservationIdsForSet( setId );
                switch( weightSettings.type_ )
                {
                    case ObservationWeightSettings::Type::default_weights:
                        return;
                    case ObservationWeightSettings::Type::constant_scalar:
                        setConstantSingleObservationScalarWeightForSet( setId, weightSettings.scalarWeight_ );
                        return;
                    case ObservationWeightSettings::Type::scalar_per_observation:
                        for( std::size_t i = 0; i < weightSettings.scalarWeights_.size( ); ++i )
                        {
                            observationWeights_.setScalarWeight( observationIds.at( i ), weightSettings.scalarWeights_.at( i ) );
                        }
                        return;
                    case ObservationWeightSettings::Type::constant_block:
                        for( const unsigned int observationId : observationIds )
                        {
                            setWeightMatrixForObservation( observationId, weightSettings.weightBlock_ );
                        }
                        return;
                    case ObservationWeightSettings::Type::block_per_observation:
                        for( std::size_t i = 0; i < weightSettings.weightBlocks_.size( ); ++i )
                        {
                            setWeightMatrixForObservation( observationIds.at( i ), weightSettings.weightBlocks_.at( i ) );
                        }
                        return;
                    case ObservationWeightSettings::Type::set_block:
                        setWeightMatrixForSet( setId, weightSettings.weightBlock_ );
                        return;
                    default:
                        throw std::runtime_error( "Error when applying observation weight settings, weight settings type is invalid." );
                }
            } );
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
    if( hasWeightMatrixForSet( setId ) && !observations.empty( ) )
    {
        throw std::runtime_error(
                "Error when appending observations to dataset, the target set has a full set-level weight block. "
                "Provide a new complete set block by replacing the set data explicitly." );
    }

    std::vector< unsigned int > updatedSourceObservationIds = observationIdsBySet_.at( setId );
    validateObservationSetData( setId, observations, times, dependentVariables, weights, residuals );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedObservations = getObservationsForSet( setId );
    std::vector< TimeType > updatedTimes = getObservationTimesForSet( setId );
    std::vector< Eigen::VectorXd > updatedDependentVariables = getDependentVariablesForSet( setId );
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > updatedWeights = getWeightsForSet( setId );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedResiduals = getResidualsForSet( setId );

    updatedObservations.insert( updatedObservations.end( ), observations.begin( ), observations.end( ) );
    updatedTimes.insert( updatedTimes.end( ), times.begin( ), times.end( ) );
    updatedSourceObservationIds.insert( updatedSourceObservationIds.end( ), observations.size( ), invalidObservationId( ) );
    updatedWeights.insert( updatedWeights.end( ), weights.begin( ), weights.end( ) );
    updatedResiduals.insert( updatedResiduals.end( ), residuals.begin( ), residuals.end( ) );

    if( !dependentVariables.empty( ) || !updatedDependentVariables.empty( ) )
    {
        if( updatedDependentVariables.empty( ) )
        {
            updatedDependentVariables.resize( getNumberOfObservationsForSet( setId ) );
        }
        if( dependentVariables.empty( ) )
        {
            updatedDependentVariables.resize( updatedObservations.size( ) );
        }
        else
        {
            updatedDependentVariables.insert( updatedDependentVariables.end( ), dependentVariables.begin( ), dependentVariables.end( ) );
        }
    }

    if( weights.empty( ) )
    {
        const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
        for( std::size_t i = 0; i < observations.size( ); ++i )
        {
            updatedWeights.push_back( Eigen::Matrix< double, Eigen::Dynamic, 1 >::Ones( observableSize, 1 ) );
        }
    }
    if( residuals.empty( ) )
    {
        const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
        for( std::size_t i = 0; i < observations.size( ); ++i )
        {
            updatedResiduals.push_back( Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( observableSize, 1 ) );
        }
    }

    if( sortObservations && updatedTimes.size( ) > 1 )
    {
        const std::vector< std::size_t > permutation = getTimeSortingPermutation( updatedTimes );
        reorderVector( updatedObservations, permutation );
        reorderVector( updatedTimes, permutation );
        reorderVector( updatedSourceObservationIds, permutation );
        reorderVector( updatedWeights, permutation );
        reorderVector( updatedResiduals, permutation );
        if( !updatedDependentVariables.empty( ) )
        {
            reorderVector( updatedDependentVariables, permutation );
        }
    }

    replaceObservationSetDataWithSourceRows( setId,
                                             updatedObservations,
                                             updatedTimes,
                                             updatedDependentVariables,
                                             updatedWeights,
                                             updatedResiduals,
                                             updatedSourceObservationIds );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::removeObservationsFromSet( const unsigned int setId,
                                                                                              std::vector< unsigned int > indicesToRemove )
{
    std::sort( indicesToRemove.begin( ), indicesToRemove.end( ) );
    indicesToRemove.erase( std::unique( indicesToRemove.begin( ), indicesToRemove.end( ) ), indicesToRemove.end( ) );

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedObservations = getObservationsForSet( setId );
    std::vector< TimeType > updatedTimes = getObservationTimesForSet( setId );
    std::vector< Eigen::VectorXd > updatedDependentVariables = getDependentVariablesForSet( setId );
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > updatedWeights = getWeightsForSet( setId );
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedResiduals = getResidualsForSet( setId );
    std::vector< unsigned int > updatedSourceObservationIds = observationIdsBySet_.at( setId );

    for( std::vector< unsigned int >::reverse_iterator indexIterator = indicesToRemove.rbegin( ); indexIterator != indicesToRemove.rend( );
         ++indexIterator )
    {
        const unsigned int indexToRemove = *indexIterator;
        if( indexToRemove >= updatedObservations.size( ) )
        {
            throw std::runtime_error( "Error when removing observations from dataset, index is out of bounds." );
        }
        updatedObservations.erase( updatedObservations.begin( ) + indexToRemove );
        updatedTimes.erase( updatedTimes.begin( ) + indexToRemove );
        updatedSourceObservationIds.erase( updatedSourceObservationIds.begin( ) + indexToRemove );
        updatedWeights.erase( updatedWeights.begin( ) + indexToRemove );
        updatedResiduals.erase( updatedResiduals.begin( ) + indexToRemove );
        if( !updatedDependentVariables.empty( ) )
        {
            updatedDependentVariables.erase( updatedDependentVariables.begin( ) + indexToRemove );
        }
    }

    replaceObservationSetDataWithSourceRows( setId,
                                             updatedObservations,
                                             updatedTimes,
                                             updatedDependentVariables,
                                             updatedWeights,
                                             updatedResiduals,
                                             updatedSourceObservationIds );
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
    Eigen::VectorXd rmsResiduals = Eigen::VectorXd::Zero( observableSize );
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = getResidualsForSet( setId );

    for( unsigned int i = 0; i < observableSize; ++i )
    {
        for( std::size_t j = 0; j < numberOfObservations; ++j )
        {
            rmsResiduals[ i ] += residuals.at( j )( i, 0 ) * residuals.at( j )( i, 0 );
        }
        rmsResiduals[ i ] = std::sqrt( rmsResiduals[ i ] / static_cast< double >( numberOfObservations ) );
    }
    return rmsResiduals;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::VectorXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getMeanResidualsForSet( const unsigned int setId ) const
{
    const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
    const std::size_t numberOfObservations = getNumberOfObservationsForSet( setId );
    Eigen::VectorXd meanResiduals = Eigen::VectorXd::Zero( observableSize );
    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals = getResidualsForSet( setId );

    for( unsigned int i = 0; i < observableSize; ++i )
    {
        for( std::size_t j = 0; j < numberOfObservations; ++j )
        {
            meanResiduals[ i ] += residuals.at( j )( i, 0 );
        }
        meanResiduals[ i ] /= static_cast< double >( numberOfObservations );
    }
    return meanResiduals;
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
    return scalarComponentRows_.size( );
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
    return observationRows_.at( observationId );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const std::vector< ObservationScalarComponentRow >& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponentRows( )
        const
{
    return scalarComponentRows_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
const ObservationScalarComponentRow& ObservationDataset< ObservationScalarType, TimeType, Dummy >::getScalarComponentRow(
        const unsigned int scalarComponentId ) const
{
    return scalarComponentRows_.at( scalarComponentId );
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
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
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
        times.push_back( observationRows_.at( observationId ).time_ );
    }
    return times;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
TimeType ObservationDataset< ObservationScalarType, TimeType, Dummy >::getObservationTime( const unsigned int observationId ) const
{
    return observationRows_.at( observationId ).time_;
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
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
    return observationWeights_.getObservationWeightVector( observationId, row.scalarSize_ );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightMatrixForObservation(
        const unsigned int observationId ) const
{
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
    return observationWeights_.getObservationWeightMatrix( observationId, row.scalarSize_ );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getWeightMatrixForSet( const unsigned int setId ) const
{
    if( observationWeights_.hasSetWeightBlock( setId ) )
    {
        return observationWeights_.getSetWeightBlock( setId );
    }

    Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Zero( getTotalScalarSizeForSet( setId ), getTotalScalarSizeForSet( setId ) );
    std::size_t currentIndex = 0;
    for( const unsigned int observationId : observationIdsBySet_.at( setId ) )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        weightMatrix.block( currentIndex, currentIndex, row.scalarSize_, row.scalarSize_ ) =
                observationWeights_.getObservationWeightMatrix( observationId, row.scalarSize_ );
        currentIndex += row.scalarSize_;
    }
    return weightMatrix;
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
    const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
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
        const Eigen::VectorXd dependentVariable = observationRows_.at( observationId ).dependentVariableValues_;
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
    return observationRows_.at( observationId ).dependentVariableValues_;
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
Eigen::MatrixXd ObservationDataset< ObservationScalarType, TimeType, Dummy >::getSingleDependentVariableForSet(
        const unsigned int setId,
        const std::pair< int, int >& dependentVariableIndexAndSize ) const
{
    const std::vector< Eigen::VectorXd > observationsDependentVariables = getDependentVariablesForSet( setId );
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
