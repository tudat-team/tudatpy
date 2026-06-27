/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_DATASET_H
#define TUDAT_OBSERVATION_DATASET_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"
#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace observation_models
{

using ObservationSetId = std::size_t;
using ObservationId = std::size_t;
using ScalarComponentId = std::size_t;
using LinkDefinitionId = std::size_t;
using AncillarySettingsId = std::size_t;
using DependentVariableLayoutId = std::size_t;

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
struct ObservationSetMetadata {
    ObservableType observableType_;
    LinkDefinitionId linkDefinitionId_;
    LinkEndType referenceLinkEnd_;
    unsigned int observableSize_;
    AncillarySettingsId ancillarySettingsId_;
    DependentVariableLayoutId dependentVariableLayoutId_;
};

template< typename TimeType = double >
struct ObservationDatasetRow {
    TimeType time_;
    ObservationSetId setId_;
    ScalarComponentId firstScalarComponent_;
    unsigned int scalarSize_;
    unsigned int indexInSet_;
    bool isActive_;
    std::string rejectionReason_;
};

struct ObservationScalarComponentRow {
    ObservationId observationId_;
    unsigned int componentIndex_;
};

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationDataset;

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class EstimationVectorProjection
{
public:
    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& getObservationVector( ) const
    {
        return observations_;
    }

    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& getResidualVector( ) const
    {
        return residuals_;
    }

    const Eigen::VectorXd& getWeightVector( ) const
    {
        return weights_;
    }

    const std::vector< TimeType >& getTimes( ) const
    {
        return times_;
    }

    const std::vector< ObservationId >& getObservationIds( ) const
    {
        return observationIds_;
    }

    const std::vector< ObservationSetId >& getSetIds( ) const
    {
        return setIds_;
    }

    const std::vector< ScalarComponentId >& getScalarComponentIds( ) const
    {
        return scalarComponentIds_;
    }

private:
    template< typename DatasetObservationScalarType,
              typename DatasetTimeType,
              typename std::enable_if< is_state_scalar_and_time_type< DatasetObservationScalarType, DatasetTimeType >::value, int >::type >
    friend class ObservationDataset;

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations_;
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals_;
    Eigen::VectorXd weights_;
    std::vector< TimeType > times_;
    std::vector< ObservationId > observationIds_;
    std::vector< ObservationSetId > setIds_;
    std::vector< ScalarComponentId > scalarComponentIds_;
};

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
class ObservationDataset
{
public:
    ObservationDataset( ) = default;

    ObservationSetId addObservationSet(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights =
                    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
            const bool sortObservations = false,
            const bool eraseDuplicateObservations = false )
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > preparedObservations = observations;
        std::vector< TimeType > preparedTimes = times;
        std::vector< Eigen::VectorXd > preparedDependentVariables = dependentVariables;
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > preparedWeights = weights;
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > preparedResiduals = residuals;

        const unsigned int observableSize = preparedObservations.empty( )
                ? static_cast< unsigned int >( getObservableSize( observableType ) )
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
                    const double currentObsValue = preparedTimes.at( i );
                    const double previousObsValue = preparedTimes.at( i - 1 );
                    if( std::abs( currentObsValue - previousObsValue ) <=
                        1.0e-12 * std::max( std::abs( currentObsValue ), std::abs( previousObsValue ) ) )
                    {
                        indicesToRemove.push_back( i );
                    }
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
                          << "duplicate observations when creating instance of SingleObservationSet" << std::endl;
            }
        }

        const LinkDefinitionId linkDefinitionId = registerLinkDefinition( linkDefinition );
        const AncillarySettingsId ancillarySettingsId = registerAncillarySettings( ancillarySettings );
        const DependentVariableLayoutId dependentVariableLayoutId = registerDependentVariableLayout( dependentVariableBookkeeping );

        const ObservationSetId setId = setMetadata_.size( );
        setMetadata_.push_back(
                { observableType, linkDefinitionId, referenceLinkEnd, observableSize, ancillarySettingsId, dependentVariableLayoutId } );

        observationIdsBySet_.push_back( std::vector< ObservationId >( ) );
        observationIdsBySet_.back( ).reserve( preparedObservations.size( ) );

        for( std::size_t i = 0; i < preparedObservations.size( ); ++i )
        {
            const ObservationId observationId = observationRows_.size( );
            const ScalarComponentId firstScalarComponent = scalarComponentRows_.size( );

            observationRows_.push_back(
                    { preparedTimes.at( i ), setId, firstScalarComponent, observableSize, static_cast< unsigned int >( i ), true, "" } );
            observationIdsBySet_.back( ).push_back( observationId );

            for( unsigned int componentIndex = 0; componentIndex < observableSize; ++componentIndex )
            {
                scalarComponentRows_.push_back( { observationId, componentIndex } );
                observedValues_.push_back( preparedObservations.at( i )( componentIndex ) );
                residualValues_.push_back( preparedResiduals.empty( ) ? ObservationScalarType( 0 )
                                                                      : preparedResiduals.at( i )( componentIndex ) );
                legacyWeights_.push_back( preparedWeights.empty( ) ? 1.0 : preparedWeights.at( i )( componentIndex ) );
            }

            if( preparedDependentVariables.empty( ) )
            {
                dependentVariableValues_.push_back( Eigen::VectorXd( ) );
            }
            else
            {
                dependentVariableValues_.push_back( preparedDependentVariables.at( i ) );
            }
        }

        return setId;
    }

    ObservationSetId addObservationSetFromDataset( const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
                                                   const ObservationSetId sourceSetId )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& sourceMetadata =
                sourceDataset.getObservationSetMetadata( sourceSetId );
        return addObservationSet( sourceMetadata.observableType_,
                                  sourceDataset.getLinkDefinition( sourceMetadata.linkDefinitionId_ ),
                                  sourceDataset.getObservationsForSet( sourceSetId ),
                                  sourceDataset.getObservationTimesForSet( sourceSetId ),
                                  sourceMetadata.referenceLinkEnd_,
                                  sourceDataset.getDependentVariablesForSet( sourceSetId ),
                                  sourceDataset.getDependentVariableBookkeeping( sourceMetadata.dependentVariableLayoutId_ ),
                                  sourceDataset.getAncillarySettings( sourceMetadata.ancillarySettingsId_ ),
                                  sourceDataset.getWeightsForSet( sourceSetId ),
                                  sourceDataset.getResidualsForSet( sourceSetId ) );
    }

    void setLinkDefinition( const ObservationSetId setId, const LinkDefinition& linkDefinition )
    {
        setMetadata_.at( setId ).linkDefinitionId_ = registerLinkDefinition( linkDefinition );
    }

    void setDependentVariableBookkeeping(
            const ObservationSetId setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping )
    {
        setMetadata_.at( setId ).dependentVariableLayoutId_ = registerDependentVariableLayout( dependentVariableBookkeeping );
    }

    void setObservationsForSet( const ObservationSetId setId,
                                const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations )
    {
        const std::vector< ObservationId >& observationIds = observationIdsBySet_.at( setId );
        if( observations.size( ) != observationIds.size( ) )
        {
            throw std::runtime_error( "Error when setting dataset observations, number of observations is inconsistent." );
        }
        for( std::size_t i = 0; i < observations.size( ); ++i )
        {
            setObservationValue( observationIds.at( i ), observations.at( i ) );
        }
    }

    void setObservationVectorForSet( const ObservationSetId setId,
                                     const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& observationVector )
    {
        const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
        const std::vector< ObservationId >& observationIds = observationIdsBySet_.at( setId );
        if( observationVector.size( ) != static_cast< int >( observationIds.size( ) * observableSize ) )
        {
            throw std::runtime_error( "Error when setting dataset observations, vector size is inconsistent." );
        }
        for( std::size_t i = 0; i < observationIds.size( ); ++i )
        {
            setObservationValue( observationIds.at( i ), observationVector.segment( i * observableSize, observableSize ) );
        }
    }

    void setResidualsForSet( const ObservationSetId setId,
                             const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals )
    {
        const std::vector< ObservationId >& observationIds = observationIdsBySet_.at( setId );
        if( residuals.size( ) != observationIds.size( ) )
        {
            throw std::runtime_error( "Error when setting dataset residuals, number of observations is inconsistent." );
        }
        for( std::size_t i = 0; i < residuals.size( ); ++i )
        {
            setResidualValue( observationIds.at( i ), residuals.at( i ) );
        }
    }

    void setResidualVectorForSet( const ObservationSetId setId,
                                  const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector )
    {
        const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
        const std::vector< ObservationId >& observationIds = observationIdsBySet_.at( setId );
        if( residualVector.size( ) != static_cast< int >( observationIds.size( ) * observableSize ) )
        {
            throw std::runtime_error( "Error when setting dataset residuals, vector size is inconsistent." );
        }
        for( std::size_t i = 0; i < observationIds.size( ); ++i )
        {
            setResidualValue( observationIds.at( i ), residualVector.segment( i * observableSize, observableSize ) );
        }
    }

    void setConstantWeightForSet( const ObservationSetId setId, const double weight )
    {
        setConstantWeightForSet( setId, weight * Eigen::VectorXd::Ones( getObservationSetMetadata( setId ).observableSize_ ) );
    }

    void setConstantWeightForSet( const ObservationSetId setId, const Eigen::VectorXd& weight )
    {
        if( weight.size( ) != static_cast< int >( getObservationSetMetadata( setId ).observableSize_ ) )
        {
            throw std::runtime_error( "Error when setting dataset weights, weight size is inconsistent." );
        }
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            setWeightValue( observationId, weight );
        }
    }

    void setWeightVectorForSet( const ObservationSetId setId, const Eigen::VectorXd& weightVector )
    {
        const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
        const std::vector< ObservationId >& observationIds = observationIdsBySet_.at( setId );
        if( weightVector.size( ) != static_cast< int >( observationIds.size( ) * observableSize ) )
        {
            throw std::runtime_error( "Error when setting dataset weights, vector size is inconsistent." );
        }
        for( std::size_t i = 0; i < observationIds.size( ); ++i )
        {
            setWeightValue( observationIds.at( i ), weightVector.segment( i * observableSize, observableSize ) );
        }
    }

    void replaceObservationSetData( const ObservationSetId setId,
                                    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                                    const std::vector< TimeType >& times,
                                    const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
                                    const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights =
                                            std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                    const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                                            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ) )
    {
        validateObservationSetData( setId, observations, times, dependentVariables, weights, residuals );

        ObservationDataset< ObservationScalarType, TimeType > rebuiltDataset;
        for( ObservationSetId currentSetId = 0; currentSetId < setMetadata_.size( ); ++currentSetId )
        {
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = setMetadata_.at( currentSetId );
            if( currentSetId == setId )
            {
                rebuiltDataset.addObservationSet( metadata.observableType_,
                                                  linkDefinitionRegistry_.at( metadata.linkDefinitionId_ ),
                                                  observations,
                                                  times,
                                                  metadata.referenceLinkEnd_,
                                                  dependentVariables,
                                                  dependentVariableLayoutRegistry_.at( metadata.dependentVariableLayoutId_ ),
                                                  ancillarySettingsRegistry_.at( metadata.ancillarySettingsId_ ),
                                                  weights,
                                                  residuals );
            }
            else
            {
                rebuiltDataset.addObservationSet( metadata.observableType_,
                                                  linkDefinitionRegistry_.at( metadata.linkDefinitionId_ ),
                                                  getObservationsForSet( currentSetId ),
                                                  getObservationTimesForSet( currentSetId ),
                                                  metadata.referenceLinkEnd_,
                                                  getDependentVariablesForSet( currentSetId ),
                                                  dependentVariableLayoutRegistry_.at( metadata.dependentVariableLayoutId_ ),
                                                  ancillarySettingsRegistry_.at( metadata.ancillarySettingsId_ ),
                                                  getWeightsForSet( currentSetId ),
                                                  getResidualsForSet( currentSetId ) );
            }
        }

        *this = rebuiltDataset;
    }

    void addObservationsToSet( const ObservationSetId setId,
                               const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
                               const std::vector< TimeType >& times,
                               const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
                               const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >& weights =
                                       std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                               const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                                       std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ),
                               const bool sortObservations = true )
    {
        validateObservationSetData( setId, observations, times, dependentVariables, weights, residuals );

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedObservations = getObservationsForSet( setId );
        std::vector< TimeType > updatedTimes = getObservationTimesForSet( setId );
        std::vector< Eigen::VectorXd > updatedDependentVariables = getDependentVariablesForSet( setId );
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > updatedWeights = getWeightsForSet( setId );
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedResiduals = getResidualsForSet( setId );

        updatedObservations.insert( updatedObservations.end( ), observations.begin( ), observations.end( ) );
        updatedTimes.insert( updatedTimes.end( ), times.begin( ), times.end( ) );
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
                updatedDependentVariables.insert(
                        updatedDependentVariables.end( ), dependentVariables.begin( ), dependentVariables.end( ) );
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
            reorderVector( updatedWeights, permutation );
            reorderVector( updatedResiduals, permutation );
            if( !updatedDependentVariables.empty( ) )
            {
                reorderVector( updatedDependentVariables, permutation );
            }
        }

        replaceObservationSetData( setId, updatedObservations, updatedTimes, updatedDependentVariables, updatedWeights, updatedResiduals );
    }

    void removeObservationsFromSet( const ObservationSetId setId, std::vector< unsigned int > indicesToRemove )
    {
        std::sort( indicesToRemove.begin( ), indicesToRemove.end( ) );
        indicesToRemove.erase( std::unique( indicesToRemove.begin( ), indicesToRemove.end( ) ), indicesToRemove.end( ) );

        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedObservations = getObservationsForSet( setId );
        std::vector< TimeType > updatedTimes = getObservationTimesForSet( setId );
        std::vector< Eigen::VectorXd > updatedDependentVariables = getDependentVariablesForSet( setId );
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > updatedWeights = getWeightsForSet( setId );
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > updatedResiduals = getResidualsForSet( setId );

        for( std::vector< unsigned int >::reverse_iterator indexIterator = indicesToRemove.rbegin( );
             indexIterator != indicesToRemove.rend( );
             ++indexIterator )
        {
            const unsigned int indexToRemove = *indexIterator;
            if( indexToRemove >= updatedObservations.size( ) )
            {
                throw std::runtime_error( "Error when removing observations from dataset, index is out of bounds." );
            }
            updatedObservations.erase( updatedObservations.begin( ) + indexToRemove );
            updatedTimes.erase( updatedTimes.begin( ) + indexToRemove );
            updatedWeights.erase( updatedWeights.begin( ) + indexToRemove );
            updatedResiduals.erase( updatedResiduals.begin( ) + indexToRemove );
            if( !updatedDependentVariables.empty( ) )
            {
                updatedDependentVariables.erase( updatedDependentVariables.begin( ) + indexToRemove );
            }
        }

        replaceObservationSetData( setId, updatedObservations, updatedTimes, updatedDependentVariables, updatedWeights, updatedResiduals );
    }

    void removeObservationFromSet( const ObservationSetId setId, const unsigned int indexToRemove )
    {
        removeObservationsFromSet( setId, std::vector< unsigned int >( { indexToRemove } ) );
    }

    std::vector< unsigned int > getFilteredObservationIndices( const ObservationSetId setId,
                                                               const std::shared_ptr< ObservationFilterBase >& observationFilter ) const
    {
        const unsigned int observableSize = getObservationSetMetadata( setId ).observableSize_;
        const unsigned int numberOfObservations = static_cast< unsigned int >( getNumberOfObservationsForSet( setId ) );
        const bool useOppositeCondition = observationFilter->useOppositeCondition( );

        std::vector< unsigned int > indicesToRemove;
        switch( observationFilter->getFilterType( ) )
        {
            case residual_filtering: {
                Eigen::VectorXd residualCutOff = Eigen::VectorXd::Zero( observableSize );
                if( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) != nullptr )
                {
                    residualCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter )->getFilterValue( ) *
                            Eigen::VectorXd::Ones( observableSize );
                }
                else if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter ) != nullptr )
                {
                    if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( ).size( ) !=
                        static_cast< int >( observableSize ) )
                    {
                        throw std::runtime_error(
                                "Error when performing residual filtering, size of the residual cut off vector inconsistent with "
                                "observable size." );
                    }
                    residualCutOff =
                            std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( );
                }

                for( unsigned int j = 0; j < numberOfObservations; ++j )
                {
                    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservationResidual =
                            getResidualValue( observationIdsBySet_.at( setId ).at( j ) );
                    bool removeObservation = false;
                    for( unsigned int k = 0; k < observableSize; ++k )
                    {
                        if( ( !useOppositeCondition && ( std::fabs( singleObservationResidual[ k ] ) > residualCutOff[ k ] ) ) ||
                            ( useOppositeCondition && ( std::fabs( singleObservationResidual[ k ] ) <= residualCutOff[ k ] ) ) )
                        {
                            removeObservation = true;
                        }
                    }
                    if( removeObservation )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case absolute_value_filtering: {
                Eigen::VectorXd absoluteValueCutOff = Eigen::VectorXd::Zero( observableSize );
                if( std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter ) != nullptr )
                {
                    absoluteValueCutOff = std::dynamic_pointer_cast< ObservationFilter< double > >( observationFilter )->getFilterValue( ) *
                            Eigen::VectorXd::Ones( observableSize );
                }
                else if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter ) != nullptr )
                {
                    if( std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( ).size( ) !=
                        static_cast< int >( observableSize ) )
                    {
                        throw std::runtime_error(
                                "Error when performing observation value filtering, size of the filter value inconsistent with observable "
                                "size." );
                    }
                    absoluteValueCutOff =
                            std::dynamic_pointer_cast< ObservationFilter< Eigen::VectorXd > >( observationFilter )->getFilterValue( );
                }

                for( unsigned int j = 0; j < numberOfObservations; ++j )
                {
                    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > singleObservation =
                            getObservationValue( observationIdsBySet_.at( setId ).at( j ) );
                    bool removeObservation = false;
                    for( unsigned int k = 0; k < observableSize; ++k )
                    {
                        if( ( !useOppositeCondition && ( singleObservation[ k ] > absoluteValueCutOff[ k ] ) ) ||
                            ( useOppositeCondition && ( singleObservation[ k ] <= absoluteValueCutOff[ k ] ) ) )
                        {
                            removeObservation = true;
                        }
                    }
                    if( removeObservation )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case epochs_filtering: {
                const std::vector< double > filterEpochs =
                        std::dynamic_pointer_cast< ObservationFilter< std::vector< double > > >( observationFilter )->getFilterValue( );
                for( unsigned int j = 0; j < numberOfObservations; ++j )
                {
                    const TimeType singleObservationTime = getObservationTime( observationIdsBySet_.at( setId ).at( j ) );
                    if( ( !useOppositeCondition &&
                          ( std::count( filterEpochs.begin( ), filterEpochs.end( ), singleObservationTime ) > 0 ) ) ||
                        ( useOppositeCondition &&
                          ( std::count( filterEpochs.begin( ), filterEpochs.end( ), singleObservationTime ) == 0 ) ) )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case time_bounds_filtering: {
                const std::pair< double, double > timeBounds =
                        std::dynamic_pointer_cast< ObservationFilter< std::pair< double, double > > >( observationFilter )
                                ->getFilterValue( );
                for( unsigned int j = 0; j < numberOfObservations; ++j )
                {
                    const TimeType singleObservationTime = getObservationTime( observationIdsBySet_.at( setId ).at( j ) );
                    if( ( !useOppositeCondition &&
                          ( ( singleObservationTime >= timeBounds.first ) && ( singleObservationTime <= timeBounds.second ) ) ) ||
                        ( useOppositeCondition &&
                          ( ( singleObservationTime < timeBounds.first ) || ( singleObservationTime > timeBounds.second ) ) ) )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            case dependent_variable_filtering: {
                if( std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter ) == nullptr )
                {
                    throw std::runtime_error(
                            "Error when performing dependent variable filtering, inconsistent filter input (should be "
                            "ObservationDependentVariableFilter object)." );
                }

                const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > settings =
                        std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter )
                                ->getDependentVariableSettings( );
                const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
                const unsigned int dependentVariableSize =
                        getObservationDependentVariableSize( settings, getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_ );

                const Eigen::VectorXd dependentVariableCutOff =
                        std::dynamic_pointer_cast< ObservationDependentVariableFilter >( observationFilter )->getFilterValue( );
                if( dependentVariableCutOff.size( ) != static_cast< int >( dependentVariableSize ) )
                {
                    throw std::runtime_error(
                            "Error when performing dependent variable filtering, size of the dependent variable cut off vector "
                            "inconsistent with dependent variable size." );
                }

                const Eigen::MatrixXd singleDependentVariableValues = getSingleDependentVariableForSet( setId, settings );
                if( ( singleDependentVariableValues.rows( ) != numberOfObservations ) ||
                    ( singleDependentVariableValues.cols( ) != dependentVariableSize ) )
                {
                    throw std::runtime_error(
                            "Error when performing dependent variable filtering, size of observation dependent variables is inconsistent "
                            "with the number of observations and presupposed dependent variable size." );
                }

                for( unsigned int j = 0; j < numberOfObservations; ++j )
                {
                    bool removeObservation = false;
                    for( unsigned int k = 0; k < dependentVariableSize; ++k )
                    {
                        if( ( !useOppositeCondition && ( singleDependentVariableValues( j, k ) ) > dependentVariableCutOff[ k ] ) ||
                            ( useOppositeCondition && ( singleDependentVariableValues( j, k ) ) <= dependentVariableCutOff[ k ] ) )
                        {
                            removeObservation = true;
                        }
                    }
                    if( removeObservation )
                    {
                        indicesToRemove.push_back( j );
                    }
                }
                break;
            }
            default:
                throw std::runtime_error( "Observation filter type not recognised." );
        }

        return indicesToRemove;
    }

    void moveObservationsToSet( const ObservationSetId sourceSetId,
                                ObservationDataset< ObservationScalarType, TimeType >& targetDataset,
                                const ObservationSetId targetSetId,
                                const std::vector< unsigned int >& indices,
                                const bool removeFromSource = true )
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
        std::vector< TimeType > times;
        std::vector< Eigen::VectorXd > dependentVariables;
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;

        const bool hasDependentVariables = !getDependentVariablesForSet( sourceSetId ).empty( );
        for( const unsigned int index : indices )
        {
            if( index >= getNumberOfObservationsForSet( sourceSetId ) )
            {
                throw std::runtime_error( "Error when moving observations in dataset, index is out of bounds." );
            }
            const ObservationId observationId = observationIdsBySet_.at( sourceSetId ).at( index );
            observations.push_back( getObservationValue( observationId ) );
            times.push_back( getObservationTime( observationId ) );
            weights.push_back( getWeightValue( observationId ) );
            residuals.push_back( getResidualValue( observationId ) );
            if( hasDependentVariables )
            {
                dependentVariables.push_back( getDependentVariables( observationId ) );
            }
        }

        targetDataset.addObservationsToSet( targetSetId, observations, times, dependentVariables, weights, residuals, true );
        if( removeFromSource )
        {
            removeObservationsFromSet( sourceSetId, indices );
        }
    }

    std::pair< TimeType, TimeType > getTimeBoundsForSet( const ObservationSetId setId ) const
    {
        const std::vector< TimeType > observationTimes = getObservationTimesForSet( setId );
        if( observationTimes.empty( ) )
        {
            return std::make_pair( TUDAT_NAN, TUDAT_NAN );
        }
        return std::make_pair( *std::min_element( observationTimes.begin( ), observationTimes.end( ) ),
                               *std::max_element( observationTimes.begin( ), observationTimes.end( ) ) );
    }

    void eraseDuplicateObservationsFromSet( const ObservationSetId setId, const bool printWarning = true )
    {
        const std::vector< TimeType > observationTimes = getObservationTimesForSet( setId );
        std::vector< unsigned int > indicesToRemove;
        for( unsigned int i = 1; i < observationTimes.size( ); ++i )
        {
            if( observationTimes.at( i ) == observationTimes.at( i - 1 ) )
            {
                const double currentObsValue = observationTimes.at( i );
                const double previousObsValue = observationTimes.at( i - 1 );
                if( std::abs( currentObsValue - previousObsValue ) <=
                    1.0e-12 * std::max( std::abs( currentObsValue ), std::abs( previousObsValue ) ) )
                {
                    indicesToRemove.push_back( i );
                }
            }
        }

        if( !indicesToRemove.empty( ) )
        {
            const std::size_t beforeCount = getNumberOfObservationsForSet( setId );
            removeObservationsFromSet( setId, indicesToRemove );
            if( printWarning )
            {
                std::cerr << "[WARNING] Detected and removed " << beforeCount - getNumberOfObservationsForSet( setId )
                          << "duplicate observations when creating instance of SingleObservationSet" << std::endl;
            }
        }
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getComputedObservationsForSet(
            const ObservationSetId setId ) const
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

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getComputedObservationVectorForSet( const ObservationSetId setId ) const
    {
        return getObservationVectorForSet( setId ) - getResidualVectorForSet( setId );
    }

    Eigen::VectorXd getRmsResidualsForSet( const ObservationSetId setId ) const
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

    Eigen::VectorXd getMeanResidualsForSet( const ObservationSetId setId ) const
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

    std::size_t getNumberOfObservationSets( ) const
    {
        return setMetadata_.size( );
    }

    std::size_t getNumberOfObservations( ) const
    {
        return observationRows_.size( );
    }

    std::size_t getTotalScalarSize( ) const
    {
        return scalarComponentRows_.size( );
    }

    const std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > >& getObservationSetMetadata( ) const
    {
        return setMetadata_;
    }

    const ObservationSetMetadata< ObservationScalarType, TimeType >& getObservationSetMetadata( const ObservationSetId setId ) const
    {
        return setMetadata_.at( setId );
    }

    const std::vector< ObservationDatasetRow< TimeType > >& getObservationRows( ) const
    {
        return observationRows_;
    }

    const ObservationDatasetRow< TimeType >& getObservationRow( const ObservationId observationId ) const
    {
        return observationRows_.at( observationId );
    }

    const std::vector< ObservationScalarComponentRow >& getScalarComponentRows( ) const
    {
        return scalarComponentRows_;
    }

    const ObservationScalarComponentRow& getScalarComponentRow( const ScalarComponentId scalarComponentId ) const
    {
        return scalarComponentRows_.at( scalarComponentId );
    }

    const std::vector< ObservationId >& getObservationIdsForSet( const ObservationSetId setId ) const
    {
        return observationIdsBySet_.at( setId );
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            observations.push_back( getObservationValue( observationId ) );
        }
        return observations;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVectorForSet( const ObservationSetId setId ) const
    {
        return createSetVector( setId, observedValues_ );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationValue( const ObservationId observationId ) const
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

    std::vector< TimeType > getObservationTimesForSet( const ObservationSetId setId ) const
    {
        std::vector< TimeType > times;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            times.push_back( observationRows_.at( observationId ).time_ );
        }
        return times;
    }

    TimeType getObservationTime( const ObservationId observationId ) const
    {
        return observationRows_.at( observationId ).time_;
    }

    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeightsForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            weights.push_back( getWeightValue( observationId ) );
        }
        return weights;
    }

    Eigen::VectorXd getWeightVectorForSet( const ObservationSetId setId ) const
    {
        Eigen::VectorXd weights = Eigen::VectorXd::Zero( getTotalScalarSizeForSet( setId ) );
        std::size_t currentIndex = 0;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            for( unsigned int i = 0; i < row.scalarSize_; ++i )
            {
                weights( currentIndex++ ) = legacyWeights_.at( row.firstScalarComponent_ + i );
            }
        }
        return weights;
    }

    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeightValue( const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        Eigen::Matrix< double, Eigen::Dynamic, 1 > value = Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( row.scalarSize_ );
        for( unsigned int i = 0; i < row.scalarSize_; ++i )
        {
            value( i ) = legacyWeights_.at( row.firstScalarComponent_ + i );
        }
        return value;
    }

    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResidualsForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            residuals.push_back( getResidualValue( observationId ) );
        }
        return residuals;
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualVectorForSet( const ObservationSetId setId ) const
    {
        return createSetVector( setId, residualValues_ );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualValue( const ObservationId observationId ) const
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

    std::vector< Eigen::VectorXd > getDependentVariablesForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::VectorXd > dependentVariables;
        bool hasNonEmptyDependentVariables = false;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            const Eigen::VectorXd dependentVariable = dependentVariableValues_.at( observationId );
            if( dependentVariable.size( ) > 0 )
            {
                hasNonEmptyDependentVariables = true;
            }
            dependentVariables.push_back( dependentVariable );
        }
        return hasNonEmptyDependentVariables ? dependentVariables : std::vector< Eigen::VectorXd >( );
    }

    Eigen::VectorXd getDependentVariables( const ObservationId observationId ) const
    {
        return dependentVariableValues_.at( observationId );
    }

    Eigen::MatrixXd getSingleDependentVariableForSet( const ObservationSetId setId,
                                                      const std::pair< int, int >& dependentVariableIndexAndSize ) const
    {
        const std::vector< Eigen::VectorXd > observationsDependentVariables = getDependentVariablesForSet( setId );
        Eigen::MatrixXd singleDependentVariable =
                Eigen::MatrixXd::Zero( getNumberOfObservationsForSet( setId ), dependentVariableIndexAndSize.second );
        for( unsigned int i = 0; i < observationsDependentVariables.size( ); ++i )
        {
            if( dependentVariableIndexAndSize.first + dependentVariableIndexAndSize.second >
                observationsDependentVariables.at( i ).size( ) )
            {
                throw std::runtime_error(
                        "Error when retrieving single observation dependent variable, required index and size incompatible with dependent "
                        "variables size." );
            }
            Eigen::VectorXd singleDependentVariableVector = observationsDependentVariables.at( i ).segment(
                    dependentVariableIndexAndSize.first, dependentVariableIndexAndSize.second );
            singleDependentVariable.block( i, 0, 1, dependentVariableIndexAndSize.second ) = singleDependentVariableVector.transpose( );
        }
        return singleDependentVariable;
    }

    Eigen::MatrixXd getSingleDependentVariableForSet(
            const ObservationSetId setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings,
            const bool returnFirstCompatibleSettings = false ) const
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
                getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
        if( dependentVariableBookkeeping == nullptr )
        {
            throw std::runtime_error(
                    "Error when getting dependent variable from observation dataset, no dependent variable bookkeeping is available." );
        }

        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                settingsIndicesAndSizes = dependentVariableBookkeeping->getSettingsIndicesAndSizes( );

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

    std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > > getCompatibleDependentVariableSettingsForSet(
            const ObservationSetId setId,
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

    std::vector< Eigen::MatrixXd > getAllCompatibleDependentVariablesForSet(
            const ObservationSetId setId,
            const std::shared_ptr< simulation_setup::ObservationDependentVariableSettings >& dependentVariableSettings ) const
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping =
                getDependentVariableBookkeeping( metadata.dependentVariableLayoutId_ );
        if( dependentVariableBookkeeping == nullptr )
        {
            return std::vector< Eigen::MatrixXd >( );
        }

        std::map< std::pair< int, int >, std::shared_ptr< simulation_setup::ObservationDependentVariableSettings > >
                settingsIndicesAndSizes = dependentVariableBookkeeping->getSettingsIndicesAndSizes( );

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

    void setDependentVariablesForSet( const ObservationSetId setId, const std::vector< Eigen::VectorXd >& dependentVariables )
    {
        const std::vector< ObservationId >& observationIds = observationIdsBySet_.at( setId );
        if( dependentVariables.size( ) != observationIds.size( ) )
        {
            throw std::runtime_error( "Error when setting dataset dependent variables, size is inconsistent." );
        }
        for( std::size_t i = 0; i < observationIds.size( ); ++i )
        {
            dependentVariableValues_.at( observationIds.at( i ) ) = dependentVariables.at( i );
        }
    }

    void clearDependentVariablesForSet( const ObservationSetId setId )
    {
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            dependentVariableValues_.at( observationId ) = Eigen::VectorXd( );
        }
    }

    std::size_t getNumberOfObservationsForSet( const ObservationSetId setId ) const
    {
        return observationIdsBySet_.at( setId ).size( );
    }

    std::size_t getTotalScalarSizeForSet( const ObservationSetId setId ) const
    {
        return observationIdsBySet_.at( setId ).size( ) * getObservationSetMetadata( setId ).observableSize_;
    }

    const LinkDefinition& getLinkDefinition( const LinkDefinitionId linkDefinitionId ) const
    {
        return linkDefinitionRegistry_.at( linkDefinitionId );
    }

    std::size_t getNumberOfLinkDefinitions( ) const
    {
        return linkDefinitionRegistry_.size( );
    }

    const std::shared_ptr< ObservationAncillarySimulationSettings >& getAncillarySettings(
            const AncillarySettingsId ancillarySettingsId ) const
    {
        return ancillarySettingsRegistry_.at( ancillarySettingsId );
    }

    const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& getDependentVariableBookkeeping(
            const DependentVariableLayoutId dependentVariableLayoutId ) const
    {
        return dependentVariableLayoutRegistry_.at( dependentVariableLayoutId );
    }

    std::vector< ObservationSetId > getObservationSetIds( const std::shared_ptr< ObservationCollectionParser >& observationParser =
                                                                  std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< ObservationSetId > setIds;
        for( ObservationSetId setId = 0; setId < setMetadata_.size( ); ++setId )
        {
            if( isObservationSetSelectedByParser( setId, observationParser ) )
            {
                setIds.push_back( setId );
            }
        }
        return setIds;
    }

    EstimationVectorProjection< ObservationScalarType, TimeType > createLegacyProjection( const bool includeInactive = true ) const
    {
        EstimationVectorProjection< ObservationScalarType, TimeType > projection;

        std::size_t projectionSize = 0;
        for( const ObservationDatasetRow< TimeType >& row : observationRows_ )
        {
            if( includeInactive || row.isActive_ )
            {
                projectionSize += row.scalarSize_;
            }
        }

        projection.observations_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( projectionSize );
        projection.residuals_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( projectionSize );
        projection.weights_ = Eigen::VectorXd::Zero( projectionSize );
        projection.times_.reserve( projectionSize );
        projection.observationIds_.reserve( projectionSize );
        projection.setIds_.reserve( projectionSize );
        projection.scalarComponentIds_.reserve( projectionSize );

        std::size_t currentIndex = 0;
        for( ObservationId observationId = 0; observationId < observationRows_.size( ); ++observationId )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            if( includeInactive || row.isActive_ )
            {
                for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
                {
                    const ScalarComponentId scalarComponentId = row.firstScalarComponent_ + componentIndex;
                    projection.observations_( currentIndex ) = observedValues_.at( scalarComponentId );
                    projection.residuals_( currentIndex ) = residualValues_.at( scalarComponentId );
                    projection.weights_( currentIndex ) = legacyWeights_.at( scalarComponentId );
                    projection.times_.push_back( row.time_ );
                    projection.observationIds_.push_back( observationId );
                    projection.setIds_.push_back( row.setId_ );
                    projection.scalarComponentIds_.push_back( scalarComponentId );
                    ++currentIndex;
                }
            }
        }

        return projection;
    }

private:
    bool isObservationSetSelectedByParser( const ObservationSetId setId,
                                           const std::shared_ptr< ObservationCollectionParser >& observationParser ) const
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = getObservationSetMetadata( setId );
        const LinkEnds linkEnds = getLinkDefinition( metadata.linkDefinitionId_ ).linkEnds_;
        const bool useOppositeCondition = observationParser->useOppositeCondition( );

        bool isSelected = false;
        switch( observationParser->getObservationParserType( ) )
        {
            case empty_parser: {
                isSelected = true;
                break;
            }
            case observable_type_parser: {
                const std::vector< ObservableType > observableTypes =
                        std::dynamic_pointer_cast< ObservationCollectionObservableTypeParser >( observationParser )->getObservableTypes( );
                isSelected = std::count( observableTypes.begin( ), observableTypes.end( ), metadata.observableType_ ) > 0;
                break;
            }
            case link_ends_parser: {
                const std::vector< LinkEnds > linkEndsVector =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndsParser >( observationParser )->getLinkEndsVector( );
                isSelected = std::count( linkEndsVector.begin( ), linkEndsVector.end( ), linkEnds ) > 0;
                break;
            }
            case link_end_string_parser: {
                const std::shared_ptr< ObservationCollectionLinkEndStringParser > linkEndStringObservationParser =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndStringParser >( observationParser );
                const std::vector< std::string > linkEndNames = linkEndStringObservationParser->getLinkEndNames( );
                const bool isReferencePoint = linkEndStringObservationParser->isReferencePoint( );

                for( const auto& linkEnd : linkEnds )
                {
                    const std::string name = isReferencePoint ? linkEnd.second.getReferencePointName( ) : linkEnd.second.bodyName_;
                    if( std::count( linkEndNames.begin( ), linkEndNames.end( ), name ) > 0 )
                    {
                        isSelected = true;
                    }
                }
                break;
            }
            case link_end_id_parser: {
                const std::vector< LinkEndId > linkEndIds =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndIdParser >( observationParser )->getLinkEndIds( );
                for( const auto& linkEnd : linkEnds )
                {
                    if( std::count( linkEndIds.begin( ), linkEndIds.end( ), linkEnd.second ) > 0 )
                    {
                        isSelected = true;
                    }
                }
                break;
            }
            case link_end_type_parser: {
                const std::vector< LinkEndType > linkEndTypes =
                        std::dynamic_pointer_cast< ObservationCollectionLinkEndTypeParser >( observationParser )->getLinkEndTypes( );
                for( const auto& linkEnd : linkEnds )
                {
                    if( std::count( linkEndTypes.begin( ), linkEndTypes.end( ), linkEnd.first ) > 0 )
                    {
                        isSelected = true;
                    }
                }
                break;
            }
            case single_link_end_parser: {
                const std::vector< std::pair< LinkEndType, LinkEndId > > singleLinkEnds =
                        std::dynamic_pointer_cast< ObservationCollectionSingleLinkEndParser >( observationParser )->getSingleLinkEnds( );
                for( const auto& singleLinkEnd : singleLinkEnds )
                {
                    if( linkEnds.count( singleLinkEnd.first ) > 0 && linkEnds.at( singleLinkEnd.first ) == singleLinkEnd.second )
                    {
                        isSelected = true;
                    }
                }
                break;
            }
            case time_bounds_parser: {
                const std::vector< std::pair< double, double > > timeBoundsVector =
                        std::dynamic_pointer_cast< ObservationCollectionTimeBoundsParser >( observationParser )->getTimeBoundsVector( );
                const std::pair< TimeType, TimeType > setTimeBounds = getTimeBoundsForSet( setId );
                for( const auto& timeBounds : timeBoundsVector )
                {
                    if( ( setTimeBounds.first >= timeBounds.first ) && ( setTimeBounds.second <= timeBounds.second ) )
                    {
                        isSelected = true;
                    }
                }
                break;
            }
            case ancillary_settings_parser: {
                const std::vector< std::shared_ptr< ObservationAncillarySimulationSettings > > ancillarySettings =
                        std::dynamic_pointer_cast< ObservationCollectionAncillarySettingsParser >( observationParser )
                                ->getAncillarySettings( );
                const std::shared_ptr< ObservationAncillarySimulationSettings >& currentSettings =
                        getAncillarySettings( metadata.ancillarySettingsId_ );
                for( const auto& settings : ancillarySettings )
                {
                    if( currentSettings == settings )
                    {
                        isSelected = true;
                    }
                    else if( currentSettings != nullptr && settings != nullptr &&
                             currentSettings->getDoubleData( ) == settings->getDoubleData( ) &&
                             currentSettings->getDoubleVectorData( ) == settings->getDoubleVectorData( ) )
                    {
                        isSelected = true;
                    }
                }
                break;
            }
            case multi_type_parser: {
                const std::shared_ptr< ObservationCollectionMultiTypeParser > multiTypeParser =
                        std::dynamic_pointer_cast< ObservationCollectionMultiTypeParser >( observationParser );
                const std::vector< std::shared_ptr< ObservationCollectionParser > > observationParsers =
                        multiTypeParser->getObservationParsers_( );

                if( multiTypeParser->areConditionsCombined( ) )
                {
                    isSelected = true;
                    for( const std::shared_ptr< ObservationCollectionParser >& parser : observationParsers )
                    {
                        if( !isObservationSetSelectedByParser( setId, parser ) )
                        {
                            isSelected = false;
                        }
                    }
                }
                else
                {
                    for( const std::shared_ptr< ObservationCollectionParser >& parser : observationParsers )
                    {
                        if( isObservationSetSelectedByParser( setId, parser ) )
                        {
                            isSelected = true;
                        }
                    }
                }
                return isSelected;
            }
            default:
                throw std::runtime_error( "Observation parser type not recognised." );
        }

        return useOppositeCondition ? !isSelected : isSelected;
    }

    void validateObservationSetData( const ObservationSetId setId,
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

    static std::vector< std::size_t > getTimeSortingPermutation( const std::vector< TimeType >& observationTimes )
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

    template< typename VectorType >
    static void reorderVector( std::vector< VectorType >& data, const std::vector< std::size_t >& permutation )
    {
        std::vector< VectorType > reorderedData( data.size( ) );
        for( std::size_t i = 0; i < data.size( ); ++i )
        {
            reorderedData.at( i ) = data.at( permutation.at( i ) );
        }
        data.swap( reorderedData );
    }

    template< typename VectorType >
    static void removeEntries( std::vector< VectorType >& data, const std::vector< unsigned int >& indicesToRemove )
    {
        for( std::vector< unsigned int >::const_reverse_iterator indexIterator = indicesToRemove.rbegin( );
             indexIterator != indicesToRemove.rend( );
             ++indexIterator )
        {
            data.erase( data.begin( ) + *indexIterator );
        }
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > createSetVector(
            const ObservationSetId setId,
            const std::vector< ObservationScalarType >& scalarValues ) const
    {
        Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > vector =
                Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( getTotalScalarSizeForSet( setId ) );
        std::size_t currentIndex = 0;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            for( unsigned int i = 0; i < row.scalarSize_; ++i )
            {
                vector( currentIndex++ ) = scalarValues.at( row.firstScalarComponent_ + i );
            }
        }
        return vector;
    }

    void setObservationValue( const ObservationId observationId,
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

    void setResidualValue( const ObservationId observationId, const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residual )
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

    void setWeightValue( const ObservationId observationId, const Eigen::VectorXd& weight )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        if( weight.size( ) != static_cast< int >( row.scalarSize_ ) )
        {
            throw std::runtime_error( "Error when setting dataset weight, scalar size is inconsistent." );
        }
        for( unsigned int i = 0; i < row.scalarSize_; ++i )
        {
            legacyWeights_.at( row.firstScalarComponent_ + i ) = weight( i );
        }
    }

    LinkDefinitionId registerLinkDefinition( const LinkDefinition& linkDefinition )
    {
        for( LinkDefinitionId i = 0; i < linkDefinitionRegistry_.size( ); ++i )
        {
            if( linkDefinitionRegistry_.at( i ) == linkDefinition )
            {
                return i;
            }
        }

        linkDefinitionRegistry_.push_back( linkDefinition );
        return linkDefinitionRegistry_.size( ) - 1;
    }

    AncillarySettingsId registerAncillarySettings( const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings )
    {
        for( AncillarySettingsId i = 0; i < ancillarySettingsRegistry_.size( ); ++i )
        {
            if( ancillarySettingsRegistry_.at( i ) == ancillarySettings )
            {
                return i;
            }
        }

        ancillarySettingsRegistry_.push_back( ancillarySettings );
        return ancillarySettingsRegistry_.size( ) - 1;
    }

    DependentVariableLayoutId registerDependentVariableLayout(
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& bookkeeping )
    {
        for( DependentVariableLayoutId i = 0; i < dependentVariableLayoutRegistry_.size( ); ++i )
        {
            if( dependentVariableLayoutRegistry_.at( i ) == bookkeeping )
            {
                return i;
            }
        }

        dependentVariableLayoutRegistry_.push_back( bookkeeping );
        return dependentVariableLayoutRegistry_.size( ) - 1;
    }

    std::vector< ObservationDatasetRow< TimeType > > observationRows_;
    std::vector< ObservationScalarComponentRow > scalarComponentRows_;
    std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > > setMetadata_;
    std::vector< std::vector< ObservationId > > observationIdsBySet_;

    std::vector< LinkDefinition > linkDefinitionRegistry_;
    std::vector< std::shared_ptr< ObservationAncillarySimulationSettings > > ancillarySettingsRegistry_;
    std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > > dependentVariableLayoutRegistry_;

    std::vector< ObservationScalarType > observedValues_;
    std::vector< ObservationScalarType > residualValues_;
    std::vector< double > legacyWeights_;
    std::vector< Eigen::VectorXd > dependentVariableValues_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_DATASET_H
