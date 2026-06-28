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
#include <Eigen/SparseCore>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/simulation/estimation_setup/estimationVectorProjection.h"
#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationOutput.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"
#include "tudat/simulation/estimation_setup/observationWeights.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"
#include "tudat/basics/tudatTypeTraits.h"

namespace tudat
{

namespace observation_models
{
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
bool isObservationSetSelectedByLegacyParser( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                                             const ObservationSetId setId,
                                             const std::shared_ptr< ObservationCollectionParser >& observationParser );

//! Queryable backend for Tudat observation data.
/*!
 * ObservationDataset stores observations as event rows plus scalar-component
 * arrays. A vector observable has one observation row and N scalar components:
 * observedValues_ is flat by scalar component, while observationRows_ records
 * the event boundary and scalarComponentRows_ provides the reverse mapping.
 * Set-level metadata is stored once in registries and referenced by id. Flat
 * estimator vectors are derived by createLegacyProjection(), not used as the
 * primary data model. ObservationDatasetViewer instances are invalidated by
 * structural mutations that add, remove or rebuild observation rows.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationDataset
{
public:
    ObservationDataset( ) = default;

    //! Add one logical observation set and register its metadata.
    /*!
     * Inputs are supplied as per-observation vectors. The method validates
     * observable dimensions, optionally sorts by time and removes duplicates,
     * registers shared metadata, then appends event rows and scalar components.
     */
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
            }

            if( preparedWeights.empty( ) )
            {
                observationWeights_.appendScalarWeight( 1.0 );
            }
            else
            {
                observationWeights_.appendDiagonalWeightVector( preparedWeights.at( i ) );
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

        ++structuralVersion_;
        return setId;
    }

    //! Copy one set from another dataset into this dataset.
    ObservationSetId addObservationSetFromDataset( const ObservationDataset< ObservationScalarType, TimeType >& sourceDataset,
                                                   const ObservationSetId sourceSetId )
    {
        const ObservationSetMetadata< ObservationScalarType, TimeType >& sourceMetadata =
                sourceDataset.getObservationSetMetadata( sourceSetId );
        const ObservationSetId newSetId =
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
            const std::vector< ObservationId >& sourceObservationIds = sourceDataset.getObservationIdsForSet( sourceSetId );
            const std::vector< ObservationId >& targetObservationIds = getObservationIdsForSet( newSetId );
            for( std::size_t i = 0; i < sourceObservationIds.size( ); ++i )
            {
                setWeightMatrixForObservation( targetObservationIds.at( i ),
                                               sourceDataset.getWeightMatrixForObservation( sourceObservationIds.at( i ) ) );
            }
        }
        return newSetId;
    }

    //! Add a set using one compact scalar weight for every observation event.
    ObservationSetId addObservationSetWithScalarWeight(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const double weight,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ) )
    {
        const ObservationSetId setId = addObservationSet( observableType,
                                                          linkDefinition,
                                                          observations,
                                                          times,
                                                          referenceLinkEnd,
                                                          dependentVariables,
                                                          dependentVariableBookkeeping,
                                                          ancillarySettings,
                                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                          residuals );
        setConstantWeightForSet( setId, weight );
        return setId;
    }

    //! Add a set using one compact scalar weight per observation event.
    ObservationSetId addObservationSetWithScalarWeights(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const std::vector< double >& weights,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ) )
    {
        if( weights.size( ) != observations.size( ) )
        {
            throw std::runtime_error( "Error when adding observation set with scalar weights, weight count is inconsistent." );
        }
        const ObservationSetId setId = addObservationSet( observableType,
                                                          linkDefinition,
                                                          observations,
                                                          times,
                                                          referenceLinkEnd,
                                                          dependentVariables,
                                                          dependentVariableBookkeeping,
                                                          ancillarySettings,
                                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                          residuals );
        const std::vector< ObservationId >& observationIds = getObservationIdsForSet( setId );
        for( std::size_t i = 0; i < weights.size( ); ++i )
        {
            observationWeights_.setScalarWeight( observationIds.at( i ), weights.at( i ) );
        }
        return setId;
    }

    //! Add a set using one observable-size N x N weight block for every observation event.
    ObservationSetId addObservationSetWithWeightBlock(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const Eigen::MatrixXd& weightBlock,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ) )
    {
        const ObservationSetId setId = addObservationSet( observableType,
                                                          linkDefinition,
                                                          observations,
                                                          times,
                                                          referenceLinkEnd,
                                                          dependentVariables,
                                                          dependentVariableBookkeeping,
                                                          ancillarySettings,
                                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                          residuals );
        for( const ObservationId observationId : getObservationIdsForSet( setId ) )
        {
            setWeightMatrixForObservation( observationId, weightBlock );
        }
        return setId;
    }

    //! Add a set using one observable-size N x N weight block per observation event.
    ObservationSetId addObservationSetWithWeightBlocks(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const std::vector< Eigen::MatrixXd >& weightBlocks,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ) )
    {
        if( weightBlocks.size( ) != observations.size( ) )
        {
            throw std::runtime_error( "Error when adding observation set with weight blocks, weight count is inconsistent." );
        }
        const ObservationSetId setId = addObservationSet( observableType,
                                                          linkDefinition,
                                                          observations,
                                                          times,
                                                          referenceLinkEnd,
                                                          dependentVariables,
                                                          dependentVariableBookkeeping,
                                                          ancillarySettings,
                                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                          residuals );
        const std::vector< ObservationId >& observationIds = getObservationIdsForSet( setId );
        for( std::size_t i = 0; i < weightBlocks.size( ); ++i )
        {
            setWeightMatrixForObservation( observationIds.at( i ), weightBlocks.at( i ) );
        }
        return setId;
    }

    //! Add a set using one full M x M set-level weight block for the new batch.
    ObservationSetId addObservationSetWithSetWeightBlock(
            const ObservableType observableType,
            const LinkDefinition& linkDefinition,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType >& times,
            const LinkEndType referenceLinkEnd,
            const Eigen::MatrixXd& setWeightBlock,
            const std::vector< Eigen::VectorXd >& dependentVariables = std::vector< Eigen::VectorXd >( ),
            const std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping >& dependentVariableBookkeeping = nullptr,
            const std::shared_ptr< ObservationAncillarySimulationSettings >& ancillarySettings = nullptr,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& residuals =
                    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ) )
    {
        const ObservationSetId setId = addObservationSet( observableType,
                                                          linkDefinition,
                                                          observations,
                                                          times,
                                                          referenceLinkEnd,
                                                          dependentVariables,
                                                          dependentVariableBookkeeping,
                                                          ancillarySettings,
                                                          std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >( ),
                                                          residuals );
        setWeightMatrixForSet( setId, setWeightBlock );
        return setId;
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

    //! Replace the vector-valued measurements for all observation events in one set.
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

    //! Replace a set's measurements from a legacy flat scalar-component vector.
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

    //! Replace the vector-valued residuals for all observation events in one set.
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

    //! Replace a set's residuals from a legacy flat scalar-component vector.
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

    //! Apply one compact scalar weight to every observation event in one set.
    void setConstantWeightForSet( const ObservationSetId setId, const double weight )
    {
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            observationWeights_.setScalarWeight( observationId, weight );
        }
    }

    //! Apply one component-wise weight vector to every observation event in one set.
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

    //! Replace a set's weights from a legacy flat scalar-component vector.
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

    //! Store one full M x M weight block for an observation set.
    void setWeightMatrixForSet( const ObservationSetId setId, const Eigen::MatrixXd& weightMatrix )
    {
        if( weightMatrix.rows( ) != static_cast< int >( getTotalScalarSizeForSet( setId ) ) ||
            weightMatrix.cols( ) != static_cast< int >( getTotalScalarSizeForSet( setId ) ) )
        {
            throw std::runtime_error( "Error when setting dataset set weight matrix, matrix size is inconsistent." );
        }
        observationWeights_.setSetWeightBlock( setId, weightMatrix );
    }

    //! Return whether a full set-level M x M weight block is stored for a set.
    bool hasWeightMatrixForSet( const ObservationSetId setId ) const
    {
        return observationWeights_.hasSetWeightBlock( setId );
    }

    //! Store one observable-size N x N weight block for an observation event.
    void setWeightMatrixForObservation( const ObservationId observationId, const Eigen::MatrixXd& weightMatrix )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        if( weightMatrix.rows( ) != static_cast< int >( row.scalarSize_ ) || weightMatrix.cols( ) != static_cast< int >( row.scalarSize_ ) )
        {
            throw std::runtime_error( "Error when setting dataset observation weight matrix, matrix size is inconsistent." );
        }
        observationWeights_.setWeightBlock( observationId, weightMatrix );
    }

    //! Return whether an observation event stores an explicit N x N weight block.
    bool hasWeightMatrixForObservation( const ObservationId observationId ) const
    {
        return observationWeights_.hasObservationWeightBlock( observationId );
    }

    //! Add an advanced off-diagonal weight block over selected scalar components.
    void addExtraWeightBlock( const ObservationWeightBlock& weightBlock )
    {
        observationWeights_.addExtraWeightBlock( weightBlock );
    }

    //! Store a dense weight block selected by observation ids.
    /*!
     * This is the public interface for correlations between unrelated
     * observations. Empty component lists select all scalar components of each
     * observation. Non-empty component lists are applied to every observation in
     * the corresponding row or column selection.
     */
    void setWeightBlock( const std::vector< ObservationId >& rowObservationIds,
                         const std::vector< ObservationId >& columnObservationIds,
                         const Eigen::MatrixXd& weightBlock,
                         const std::vector< unsigned int >& rowComponents = std::vector< unsigned int >( ),
                         const std::vector< unsigned int >& columnComponents = std::vector< unsigned int >( ),
                         const bool makeSymmetric = false )
    {
        const std::vector< ScalarComponentId > rowScalarComponentIds =
                getScalarComponentIdsForObservationSelection( rowObservationIds, rowComponents );
        const std::vector< ScalarComponentId > columnScalarComponentIds =
                getScalarComponentIdsForObservationSelection( columnObservationIds, columnComponents );

        if( weightBlock.rows( ) != static_cast< int >( rowScalarComponentIds.size( ) ) ||
            weightBlock.cols( ) != static_cast< int >( columnScalarComponentIds.size( ) ) )
        {
            throw std::runtime_error( "Error when setting dataset weight block, matrix size is inconsistent with selected observations." );
        }

        ObservationWeightBlock datasetWeightBlock;
        datasetWeightBlock.rowScalarComponentIds_ = rowScalarComponentIds;
        datasetWeightBlock.columnScalarComponentIds_ = columnScalarComponentIds;
        datasetWeightBlock.weightBlock_ = weightBlock;
        addExtraWeightBlock( datasetWeightBlock );

        if( makeSymmetric )
        {
            if( rowScalarComponentIds == columnScalarComponentIds )
            {
                if( weightBlock.rows( ) != weightBlock.cols( ) || !weightBlock.isApprox( weightBlock.transpose( ) ) )
                {
                    throw std::runtime_error(
                            "Error when setting symmetric dataset weight block, block with identical row and column selection is not "
                            "symmetric." );
                }
            }
            else
            {
                ObservationWeightBlock transposedWeightBlock;
                transposedWeightBlock.rowScalarComponentIds_ = columnScalarComponentIds;
                transposedWeightBlock.columnScalarComponentIds_ = rowScalarComponentIds;
                transposedWeightBlock.weightBlock_ = weightBlock.transpose( );
                addExtraWeightBlock( transposedWeightBlock );
            }
        }
    }

    //! Return the advanced scalar-component weight blocks stored on the dataset.
    const std::vector< ObservationWeightBlock >& getExtraWeightBlocks( ) const
    {
        return observationWeights_.getExtraWeightBlocks( );
    }

    //! Return whether advanced scalar-component weight blocks are stored on the dataset.
    bool hasExtraWeightBlocks( ) const
    {
        return observationWeights_.hasExtraWeightBlocks( );
    }

    //! Set a scalar constant weight on all observations selected by a parser.
    void setConstantWeight(
            const double weight = 1.0,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< ObservationSetId > setIds = getObservationSetIds( observationParser );
        if( setIds.empty( ) )
        {
            std::cerr << "Warning when setting constant weights, no observation dataset set found for specified observation parser. "
                         "Weights not set";
        }
        for( const ObservationSetId setId : setIds )
        {
            setConstantWeightForSet( setId, weight );
        }
    }

    //! Set a vector-valued constant weight on all observations selected by a parser.
    void setConstantWeight(
            const Eigen::VectorXd weight,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< ObservationSetId > setIds = getObservationSetIds( observationParser );
        if( setIds.empty( ) )
        {
            std::cerr << "Warning when setting constant weights, no observation dataset set found for specified observation parser. "
                         "Weights not set";
        }
        for( const ObservationSetId setId : setIds )
        {
            setConstantWeightForSet( setId, weight );
        }
    }

    //! Apply scalar constant weights for multiple parser-defined selections.
    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, double > weightsPerObservationParser )
    {
        for( auto parserIt : weightsPerObservationParser )
        {
            setConstantWeight( parserIt.second, parserIt.first );
        }
    }

    //! Apply vector constant weights for multiple parser-defined selections.
    void setConstantWeightPerObservable(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
    {
        for( auto parserIt : weightsPerObservationParser )
        {
            setConstantWeight( parserIt.second, parserIt.first );
        }
    }

    //! Set scalar-component weights from a flat vector for parser-selected sets.
    void setTabulatedWeights(
            const Eigen::VectorXd tabulatedWeights,
            const std::shared_ptr< ObservationCollectionParser > observationParser = std::make_shared< ObservationCollectionParser >( ) )
    {
        std::vector< ObservationSetId > setIds = getObservationSetIds( observationParser );
        if( setIds.empty( ) )
        {
            std::cerr << "Warning when setting tabulated weights, no observation dataset set found for specified observation parser. "
                         "Weights not set";
            return;
        }

        bool areSetsSameSize = true;
        int totalSizeAllSets = static_cast< int >( getTotalScalarSizeForSet( setIds.at( 0 ) ) );
        for( unsigned int i = 1; i < setIds.size( ); ++i )
        {
            const int currentSetSize = static_cast< int >( getTotalScalarSizeForSet( setIds.at( i ) ) );
            totalSizeAllSets += currentSetSize;
            if( currentSetSize != static_cast< int >( getTotalScalarSizeForSet( setIds.at( 0 ) ) ) )
            {
                areSetsSameSize = false;
            }
        }

        int startSet = 0;
        for( const ObservationSetId setId : setIds )
        {
            const int currentSetSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
            if( tabulatedWeights.size( ) == totalSizeAllSets )
            {
                setWeightVectorForSet( setId, tabulatedWeights.segment( startSet, currentSetSize ) );
                startSet += currentSetSize;
            }
            else if( areSetsSameSize && tabulatedWeights.size( ) == static_cast< int >( getTotalScalarSizeForSet( setIds.at( 0 ) ) ) )
            {
                setWeightVectorForSet( setId, tabulatedWeights );
            }
            else
            {
                throw std::runtime_error(
                        "Error when setting tabulated weights, the size of the input weight vector should be consistent with either the "
                        "size of each individual observation set, or the combined size of all selected observation sets." );
            }
        }
    }

    //! Apply tabulated scalar-component weights for multiple parser-defined selections.
    void setTabulatedWeights(
            const std::map< std::shared_ptr< ObservationCollectionParser >, Eigen::VectorXd > weightsPerObservationParser )
    {
        for( auto parserIt : weightsPerObservationParser )
        {
            setTabulatedWeights( parserIt.second, parserIt.first );
        }
    }

    //! Replace all rows and scalar components of one set while preserving its set id.
    /*!
     * This operation rebuilds the dataset contents in-place so existing shared
     * dataset pointers remain valid. Other set ids are preserved in order.
     */
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
                const ObservationSetId newSetId =
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
                if( hasWeightMatrixForSet( currentSetId ) )
                {
                    rebuiltDataset.setWeightMatrixForSet( newSetId, getWeightMatrixForSet( currentSetId ) );
                }
            }
        }

        const std::size_t newStructuralVersion = structuralVersion_ + 1;
        *this = rebuiltDataset;
        structuralVersion_ = newStructuralVersion;
    }

    //! Append per-observation data to an existing set.
    /*!
     * Missing weights and residuals are filled with unit weights and zero
     * residuals. If requested, the complete set is sorted by time after append.
     */
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

    //! Remove observation events from a set by per-set observation index.
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

    //! Evaluate an ObservationFilter and return per-set observation indices to remove.
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

    //! Copy or move selected observations between sets, preserving scalar data.
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

    //! Remove duplicate epochs from one set using the legacy duplicate criterion.
    void eraseDuplicateObservationsFromSet( const ObservationSetId setId, const bool printWarning = true )
    {
        const std::vector< TimeType > observationTimes = getObservationTimesForSet( setId );
        std::vector< unsigned int > indicesToRemove;
        for( unsigned int i = 1; i < observationTimes.size( ); ++i )
        {
            if( observationTimes.at( i ) == observationTimes.at( i - 1 ) )
            {
                indicesToRemove.push_back( i );
            }
        }

        if( !indicesToRemove.empty( ) )
        {
            const std::size_t beforeCount = getNumberOfObservationsForSet( setId );
            removeObservationsFromSet( setId, indicesToRemove );
            if( printWarning )
            {
                std::cerr << "[WARNING] Detected and removed " << beforeCount - getNumberOfObservationsForSet( setId )
                          << "duplicate observations when creating observation dataset" << std::endl;
            }
        }
    }

    //! Return computed observations as observed-minus-residual values for one set.
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

    //! Return computed observations as one flat scalar vector for one set.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getComputedObservationVectorForSet( const ObservationSetId setId ) const
    {
        return getObservationVectorForSet( setId ) - getResidualVectorForSet( setId );
    }

    //! Compute component-wise RMS residuals over all observations in one set.
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

    //! Compute component-wise mean residuals over all observations in one set.
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

    //! Return one vector-valued measurement per observation event in a set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getObservationsForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > observations;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            observations.push_back( getObservationValue( observationId ) );
        }
        return observations;
    }

    //! Return all scalar components of a set as one legacy flat vector.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationVectorForSet( const ObservationSetId setId ) const
    {
        return createSetVector( setId, observedValues_ );
    }

    //! Reconstruct the vector-valued measurement for one observation event.
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

    //! Return one reference-link-end time per observation event in a set.
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

    //! Return one vector of scalar-component weights per observation event.
    std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > getWeightsForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > > weights;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            weights.push_back( getWeightValue( observationId ) );
        }
        return weights;
    }

    //! Return all scalar-component weights of a set as one legacy flat vector.
    Eigen::VectorXd getWeightVectorForSet( const ObservationSetId setId ) const
    {
        Eigen::VectorXd weights = Eigen::VectorXd::Zero( getTotalScalarSizeForSet( setId ) );
        std::size_t currentIndex = 0;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            weights.segment( currentIndex, row.scalarSize_ ) =
                    observationWeights_.getObservationWeightVector( observationId, row.scalarSize_ );
            currentIndex += row.scalarSize_;
        }
        if( observationWeights_.hasSetWeightBlock( setId ) )
        {
            weights = observationWeights_.getSetWeightBlock( setId ).diagonal( );
        }
        return weights;
    }

    //! Reconstruct the scalar-component weight vector for one observation event.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > getWeightValue( const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        return observationWeights_.getObservationWeightVector( observationId, row.scalarSize_ );
    }

    //! Reconstruct the weight matrix for one observation event.
    Eigen::MatrixXd getWeightMatrixForObservation( const ObservationId observationId ) const
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        return observationWeights_.getObservationWeightMatrix( observationId, row.scalarSize_ );
    }

    //! Return a set's full weight matrix, materializing compact weights if needed.
    Eigen::MatrixXd getWeightMatrixForSet( const ObservationSetId setId ) const
    {
        if( observationWeights_.hasSetWeightBlock( setId ) )
        {
            return observationWeights_.getSetWeightBlock( setId );
        }

        Eigen::MatrixXd weightMatrix = Eigen::MatrixXd::Zero( getTotalScalarSizeForSet( setId ), getTotalScalarSizeForSet( setId ) );
        std::size_t currentIndex = 0;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            weightMatrix.block( currentIndex, currentIndex, row.scalarSize_, row.scalarSize_ ) =
                    observationWeights_.getObservationWeightMatrix( observationId, row.scalarSize_ );
            currentIndex += row.scalarSize_;
        }
        return weightMatrix;
    }

    //! Return one vector-valued residual per observation event in a set.
    std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > getResidualsForSet( const ObservationSetId setId ) const
    {
        std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > residuals;
        for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
        {
            residuals.push_back( getResidualValue( observationId ) );
        }
        return residuals;
    }

    //! Return all scalar residual components of a set as one legacy flat vector.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getResidualVectorForSet( const ObservationSetId setId ) const
    {
        return createSetVector( setId, residualValues_ );
    }

    //! Reconstruct the vector-valued residual for one observation event.
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

    //! Return per-observation dependent-variable vectors for one set.
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

    //! Extract one dependent-variable block by column start and size.
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

    //! Extract one dependent-variable block compatible with the requested settings.
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

    //! List dependent-variable settings in a set that match the requested settings.
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

    //! Extract every dependent-variable block compatible with the requested settings.
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

    //! Replace per-observation dependent-variable vectors for one set.
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

    //! Clear dependent-variable values for every observation event in one set.
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

    //! Translate a legacy ObservationCollectionParser into dataset set ids.
    std::vector< ObservationSetId > getObservationSetIds( const std::shared_ptr< ObservationCollectionParser >& observationParser =
                                                                  std::make_shared< ObservationCollectionParser >( ) ) const
    {
        std::vector< ObservationSetId > setIds;
        for( const ObservationSetId setId : getSetIdsInLegacyOrder( ) )
        {
            if( isObservationSetSelectedByLegacyParser( *this, setId, observationParser ) )
            {
                setIds.push_back( setId );
            }
        }
        return setIds;
    }

    //! Return observation ids for rows that satisfy a new row-level condition.
    std::vector< ObservationId > getObservationIdsMatchingCondition(
            const ObservationCondition< ObservationScalarType, TimeType >& condition ) const
    {
        std::vector< ObservationId > observationIds;
        for( ObservationId observationId = 0; observationId < observationRows_.size( ); ++observationId )
        {
            if( condition( *this, observationId ) )
            {
                observationIds.push_back( observationId );
            }
        }
        return observationIds;
    }

    //! Create a read-only view over observations satisfying a row-level condition.
    ObservationDatasetViewer< ObservationScalarType, TimeType > createViewer(
            const ObservationCondition< ObservationScalarType, TimeType >& condition ) const;

    //! Create a new dataset containing only observations satisfying a condition.
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createNewAndKeep(
            const ObservationCondition< ObservationScalarType, TimeType >& condition ) const
    {
        std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > reducedDataset =
                std::make_shared< ObservationDataset< ObservationScalarType, TimeType > >( );

        for( ObservationSetId setId = 0; setId < setMetadata_.size( ); ++setId )
        {
            std::vector< ObservationId > selectedObservationIds;
            for( const ObservationId observationId : observationIdsBySet_.at( setId ) )
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

                for( const ObservationId observationId : selectedObservationIds )
                {
                    observations.push_back( getObservationValue( observationId ) );
                    times.push_back( getObservationTime( observationId ) );
                    weights.push_back( getWeightValue( observationId ) );
                    residuals.push_back( getResidualValue( observationId ) );
                    if( hasDependentVariables )
                    {
                        dependentVariables.push_back( getDependentVariables( observationId ) );
                    }
                }

                const ObservationSetId newSetId =
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

                if( hasWeightMatrixForSet( setId ) )
                {
                    std::vector< std::size_t > selectedSetScalarIndices;
                    for( const ObservationId observationId : selectedObservationIds )
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
                else
                {
                    const std::vector< ObservationId >& newObservationIds = reducedDataset->getObservationIdsForSet( newSetId );
                    for( std::size_t i = 0; i < selectedObservationIds.size( ); ++i )
                    {
                        reducedDataset->setWeightMatrixForObservation( newObservationIds.at( i ),
                                                                       getWeightMatrixForObservation( selectedObservationIds.at( i ) ) );
                    }
                }

                for( std::size_t i = 0; i < selectedObservationIds.size( ); ++i )
                {
                    const ObservationDatasetRow< TimeType >& sourceRow = observationRows_.at( selectedObservationIds.at( i ) );
                    ObservationDatasetRow< TimeType >& targetRow =
                            reducedDataset->observationRows_.at( reducedDataset->observationIdsBySet_.at( newSetId ).at( i ) );
                    targetRow.isActive_ = sourceRow.isActive_;
                    targetRow.rejectionReason_ = sourceRow.rejectionReason_;
                }
            }
        }

        return reducedDataset;
    }

    //! Create a new dataset containing all observations except those satisfying a condition.
    std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > > createNewAndDrop(
            const ObservationCondition< ObservationScalarType, TimeType >& condition ) const
    {
        return createNewAndKeep( !condition );
    }

    //! Mark matching observation rows as rejected without physically deleting data.
    void rejectObservations( const ObservationCondition< ObservationScalarType, TimeType >& condition, const std::string& reason = "" )
    {
        for( ObservationId observationId = 0; observationId < observationRows_.size( ); ++observationId )
        {
            if( condition( *this, observationId ) )
            {
                observationRows_.at( observationId ).isActive_ = false;
                observationRows_.at( observationId ).rejectionReason_ = reason;
            }
        }
    }

    //! Restore matching observation rows to active status.
    void restoreObservations( const ObservationCondition< ObservationScalarType, TimeType >& condition )
    {
        for( ObservationId observationId = 0; observationId < observationRows_.size( ); ++observationId )
        {
            if( condition( *this, observationId ) )
            {
                observationRows_.at( observationId ).isActive_ = true;
                observationRows_.at( observationId ).rejectionReason_.clear( );
            }
        }
    }

    //! Materialize the active-by-default projection used by estimation routines.
    EstimationVectorProjection< ObservationScalarType, TimeType > createEstimationProjection( const bool includeRejected = false ) const
    {
        return createProjectionFromObservationIds( getAllObservationIds( ), includeRejected );
    }

    //! Materialize the flat legacy/estimation vector view of the dataset.
    /*!
     * Projection order follows dataset observation-row order and scalar
     * component order within each row. Times and ids are repeated per scalar
     * component, matching legacy concatenated-vector conventions.
     */
    EstimationVectorProjection< ObservationScalarType, TimeType > createLegacyProjection( const bool includeInactive = true ) const
    {
        return createProjectionFromObservationIds( getObservationIdsInLegacyOrder( ), includeInactive );
    }

    //! Return set ids in the legacy observable-type/link-ends/index ordering.
    std::vector< ObservationSetId > getSetIdsInLegacyOrder( ) const
    {
        std::map< ObservableType, std::map< LinkEnds, std::vector< ObservationSetId > > > setIdsByObservableAndLinkEnds;
        for( ObservationSetId setId = 0; setId < setMetadata_.size( ); ++setId )
        {
            const ObservationSetMetadata< ObservationScalarType, TimeType >& metadata = setMetadata_.at( setId );
            setIdsByObservableAndLinkEnds[ metadata.observableType_ ][ linkDefinitionRegistry_.at( metadata.linkDefinitionId_ ).linkEnds_ ]
                    .push_back( setId );
        }

        std::vector< ObservationSetId > setIds;
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

    std::size_t getStructuralVersion( ) const
    {
        return structuralVersion_;
    }

private:
    template< typename DatasetObservationScalarType,
              typename DatasetTimeType,
              typename std::enable_if< is_state_scalar_and_time_type< DatasetObservationScalarType, DatasetTimeType >::value, int >::type >
    friend class ObservationDatasetViewer;

    //! Return all observation ids in dataset row order.
    std::vector< ObservationId > getAllObservationIds( ) const
    {
        std::vector< ObservationId > observationIds;
        observationIds.reserve( observationRows_.size( ) );
        for( ObservationId observationId = 0; observationId < observationRows_.size( ); ++observationId )
        {
            observationIds.push_back( observationId );
        }
        return observationIds;
    }

    //! Return observation ids in legacy set order and row order within each set.
    std::vector< ObservationId > getObservationIdsInLegacyOrder( ) const
    {
        std::vector< ObservationId > observationIds;
        observationIds.reserve( observationRows_.size( ) );
        for( const ObservationSetId setId : getSetIdsInLegacyOrder( ) )
        {
            const std::vector< ObservationId >& setObservationIds = observationIdsBySet_.at( setId );
            observationIds.insert( observationIds.end( ), setObservationIds.begin( ), setObservationIds.end( ) );
        }
        return observationIds;
    }

    //! Materialize a projection for selected observation ids.
    EstimationVectorProjection< ObservationScalarType, TimeType > createProjectionFromObservationIds(
            const std::vector< ObservationId >& selectedObservationIds,
            const bool includeInactive = true ) const
    {
        EstimationVectorProjection< ObservationScalarType, TimeType > projection;

        std::size_t projectionSize = 0;
        bool materializeWeightMatrix = false;
        std::map< ScalarComponentId, bool > selectedScalarComponentIds;
        std::map< ObservationSetId, bool > selectedSetIds;
        for( const ObservationId observationId : selectedObservationIds )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            if( includeInactive || row.isActive_ )
            {
                projectionSize += row.scalarSize_;
                selectedSetIds[ row.setId_ ] = true;
                if( !observationWeights_.isObservationWeightDiagonalOnly( observationId, row.scalarSize_ ) )
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
            if( observationWeights_.hasSetWeightBlock( selectedSet.first ) &&
                !observationWeights_.isSetWeightBlockDiagonalOnly( selectedSet.first ) )
            {
                materializeWeightMatrix = true;
            }
        }

        for( const ObservationWeightBlock& extraWeightBlock : observationWeights_.getExtraWeightBlocks( ) )
        {
            for( std::size_t i = 0; i < extraWeightBlock.rowScalarComponentIds_.size( ); ++i )
            {
                const ScalarComponentId rowScalarComponentId = extraWeightBlock.rowScalarComponentIds_.at( i );
                if( selectedScalarComponentIds.count( rowScalarComponentId ) == 0 )
                {
                    continue;
                }
                for( std::size_t j = 0; j < extraWeightBlock.columnScalarComponentIds_.size( ); ++j )
                {
                    const ScalarComponentId columnScalarComponentId = extraWeightBlock.columnScalarComponentIds_.at( j );
                    if( selectedScalarComponentIds.count( columnScalarComponentId ) > 0 &&
                        rowScalarComponentId != columnScalarComponentId && extraWeightBlock.weightBlock_( i, j ) != 0.0 )
                    {
                        materializeWeightMatrix = true;
                    }
                }
            }
        }

        projection.observations_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( projectionSize );
        projection.residuals_ = Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( projectionSize );
        projection.weights_ = Eigen::VectorXd::Zero( projectionSize );
        if( materializeWeightMatrix )
        {
            projection.weightMatrix_.resize( projectionSize, projectionSize );
        }
        projection.times_.reserve( projectionSize );
        projection.observationIds_.reserve( projectionSize );
        projection.setIds_.reserve( projectionSize );
        projection.scalarComponentIds_.reserve( projectionSize );

        std::size_t currentIndex = 0;
        std::map< ScalarComponentId, std::size_t > projectionIndexByScalarComponent;
        std::map< std::pair< std::size_t, std::size_t >, double > sparseWeightEntries;
        auto setProjectionWeightEntry = [ &sparseWeightEntries ](
                                                const std::size_t rowIndex, const std::size_t columnIndex, const double weight ) {
            const std::pair< std::size_t, std::size_t > indexPair = std::make_pair( rowIndex, columnIndex );
            if( weight != 0.0 )
            {
                sparseWeightEntries[ indexPair ] = weight;
            }
            else
            {
                sparseWeightEntries.erase( indexPair );
            }
        };
        for( const ObservationId observationId : selectedObservationIds )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            if( includeInactive || row.isActive_ )
            {
                const Eigen::VectorXd observationWeightVector =
                        observationWeights_.getObservationWeightVector( observationId, row.scalarSize_ );
                for( unsigned int componentIndex = 0; componentIndex < row.scalarSize_; ++componentIndex )
                {
                    const ScalarComponentId scalarComponentId = row.firstScalarComponent_ + componentIndex;
                    projection.observations_( currentIndex ) = observedValues_.at( scalarComponentId );
                    projection.residuals_( currentIndex ) = residualValues_.at( scalarComponentId );
                    projection.weights_( currentIndex ) = observationWeightVector( componentIndex );
                    projection.times_.push_back( row.time_ );
                    projection.observationIds_.push_back( observationId );
                    projection.setIds_.push_back( row.setId_ );
                    projection.scalarComponentIds_.push_back( scalarComponentId );
                    projectionIndexByScalarComponent[ scalarComponentId ] = currentIndex;
                    ++currentIndex;
                }
                if( materializeWeightMatrix )
                {
                    if( !observationWeights_.hasSetWeightBlock( row.setId_ ) )
                    {
                        const Eigen::MatrixXd observationWeightMatrix =
                                observationWeights_.getObservationWeightMatrix( observationId, row.scalarSize_ );
                        const std::size_t observationStartIndex = currentIndex - row.scalarSize_;
                        for( unsigned int rowComponentIndex = 0; rowComponentIndex < row.scalarSize_; ++rowComponentIndex )
                        {
                            for( unsigned int columnComponentIndex = 0; columnComponentIndex < row.scalarSize_; ++columnComponentIndex )
                            {
                                setProjectionWeightEntry( observationStartIndex + rowComponentIndex,
                                                          observationStartIndex + columnComponentIndex,
                                                          observationWeightMatrix( rowComponentIndex, columnComponentIndex ) );
                            }
                        }
                    }
                }
            }
        }

        for( ObservationSetId setId = 0; setId < setMetadata_.size( ); ++setId )
        {
            if( observationWeights_.hasSetWeightBlock( setId ) )
            {
                const Eigen::MatrixXd& setWeightBlock = observationWeights_.getSetWeightBlock( setId );
                for( const ObservationId rowObservationId : observationIdsBySet_.at( setId ) )
                {
                    const ObservationDatasetRow< TimeType >& row = observationRows_.at( rowObservationId );
                    if( includeInactive || row.isActive_ )
                    {
                        for( unsigned int rowComponentIndex = 0; rowComponentIndex < row.scalarSize_; ++rowComponentIndex )
                        {
                            const ScalarComponentId rowScalarComponentId = row.firstScalarComponent_ + rowComponentIndex;
                            if( projectionIndexByScalarComponent.count( rowScalarComponentId ) == 0 )
                            {
                                continue;
                            }
                            const std::size_t projectionRow = projectionIndexByScalarComponent.at( rowScalarComponentId );
                            const std::size_t setLocalRow = row.indexInSet_ * setMetadata_.at( setId ).observableSize_ + rowComponentIndex;
                            projection.weights_( projectionRow ) = setWeightBlock( setLocalRow, setLocalRow );

                            if( materializeWeightMatrix )
                            {
                                for( const ObservationId columnObservationId : observationIdsBySet_.at( setId ) )
                                {
                                    const ObservationDatasetRow< TimeType >& columnRow = observationRows_.at( columnObservationId );
                                    if( includeInactive || columnRow.isActive_ )
                                    {
                                        for( unsigned int columnComponentIndex = 0; columnComponentIndex < columnRow.scalarSize_;
                                             ++columnComponentIndex )
                                        {
                                            const ScalarComponentId columnScalarComponentId =
                                                    columnRow.firstScalarComponent_ + columnComponentIndex;
                                            if( projectionIndexByScalarComponent.count( columnScalarComponentId ) == 0 )
                                            {
                                                continue;
                                            }
                                            const std::size_t projectionColumn =
                                                    projectionIndexByScalarComponent.at( columnScalarComponentId );
                                            const std::size_t setLocalColumn =
                                                    columnRow.indexInSet_ * setMetadata_.at( setId ).observableSize_ + columnComponentIndex;
                                            setProjectionWeightEntry(
                                                    projectionRow, projectionColumn, setWeightBlock( setLocalRow, setLocalColumn ) );
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }

        for( const ObservationWeightBlock& extraWeightBlock : observationWeights_.getExtraWeightBlocks( ) )
        {
            for( std::size_t i = 0; i < extraWeightBlock.rowScalarComponentIds_.size( ); ++i )
            {
                const ScalarComponentId rowScalarComponentId = extraWeightBlock.rowScalarComponentIds_.at( i );
                if( projectionIndexByScalarComponent.count( rowScalarComponentId ) == 0 )
                {
                    continue;
                }
                for( std::size_t j = 0; j < extraWeightBlock.columnScalarComponentIds_.size( ); ++j )
                {
                    const ScalarComponentId columnScalarComponentId = extraWeightBlock.columnScalarComponentIds_.at( j );
                    if( projectionIndexByScalarComponent.count( columnScalarComponentId ) == 0 )
                    {
                        continue;
                    }
                    if( rowScalarComponentId == columnScalarComponentId )
                    {
                        projection.weights_( projectionIndexByScalarComponent.at( rowScalarComponentId ) ) =
                                extraWeightBlock.weightBlock_( i, j );
                    }
                    if( materializeWeightMatrix )
                    {
                        setProjectionWeightEntry( projectionIndexByScalarComponent.at( rowScalarComponentId ),
                                                  projectionIndexByScalarComponent.at( columnScalarComponentId ),
                                                  extraWeightBlock.weightBlock_( i, j ) );
                    }
                }
            }
        }

        if( materializeWeightMatrix )
        {
            std::vector< Eigen::Triplet< double > > sparseWeightTriplets;
            sparseWeightTriplets.reserve( sparseWeightEntries.size( ) );
            for( const auto& weightEntry : sparseWeightEntries )
            {
                sparseWeightTriplets.emplace_back( weightEntry.first.first, weightEntry.first.second, weightEntry.second );
                if( weightEntry.first.first == weightEntry.first.second )
                {
                    projection.weights_( weightEntry.first.first ) = weightEntry.second;
                }
            }
            projection.weightMatrix_.setFromTriplets( sparseWeightTriplets.begin( ), sparseWeightTriplets.end( ) );
            projection.weightMatrix_.makeCompressed( );
            projection.isDiagonalWeightOnly_ = false;
        }
        else
        {
            projection.isDiagonalWeightOnly_ = true;
        }

        return projection;
    }

public:
    //! Return flat scalar-vector start and size for each set in dataset insertion order.
    std::vector< std::pair< int, int > > getObservationSetStartAndSizeInDatasetOrder( ) const
    {
        std::vector< std::pair< int, int > > startAndSize;
        startAndSize.reserve( getNumberOfObservationSets( ) );

        int currentIndex = 0;
        for( ObservationSetId setId = 0; setId < getNumberOfObservationSets( ); ++setId )
        {
            const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
            startAndSize.push_back( std::make_pair( currentIndex, currentSize ) );
            currentIndex += currentSize;
        }
        return startAndSize;
    }

    //! Return flat scalar-vector start and size for each set in legacy projection order.
    /*!
     * This method preserves the old ObservationCollection vector-order contract:
     * observable type, then link ends, then set index within that group. Use
     * getObservationSetStartAndSizeInDatasetOrder() when indices must follow
     * dataset set ids instead.
     */
    std::vector< std::pair< int, int > > getObservationSetStartAndSize( ) const
    {
        std::vector< std::pair< int, int > > startAndSize;
        startAndSize.reserve( getNumberOfObservationSets( ) );

        int currentIndex = 0;
        for( const ObservationSetId setId : getSetIdsInLegacyOrder( ) )
        {
            const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
            startAndSize.push_back( std::make_pair( currentIndex, currentSize ) );
            currentIndex += currentSize;
        }
        return startAndSize;
    }

    //! Return flat scalar-vector start and size grouped by observable type.
    std::map< ObservableType, std::pair< int, int > > getObservableTypeStartAndSize( ) const
    {
        std::map< ObservableType, std::pair< int, int > > startAndSize;
        int currentIndex = 0;

        for( const ObservationSetId setId : getSetIdsInLegacyOrder( ) )
        {
            const ObservableType observableType = getObservationSetMetadata( setId ).observableType_;
            const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
            if( startAndSize.count( observableType ) == 0 )
            {
                startAndSize[ observableType ] = std::make_pair( currentIndex, currentSize );
            }
            else
            {
                startAndSize[ observableType ].second += currentSize;
            }
            currentIndex += currentSize;
        }
        return startAndSize;
    }

    //! Set residuals from one flat scalar vector in dataset projection order.
    void setResidualVector( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& residualVector )
    {
        if( residualVector.size( ) != static_cast< int >( getTotalScalarSize( ) ) )
        {
            throw std::runtime_error(
                    "Error when setting dataset residual vector, input size is inconsistent with total scalar observation size." );
        }

        int currentIndex = 0;
        for( ObservationSetId setId = 0; setId < getNumberOfObservationSets( ); ++setId )
        {
            const int currentSize = static_cast< int >( getTotalScalarSizeForSet( setId ) );
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > setResiduals = residualVector.segment( currentIndex, currentSize );
            setResidualVectorForSet( setId, setResiduals );
            currentIndex += currentSize;
        }
    }

private:  //! Validate per-observation vectors before replacing/appending set data.
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

    //! Return a permutation that sorts observation events by their time tags.
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

    //! Apply a precomputed observation-event permutation to one parallel vector.
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

    //! Remove observation-event entries from a parallel vector by index.
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

    //! Gather scalar values for one set into a legacy flat vector.
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

    //! Overwrite all scalar components belonging to one observation event.
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

    //! Overwrite all residual scalar components belonging to one observation event.
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

    //! Overwrite all weight scalar components belonging to one observation event.
    void setWeightValue( const ObservationId observationId, const Eigen::VectorXd& weight )
    {
        const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
        if( weight.size( ) != static_cast< int >( row.scalarSize_ ) )
        {
            throw std::runtime_error( "Error when setting dataset weight, scalar size is inconsistent." );
        }
        observationWeights_.setDiagonalWeightVector( observationId, weight );
    }

    //! Expand observation ids and optional component indices to scalar-component ids.
    std::vector< ScalarComponentId > getScalarComponentIdsForObservationSelection( const std::vector< ObservationId >& observationIds,
                                                                                   const std::vector< unsigned int >& components ) const
    {
        std::vector< ScalarComponentId > scalarComponentIds;
        for( const ObservationId observationId : observationIds )
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

    //! Return the registry id for a link definition, inserting it if needed.
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

    //! Return the registry id for ancillary settings, inserting them if needed.
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

    //! Return the registry id for dependent-variable bookkeeping, inserting it if needed.
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

    //! One row per observation event; vector observables occupy one row, not N rows.
    std::vector< ObservationDatasetRow< TimeType > > observationRows_;
    //! One row per scalar component; maps scalar storage entries back to observation rows.
    std::vector< ObservationScalarComponentRow > scalarComponentRows_;
    //! One metadata record per observation set.
    std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > > setMetadata_;
    //! For each set id, the ordered observation ids belonging to that set.
    std::vector< std::vector< ObservationId > > observationIdsBySet_;

    //! Registry of full link definitions shared by observation sets.
    std::vector< LinkDefinition > linkDefinitionRegistry_;
    //! Registry of ancillary settings pointers shared by observation sets.
    std::vector< std::shared_ptr< ObservationAncillarySimulationSettings > > ancillarySettingsRegistry_;
    //! Registry of dependent-variable bookkeeping/layout pointers shared by observation sets.
    std::vector< std::shared_ptr< simulation_setup::ObservationDependentVariableBookkeeping > > dependentVariableLayoutRegistry_;

    //! Scalar observed values; vector observables contribute one entry per component.
    std::vector< ObservationScalarType > observedValues_;
    //! Scalar residual values aligned one-to-one with observedValues_.
    std::vector< ObservationScalarType > residualValues_;
    //! Compact observation weight storage; materialized into vectors/matrices only on request.
    ObservationWeights observationWeights_;
    //! Per-observation dependent-variable vectors aligned one-to-one with observationRows_.
    std::vector< Eigen::VectorXd > dependentVariableValues_;
    //! Monotonic counter used to invalidate viewers after structural mutations.
    std::size_t structuralVersion_ = 0;
};
}  // namespace observation_models

}  // namespace tudat

#include "tudat/simulation/estimation_setup/observationLegacyParserAdapter.h"
#include "tudat/simulation/estimation_setup/observationDatasetViewer.h"

#endif  // TUDAT_OBSERVATION_DATASET_H
