/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_VECTOR_DATA_H
#define TUDAT_OBSERVATION_VECTOR_DATA_H

#include <stdexcept>
#include <vector>
#include <memory>
#include <unordered_map>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/simulation/estimation_setup/observationDatasetRows.h"

namespace tudat
{

namespace observation_models
{

//! Scalar-aligned vector representation of an ObservationDataset selection.
/*!
 * This object is not primary storage. It stores selected observation rows as
 * observation, residual and diagonal weight vectors in a fixed scalar order.
 * When block/correlation weights are present, it also stores the sparse weight
 * matrix. The id vectors map every scalar entry back to the dataset row, set
 * and scalar component from which it came.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationVectorData
{
public:
    //! Return the scalar observation vector represented by this snapshot.
    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& getObservationVector( ) const
    {
        return observations_;
    }

    //! Return residuals in the same scalar order as the observation vector.
    const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& getResidualVector( ) const
    {
        return residuals_;
    }

    //! Return diagonal weights in the same scalar order as the observation vector.
    const Eigen::VectorXd& getWeightVector( ) const
    {
        return weights_;
    }

    //! Return the full sparse weight matrix, materializing diagonal storage on demand.
    const Eigen::SparseMatrix< double >& getSparseWeightMatrix( ) const
    {
        // This lazy materialization mutates cached storage and is not safe for concurrent const access.
        if( weightMatrix_.rows( ) == 0 && weights_.size( ) > 0 )
        {
            // Diagonal-only data keep only the vector form until a matrix is explicitly requested.
            weightMatrix_.resize( weights_.size( ), weights_.size( ) );
            weightMatrix_.reserve( weights_.size( ) );
            for( int i = 0; i < weights_.size( ); ++i )
            {
                if( weights_( i ) != 0.0 )
                {
                    weightMatrix_.insert( i, i ) = weights_( i );
                }
            }
            weightMatrix_.makeCompressed( );
        }
        return weightMatrix_;
    }

    //! Return whether this snapshot contains diagonal weights only.
    bool isDiagonalWeightOnly( ) const
    {
        return isDiagonalWeightOnly_;
    }

    //! Return whether this snapshot contains at least one off-diagonal weight.
    bool hasOffDiagonalWeights( ) const
    {
        return !isDiagonalWeightOnly_;
    }

    //! Return one reference-link-end time for every scalar component.
    const std::vector< TimeType >& getTimes( ) const
    {
        return times_;
    }

    //! Return the observation identity associated with every scalar component.
    const std::vector< unsigned int >& getObservationIds( ) const
    {
        return observationIds_;
    }

    //! Return the observation-set identity associated with every scalar component.
    const std::vector< unsigned int >& getSetIds( ) const
    {
        return setIds_;
    }

    //! Return the link-definition identity associated with every scalar component.
    std::vector< unsigned int > getLinkDefinitionIds( ) const
    {
        std::vector< unsigned int > linkDefinitionIds;
        linkDefinitionIds.reserve( setIds_.size( ) );
        for( const unsigned int setId : setIds_ )
        {
            linkDefinitionIds.push_back( metadataBySet_.at( setId ).linkDefinitionId_ );
        }
        return linkDefinitionIds;
    }

    //! Return each component's index within its vector-valued observation.
    const std::vector< unsigned int >& getScalarComponentIds( ) const
    {
        return scalarComponentIds_;
    }

    //! Resolve an observation/component pair to its scalar row in this snapshot.
    int getVectorRow( const unsigned int observationId, const unsigned int componentIndex ) const
    {
        const auto row = rowMapping_.find( observationId );
        if( row == rowMapping_.end( ) || componentIndex >= row->second.second )
        {
            throw std::runtime_error( "Observation/component pair is not present in the observation vector data." );
        }
        return row->second.first + componentIndex;
    }

    //! Return the first scalar row of an observation, or -1 when absent from this snapshot.
    int getFirstVectorRowForObservation( const unsigned int observationId ) const
    {
        const auto row = rowMapping_.find( observationId );
        return row == rowMapping_.end( ) ? -1 : static_cast< int >( row->second.first );
    }

    //! Return the number of scalar rows for an observation, or zero when absent.
    unsigned int getScalarSizeForObservation( const unsigned int observationId ) const
    {
        const auto row = rowMapping_.find( observationId );
        return row == rowMapping_.end( ) ? 0 : row->second.second;
    }

    //! Return set identities in their first-appearance order in this snapshot.
    const std::vector< unsigned int >& getSetIdsInRowOrder( ) const
    {
        return setIdsInRowOrder_;
    }

    //! Return the unique observation identities represented for one set, in row order.
    const std::vector< unsigned int >& getUniqueObservationIdsForSetInRowOrder( const unsigned int setId ) const
    {
        if( setId >= uniqueObservationIdsBySet_.size( ) || uniqueObservationIdsBySet_.at( setId ).empty( ) )
        {
            throw std::runtime_error( "Error when retrieving observation vector rows, requested set is not present." );
        }
        return uniqueObservationIdsBySet_.at( setId );
    }

    //! Return the captured metadata for one represented set.
    const ObservationSetMetadata< ObservationScalarType, TimeType >& getSetMetadata( const unsigned int setId ) const
    {
        return metadataBySet_.at( setId );
    }
    //! Return the captured link definition for one represented set.
    const LinkDefinition& getLinkDefinitionForSet( const unsigned int setId ) const
    {
        return linksBySet_.at( setId );
    }
    //! Return an independent copy of captured ancillary settings for one set.
    std::shared_ptr< ObservationAncillarySimulationSettings > getAncillarySettingsForSet( const unsigned int setId ) const
    {
        const auto& settings = ancillaryBySet_.at( setId );
        return settings ? std::make_shared< ObservationAncillarySimulationSettings >( *settings ) : nullptr;
    }
    //! Return captured dependent-variable values for one observation identity.
    const Eigen::VectorXd& getDependentVariables( const unsigned int observationId ) const
    {
        return dependentVariables_.at( observationId );
    }

private:
    template< typename DatasetObservationScalarType,
              typename DatasetTimeType,
              typename std::enable_if< is_state_scalar_and_time_type< DatasetObservationScalarType, DatasetTimeType >::value, int >::type >
    friend class ObservationDataset;

    //! Observation scalar values in this object's row order.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations_;

    //! Residual scalar values in the same order as observations_.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals_;

    //! Diagonal scalar weights in the same order as observations_.
    Eigen::VectorXd weights_;

    //! Full weight matrix; left empty for diagonal-only data until requested.
    mutable Eigen::SparseMatrix< double > weightMatrix_;

    //! True when the full weight matrix has no off-diagonal entries.
    bool isDiagonalWeightOnly_ = true;

    //! Reference-link-end time for each scalar entry.
    std::vector< TimeType > times_;

    //! Observation id for each scalar entry.
    std::vector< unsigned int > observationIds_;

    //! Observation set id for each scalar entry.
    std::vector< unsigned int > setIds_;

    //! Scalar component id for each scalar entry.
    std::vector< unsigned int > scalarComponentIds_;

    //! Mapping from observation identity to its first scalar row and component count.
    std::unordered_map< unsigned int, std::pair< unsigned int, unsigned int > > rowMapping_;
    //! Weak identity token for validating writeback to the originating dataset.
    std::weak_ptr< const int > source_;
    //! Structural revision captured when this snapshot was created.
    std::size_t structuralVersion_ = 0;
    //! Value/selection revision captured when this snapshot was created.
    std::size_t vectorDataVersion_ = 0;
    //! Detached set metadata keyed by stable set identity.
    std::unordered_map< unsigned int, ObservationSetMetadata< ObservationScalarType, TimeType > > metadataBySet_;
    //! Detached link definitions keyed by stable set identity.
    std::unordered_map< unsigned int, LinkDefinition > linksBySet_;
    //! Detached ancillary settings keyed by stable set identity.
    std::unordered_map< unsigned int, std::shared_ptr< ObservationAncillarySimulationSettings > > ancillaryBySet_;
    //! Dependent-variable values keyed by stable observation identity.
    std::unordered_map< unsigned int, Eigen::VectorXd > dependentVariables_;

    //! Unique observation ids grouped by set, preserving this object's row order.
    std::vector< std::vector< unsigned int > > uniqueObservationIdsBySet_;

    //! Set ids in the order in which each set first appears in this object's rows.
    std::vector< unsigned int > setIdsInRowOrder_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_VECTOR_DATA_H
