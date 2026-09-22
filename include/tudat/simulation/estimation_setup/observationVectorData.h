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
#include <Eigen/Cholesky>
#include <Eigen/SparseCore>

#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationWeights.h"

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

    //! Return one observation's measurement covariance, using the complete weights captured from the dataset.
    /*!
     * observationId identifies any observation stored in the dataset when this snapshot was created,
     * including rejected observations and observations outside this object's vector selection.
     * The returned matrix has one row and column per observable component (e.g. 2 by 2 for angular position).
     * For correlated observations, invert the complete set, including rejected observations, and return
     * the requested observation's diagonal block of that inverse. Do not invert only its own weight block.
     * Weights connecting different sets, and singular or non-positive-definite weights, are unsupported.
     * Results are cached on first access; concurrent calls on the same object are not supported.
     */
    const Eigen::MatrixXd& getInverseWeightMatrixForObservation( const unsigned int observationId ) const
    {
        const auto observation = completeWeightObservationMapping_.find( observationId );
        if( observation == completeWeightObservationMapping_.end( ) )
        {
            throw std::runtime_error( "Observation is not present in the complete weight snapshot." );
        }
        const unsigned int setId = observation->second.second;
        const auto& metadata = completeWeightSetMetadata_.at( setId );
        if( metadata.weightStructure_ == ObservationWeightStructure::inter_set_weights )
        {
            throw std::runtime_error( "Per-observation covariance is unsupported for weights connecting different observation sets." );
        }
        const auto cached = inverseWeightsByObservation_.find( observationId );
        if( cached != inverseWeightsByObservation_.end( ) )
        {
            return cached->second;
        }
        const unsigned int dimension = metadata.observableSize_;
        const std::vector< unsigned int > ids = metadata.weightStructure_ == ObservationWeightStructure::per_set
                ? completeWeightObservationIdsBySet_.at( setId )
                : std::vector< unsigned int >{ observationId };
        std::vector< unsigned int > indices;
        indices.reserve( ids.size( ) * dimension );
        for( const unsigned int id : ids )
        {
            const unsigned int first = completeWeightObservationMapping_.at( id ).first;
            for( unsigned int component = 0; component < dimension; ++component )
            {
                indices.push_back( first + component );
            }
        }
        const Eigen::VectorXd diagonal = completeWeights_.getDiagonal( indices );
        if( ( diagonal.array( ) <= 0.0 ).any( ) )
        {
            throw std::runtime_error( "Measurement covariance requires strictly positive observation weights." );
        }
        Eigen::MatrixXd covariance;
        if( metadata.weightStructure_ == ObservationWeightStructure::diagonal )
        {
            covariance = diagonal.cwiseInverse( ).asDiagonal( );
        }
        else
        {
            Eigen::MatrixXd block = diagonal.asDiagonal( );
            const auto& entries = completeWeights_.getOffDiagonalEntries( );
            for( std::size_t i = 0; i < indices.size( ); ++i )
            {
                for( std::size_t j = i + 1; j < indices.size( ); ++j )
                {
                    const auto entry = entries.find( std::minmax( indices.at( i ), indices.at( j ) ) );
                    if( entry != entries.end( ) )
                    {
                        block( i, j ) = block( j, i ) = entry->second;
                    }
                }
            }
            const Eigen::LLT< Eigen::MatrixXd > factorization( block );
            if( factorization.info( ) != Eigen::Success )
            {
                throw std::runtime_error( "Measurement covariance requires a nonsingular, positive-definite weight matrix." );
            }
            covariance = factorization.solve( Eigen::MatrixXd::Identity( block.rows( ), block.cols( ) ) );
        }
        if( !covariance.allFinite( ) )
        {
            throw std::runtime_error( "Observation weight inversion produced a non-finite measurement covariance." );
        }
        // Invert each complete set only once, retaining only the small per-observation covariance blocks.
        for( std::size_t i = 0; i < ids.size( ); ++i )
        {
            inverseWeightsByObservation_.emplace( ids.at( i ), covariance.block( i * dimension, i * dimension, dimension, dimension ) );
        }
        return inverseWeightsByObservation_.at( observationId );
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

    //! Complete dataset weights at snapshot creation, before selection or rejection removes any observations from the vectors.
    ObservationWeights completeWeights_;
    //! Set metadata for the complete weight snapshot; indexed by dataset set id.
    std::vector< ObservationSetMetadata< ObservationScalarType, TimeType > > completeWeightSetMetadata_;
    //! All observation ids in each set, including rejected ones; the outer index is the dataset set id.
    std::vector< std::vector< unsigned int > > completeWeightObservationIdsBySet_;
    //! For each observation id, its first scalar index in completeWeights_ and its dataset set id, respectively.
    std::unordered_map< unsigned int, std::pair< unsigned int, unsigned int > > completeWeightObservationMapping_;
    //! Cached diagonal blocks of the complete inverse weights, keyed by observation id, not vector row.
    mutable std::unordered_map< unsigned int, Eigen::MatrixXd > inverseWeightsByObservation_;

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
