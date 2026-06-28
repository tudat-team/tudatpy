/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_FLATTENED_OBSERVATION_DATA_H
#define TUDAT_FLATTENED_OBSERVATION_DATA_H

#include <stdexcept>
#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "tudat/simulation/estimation_setup/observationDatasetRows.h"

namespace tudat
{

namespace observation_models
{

//! Flattened vector representation of an ObservationDataset selection.
/*!
 * This object is not primary storage. It stores selected observation rows as
 * flat observation, residual and diagonal weight vectors in a fixed row order.
 * When block/correlation weights are present, it also stores the sparse weight
 * matrix. The id vectors map every scalar entry back to the dataset row, set
 * and scalar component from which it came.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class FlattenedObservationData
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

    const Eigen::SparseMatrix< double >& getSparseWeightMatrix( ) const
    {
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

    const Eigen::SparseMatrix< double >& getWeightMatrix( ) const
    {
        return getSparseWeightMatrix( );
    }

    bool isDiagonalWeightOnly( ) const
    {
        return isDiagonalWeightOnly_;
    }

    bool hasOffDiagonalWeights( ) const
    {
        return !isDiagonalWeightOnly_;
    }

    const std::vector< TimeType >& getTimes( ) const
    {
        return times_;
    }

    const std::vector< unsigned int >& getObservationIds( ) const
    {
        return observationIds_;
    }

    const std::vector< unsigned int >& getSetIds( ) const
    {
        return setIds_;
    }

    const std::vector< unsigned int >& getScalarComponentIds( ) const
    {
        return scalarComponentIds_;
    }

    int getFlattenedRow( const unsigned int observationId, const unsigned int componentIndex ) const
    {
        if( observationId >= firstFlattenedRowByObservation_.size( ) || firstFlattenedRowByObservation_.at( observationId ) < 0 ||
            componentIndex >= scalarSizeByObservation_.at( observationId ) )
        {
            throw std::runtime_error( "Error when retrieving flattened observation row, observation/component pair is not present." );
        }
        return firstFlattenedRowByObservation_.at( observationId ) + static_cast< int >( componentIndex );
    }

    const std::vector< unsigned int >& getSetIdsInRowOrder( ) const
    {
        return setIdsInRowOrder_;
    }

    const std::vector< unsigned int >& getUniqueObservationIdsForSetInRowOrder( const unsigned int setId ) const
    {
        if( setId >= uniqueObservationIdsBySet_.size( ) || uniqueObservationIdsBySet_.at( setId ).empty( ) )
        {
            throw std::runtime_error( "Error when retrieving flattened observation rows, requested set is not present." );
        }
        return uniqueObservationIdsBySet_.at( setId );
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
    //! First row in the flat vectors for each selected observation; -1 for unselected observations.
    std::vector< int > firstFlattenedRowByObservation_;
    //! Scalar size for each selected observation; 0 for unselected observations.
    std::vector< unsigned int > scalarSizeByObservation_;
    //! Unique observation ids grouped by set, preserving this object's row order.
    std::vector< std::vector< unsigned int > > uniqueObservationIdsBySet_;
    //! Set ids in the order in which each set first appears in this object's rows.
    std::vector< unsigned int > setIdsInRowOrder_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_FLATTENED_OBSERVATION_DATA_H
