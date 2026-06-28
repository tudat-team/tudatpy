/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_ESTIMATION_VECTOR_PROJECTION_H
#define TUDAT_ESTIMATION_VECTOR_PROJECTION_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/SparseCore>

#include "tudat/simulation/estimation_setup/observationDatasetRows.h"

namespace tudat
{

namespace observation_models
{

//! Derived estimator-vector view of an ObservationDataset selection.
/*!
 * This object is a materialized projection, not primary storage. It stores the
 * flat observation, residual and diagonal weight vectors in estimator order,
 * the full weight matrix when block/correlation weights are present, and the
 * mappings needed to trace every scalar entry back to the dataset.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
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

    const Eigen::SparseMatrix< double >& getSparseWeightMatrix( ) const
    {
        if( weightMatrix_.rows( ) == 0 && weights_.size( ) > 0 )
        {
            // Diagonal-only projections keep only the vector form until a matrix is explicitly requested.
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

    //! Observation scalar values in projection/estimator order.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > observations_;
    //! Residual scalar values in the same order as observations_.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > residuals_;
    //! Diagonal scalar weights in the same order as observations_.
    Eigen::VectorXd weights_;
    //! Full projection weight matrix; left empty for diagonal-only projections until requested.
    mutable Eigen::SparseMatrix< double > weightMatrix_;
    //! True when the full projection weight matrix has no off-diagonal entries.
    bool isDiagonalWeightOnly_ = true;
    //! Reference-link-end time for each scalar entry in the projection.
    std::vector< TimeType > times_;
    //! Observation id for each scalar entry in the projection.
    std::vector< ObservationId > observationIds_;
    //! Observation set id for each scalar entry in the projection.
    std::vector< ObservationSetId > setIds_;
    //! Scalar component id for each scalar entry in the projection.
    std::vector< ScalarComponentId > scalarComponentIds_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_ESTIMATION_VECTOR_PROJECTION_H
