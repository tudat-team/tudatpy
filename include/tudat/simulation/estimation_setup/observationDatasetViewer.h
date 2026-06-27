/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATION_DATASET_VIEWER_H
#define TUDAT_OBSERVATION_DATASET_VIEWER_H

#include <cstddef>
#include <stdexcept>
#include <vector>

#include <Eigen/Core>

#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/estimationVectorProjection.h"
#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"

namespace tudat
{

namespace observation_models
{

//! Read-only view on a selected subset of an ObservationDataset.
/*!
 * A viewer references its parent dataset and a fixed list of observation ids.
 * It exposes inspection and projection operations only: no mutation, rejection,
 * restoration or metadata editing. Viewers are invalidated when the parent
 * dataset structurally changes by adding, removing or rebuilding observation
 * rows; non-structural value/status updates keep ids valid.
 */
template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type >
class ObservationDatasetViewer
{
public:
    ObservationDatasetViewer( const ObservationDataset< ObservationScalarType, TimeType >& dataset,
                              const std::vector< ObservationId >& observationIds,
                              const std::size_t structuralVersion ):
        dataset_( dataset ), observationIds_( observationIds ), structuralVersion_( structuralVersion )
    {}

    std::size_t getNumberOfObservations( ) const
    {
        checkValidity( );
        return observationIds_.size( );
    }

    const std::vector< ObservationId >& getObservationIds( ) const
    {
        checkValidity( );
        return observationIds_;
    }

    const ObservationDatasetRow< TimeType >& getObservationRow( const std::size_t viewerIndex ) const
    {
        checkValidity( );
        return dataset_.getObservationRow( observationIds_.at( viewerIndex ) );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationValue( const std::size_t viewerIndex ) const
    {
        checkValidity( );
        return dataset_.getObservationValue( observationIds_.at( viewerIndex ) );
    }

    TimeType getObservationTime( const std::size_t viewerIndex ) const
    {
        checkValidity( );
        return dataset_.getObservationTime( observationIds_.at( viewerIndex ) );
    }

    ObservationDatasetViewer createViewer( const ObservationCondition< ObservationScalarType, TimeType >& condition ) const
    {
        checkValidity( );
        std::vector< ObservationId > narrowedObservationIds;
        for( const ObservationId observationId : observationIds_ )
        {
            if( condition( dataset_, observationId ) )
            {
                narrowedObservationIds.push_back( observationId );
            }
        }
        return ObservationDatasetViewer( dataset_, narrowedObservationIds, structuralVersion_ );
    }

    EstimationVectorProjection< ObservationScalarType, TimeType > createEstimationProjection( const bool includeRejected = false ) const
    {
        checkValidity( );
        return dataset_.createProjectionFromObservationIds( observationIds_, includeRejected );
    }

    EstimationVectorProjection< ObservationScalarType, TimeType > createLegacyProjection( const bool includeInactive = true ) const
    {
        checkValidity( );
        return dataset_.createProjectionFromObservationIds( observationIds_, includeInactive );
    }

private:
    //! Throw if the parent dataset has structurally changed since viewer creation.
    void checkValidity( ) const
    {
        if( structuralVersion_ != dataset_.getStructuralVersion( ) )
        {
            throw std::runtime_error(
                    "Error when using observation dataset viewer, parent dataset has been structurally modified since viewer creation." );
        }
    }

    //! Non-owning parent dataset reference; caller must keep the dataset alive.
    const ObservationDataset< ObservationScalarType, TimeType >& dataset_;
    //! Observation ids selected by this viewer, in projection order.
    std::vector< ObservationId > observationIds_;
    //! Parent structural version captured when the viewer was created.
    std::size_t structuralVersion_;
};

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Enable >
ObservationDatasetViewer< ObservationScalarType, TimeType > ObservationDataset< ObservationScalarType, TimeType, Enable >::createViewer(
        const ObservationCondition< ObservationScalarType, TimeType >& condition ) const
{
    return ObservationDatasetViewer< ObservationScalarType, TimeType >(
            *this, getObservationIdsMatchingCondition( condition ), structuralVersion_ );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_DATASET_VIEWER_H
