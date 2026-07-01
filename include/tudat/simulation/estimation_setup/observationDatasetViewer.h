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
#include <memory>
#include <stdexcept>
#include <unordered_set>
#include <vector>

#include <Eigen/Core>

#include "tudat/simulation/estimation_setup/observationDataset.h"
#include "tudat/simulation/estimation_setup/flattenedObservationData.h"
#include "tudat/simulation/estimation_setup/observationDatasetRows.h"
#include "tudat/simulation/estimation_setup/observationCondition.h"

namespace tudat
{

namespace observation_models
{

//! Read-only view on a selected subset of an ObservationDataset.
/*!
 * A viewer references its parent dataset and a fixed list of observation ids.
 * It exposes inspection and flattened-observation-data operations only: no mutation, rejection,
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
                              const std::vector< unsigned int >& observationIds,
                              const std::size_t structuralVersion ):
        dataset_( &dataset ), observationIds_( observationIds ), structuralVersion_( structuralVersion )
    {}

    ObservationDatasetViewer( const std::shared_ptr< const ObservationDataset< ObservationScalarType, TimeType > >& dataset,
                              const std::vector< unsigned int >& observationIds,
                              const std::size_t structuralVersion ):
        ownedDataset_( dataset ), dataset_( dataset.get( ) ), observationIds_( observationIds ), structuralVersion_( structuralVersion )
    {}

    std::size_t getNumberOfObservations( ) const
    {
        checkValidity( );
        return observationIds_.size( );
    }

    const std::vector< unsigned int >& getObservationIds( ) const
    {
        checkValidity( );
        return observationIds_;
    }

    const ObservationDatasetRow< TimeType >& getObservationRow( const std::size_t viewerIndex ) const
    {
        checkValidity( );
        return dataset( ).getObservationRow( observationIds_.at( viewerIndex ) );
    }

    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getObservationValue( const std::size_t viewerIndex ) const
    {
        checkValidity( );
        return dataset( ).getObservationValue( observationIds_.at( viewerIndex ) );
    }

    TimeType getObservationTime( const std::size_t viewerIndex ) const
    {
        checkValidity( );
        return dataset( ).getObservationTime( observationIds_.at( viewerIndex ) );
    }

    ObservationDatasetViewer createViewer( const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
    {
        checkValidity( );
        std::vector< unsigned int > narrowedObservationIds;
        for( const unsigned int observationId : observationIds_ )
        {
            if( condition( dataset( ), observationId ) )
            {
                narrowedObservationIds.push_back( observationId );
            }
        }
        if( ownedDataset_ != nullptr )
        {
            return ObservationDatasetViewer( ownedDataset_, narrowedObservationIds, structuralVersion_ );
        }
        return ObservationDatasetViewer( dataset( ), narrowedObservationIds, structuralVersion_ );
    }

    FlattenedObservationData< ObservationScalarType, TimeType > createEstimationFlattenedObservationData(
            const bool includeRejected = false ) const
    {
        checkValidity( );
        return dataset( ).createFlattenedObservationDataFromObservationIds( observationIds_, includeRejected );
    }

    FlattenedObservationData< ObservationScalarType, TimeType > createOrderedFlattenedObservationData(
            const bool includeInactive = true ) const
    {
        checkValidity( );
        const std::unordered_set< unsigned int > selectedObservationIds( observationIds_.begin( ), observationIds_.end( ) );
        std::vector< unsigned int > orderedSelectedObservationIds;
        orderedSelectedObservationIds.reserve( observationIds_.size( ) );
        for( const unsigned int observationId : dataset( ).getObservationIdsInOrderedFlattenedDataOrder( ) )
        {
            if( selectedObservationIds.count( observationId ) > 0 )
            {
                orderedSelectedObservationIds.push_back( observationId );
            }
        }
        return dataset( ).createFlattenedObservationDataFromObservationIds( orderedSelectedObservationIds, includeInactive );
    }

private:
    const ObservationDataset< ObservationScalarType, TimeType >& dataset( ) const
    {
        if( dataset_ == nullptr )
        {
            throw std::runtime_error( "Error when using observation dataset viewer, parent dataset is null." );
        }
        return *dataset_;
    }

    //! Throw if the parent dataset has structurally changed since viewer creation.
    void checkValidity( ) const
    {
        if( structuralVersion_ != dataset( ).getStructuralVersion( ) )
        {
            throw std::runtime_error(
                    "Error when using observation dataset viewer, parent dataset has been structurally modified since viewer creation." );
        }
    }

    //! Optional owning parent dataset pointer used when the dataset was created by shared ownership.
    std::shared_ptr< const ObservationDataset< ObservationScalarType, TimeType > > ownedDataset_;

    //! Parent dataset pointer; backed by ownedDataset_ when shared ownership is available.
    const ObservationDataset< ObservationScalarType, TimeType >* dataset_;

    //! Observation ids selected by this viewer, in flattened-data row order.
    std::vector< unsigned int > observationIds_;

    //! Parent structural version captured when the viewer was created.
    std::size_t structuralVersion_;
};

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Enable >
ObservationDatasetViewer< ObservationScalarType, TimeType > ObservationDataset< ObservationScalarType, TimeType, Enable >::createViewer(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
{
    try
    {
        return ObservationDatasetViewer< ObservationScalarType, TimeType >(
                this->shared_from_this( ), getObservationIdsMatchingCondition( condition ), structuralVersion_ );
    }
    catch( const std::bad_weak_ptr& )
    {}

    return ObservationDatasetViewer< ObservationScalarType, TimeType >(
            *this, getObservationIdsMatchingCondition( condition ), structuralVersion_ );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATION_DATASET_VIEWER_H
