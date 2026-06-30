/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONDATASETSELECTIONIMPLEMENTATION_H
#define TUDAT_OBSERVATIONDATASETSELECTIONIMPLEMENTATION_H

#include "tudat/simulation/estimation_setup/observationDataset.h"

namespace tudat
{

namespace observation_models
{

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
std::shared_ptr< ObservationDataset< ObservationScalarType, TimeType > >
ObservationDataset< ObservationScalarType, TimeType, Dummy >::createNewAndDrop(
        const ObservationCondition< ObservationScalarType, TimeType >& condition ) const
{
    return createNewAndKeep( !condition );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::removeObservations(
        const ObservationCondition< ObservationScalarType, TimeType >& condition )
{
    std::vector< std::vector< unsigned int > > indicesToRemoveBySet( getNumberOfObservationSets( ) );
    for( unsigned int observationId = 0; observationId < observationRows_.size( ); ++observationId )
    {
        if( condition( *this, observationId ) )
        {
            const ObservationDatasetRow< TimeType >& row = observationRows_.at( observationId );
            indicesToRemoveBySet.at( row.setId_ ).push_back( row.indexInSet_ );
        }
    }

    for( unsigned int setId = 0; setId < indicesToRemoveBySet.size( ); ++setId )
    {
        if( !indicesToRemoveBySet.at( setId ).empty( ) )
        {
            removeObservationsFromSet( setId, indicesToRemoveBySet.at( setId ) );
        }
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::deleteRejectedObservations( )
{
    removeObservations( ObservationCondition< ObservationScalarType, TimeType >::rejected( ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::rejectObservations(
        const ObservationCondition< ObservationScalarType, TimeType >& condition,
        const std::string& reason )
{
    for( unsigned int observationId = 0; observationId < observationRows_.size( ); ++observationId )
    {
        if( condition( *this, observationId ) )
        {
            observationRows_.at( observationId ).isActive_ = false;
            observationRows_.at( observationId ).rejectionReason_ = reason;
        }
    }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::restoreObservations(
        const ObservationCondition< ObservationScalarType, TimeType >& condition )
{
    for( unsigned int observationId = 0; observationId < observationRows_.size( ); ++observationId )
    {
        if( condition( *this, observationId ) )
        {
            observationRows_.at( observationId ).isActive_ = true;
            observationRows_.at( observationId ).rejectionReason_.clear( );
        }
    }
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETSELECTIONIMPLEMENTATION_H
