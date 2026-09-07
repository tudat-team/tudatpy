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
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition ) const
{
    return createNewAndKeep( !condition );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::removeObservations(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    retainObservationRows( getObservationIdsMatchingCondition( !condition ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::removeRejectedObservations( )
{
    removeObservations( ObservationSelectionCondition< ObservationScalarType, TimeType >::rejected( ) );
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::rejectObservations(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition,
        const std::string& reason )
{
    const auto selected = getObservationIdsMatchingCondition( condition );
    for( const unsigned int id : selected )
    {
        auto& row = mutableObservationRow( id );
        row.isActive_ = false;
        if( !reason.empty( ) ) { row.rejectionReason_ = reason; }
    }
    if( !selected.empty( ) ) { ++selectionVersion_; }
}

template< typename ObservationScalarType,
          typename TimeType,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type Dummy >
void ObservationDataset< ObservationScalarType, TimeType, Dummy >::restoreObservations(
        const ObservationSelectionCondition< ObservationScalarType, TimeType >& condition )
{
    const auto selected = getObservationIdsMatchingCondition( condition );
    for( const unsigned int id : selected )
    {
        auto& row = mutableObservationRow( id );
        row.isActive_ = true;
    }
    if( !selected.empty( ) ) { ++selectionVersion_; }
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_OBSERVATIONDATASETSELECTIONIMPLEMENTATION_H
