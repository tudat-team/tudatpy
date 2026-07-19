/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and binary
 *    forms, with or without modification, are permitted exclusively under the
 *    terms of the Modified BSD license.
 */

#ifndef TUDAT_CREATE_OBSERVATION_COLLECTION_H
#define TUDAT_CREATE_OBSERVATION_COLLECTION_H

#include "tudat/simulation/estimation_setup/createObservationDataset.h"
#include "tudat/simulation/estimation_setup/observationCollection.h"

namespace tudat
{
namespace observation_models
{

//! Backwards-compatible adapter for the legacy collection representation.
template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
std::shared_ptr< ObservationCollection< ObservationScalarType, TimeType > > createObservationCollection(
        const std::vector< std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > > trackingDataList,
        simulation_setup::SystemOfBodies& bodies,
        const bool applyCorrections = false )
{
    return observation_models::createObservationCollection< ObservationScalarType, TimeType >(
            createObservationDatasetFromTrackingData< ObservationScalarType, TimeType >( trackingDataList, bodies, applyCorrections ) );
}

}  // namespace observation_models
}  // namespace tudat

#endif  // TUDAT_CREATE_OBSERVATION_COLLECTION_H
