/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_WEIGHTING_H
#define TUDAT_TRACKING_DATA_WEIGHTING_H

#include <Eigen/Core>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/basics/tudatTypeTraits.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingDataWeights.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace observation_models
{

Eigen::Vector3d getGroundStationPositionForTrackingDataLinkEnd( const simulation_setup::SystemOfBodies& bodies,
                                                                const LinkEndId& linkEndId );

template< typename ObservationScalarType = double,
          typename TimeType = double,
          typename std::enable_if< is_state_scalar_and_time_type< ObservationScalarType, TimeType >::value, int >::type = 0 >
void setObservationWeightsFromTrackingDataScheme(
        const std::shared_ptr< data::TrackingData< ObservationScalarType, TimeType > > trackingData,
        const simulation_setup::SystemOfBodies& bodies,
        const observation_models::LinkEnds& rawLinkEnds,
        const observation_models::LinkEndType referenceLinkEnd )
{
    if( trackingData->getWeighingScheme( ).empty( ) )
    {
        return;
    }

    if( !trackingData->getObservationWeights( ).empty( ) )
    {
        return;
    }

    if( trackingData->getWeighingScheme( ) != "VFCC17" )
    {
        throw std::runtime_error( "Error when setting observation weights from TrackingData, weighing scheme '" +
                                  trackingData->getWeighingScheme( ) + "' is not recognized." );
    }

    const Eigen::Vector3d stationPosition = getGroundStationPositionForTrackingDataLinkEnd( bodies, rawLinkEnds.at( referenceLinkEnd ) );
    const std::string stationCode = rawLinkEnds.at( referenceLinkEnd ).getReferencePointName( );
    std::vector< Eigen::Vector3d > stationPositions( trackingData->getNumberOfObservations( ), stationPosition );

    std::vector< std::map< std::string, std::string > > stringMetadata( trackingData->getNumberOfObservations( ) );
    for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
    {
        stringMetadata.at( i )[ "observatory" ] = stationCode;
    }

    const std::map< std::string, std::vector< std::string > > ancillaryStringVectorSettings =
            trackingData->getAncillarySettingsStringVector( );

    for( const std::string metadataKey : { "note2", "catalog" } )
    {
        auto metadataIterator = ancillaryStringVectorSettings.find( metadataKey );
        if( metadataIterator == ancillaryStringVectorSettings.end( ) )
        {
            continue;
        }
        if( metadataIterator->second.size( ) != trackingData->getNumberOfObservations( ) )
        {
            throw std::runtime_error( "Error when setting VFCC17 weights from TrackingData, metadata key '" + metadataKey +
                                      "' has an inconsistent number of entries." );
        }
        for( unsigned int i = 0; i < trackingData->getNumberOfObservations( ); ++i )
        {
            stringMetadata.at( i )[ metadataKey ] = metadataIterator->second.at( i );
        }
    }

    data::setVFCC17Weights< ObservationScalarType, TimeType >( trackingData, stationPositions, stringMetadata );
}

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_TRACKING_DATA_WEIGHTING_H
