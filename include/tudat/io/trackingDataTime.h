/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved.
 *    This file is part of Tudat. Redistribution and use in source and binary
 *    forms are permitted under the terms of the Modified BSD license.
 *    See https://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_TRACKING_DATA_TIME_H
#define TUDAT_TRACKING_DATA_TIME_H

#include <algorithm>
#include <optional>

#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"
#include "tudat/io/trackingData.h"

namespace tudat
{
namespace data
{

//! Return the earliest and latest tracking epochs in the requested time scale.
//! Conversions use the geocentre and require no dynamics or estimation environment.
template< typename ObservationScalarType = double, typename TimeType = double >
std::pair< TimeType, TimeType > getTrackingDataEpochBounds(
        const std::vector< std::shared_ptr< TrackingData< ObservationScalarType, TimeType > > >& trackingData,
        const basic_astrodynamics::TimeScales outputTimeScale = basic_astrodynamics::tdb_scale )
{
    if( outputTimeScale < basic_astrodynamics::tai_scale || outputTimeScale > basic_astrodynamics::ut1_scale )
    {
        throw std::invalid_argument( "Tracking epoch bounds support TAI, TT, TDB, UTC and UT1 output scales." );
    }
    std::optional< std::pair< TimeType, TimeType > > bounds;
    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > converter;
    for( const auto& data : trackingData )
    {
        if( data == nullptr )
        {
            throw std::invalid_argument( "Cannot retrieve epoch bounds from null tracking data." );
        }
        const auto& epochs = data->getObservationEpochs( );
        if( epochs.empty( ) )
        {
            continue;
        }
        const auto inputTimeScale = basic_astrodynamics::timeScaleFromString( data->getTimeScale( ) );
        const auto limits = std::minmax_element( epochs.begin( ), epochs.end( ) );
        TimeType first = *limits.first;
        TimeType last = *limits.second;
        if( inputTimeScale != outputTimeScale )
        {
            if( converter == nullptr )
            {
                converter = earth_orientation::createDefaultTimeConverter( );
            }
            first = converter->getCurrentTime< TimeType >( inputTimeScale, outputTimeScale, first );
            last = converter->getCurrentTime< TimeType >( inputTimeScale, outputTimeScale, last );
        }
        if( !bounds )
        {
            bounds = std::make_pair( first, last );
        }
        else
        {
            bounds->first = std::min( bounds->first, first );
            bounds->second = std::max( bounds->second, last );
        }
    }
    if( !bounds )
    {
        throw std::invalid_argument( "Cannot retrieve epoch bounds: tracking data contain no observations." );
    }
    return *bounds;
}

}  // namespace data
}  // namespace tudat

#endif
