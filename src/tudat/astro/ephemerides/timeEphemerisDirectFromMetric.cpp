/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/ephemerides/timeEphemerisDirectFromMetric.h"

#include <stdexcept>

#include "tudat/astro/relativity/relativisticTimeConversion.h"

namespace tudat
{

using namespace basic_astrodynamics;

std::function< Time( const Time ) > TimeEphemerisDirectFromMetric::getTimeDifferenceFunctionFromExtendedTime(
        const TimeScales inputScale, const TimeScales outputScale, const std::string& pointIdentifier )
{
    std::function< Time( const Time ) > timeDifferenceFunction;

    if( !isTimeScaleRelativistic( inputScale ) || !isTimeScaleRelativistic( outputScale ) )
    {
        throw std::runtime_error( "Error when getting relatvistic time conversion, scales are not both relativistic" );
    }
    else
    {
        if( inputScale == body_centered_coordinate_time_scale )
        {
            throw std::runtime_error( "Cannot do body-centered input time scale for direct-from-metric time ephemeris" );
        }

        if( outputScale == body_centered_coordinate_time_scale )
        {
            throw std::runtime_error( "Cannot do body-centered output time scale for direct-from-metric time ephemeris" );
        }

        if( inputScale == barycentric_coordinate_time_scale )
        {
            if( outputScale == local_proper_time_scale )
            {
                if( globalCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 &&
                    globalCoordinateToProperTimeInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    throw std::runtime_error(
                                "Error, body-point " + pointIdentifier +
                                " not found in direct-from-metric time ephemeris" );
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        if( globalCoordinateToProperTimeInterpolatorsExtended_.count( pointIdentifier ) > 0 )
                        {
                            return Time( globalCoordinateToProperTimeInterpolatorsExtended_.at( pointIdentifier )->interpolate( currentTime ) );
                        }
                        else
                        {
                            return Time( getGlobalCoordinateToProperTimeDifference( pointIdentifier, currentTime.getSeconds< double >( ) ) );
                        }
                    };
                }
            }
            else
            {
                throw std::runtime_error( "Error A when using direct from metric time ephemeris" );
            }
        }
        else if( inputScale == local_proper_time_scale )
        {
            if( outputScale == barycentric_coordinate_time_scale )
            {
                if( properTimeToGlobalCoordinateInterpolators_.count( pointIdentifier ) == 0 &&
                    properTimeToGlobalCoordinateInterpolatorsExtended_.count( pointIdentifier ) == 0 )
                {
                    throw std::runtime_error(
                                "Error, body-point " + pointIdentifier +
                                " not found in direct-from-metric time ephemeris" );
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        if( properTimeToGlobalCoordinateInterpolatorsExtended_.count( pointIdentifier ) > 0 )
                        {
                            return Time( properTimeToGlobalCoordinateInterpolatorsExtended_.at( pointIdentifier )->interpolate( currentTime ) );
                        }
                        else
                        {
                            return Time( getProperToGlobalCoordinateTimeDifference( pointIdentifier, currentTime.getSeconds< double >( ) ) );
                        }
                    };
                }
            }
            else
            {
                throw std::runtime_error( "Error B when using direct from metric time ephemeris" );
            }
        }
    }

    return timeDifferenceFunction;
}

} // namespace tudat
