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

#ifndef TUDAT_TIME_EPHEMERIS_DIRECT_FROM_METRIC_H
#define TUDAT_TIME_EPHEMERIS_DIRECT_FROM_METRIC_H

#include <functional>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>

#include "tudat/astro/ephemerides/timeEphemeris.h"
#include "tudat/math/interpolators/oneDimensionalInterpolator.h"

namespace tudat
{

//! Time ephemeris using direct metric-based coordinate-time/proper-time interpolation.
template< typename TimeType = double, typename StateScalarType = double >
class TimeEphemerisDirectFromMetric : public TimeEphemeris
{
public:
    using TimeDifferenceInterpolator = std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, StateScalarType > >;

    TimeEphemerisDirectFromMetric(
            const std::string& centralBodyName,
            const std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ),
            const std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators =
            ( std::map< std::string, TimeDifferenceInterpolator >( ) ) ) :
        TimeEphemeris( centralBodyName ),
        globalCoordinateToProperTimeInterpolators_( globalCoordinateToProperTimeInterpolators ),
        properTimeToGlobalCoordinateInterpolators_( properTimeToGlobalCoordinateInterpolators )
    { }

    void resetGlobalToProperTimeInterpolators(
            const TimeDifferenceInterpolator globalCoordinateToProperTimeInterpolator,
            const TimeDifferenceInterpolator properTimeToGlobalCoordinateInterpolator,
            const std::string& referencePoint )
    {
        globalCoordinateToProperTimeInterpolators_[ referencePoint ] = globalCoordinateToProperTimeInterpolator;
        properTimeToGlobalCoordinateInterpolators_[ referencePoint ] = properTimeToGlobalCoordinateInterpolator;
    }

    std::function< Time( const Time ) > getTimeDifferenceFunctionFromExtendedTime(
            const basic_astrodynamics::TimeScales inputScale,
            const basic_astrodynamics::TimeScales outputScale,
            const std::string& pointIdentifier = "" ) override;

    StateScalarType getGlobalCoordinateToProperTimeDifference( const std::string& pointId, const TimeType evaluationTime )
    {
        return globalCoordinateToProperTimeInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

    StateScalarType getProperToGlobalCoordinateTimeDifference( const std::string& pointId, const TimeType evaluationTime )
    {
        return properTimeToGlobalCoordinateInterpolators_.at( pointId )->interpolate( evaluationTime );
    }

protected:
    std::map< std::string, TimeDifferenceInterpolator > globalCoordinateToProperTimeInterpolators_;
    std::map< std::string, TimeDifferenceInterpolator > properTimeToGlobalCoordinateInterpolators_;
};

template< typename TimeType, typename StateScalarType >
std::function< Time( const Time ) > TimeEphemerisDirectFromMetric< TimeType, StateScalarType >::getTimeDifferenceFunctionFromExtendedTime(
        const basic_astrodynamics::TimeScales inputScale,
        const basic_astrodynamics::TimeScales outputScale,
        const std::string& pointIdentifier )
{
    std::function< Time( const Time ) > timeDifferenceFunction;

    if( !basic_astrodynamics::isTimeScaleRelativistic( inputScale ) ||
            !basic_astrodynamics::isTimeScaleRelativistic( outputScale ) )
    {
        throw std::runtime_error( "Error when getting relatvistic time conversion, scales are not both relativistic" );
    }
    else
    {
        if( inputScale == basic_astrodynamics::body_centered_coordinate_time_scale )
        {
            throw std::runtime_error( "Cannot do body-centered input time scale for direct-from-metric time ephemeris" );
        }

        if( outputScale == basic_astrodynamics::body_centered_coordinate_time_scale )
        {
            throw std::runtime_error( "Cannot do body-centered output time scale for direct-from-metric time ephemeris" );
        }

        if( inputScale == basic_astrodynamics::barycentric_coordinate_time_scale )
        {
            if( outputScale == basic_astrodynamics::local_proper_time_scale )
            {
                if( globalCoordinateToProperTimeInterpolators_.count( pointIdentifier ) == 0 )
                {
                    throw std::runtime_error(
                                "Error, body-point " + pointIdentifier +
                                " not found in direct-from-metric time ephemeris" );
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        const TimeType interpolationTime =
                                detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime );
                        return Time( static_cast< long double >(
                                         getGlobalCoordinateToProperTimeDifference( pointIdentifier, interpolationTime ) ) );
                    };
                }
            }
            else
            {
                throw std::runtime_error( "Error A when using direct from metric time ephemeris" );
            }
        }
        else if( inputScale == basic_astrodynamics::local_proper_time_scale )
        {
            if( outputScale == basic_astrodynamics::barycentric_coordinate_time_scale )
            {
                if( properTimeToGlobalCoordinateInterpolators_.count( pointIdentifier ) == 0 )
                {
                    throw std::runtime_error(
                                "Error, body-point " + pointIdentifier +
                                " not found in direct-from-metric time ephemeris" );
                }
                else
                {
                    timeDifferenceFunction = [=]( const Time currentTime )
                    {
                        const TimeType interpolationTime =
                                detail::convertInterpolatorTimeFromExtendedTime< TimeType >( currentTime );
                        return Time( static_cast< long double >(
                                         getProperToGlobalCoordinateTimeDifference( pointIdentifier, interpolationTime ) ) );
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

}  // namespace tudat

#endif  // TUDAT_TIME_EPHEMERIS_DIRECT_FROM_METRIC_H
