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

#include "tudat/astro/ephemerides/timeEphemerisWithFirstOrderDirectConversion.h"

#include "tudat/astro/relativity/relativisticTimeConversion.h"

namespace tudat
{

TimeEphemerisWithFirstOrderDirectConversion::TimeEphemerisWithFirstOrderDirectConversion(
        const TimeDifferenceInterpolator barycenterToPlanetCenterCoordinateTimeInterpolator,
        const TimeDifferenceInterpolator planetCenterToBarycenterCoordinateTimeInterpolator,
        const std::string& centralBodyName,
        const std::function< Eigen::Vector6d( const double ) > centralBodyStateFunction,
        const std::map< std::string, std::function< Eigen::Vector3d( const double ) > >& groundStationPositionFunctions,
        const std::map< std::string, TimeDifferenceInterpolator > planetCoordinateToProperTimeInterpolators,
        const std::map< std::string, TimeDifferenceInterpolator > properTimeToPlanetCoordinateInterpolators ) :
    TimeEphemerisFromPostNewtonianExpansion(
        barycenterToPlanetCenterCoordinateTimeInterpolator,
        planetCenterToBarycenterCoordinateTimeInterpolator,
        centralBodyName,
        groundStationPositionFunctions,
        planetCoordinateToProperTimeInterpolators,
        properTimeToPlanetCoordinateInterpolators ),
    centralBodyStateFunction_( centralBodyStateFunction )
{ }

double TimeEphemerisWithFirstOrderDirectConversion::calculateDirectTimeDifferenceTerm(
        const Eigen::Vector3d positionVectorFromReferencePoint,
        const double currentTime )
{
    return relativity::calculateFirstOrderTcbToTcgDirectCorrection(
                positionVectorFromReferencePoint, centralBodyStateFunction_( currentTime ).segment( 3, 3 ) );
}

} // namespace tudat
