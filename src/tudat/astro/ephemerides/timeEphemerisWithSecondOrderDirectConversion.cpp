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

#include "tudat/astro/ephemerides/timeEphemerisWithSecondOrderDirectConversion.h"

#include <stdexcept>

namespace tudat
{

TimeEphemerisWithSecondOrderDirectConversion::TimeEphemerisWithSecondOrderDirectConversion(
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

double TimeEphemerisWithSecondOrderDirectConversion::calculateDirectTimeDifferenceTerm(
        const Eigen::Vector3d positionVectorFromReferencePoint,
        const double currentTime )
{
    throw std::runtime_error( "Error, second order time ephemeris direct difference term not yet implemented" );
    return 0.0;
}

} // namespace tudat
