/*! \file unitConversions.h
 *    This file contains a namespace with selected unit conversions
 *    commonly used in astrodynamics.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 1
 *    Check status      : Unchecked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Date created      : 6 September, 2010
 *    Last modified     : 6 September, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    author        comment
 *      100906    J. Melman     First creation of code.
 */

#ifndef UNIT_CONVERSIONS_H
#define UNIT_CONVERSIONS_H

#include <cmath>
#include "physicalConstants.h"

namespace unit_conversions
{

//! Convert from radians to degrees.
/*!
 *  Convert an angle from radians to degrees.
 *  \param angleInRadians angle in radians.
 *  \return angleInDegrees angle in degrees.
 */
template < typename T >
T convertRadiansToDegrees( T angleInRadians )
{
    return angleInRadians / M_PI * 180.0;
}

//! Convert from degrees to radians.
/*!
 *  Convert an angle from degrees to radians.
 *  \param angleInDegrees angle in degrees.
 *  \return angleInRadians angle in radians.
 */
template < typename T >
T convertDegreesToRadians( T angleInDegrees )
{
    return angleInDegrees / 180.0 * M_PI;
}

//! Convert from meters to kilometers.
/*!
 *  Convert a distance from meters to kilometers.
 *  \param distanceInMeters distance in meters.
 *  \return distanceInKilometers distance in kilometers.
 */
template < typename T >
T convertMetersToKilometers( T distanceInMeters )
{
    return distanceInMeters / 1000.;
}

//! Convert from kilometers to meters.
/*!
 *  Convert a distance from kilometers to meters.
 *  \param distanceInKilometers distance in kilometers.
 *  \return distanceInMeters distance in meters.
 */
template < typename T >
T convertKilometersToMeters( T distanceInKilometers )
{
    return distanceInKilometers * 1000.;
}

//! Convert from meters to astronomical units.
/*!
 *  Convert a distance from meters to astronomical units.
 *  \param distanceInMeters distance in meters.
 *  \return distanceInAstronomicalUnits distance in astronomical units.
 */
template < typename T >
T convertMetersToAstronomicalUnits( T distanceInMeters )
{
    return distanceInMeters / PhysicalConstants::ASTRONOMICAL_UNIT;
}

//! Convert from astronomical units to meters.
/*!
 *  Convert a distance from astronomical units to meters.
 *  \param distanceInAstronomicalUnits distance in astronomical units.
 *  \return distanceInMeters distance in meters.
 */
template < typename T >
T convertAstronomicalUnitsToMeters( T distanceInAstronomicalUnits )
{
    return distanceInAstronomicalUnits * PhysicalConstants::ASTRONOMICAL_UNIT;
}

}

#endif // UNIT_CONVERSIONS_H

// End of file.
