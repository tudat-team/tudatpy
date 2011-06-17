/*! \file unitConversions.h
 *    This file contains a namespace with selected unit conversions
 *    commonly used in astrodynamics.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 3
 *    Check status      : Checked
 *
 *    Author            : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Author            : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : D.Dirkx@student.tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 6 September, 2010
 *    Last modified     : 11 April, 2011
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
 *      YYMMDD    Author            Comment
 *      100906    J. Melman         First creation of code.
 *      110411    K. Kumar          Added convertDegreesToArcminutes() and
 *                                  convertArcminutesToArcseconds().
 *      110609    F.M. Engelen      Added Rankine, feet, and pound per square feet
 *                                  to SI conversion.
 */

#ifndef UNIT_CONVERSIONS_H
#define UNIT_CONVERSIONS_H

#include <cmath>
#include "physicalConstants.h"

//! Unit conversions namespace.
/*!
 *  Unit conversions namespace.
 */
namespace unit_conversions
{

//! Convert from radians to degrees.
/*!
 *  Convert an angle from radians to degrees.
 *  \param angleInRadians Angle in radians.
 *  \return Angle in degrees.
 */
template < typename T >
T convertRadiansToDegrees( T angleInRadians )
{
    return angleInRadians / M_PI * 180.0;
}

//! Convert from degrees to radians.
/*!
 *  Convert an angle from degrees to radians.
 *  \param angleInDegrees Angle in degrees.
 *  \return Angle in radians.
 */
template < typename T >
T convertDegreesToRadians( T angleInDegrees )
{
    return angleInDegrees / 180.0 * M_PI;
}

//! Convert from degrees to arcminutes.
/*!
 *  Convert an angle from degrees to arcminutes.
 *  \param angleInDegrees Angle in degrees.
 *  \return Angle in arcminutes.
 */
template < typename T >
T convertDegreesToArcminutes( T angleInDegrees )
{
    return angleInDegrees * 60.0;
}

//! Convert from arcminutes to arcseconds.
/*!
 *  Convert an angle from arcminutes to arcseconds.
 *  \param angleInArcminutes Angle in arcminutes.
 *  \return Angle in arcseconds.
 */
template < typename T >
T convertArcminutesToArcseconds( T angleInArcminutes )
{
    return angleInArcminutes * 60.0;
}

//! Convert from meters to kilometers.
/*!
 *  Convert a distance from meters to kilometers.
 *  \param distanceInMeters Distance in meters.
 *  \return Distance in kilometers.
 */
template < typename T >
T convertMetersToKilometers( T distanceInMeters )
{
    return distanceInMeters / 1000.;
}

//! Convert from kilometers to meters.
/*!
 *  Convert a distance from kilometers to meters.
 *  \param distanceInKilometers Distance in kilometers.
 *  \return Distance in meters.
 */
template < typename T >
T convertKilometersToMeters( T distanceInKilometers )
{
    return distanceInKilometers * 1000.;
}

//! Convert from meters to astronomical units.
/*!
 *  Convert a distance from meters to astronomical units.
 *  \param distanceInMeters Distance in meters.
 *  \return Distance in astronomical units.
 */
template < typename T >
T convertMetersToAstronomicalUnits( T distanceInMeters )
{
    return distanceInMeters / PhysicalConstants::ASTRONOMICAL_UNIT;
}

//! Convert from astronomical units to meters.
/*!
 *  Convert a distance from astronomical units to meters.
 *  \param distanceInAstronomicalUnits Distance in astronomical units.
 *  \return Distance in meters.
 */
template < typename T >
T convertAstronomicalUnitsToMeters( T distanceInAstronomicalUnits )
{
    return distanceInAstronomicalUnits * PhysicalConstants::ASTRONOMICAL_UNIT;
}

//! Convert from Rankine to Kelvin.
/*!
 *  Convert a temperature in Rankine to Kelvin.
 *  \param temperatureInRankine Temperature in Rankine.
 *  \return Temperature in Kelvin.
 */
template < typename T >
T convertRankineToKelvin( T temperatureInRankine )
{
   return temperatureInRankine * 5.0 / 9.0;
}

//! Convert from feet to meters.
/*!
 *  Convert a distance in feet to meters.
 *  \param distanceInFeet Distance in feet.
 *  \return distance in meters.
 */
template < typename T >
T convertFeetToMeter( T distanceInFeet )
{
    return distanceInFeet * 0.3048 ;
}

//! Convert from pound per square feet to Newton per square meter (Pascal)
/*!
 * Convert a pressure in pound per square feet to Newton per square meter.
 * \param pressureInPoundPerSquareFeet pressure in pound per square feet.
 * \return pressure in Newton per square meter (Pascal)
 */
template < typename T >
T convertPoundPerSquareFeetToPascal( T pressureInPoundPerSquareFeet )
{
    return pressureInPoundPerSquareFeet * 47.880259;
}

}

#endif // UNIT_CONVERSIONS_H

// End of file.
