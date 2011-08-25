/*! \file unitConversions.h
 *    This file contains a namespace with selected unit conversions
 *    commonly used in astrodynamics.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 6
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
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : J. Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 6 September, 2010
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
 *      The behaviour of the template conversion functions has not been tested
 *      for integer data types.
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
 *      110808    J. Melman         Added time conversion templates.
 *      110809    K. Kumar          Minor comment changes; added note over ints.
 *      110810    K. Kumar          Doxygen comment corrections.
 */

#ifndef UNIT_CONVERSIONS_H
#define UNIT_CONVERSIONS_H

// Include statements.
#include <cmath>
#include "physicalConstants.h"

//! Unit conversions namespace.
/*!
 *  Unit conversions namespace.
 */
namespace unit_conversions
{

//! Convert angle in radians to degrees.
/*!
 * Converts angle given in radians to degrees.
 * \param angleInRadians Angle in radians.
 * \return Angle in degrees.
 */
template < typename T >
T convertRadiansToDegrees( T angleInRadians )
{
    return angleInRadians / M_PI * 180.0;
}

//! Convert angle in degrees to radians.
/*!
 * Converts angle given in degrees to radians.
 * \param angleInDegrees Angle in degrees.
 * \return Angle in radians.
 */
template < typename T >
T convertDegreesToRadians( T angleInDegrees )
{
    return angleInDegrees / 180.0 * M_PI;
}

//! Convert angle in degrees to arcminutes.
/*!
 * Converts angle given in degrees to arcminutes.
 * \param angleInDegrees Angle in degrees.
 * \return Angle in arcminutes.
 */
template < typename T >
T convertDegreesToArcminutes( T angleInDegrees )
{
    return angleInDegrees * 60.0;
}

//! Convert angle in arcminutes to arcseconds.
/*!
 * Converts angle given in arcminutes to arcseconds.
 * \param angleInArcminutes Angle in arcminutes.
 * \return Angle in arcseconds.
 */
template < typename T >
T convertArcminutesToArcseconds( T angleInArcminutes )
{
    return angleInArcminutes * 60.0;
}

//! Convert distance in meters to kilometers.
/*!
 * Converts distance given in meters to kilometers.
 * \param distanceInMeters Distance in meters.
 * \return Distance in kilometers.
 */
template < typename T >
T convertMetersToKilometers( T distanceInMeters )
{
    return distanceInMeters / 1000.0;
}

//! Convert distance in kilometers to meters.
/*!
 * Converts distance given in kilometers to meters.
 * \param distanceInKilometers Distance in kilometers.
 * \return Distance in meters.
 */
template < typename T >
T convertKilometersToMeters( T distanceInKilometers )
{
    return distanceInKilometers * 1000.0;
}

//! Convert distance in meters to astronomical units.
/*!
 * Converts distance given in meters to astronomical units.
 * \param distanceInMeters Distance in meters.
 * \return Distance in astronomical units.
 */
template < typename T >
T convertMetersToAstronomicalUnits( T distanceInMeters )
{
    return distanceInMeters / PhysicalConstants::ASTRONOMICAL_UNIT;
}

//! Convert distance in astronomical units to meters.
/*!
 * Converts distance given in astronomical units to meters.
 * \param distanceInAstronomicalUnits Distance in astronomical units.
 * \return Distance in meters.
 */
template < typename T >
T convertAstronomicalUnitsToMeters( T distanceInAstronomicalUnits )
{
    return distanceInAstronomicalUnits * PhysicalConstants::ASTRONOMICAL_UNIT;
}

//! Convert time in seconds to minutes.
/*!
 * Converts time given in seconds to minutes.
 * \param timeInSeconds Time in seconds.
 * \return Time in minutes.
 */
template < typename T >
T convertSecondsToMinutes( T timeInSeconds )
{
    return timeInSeconds / 60.0;
}

//! Convert time in minutes to seconds.
/*!
 * Converts time given in minutes to seconds.
 * \param timeInMinutes Time in minutes.
 * \return Time in seconds.
 */
template < typename T >
T convertMinutesToSeconds( T timeInMinutes )
{
    return timeInMinutes * 60.0;
}

//! Convert time in seconds to hours.
/*!
 * Converts time given in seconds to hours.
 * \param timeInSeconds Time in seconds.
 * \return Time in hours.
 */
template < typename T >
T convertSecondsToHours( T timeInSeconds )
{
    return timeInSeconds / 3600.0;
}

//! Convert time in hours to seconds.
/*!
 * Converts time given in hours to seconds.
 * \param timeInHours Time in hours.
 * \return Time in seconds.
 */
template < typename T >
T convertHoursToSeconds( T timeInHours )
{
    return timeInHours * 3600.0;
}

//! Convert time in seconds to Julian days.
/*!
 * Converts time given in seconds to Julian days.
 * \param timeInSeconds Time in seconds.
 * \return Time in Julian days.
 */
template < typename T >
T convertSecondsToJulianDays( T timeInSeconds )
{
    return timeInSeconds / PhysicalConstants::JULIAN_DAY;
}

//! Convert time in Julian days to seconds.
/*!
 * Converts time given in Julian days to seconds.
 * \param timeInJulianDays Time in Julian days.
 * \return Time in seconds.
 */
template < typename T >
T convertJulianDaysToSeconds( T timeInJulianDays )
{
    return timeInJulianDays * PhysicalConstants::JULIAN_DAY;
}

//! Convert time in seconds to sidereal days.
/*!
 * Converts time given in seconds to sidereal days.
 * \param timeInSeconds Time in seconds.
 * \return Time in sidereal days.
 */
template < typename T >
T convertSecondsToSiderealDays( T timeInSeconds )
{
    return timeInSeconds / PhysicalConstants::SIDEREAL_DAY;
}

//! Convert time in sidereal days to seconds.
/*!
 * Converts time given in sidereal days to seconds.
 * \param timeInSiderealDays Time in sidereal days.
 * \return Time in seconds.
 */
template < typename T >
T convertSiderealDaysToSeconds( T timeInSiderealDays )
{
    return timeInSiderealDays * PhysicalConstants::SIDEREAL_DAY;
}

//! Convert time in Julian days to Julian years.
/*!
 * Converts time given in Julian days to Julian years.
 * \param timeInJulianDays Time in Julian days.
 * \return Time in Julian years.
 */
template < typename T >
T convertJulianDaysToJulianYears( T timeInJulianDays )
{
    return timeInJulianDays / PhysicalConstants::JULIAN_YEAR_IN_DAYS;
}

//! Convert time in Julian years to Julian days.
/*!
 * Converts time given in Julian years to Julian days.
 * \param timeInJulianYears Time in Julian years.
 * \return Time in Julian days.
 */
template < typename T >
T convertJulianYearsToJulianDays( T timeInJulianYears )
{
    return timeInJulianYears * PhysicalConstants::JULIAN_YEAR_IN_DAYS;
}

//! Convert temperature in Rankine to Kelvin.
/*!
 * Converts temperature given in Rankine to Kelvin.
 * \param temperatureInRankine Temperature in Rankine.
 * \return Temperature in Kelvin.
 */
template < typename T >
T convertRankineToKelvin( T temperatureInRankine )
{
   return temperatureInRankine * 5.0 / 9.0;
}

//! Convert distance in feet to meters.
/*!
 * Converts distance given in feet to meters.
 * \param distanceInFeet Distance in feet.
 * \return Distance in meters.
 */
template < typename T >
T convertFeetToMeter( T distanceInFeet )
{
    return distanceInFeet * 0.3048;
}

//! Convert pressure in pound per square feet to Newton per square meter.
/*!
 * Converts pressure given in pound per square feet to Newton per square meter (Pascal).
 * \param pressureInPoundPerSquareFeet Pressure in pound per square feet.
 * \return Pressure in Newton per square meter (Pascal).
 */
template < typename T >
T convertPoundPerSquareFeetToPascal( T pressureInPoundPerSquareFeet )
{
    return pressureInPoundPerSquareFeet * 47.880259;
}

}

#endif // UNIT_CONVERSIONS_H

// End of file.
