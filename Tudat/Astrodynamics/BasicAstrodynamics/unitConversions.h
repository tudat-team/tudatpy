/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      100906    J. Melman         First creation of code.
 *      110411    K. Kumar          Added convertDegreesToArcminutes() and
 *                                  convertArcminutesToArcseconds().
 *      110609    F.M. Engelen      Added Rankine, feet, pound-per-square-foot to SI conversion.
 *      110808    J. Melman         Added time conversion templates.
 *      110809    K. Kumar          Minor comment changes; added note over ints.
 *      110810    K. Kumar          Doxygen comment corrections.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120127    K. Kumar          Added references in preparation for commit to Tudat Core.
 *      120127    D. Dirkx          Minor textual changes and made consistent with rewrite of
 *                                  physical constants.
 *      121212    K. Kumar          Migrated namespace to directory-based protocol and added
 *                                  backwards compatibility.
 *
 *    References
 *      Wikipedia. http://en.wikipedia.org/wiki/Temperature_conversion_formulas,
 *          last accessed: 27, January 2012(a).
 *      Wikipedia. http://en.wikipedia.org/wiki/Conversion_of_units#Length,
 *          last accessed: 27, January 2012(b).
 *       Wikipedia. http://en.wikipedia.org/wiki/Conversion_of_units#Pressure_or_mechanical_stress,
 *          last accessed: 27, January 2012(c).
 *
 *    Notes
 *      The behaviour of the template conversion functions has not been tested for integer
 *      datatypes.
 *
 *      Backwards compatibility of namespaces is implemented for Tudat Core 2 in this file. The
 *      code block marked "DEPRECATED!" at the end of the file should be removed in Tudat Core 3.
 *
 */

#ifndef TUDAT_UNIT_CONVERSIONS_H
#define TUDAT_UNIT_CONVERSIONS_H

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

namespace tudat
{

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
{ return angleInRadians / mathematical_constants::PI * 180.0; }

//! Convert angle in degrees to radians.
/*!
 * Converts angle given in degrees to radians.
 * \param angleInDegrees Angle in degrees.
 * \return Angle in radians.
 */
template < typename T >
T convertDegreesToRadians( T angleInDegrees )
{ return angleInDegrees / 180.0 * mathematical_constants::PI; }

//! Convert angle in degrees to arcminutes.
/*!
 * Converts angle given in degrees to arcminutes.
 * \param angleInDegrees Angle in degrees.
 * \return Angle in arcminutes.
 */
template < typename T >
T convertDegreesToArcminutes( T angleInDegrees )
{ return angleInDegrees * 60.0; }

//! Convert angle in arcminutes to arcseconds.
/*!
 * Converts angle given in arcminutes to arcseconds.
 * \param angleInArcminutes Angle in arcminutes.
 * \return Angle in arcseconds.
 */
template < typename T >
T convertArcminutesToArcseconds( T angleInArcminutes )
{ return angleInArcminutes * 60.0; }

//! Convert distance in meters to kilometers.
/*!
 * Converts distance given in meters to kilometers.
 * \param distanceInMeters Distance in meters.
 * \return Distance in kilometers.
 */
template < typename T >
T convertMetersToKilometers( T distanceInMeters )
{ return distanceInMeters / 1000.0; }

//! Convert distance in kilometers to meters.
/*!
 * Converts distance given in kilometers to meters.
 * \param distanceInKilometers Distance in kilometers.
 * \return Distance in meters.
 */
template < typename T >
T convertKilometersToMeters( T distanceInKilometers )
{ return distanceInKilometers * 1000.0; }

//! Convert distance in meters to astronomical units.
/*!
 * Converts distance given in meters to astronomical units.
 * \param distanceInMeters Distance in meters.
 * \return Distance in astronomical units.
 */
template < typename T >
T convertMetersToAstronomicalUnits( T distanceInMeters )
{ return distanceInMeters / physical_constants::ASTRONOMICAL_UNIT; }

//! Convert distance in astronomical units to meters.
/*!
 * Converts distance given in astronomical units to meters.
 * \param distanceInAstronomicalUnits Distance in astronomical units.
 * \return Distance in meters.
 */
template < typename T >
T convertAstronomicalUnitsToMeters( T distanceInAstronomicalUnits )
{ return distanceInAstronomicalUnits * physical_constants::ASTRONOMICAL_UNIT; }

//! Convert time in seconds to minutes.
/*!
 * Converts time given in seconds to minutes.
 * \param timeInSeconds Time in seconds.
 * \return Time in minutes.
 */
template < typename T >
T convertSecondsToMinutes( T timeInSeconds )
{ return timeInSeconds / 60.0; }

//! Convert time in minutes to seconds.
/*!
 * Converts time given in minutes to seconds.
 * \param timeInMinutes Time in minutes.
 * \return Time in seconds.
 */
template < typename T >
T convertMinutesToSeconds( T timeInMinutes )
{ return timeInMinutes * 60.0; }

//! Convert time in seconds to hours.
/*!
 * Converts time given in seconds to hours.
 * \param timeInSeconds Time in seconds.
 * \return Time in hours.
 */
template < typename T >
T convertSecondsToHours( T timeInSeconds )
{ return timeInSeconds / 3600.0; }

//! Convert time in hours to seconds.
/*!
 * Converts time given in hours to seconds.
 * \param timeInHours Time in hours.
 * \return Time in seconds.
 */
template < typename T >
T convertHoursToSeconds( T timeInHours )
{ return timeInHours * 3600.0; }

//! Convert time in seconds to Julian days.
/*!
 * Converts time given in seconds to Julian days.
 * \param timeInSeconds Time in seconds.
 * \return Time in Julian days.
 */
template < typename T >
T convertSecondsToJulianDays( T timeInSeconds )
{ return timeInSeconds / physical_constants::JULIAN_DAY; }

//! Convert time in Julian days to seconds.
/*!
 * Converts time given in Julian days to seconds.
 * \param timeInJulianDays Time in Julian days.
 * \return Time in seconds.
 */
template < typename T >
T convertJulianDaysToSeconds( T timeInJulianDays )
{ return timeInJulianDays * physical_constants::JULIAN_DAY; }

//! Convert time in seconds to sidereal days.
/*!
 * Converts time given in seconds to sidereal days.
 * \param timeInSeconds Time in seconds.
 * \return Time in sidereal days.
 */
template < typename T >
T convertSecondsToSiderealDays( T timeInSeconds )
{ return timeInSeconds / physical_constants::SIDEREAL_DAY; }

//! Convert time in sidereal days to seconds.
/*!
 * Converts time given in sidereal days to seconds.
 * \param timeInSiderealDays Time in sidereal days.
 * \return Time in seconds.
 */
template < typename T >
T convertSiderealDaysToSeconds( T timeInSiderealDays )
{ return timeInSiderealDays * physical_constants::SIDEREAL_DAY; }

//! Convert time in Julian days to Julian years.
/*!
 * Converts time given in Julian days to Julian years.
 * \param timeInJulianDays Time in Julian days.
 * \return Time in Julian years.
 */
template < typename T >
T convertJulianDaysToJulianYears( T timeInJulianDays )
{ return timeInJulianDays / physical_constants::JULIAN_YEAR_IN_DAYS; }

//! Convert time in Julian years to Julian days.
/*!
 * Converts time given in Julian years to Julian days.
 * \param timeInJulianYears Time in Julian years.
 * \return Time in Julian days.
 */
template < typename T >
T convertJulianYearsToJulianDays( T timeInJulianYears )
{ return timeInJulianYears * physical_constants::JULIAN_YEAR_IN_DAYS; }

//! Convert temperature in Rankine to Kelvin.
/*!
 * Converts temperature given in Rankine to Kelvin (Wikipedia, 2012a).
 * \param temperatureInRankine Temperature in Rankine.
 * \return Temperature in Kelvin.
 */
template < typename T >
T convertRankineToKelvin( T temperatureInRankine )
{ return temperatureInRankine * 5.0 / 9.0; }

//! Convert distance in feet to meters.
/*!
 * Converts distance given in feet to meters (Wikipedia, 2012b).
 * \param distanceInFeet Distance in feet.
 * \return Distance in meters.
 */
template < typename T >
T convertFeetToMeter( T distanceInFeet )
{ return distanceInFeet * 0.3048; }

//! Convert pressure in pound per square feet to Newton per square meter.
/*!
 * Converts pressure given in pound per square feet to Newton per square meter (Pascal)
 * (Wikipedia, 2012c).
 * \param pressureInPoundPerSquareFeet Pressure in pound per square feet.
 * \return Pressure in Newton per square meter (Pascal).
 */
template < typename T >
T convertPoundPerSquareFeetToPascal( T pressureInPoundPerSquareFeet )
{ return pressureInPoundPerSquareFeet * 47.880259; }

} // namespace unit_conversions

} // namespace tudat

#endif // TUDAT_UNIT_CONVERSIONS_H
