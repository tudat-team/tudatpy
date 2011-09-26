/*! \file unitTestUnitConversions.cpp
 *    This unit test will test the unit conversions that are
 *    defined in unitConversions.h.
 *
 *    Path              : /Astrodynamics/
 *    Version           : 7
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
 *    Date created      : 10 September, 2010
 *    Last modified     : 9 August, 2011
 *
 *    References
 *
 *    Notes
 *      At the moment, not all conversion routines are test both ways. This
 *      should be modified in a next version.
 *
 *      Test runs code and verifies result against expected value.
 *      If the tested code is erroneous, the test function returns a boolean
 *      true; if the code is correct, the function returns a boolean false.
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110111    J. Melman         First creation of code.
 *      110124    J. Melman         Adapted to the offical Tudat standards.
 *      110411    K. Kumar          Added unit tests for
 *                                  convertDegreesToArcminutes() and
 *                                  convertArcminutesToArcseconds().
 *      110615    F.M. Engelen      Added Rankine, feet, and pound/m^2
 *                                  conversion. Solved bug with respect to
 *                                  absolute precision and relative precision.
 *      110808    J. Melman         Added time conversion unit tests.
 *      110809    K. Kumar          Minor corrections.
 *      110810    J. Melman         Added 3 more time conversion unit tests.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 */

// Include statements.
#include <cmath>
#include "Mathematics/unitTestUnitConversions.h"

// Using directives.
using std::cerr;
using std::endl;
using std::fabs;
using mathematics::MACHINE_PRECISION_DOUBLES;

//! Namespace for all unit tests.
namespace unit_tests
{

//! Test unit conversions.
bool testUnitConversions( )
{
    // Test result initialised to false.
    bool isUnitConversionsErroneous = false;

    // Test conversion from kilometers to meters.
    if ( fabs( unit_conversions::convertKilometersToMeters( 1.0e6 ) -
               1.0e6 * 1.0e3 ) / 1.0e9  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from kilometers to meters does not "
             << "function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from degrees to radians.
    // The tested constants are only given to the 8th significant number.
    if ( fabs( unit_conversions::convertDegreesToRadians(
                   PhysicalConstants::OBLIQUITY_ECLIPTIC_IN_DEGREES ) -
               PhysicalConstants::OBLIQUITY_ECLIPTIC ) > 1.0e-8 )
    {
        cerr << "The conversion from degrees to radians does not "
             << "function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from degrees to arcminutes.
    if ( fabs( unit_conversions::convertDegreesToArcminutes( 43.2 ) -
               43.2 * 60.0 ) / ( 43.2 * 60.0 )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from degrees to arcminutes does not "
             << "function correctly." << endl;
        cerr << "The calculated value is: "
             << unit_conversions::convertDegreesToArcminutes( 43.2 );
        cerr << ", which has an error of "
             << unit_conversions::convertDegreesToArcminutes( 43.2 ) - 43.2 * 60.0 << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from arcminutes to arcseconds.
    if ( fabs( unit_conversions::convertArcminutesToArcseconds( 125.9 ) -
               125.9 * 60.0 ) / ( 125.9 * 60.0 )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from arcminutes to arcseconds does not function correctly."
             << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from astronomical units to meters.
    // Case of Neptune's semi-major axis used (source: Wikipedia).
    if ( fabs( unit_conversions::convertAstronomicalUnitsToMeters(
                   30.10366151 ) - 4.503443661e+12 ) > 1.0e3 )
    {
        cerr << "The conversion from astronomical units to meters does not function correctly."
             << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from minutes to seconds.
    if ( fabs( unit_conversions::convertMinutesToSeconds( 12.0 ) -
               12.0 * 60.0 ) / ( 12.0 * 60.0 )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from minutes to seconds does not function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from seconds to minutes.
    if ( fabs( unit_conversions::convertSecondsToMinutes( 12.0 ) - 0.2 ) /
         0.2 > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from seconds to minutes does not function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from hours to Julian years.
    // Three tests in one.
    if ( fabs( unit_conversions::convertJulianDaysToJulianYears(
                   unit_conversions::convertSecondsToJulianDays(
                       unit_conversions::convertHoursToSeconds( 24.0 ) ) ) -
               1.0 / 365.25 ) / ( 1.0 / 365.25 )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from hours to Julian years does not function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from Julian years to hours.
    // Three tests in one.
    if ( fabs( unit_conversions::convertSecondsToHours(
                   unit_conversions::convertJulianDaysToSeconds(
                       unit_conversions::convertJulianYearsToJulianDays( 1.0 / 365.25 ) ) ) -
               24.0 ) / ( 24.0 )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from Julian years to hours does not "
             << "function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from sidereal days to seconds.
    if ( fabs( unit_conversions::convertSiderealDaysToSeconds( 7.0 ) -
               7.0 * PhysicalConstants::SIDEREAL_DAY )
         / ( 7.0 * PhysicalConstants::SIDEREAL_DAY )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from sidereal days to seconds does not "
             << "function correctly." << endl;
        cerr << "The calculated value is: "
             << unit_conversions::convertSiderealDaysToSeconds( 7.0 );
        cerr << ", which has an error of "
             << unit_conversions::convertSiderealDaysToSeconds( 7.0 )
                - ( 7.0 * PhysicalConstants::SIDEREAL_DAY ) << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from seconds to sidereal days.
    if ( fabs( unit_conversions::convertSecondsToSiderealDays( 100.0 ) -
               100.0 / PhysicalConstants::SIDEREAL_DAY )
         / ( 100.0 / PhysicalConstants::SIDEREAL_DAY )  > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from seconds to sidereal days does not "
             << "function correctly." << endl;
        cerr << "The calculated value is: "
             << unit_conversions::convertSecondsToSiderealDays( 100.0 );
        cerr << ", which has an error of "
             << unit_conversions::convertSecondsToSiderealDays( 100.0 )
                - ( 100.0 / PhysicalConstants::SIDEREAL_DAY ) << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from temperature in Rankine to Kelvin.
    // meltingtemperature ice (source wikipedia).
    if ( fabs( unit_conversions::convertRankineToKelvin( 491.67 ) - 273.15 ) / 273.15
         > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from temperature in Rankine to Kelvin does not "
             << "function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from distance in feet to meters.
    // Case: length of a statute mile (source wikipedia).
    if ( fabs( unit_conversions::convertFeetToMeter( 5280.0 ) - 1609.344 ) / 1609.344
         > MACHINE_PRECISION_DOUBLES )
    {
        cerr << "The conversion from distance in feet to meters does not "
             << "function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Test conversion from pressure in pound per square feet to Pascal.
    // Case: atmospheric pressure at sea level. (source wikipedia).
    if ( fabs( unit_conversions::convertPoundPerSquareFeetToPascal( 2116.21662367394 ) - 101325.0 )
         > 1.0e-4 )
    {
        cerr << "The conversion from pressure in pound per square feet to Pascal does not "
             << "function correctly." << endl;
        isUnitConversionsErroneous = true;
    }

    // Return test result.
    // If test is successful return false; if test fails, return true.
    return isUnitConversionsErroneous;
}

}


// End of file.
