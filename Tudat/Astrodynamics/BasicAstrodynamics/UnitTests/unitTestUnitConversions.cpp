/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wikipedia. http://en.wikipedia.org/wiki/Neptune, last accessed: 27 January, 2012(a).
 *      Wikipedia. http://en.wikipedia.org/wiki/Mile, last accessed: 27 January, 2012(b).
 *      Wikipedia. http://en.wikipedia.org/wiki/Atmospheric_pressure,
 *          last accessed: 27 January, 2012(c).
 *      Wikipedia. http://en.wikipedia.org/wiki/Temperature_conversion_formulas,
 *          last accessed: 27, January 2012(d).
 *
 *    Notes
 *      At the moment, not all conversion routines are tested both ways. This should be corrected
 *      in an update version.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

namespace tudat
{
namespace unit_tests
{

//! Test suit for unit conversions.
BOOST_AUTO_TEST_SUITE( test_unit_conversions )

//! Test conversion from kilometers to meters.
BOOST_AUTO_TEST_CASE( testConversionFromKilometersToMeters )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertKilometersToMeters( 1.0e6 ),
                                1.0e6 * 1.0e3, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from degrees to radians.
BOOST_AUTO_TEST_CASE( testConversionFromDegreesToRadians )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertDegreesToRadians( 45.0 ),
                                mathematical_constants::PI / 4.0,
                                std::numeric_limits< double >::epsilon( ) );
}


//! Test conversion from degrees to arcminutes.
BOOST_AUTO_TEST_CASE( testConversionFromDegreesToArcminutes )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertDegreesToArcminutes( 43.2 ),
                                43.2 * 60.0, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from arcminutes to arcseconds.
BOOST_AUTO_TEST_CASE( testConversionFromArcminutesToArcSeconds )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertArcminutesToArcseconds( 125.9 ),
                                125.9 * 60.0, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from astronomical units to meters.
BOOST_AUTO_TEST_CASE( testConversionFromAstronomicalUnitsToMeters )
{
    // Case: Neptune's semi-major axis (Wikipedia, 2012a).
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertAstronomicalUnitsToMeters( 30.10366151 ),
                                4.503443661e+12, 1.0e-9 );
}

//! Test conversion from minutes to seconds.
BOOST_AUTO_TEST_CASE( testConversionFromMinutesToSeconds )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertMinutesToSeconds( 12.0 ),
                                12.0 * 60.0, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from seconds to minutes.
BOOST_AUTO_TEST_CASE( testConversionFromSecondsToMinutes )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertSecondsToMinutes( 12.0 ),
                                0.2, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from hours to Julian years.
BOOST_AUTO_TEST_CASE( testConversionFromHoursToJulianYears )
{
    using namespace unit_conversions;

    BOOST_CHECK_CLOSE_FRACTION( convertJulianDaysToJulianYears(
                                    convertSecondsToJulianDays(
                                        convertHoursToSeconds( 24.0 ) ) ),
                                1.0 / 365.25, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from Julian years to hours.
BOOST_AUTO_TEST_CASE( testConversionFromJulianYearsToHours )
{
    using namespace unit_conversions;

    BOOST_CHECK_CLOSE_FRACTION( convertSecondsToHours(
                                    convertJulianDaysToSeconds(
                                        convertJulianYearsToJulianDays(
                                            1.0 / 365.25 ) ) ),
                                24.0, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from sidereal days to seconds.
BOOST_AUTO_TEST_CASE( testConversionFromSiderealDaysToSeconds )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertSiderealDaysToSeconds( 7.0 ),
                                7.0 * physical_constants::SIDEREAL_DAY,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from seconds to sidereal days.
BOOST_AUTO_TEST_CASE( testConversionFromSecondsToSiderealDays )
{
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertSecondsToSiderealDays( 100.0 ),
                                100.0
                                / physical_constants::SIDEREAL_DAY,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion of temperature in Rankine to Kelvin.
BOOST_AUTO_TEST_CASE( testConversionFromRankineToKelvin )
{
    // Case: 0 deg Celcius (Wikipedia, 2011d).
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertRankineToKelvin( 491.67 ),
                                273.15, std::numeric_limits< double >::epsilon( ) );
}

//! Test conversion from distance in feet to meters.
BOOST_AUTO_TEST_CASE( testConversionFromFeetToMeters )
{
    // Case: length of a statute mile (Wikipedia, 2011b).
    BOOST_CHECK_CLOSE_FRACTION( unit_conversions::
                                convertFeetToMeter( 5280.0 ),
                                1609.344, std::numeric_limits< double >::epsilon( ) );
}


//! Test conversion from pounds-per-square-feet to Pascal.
BOOST_AUTO_TEST_CASE( testConversionFromPoundsPerSquareFeetToPascal )
{
    // Case: atmospheric pressure at sea level (Wikipedia, 2011c).
    BOOST_CHECK_CLOSE_FRACTION(
                unit_conversions::
                convertPoundPerSquareFeetToPascal( 2116.21662367394 ),
                101325.0, 1.0e-9 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
