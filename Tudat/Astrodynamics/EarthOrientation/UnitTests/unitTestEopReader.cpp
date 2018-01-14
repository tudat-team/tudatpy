/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/Interpolators/cubicSplineInterpolator.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::earth_orientation;
using namespace tudat::interpolators;

BOOST_AUTO_TEST_SUITE( test_eop_reader )

//! Tests whether EOP data is properly read and processed, and whether data in correctly interpolated (linearly)
BOOST_AUTO_TEST_CASE( testEopReaderData )
{
    boost::shared_ptr< EOPReader > eopReader = boost::make_shared< EOPReader >(
                tudat::input_output::getEarthOrientationDataFilesPath( ) + "eopc04_08_IAU2000.62-now.txt",
                "C04", basic_astrodynamics::iau_2000_a );

    // Define current time/
    int year = 2007;
    int month = 4;
    int day = 5;
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;
    double utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 )
            * physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 5-4-2007
    double arcSecondToRadian = 4.848136811095359935899141E-6;
    double expectedXp = 0.033227 * arcSecondToRadian;
    double expectedYp = 0.483135 * arcSecondToRadian;
    double expectedUT1COffset = -0.0714209;
    double expecteddX = 0.000247 * arcSecondToRadian;
    double expecteddY = -0.000280 * arcSecondToRadian;

    // Create interpolator for x_{p} and y_{p} (polar motion) and retrieve current value
    boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            tudat::earth_orientation::createStandardEarthOrientationCalculator( eopReader )->getPolarMotionCalculator( )->
            getDailyIersValueInterpolator( );
    Eigen::Vector2d cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Create interpolator for dX and dY (nutation correction) and retrieve current value
    boost::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            tudat::earth_orientation::createStandardEarthOrientationCalculator( eopReader )->getPrecessionNutationCalculator( )->
            getDailyCorrectionInterpolator( );
    Eigen::Vector2d cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Create interpolator for d(UT1-UTC) correction  and retrieve current value
    boost::shared_ptr< OneDimensionalInterpolator< double, double > > ut1MinusUtcInterpolator =
            createDefaultTimeConverter( eopReader )->getDailyUtcUt1CorrectionInterpolator( );
    double utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Check interpolated values against test data.
    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY, std::numeric_limits< double >::epsilon( ) );

    // Reset time
    int day2 = 6;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day2, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 )
            * physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 6-4-2007
    double expectedXp2 = 0.035739  * arcSecondToRadian;
    double expectedYp2 = 0.484209 * arcSecondToRadian;
    double expectedUT1COffset2 = -0.0727562;
    double expecteddX2 = 0.000252 * arcSecondToRadian;
    double expecteddY2 = -0.000294 * arcSecondToRadian;


    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Check interpolated values against test data.
    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX2, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY2, std::numeric_limits< double >::epsilon( ) );

    // Reset time to halfway between previous times
    int hour2 = 12;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day, hour2, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 )
            * physical_constants::JULIAN_DAY;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, ( expectedUT1COffset + expectedUT1COffset2 ) / 2.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), ( expectedXp + expectedXp2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), ( expectedYp + expectedYp2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), ( expecteddX + expecteddX2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), ( expecteddY + expecteddY2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );

    // Reset time to halfway between previous times
    int hour3 = 1;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                year, month, day, hour3, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 )
            * physical_constants::JULIAN_DAY;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1, expectedUT1COffset + ( expectedUT1COffset2 - expectedUT1COffset ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), expectedXp + ( expectedXp2 - expectedXp ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), expectedYp + ( expectedYp2 - expectedYp ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), expecteddX + ( expecteddX2 - expecteddX ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), expecteddY + ( expecteddY2 - expecteddY ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Check whether UT1-UTC interpolator correctly handles leap secons
BOOST_AUTO_TEST_CASE( testLeapSecondIdentification )
{
    // Read EOP file and get UT1-UTC interpolator
    boost::shared_ptr< EOPReader > eopReader = boost::make_shared< EOPReader >(
                tudat::input_output::getEarthOrientationDataFilesPath( ) + "eopc04_08_IAU2000.62-now.txt",
                "C04", basic_astrodynamics::iau_2000_a );
    boost::shared_ptr< OneDimensionalInterpolator< double, double > > ut1MinusUtcInterpolator =
            createDefaultTimeConverter( eopReader )->getDailyUtcUt1CorrectionInterpolator( );

    // Define list of leap seconds
    Eigen::Matrix< int, Eigen::Dynamic, 3 > leapSecondDays;
    leapSecondDays.resize( 27, 3 );
    leapSecondDays << 1, 7, 1972,
            1, 1, 1973,
            1, 1, 1974,
            1, 1, 1975,
            1, 1, 1976,
            1, 1, 1977,
            1, 1, 1978,
            1, 1, 1979,
            1, 1, 1980,
            1, 7, 1981,
            1, 7, 1982,
            1, 7, 1983,
            1, 7, 1985,
            1, 1, 1988,
            1, 1, 1990,
            1, 1, 1991,
            1, 7, 1992,
            1, 7, 1993,
            1, 7, 1994,
            1, 1, 1996,
            1, 7, 1997,
            1, 1, 1999,
            1, 1, 2006,
            1, 1, 2009,
            1, 7, 2012,
            1, 7, 2015,
            1, 1, 2017;

    // Declare test variables.
    double differenceMinusOneDay, differenceMinusOneMilliSecond,
            differenceAtEpochPlusOneMilliSecond, differencePlusOneDay;

    // Set check time: at midnight
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;

    // Iterate over all leap seconds
    for( int i = 0; i < leapSecondDays.rows( ); i++ )
    {
        // Get UTC time of current leap second
        int year, month, day;
        year = leapSecondDays( i, 2 );
        month = leapSecondDays( i, 1 );
        day = leapSecondDays( i, 0 );
        double utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                    year, month, day, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 )
                * physical_constants::JULIAN_DAY;

        // Define test times
        differenceMinusOneDay = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 - physical_constants::JULIAN_DAY );
        differenceMinusOneMilliSecond = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 - 1.0E-3 );
        differenceAtEpochPlusOneMilliSecond = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 + 1.0E-3 );
        differencePlusOneDay = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 + physical_constants::JULIAN_DAY );

        // Check if leap second is properly handled
        BOOST_CHECK_SMALL( std::fabs( differenceMinusOneDay - differenceMinusOneMilliSecond ), 0.005 );
        BOOST_CHECK_SMALL( std::fabs( differencePlusOneDay - differenceAtEpochPlusOneMilliSecond ), 0.005 );

        BOOST_CHECK_SMALL( std::fabs( differenceMinusOneMilliSecond - differenceAtEpochPlusOneMilliSecond + 1.0 ), 1.0E-10 );
        BOOST_CHECK_SMALL( std::fabs( differenceMinusOneDay - differencePlusOneDay + 1.0 ), 0.01 );

    }
}



BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




