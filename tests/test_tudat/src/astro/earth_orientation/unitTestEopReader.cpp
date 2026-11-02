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

#define BOOST_TEST_MAIN

#include <limits>
#include <sstream>
#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/io/basicInputOutput.h"

#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/astro/earth_orientation/readEarthOrientationData.h"
#include "tudat/math/interpolators/jumpDataLinearInterpolator.h"

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
    std::shared_ptr< EOPReader > eopReader = std::make_shared< EOPReader >(
            tudat::paths::getEarthOrientationDataFilesPath( ) + "/eopc04_14_IAU2000.62-now.txt", "C04", basic_astrodynamics::iau_2000_a );

    // Define current time/
    int year = 2007;
    int month = 4;
    int day = 5;
    int hour = 0;
    int minute = 0;
    double seconds = 0.0;
    double utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                          year, month, day, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 5-4-2007
    double arcSecondToRadian = 4.848136811095359935899141E-6;
    double expectedXp = 0.033194 * arcSecondToRadian;
    double expectedYp = 0.483144 * arcSecondToRadian;
    double expectedUT1COffset = -0.0714163;
    double expecteddX = 0.000250 * arcSecondToRadian;
    double expecteddY = -0.000302 * arcSecondToRadian;

    // Create interpolator for x_{p} and y_{p} (polar motion) and retrieve current value
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector2d > > cipInItrsInterpolator =
            tudat::earth_orientation::createStandardEarthOrientationCalculator( eopReader )
                    ->getPolarMotionCalculator( )
                    ->getDailyIersValueInterpolator( );
    Eigen::Vector2d cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Create interpolator for dX and dY (nutation correction) and retrieve current value
    std::shared_ptr< OneDimensionalInterpolator< double, Eigen::Vector2d > > cipInGcrsCorrectionInterpolator =
            tudat::earth_orientation::createStandardEarthOrientationCalculator( eopReader )
                    ->getPrecessionNutationCalculator( )
                    ->getDailyCorrectionInterpolator( );
    Eigen::Vector2d cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );

    // Create interpolator for d(UT1-UTC) correction  and retrieve current value
    std::shared_ptr< OneDimensionalInterpolator< double, double > > ut1MinusUtcInterpolator =
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
                                   year, month, day2, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Set EOP corrections read manually from file for 6-4-2007
    double expectedXp2 = 0.035736 * arcSecondToRadian;
    double expectedYp2 = 0.484204 * arcSecondToRadian;
    double expectedUT1COffset2 = -0.0728009;
    double expecteddX2 = 0.000254 * arcSecondToRadian;
    double expecteddY2 = -0.000309 * arcSecondToRadian;

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
                                   year, month, day, hour2, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION(
            utcMinusUt1, ( expectedUT1COffset + expectedUT1COffset2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.x( ), ( expectedXp + expectedXp2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInItrs.y( ), ( expectedYp + expectedYp2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.x( ), ( expecteddX + expecteddX2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( cipInGcrs.y( ), ( expecteddY + expecteddY2 ) / 2.0, std::numeric_limits< double >::epsilon( ) );

    // Reset time to halfway between previous times
    int hour3 = 1;
    utcSecondsSinceJ2000 = tudat::basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                                   year, month, day, hour3, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY;

    // Read correction interpolators
    cipInItrs = cipInItrsInterpolator->interpolate( utcSecondsSinceJ2000 );
    cipInGcrs = cipInGcrsCorrectionInterpolator->interpolate( utcSecondsSinceJ2000 );
    utcMinusUt1 = ut1MinusUtcInterpolator->interpolate( utcSecondsSinceJ2000 );

    BOOST_CHECK_CLOSE_FRACTION( utcMinusUt1,
                                expectedUT1COffset + ( expectedUT1COffset2 - expectedUT1COffset ) / 24.0,
                                std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInItrs.x( ), expectedXp + ( expectedXp2 - expectedXp ) / 24.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInItrs.y( ), expectedYp + ( expectedYp2 - expectedYp ) / 24.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInGcrs.x( ), expecteddX + ( expecteddX2 - expecteddX ) / 24.0, std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            cipInGcrs.y( ), expecteddY + ( expecteddY2 - expecteddY ) / 24.0, std::numeric_limits< double >::epsilon( ) );
}

//! Check whether UT1-UTC interpolator correctly handles leap secons
BOOST_AUTO_TEST_CASE( testLeapSecondIdentification )
{
    // Read EOP file and get UT1-UTC interpolator
    std::shared_ptr< EOPReader > eopReader = std::make_shared< EOPReader >(
            tudat::paths::getEarthOrientationDataFilesPath( ) + "/eopc04_14_IAU2000.62-now.txt", "C04", basic_astrodynamics::iau_2000_a );
    std::shared_ptr< OneDimensionalInterpolator< double, double > > ut1MinusUtcInterpolator =
            createDefaultTimeConverter( eopReader )->getDailyUtcUt1CorrectionInterpolator( );

    // Define list of leap seconds
    Eigen::Matrix< int, Eigen::Dynamic, 3 > leapSecondDays;
    leapSecondDays.resize( 27, 3 );
    leapSecondDays << 1, 7, 1972, 1, 1, 1973, 1, 1, 1974, 1, 1, 1975, 1, 1, 1976, 1, 1, 1977, 1, 1, 1978, 1, 1, 1979, 1, 1, 1980, 1, 7,
            1981, 1, 7, 1982, 1, 7, 1983, 1, 7, 1985, 1, 1, 1988, 1, 1, 1990, 1, 1, 1991, 1, 7, 1992, 1, 7, 1993, 1, 7, 1994, 1, 1, 1996, 1,
            7, 1997, 1, 1, 1999, 1, 1, 2006, 1, 1, 2009, 1, 7, 2012, 1, 7, 2015, 1, 1, 2017;

    // Declare test variables.
    double differenceMinusOneDay, differenceMinusOneMilliSecond, differenceAtEpochPlusOneMilliSecond, differencePlusOneDay;

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
                                              year, month, day, hour, minute, seconds, basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
                physical_constants::JULIAN_DAY;

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

//! Check reading of modern (20 C04) C04 data, using an extract of the official IERS file (header and first three rows)
BOOST_AUTO_TEST_CASE( testModernC04Reading )
{
    std::istringstream stream(
            R"(# EARTH ORIENTATION PARAMETER (EOP) PRODUCT CENTER CENTER (PARIS OBSERVATORY) - INTERNATIONAL EARTH ROTATION AND REFERENCE SYSTEMS SERVICE  
# EOP (IERS) 20 C04 TIME SERIES consistent with ITRF 2020 - sampled at 12h UTC       
# Description: https://hpiers.obspm.fr/eoppc/eop/eopc04/readme        Contact: christian.bizouard@obspm.fr                              
# Reference Precession-Nutation Model: IAU 2000
# format(4(i4),f10.2,2(f12.6),f12.7,2(f12.6),2(f12.6),f12.7,2(f12.6),f12.7,2(f12.6),2(f12.6),f12.7)
# YR  MM  DD  HH       MJD        x(")        y(")  UT1-UTC(s)       dX(")       dY(")  xrt("/day)  yrt("/day)      LOD(s)        x Er        y Er  UT1-UTC Er       dX Er       dY Er      xrt Er      yrt Er      LOD Er
1984   1   1  12  45700.50   -0.134064    0.093057   0.3968898    0.002255   -0.002791   -0.002644    0.003317   0.0014139    0.001250    0.000992   0.0002417    0.000349    0.000351    0.000873    0.000771   0.0000730
1984   1   2  12  45701.50   -0.136239    0.096494   0.3955314    0.001470   -0.001946   -0.002689    0.003144   0.0013168    0.001201    0.000990   0.0002417    0.000349    0.000351    0.000841    0.000770   0.0000730
1984   1   3  12  45702.50   -0.139537    0.099240   0.3942317    0.000860   -0.001281   -0.003560    0.002270   0.0012952    0.001156    0.000987   0.0002417    0.000349    0.000351    0.000804    0.000786   0.0000730
)" );
    std::map< double, EarthOrientationRecord > records = readEarthOrientationData( stream, EarthOrientationFormat::c04 );

    const double arcsecond = mathematical_constants::PI / 648000.0;
    BOOST_REQUIRE_EQUAL( records.size( ), 3 );
    BOOST_CHECK_EQUAL( records.begin( )->first, 45700.5 );
    BOOST_CHECK_EQUAL( records.rbegin( )->first, 45702.5 );

    const EarthOrientationRecord& firstRecord = records.begin( )->second;
    BOOST_CHECK_SMALL( firstRecord.polarMotion[ 0 ] + 0.134064 * arcsecond, 1.0E-20 );
    BOOST_CHECK_SMALL( firstRecord.polarMotion[ 1 ] - 0.093057 * arcsecond, 1.0E-20 );
    BOOST_CHECK_EQUAL( firstRecord.ut1MinusUtc, 0.3968898 );
    BOOST_REQUIRE( firstRecord.lengthOfDay.has_value( ) );
    BOOST_CHECK_EQUAL( *firstRecord.lengthOfDay, 0.0014139 );
    BOOST_REQUIRE( firstRecord.celestialPoleOffsets.has_value( ) );
    BOOST_CHECK_SMALL( ( *firstRecord.celestialPoleOffsets )[ 0 ] - 0.002255 * arcsecond, 1.0E-22 );
    BOOST_CHECK_SMALL( ( *firstRecord.celestialPoleOffsets )[ 1 ] + 0.002791 * arcsecond, 1.0E-22 );
}

//! Check reading of USNO finals2000A data, using an extract of the official file (three observed rows, one prediction without
//! LOD and celestial pole offsets, and one date-only future placeholder)
BOOST_AUTO_TEST_CASE( testFinals2000aReading )
{
    std::istringstream stream(
            R"(73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916  P    -0.766    0.199    -0.720    0.300   .143000   .137000   .8075000   -18.637    -3.667  
73 1 3 41685.00 I  0.118980 0.011039  0.135656 0.013616  I 0.8056163 0.0002710  3.5563 0.1916  P    -0.751    0.199    -0.701    0.300   .141000   .134000   .8044000   -18.636    -3.571  
73 1 4 41686.00 I  0.117227 0.011039  0.134348 0.013616  I 0.8027895 0.0002710  2.6599 0.1916  P    -0.738    0.199    -0.662    0.300   .139000   .131000   .8012000   -18.669    -3.621  
261222 61396.00 P  0.091510 0.007564  0.349919 0.010511  P-0.1142750 0.0072970                                                                                                             
2710 5 61683.00                                                                                                                                                                            
)" );
    std::map< double, EarthOrientationRecord > records = readEarthOrientationData( stream, EarthOrientationFormat::finals2000a );

    const double arcsecond = mathematical_constants::PI / 648000.0;
    BOOST_REQUIRE_EQUAL( records.size( ), 4 );

    const EarthOrientationRecord& firstRecord = records.begin( )->second;
    BOOST_CHECK_SMALL( firstRecord.polarMotion[ 0 ] - 0.120733 * arcsecond, 1.0E-20 );
    BOOST_CHECK_EQUAL( firstRecord.ut1MinusUtc, 0.8084178 );

    BOOST_REQUIRE_EQUAL( records.count( 61396.0 ), 1 );
    BOOST_CHECK( !records.at( 61396.0 ).lengthOfDay.has_value( ) );
    BOOST_CHECK( !records.at( 61396.0 ).celestialPoleOffsets.has_value( ) );
    BOOST_CHECK_EQUAL( records.count( 61683.0 ), 0 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
