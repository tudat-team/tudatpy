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
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Basics/utilityMacros.h"

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{
namespace unit_tests
{

struct SofaTimeOutput
{
    double expectedUtcFractionalDays;
    double expectedUtcDays;
    double expectedUt1Seconds;
    double expectedUt1Days;
    double expectedTaiFractionalDays;
    double expectedTaiDays;
    double expectedTtFractionalDays;
    double expectedTtDays;
    double expectedTdbFractionalDays;
    double expectedTdbDays;
};

SofaTimeOutput getSofaDirectTimes( )
{
    SofaTimeOutput sofaTimes;

    int latnd, latnm, lonwd, lonwm, j, iy, mo, id, ih, im;
    double slatn, slonw, hm, elon, phi, xyz[3], u, v, sec,
            utc1, utc2, dut, ut11, ut12, ut, tai1, tai2, tt1, tt2,
            tcg1, tcg2, dtr, tdb1, tdb2, tcb1, tcb2;

    TUDAT_UNUSED_PARAMETER( j );

    /* Site terrestrial coordinates (WGS84). */
    latnd = 19;
    latnm = 28;
    slatn = 52.5;
    lonwd = 155;
    lonwm = 55;
    slonw = 59.6;
    hm = 0.0;

    /* Transform to geocentric. */
    j = iauAf2a ( '+', latnd, latnm, slatn, &phi );
    j = iauAf2a ( '-', lonwd, lonwm, slonw, &elon );
    j = iauGd2gc ( 1, elon, phi, hm, xyz );
    u = sqrt ( xyz[0]*xyz[0] + xyz[1]*xyz[1] );
    v = xyz[2];

    /* UTC date and time. */
    iy = 2006;
    mo = 1;
    id = 15;
    ih = 21;
    im = 24;
    sec = 37.5;

    /* Transform into internal format. */
    j = iauDtf2d ( "UTC", iy, mo, id, ih, im, sec, &utc1, &utc2 );

    /* UT1-UTC (s, from IERS). */
    dut = 0.3340960443019867; // Value modified to coincide with value in code (difference is of order microsecond, only influences utc<->ut1).

    /* UTC -> UT1. */

    j = iauUtcut1 ( utc1, utc2, dut, &ut11, &ut12 );

    /* Extract fraction for TDB-TT calculation, later. */
    ut = fmod ( fmod(ut11,1.0) + fmod(ut12,1.0), 1.0 ) + 0.5;

    /* UTC -> TAI -> TT -> TCG. */
    j = iauUtctai ( utc1, utc2, &tai1, &tai2 );
    j = iauTaitt ( tai1, tai2, &tt1, &tt2 );
    j = iauTttcg ( tt1, tt2, &tcg1, &tcg2 );

    /* TDB-TT (using TT as a substitute for TDB). */
    dtr = iauDtdb ( tt1, tt2, ut, elon, u / 1000.0, v / 1000.0 );

    /* TT -> TDB -> TCB. */
    j = iauTttdb ( tt1, tt2, dtr, &tdb1, &tdb2 );
    j = iauTdbtcb ( tdb1, tdb2, &tcb1, &tcb2 );

    /* Report. */
    sofaTimes.expectedUtcDays = utc1;
    sofaTimes.expectedUtcFractionalDays = utc2;

    sofaTimes.expectedUt1Days = ut11;
    sofaTimes.expectedUt1Seconds = ut12;

    sofaTimes.expectedTaiDays = tai1;
    sofaTimes.expectedTaiFractionalDays = tai2;

    sofaTimes.expectedTtDays = tt1;
    sofaTimes.expectedTtFractionalDays = tt2;

    sofaTimes.expectedTdbDays = tdb1;
    sofaTimes.expectedTdbFractionalDays = tdb2;
    return sofaTimes;
}


double convertSofaOutputToSecondsSinceJ2000(
        const double fullJulianDays, const double fractionalJulianDays )
{
    return ( fullJulianDays - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) *
            physical_constants::JULIAN_DAY + fractionalJulianDays * physical_constants::JULIAN_DAY ;
}

using namespace tudat::earth_orientation;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_time_scale_converters )

//! Compare time conversions to example from SOFA Cookbook. Only change w.r.t. cookbook is UT1-UTC, which is set to correspond to
//! value used here
BOOST_AUTO_TEST_CASE( testDifferentTimeScaleConversions )
{
    // Get times from SOFA cookbook
    SofaTimeOutput sofaTimes = getSofaDirectTimes( );

    // Create default time converter.
    std::shared_ptr< TerrestrialTimeScaleConverter > timeScaleConverter =
            createStandardEarthOrientationCalculator( )->getTerrestrialTimeScaleConverter( );

    // Set station position
    Eigen::Vector3d stationCartesianPosition;
    stationCartesianPosition << -5492333.306498738, -2453018.508911721, 2113645.653406073;

    // Retrieve SOFA times in double and Time precision.
    std::map< TimeScales, double > sofaSecondsSinceJ2000;
    sofaSecondsSinceJ2000[ tt_scale ] = convertSofaOutputToSecondsSinceJ2000(
                sofaTimes.expectedTtDays, sofaTimes.expectedTtFractionalDays );
    sofaSecondsSinceJ2000[ utc_scale ] = convertSofaOutputToSecondsSinceJ2000(
                sofaTimes.expectedUtcDays, sofaTimes.expectedUtcFractionalDays );
    sofaSecondsSinceJ2000[ ut1_scale ] = convertSofaOutputToSecondsSinceJ2000(
                sofaTimes.expectedUt1Days, sofaTimes.expectedUt1Seconds );
    sofaSecondsSinceJ2000[ tai_scale ] = convertSofaOutputToSecondsSinceJ2000(
                sofaTimes.expectedTaiDays, sofaTimes.expectedTaiFractionalDays );
    sofaSecondsSinceJ2000[ tdb_scale ] = convertSofaOutputToSecondsSinceJ2000(
                sofaTimes.expectedTdbDays, sofaTimes.expectedTdbFractionalDays );

    // Define list of time scales
    std::vector< TimeScales > originScales;
    originScales.push_back( tt_scale );
    originScales.push_back( utc_scale );
    originScales.push_back( ut1_scale );
    originScales.push_back( tai_scale );
    originScales.push_back( tdb_scale );

    // Test conversion between time scales, single update function at each iteration
    for( unsigned int i = 0; i < originScales.size( ); i++ )
    {
        timeScaleConverter->resetTimes< double >( );
        timeScaleConverter->updateTimes( originScales.at( i ), sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );

        double ut1 = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), ut1_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( ut1 - sofaSecondsSinceJ2000[ ut1_scale ], std::numeric_limits< double >::epsilon( ) );

        double utc = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), utc_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( utc - sofaSecondsSinceJ2000[ utc_scale ], std::numeric_limits< double >::epsilon( ) );

        double tdb = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), tdb_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tdb - sofaSecondsSinceJ2000[ tdb_scale ], std::numeric_limits< double >::epsilon( ) );

        double tai = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), tai_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tai -sofaSecondsSinceJ2000[ tai_scale ], std::numeric_limits< double >::epsilon( ) );

        double tt = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), tt_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tt -sofaSecondsSinceJ2000[ tt_scale ], std::numeric_limits< double >::epsilon( ) );

    }

    // Test conversion between time scales, without caching of values (reset function called after each iteration).
    for( unsigned int i = 0; i < originScales.size( ); i++ )
    {
        timeScaleConverter->resetTimes< double >( );
        timeScaleConverter->updateTimes(
                    originScales.at( i ), sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );

        double ut1 = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), ut1_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( ut1 - sofaSecondsSinceJ2000[ ut1_scale ], std::numeric_limits< double >::epsilon( ) );

        timeScaleConverter->resetTimes< double >( );
        double utc = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), utc_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( utc - sofaSecondsSinceJ2000[ utc_scale ], std::numeric_limits< double >::epsilon( ) );

        timeScaleConverter->resetTimes< double >( );
        double tdb = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), tdb_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tdb - sofaSecondsSinceJ2000[ tdb_scale ], std::numeric_limits< double >::epsilon( ) );

        timeScaleConverter->resetTimes< double >( );
        double tai = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), tai_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tai -sofaSecondsSinceJ2000[ tai_scale ], std::numeric_limits< double >::epsilon( ) );

        timeScaleConverter->resetTimes< double >( );
        double tt = timeScaleConverter->getCurrentTime(
                    originScales.at( i ), tt_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tt -sofaSecondsSinceJ2000[ tt_scale ], std::numeric_limits< double >::epsilon( ) );

    }

    // Test conversion between time scales, single update function at each iteration, with Time input.
    for( unsigned int i = 0; i < originScales.size( ); i++ )
    {
        timeScaleConverter->resetTimes< Time >( );
        timeScaleConverter->updateTimes< Time >( originScales.at( i ), sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );

        double ut1 = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), ut1_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( ut1 - sofaSecondsSinceJ2000[ ut1_scale ], std::numeric_limits< double >::epsilon( ) );

        double utc = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), utc_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( utc - sofaSecondsSinceJ2000[ utc_scale ], std::numeric_limits< double >::epsilon( ) );

        double tdb = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), tdb_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tdb - sofaSecondsSinceJ2000[ tdb_scale ], std::numeric_limits< double >::epsilon( ) );

        double tai = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), tai_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tai -sofaSecondsSinceJ2000[ tai_scale ], std::numeric_limits< double >::epsilon( ) );

        double tt = timeScaleConverter->getCurrentTime< Time >( originScales.at( i ), tt_scale, sofaSecondsSinceJ2000[ originScales.at( i ) ], stationCartesianPosition );
        BOOST_CHECK_SMALL( tt -sofaSecondsSinceJ2000[ tt_scale ], std::numeric_limits< double >::epsilon( ) );

    }
}


//! Test if converting back and forth is done correctly
BOOST_AUTO_TEST_CASE( testTimeScaleConversionPrecisionWithTimeType )
{
    // Create time scale converters
    std::shared_ptr< TerrestrialTimeScaleConverter > timeScaleConverter =
            createStandardEarthOrientationCalculator( )->getTerrestrialTimeScaleConverter( );
    std::shared_ptr< TerrestrialTimeScaleConverter > comparisonTimeScaleConverter =
            createStandardEarthOrientationCalculator( )->getTerrestrialTimeScaleConverter( );

    // Define time scales
    std::vector< TimeScales > originScales;
    originScales.push_back( tt_scale );
    originScales.push_back( utc_scale );
    originScales.push_back( ut1_scale );
    originScales.push_back( tai_scale );
    originScales.push_back( tdb_scale );

    // Define maps of times
    std::map< TimeScales, Time > currentTimes;
    std::map< TimeScales, Time > comparisonCurrentTimes;

    // Define current time/position
    Time baseTime = Time( -20.0 * physical_constants::JULIAN_YEAR );
    Eigen::Vector3d stationCartesianPosition;
    stationCartesianPosition << -5492333.306498738, -2453018.508911721, 2113645.653406073;

    for( unsigned int i = 0; i < originScales.size( ); i++ )
    {
        currentTimes.clear( );

        // Convert time in scale( i ) to all other times.
        timeScaleConverter->resetTimes< Time >( );
        timeScaleConverter->updateTimes< Time >( originScales.at( i ), baseTime, stationCartesianPosition );
        for( unsigned int j = 0; j < originScales.size( ); j++ )
        {
            currentTimes[ originScales.at( j ) ] =
                    timeScaleConverter->getCurrentTime< Time >(
                        originScales.at( i ), originScales.at( j ), baseTime, stationCartesianPosition );
        }

        // Convert back and compare results. Tolerance is set at only ps level, since TDB-TT, and UT1-UTC computations, are still
        // done at only double precision.
        for( unsigned int j = 0; j < originScales.size( ); j++ )
        {
            Time currentBackConvertedTime =
                    comparisonTimeScaleConverter->getCurrentTime< Time >(
                        originScales.at( j ), originScales.at( i ), currentTimes[ originScales.at( j ) ],
                    stationCartesianPosition );
            BOOST_CHECK_SMALL( std::fabs( static_cast< long double >( currentBackConvertedTime - baseTime ) ), 1.0E-12L );
        }
    }
}

//! Test validity of time scale converter around leap seconds
BOOST_AUTO_TEST_CASE( testTimeScaleConversionDuringLeapSeconds )
{
    std::shared_ptr< TerrestrialTimeScaleConverter > timeScaleConverter =
            createStandardEarthOrientationCalculator( )->getTerrestrialTimeScaleConverter( );

    // Define leap second list (day, month, year)
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

    // Convert UTC to TAI, and back, in microsecond before and after leap second, and check results in double precision.
    for( int i = 0; i < leapSecondDays.rows( ); i++ )
    {
        double utcTimeOfLeapSeconds = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                    leapSecondDays( i, 2 ), leapSecondDays( i, 1 ), leapSecondDays( i, 0 ), 0, 0, 0.0,
                    basic_astrodynamics::JULIAN_DAY_ON_J2000 );

        double taiPreLeap = timeScaleConverter->getCurrentTime< double >(
                    utc_scale, tai_scale, utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY - 1.0E-6,
                    Eigen::Vector3d::Zero( ) );
        double taiPostLeap = timeScaleConverter->getCurrentTime< double >(
                    utc_scale, tai_scale, utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY + 1.0E-6,
                    Eigen::Vector3d::Zero( ) );
        BOOST_CHECK_SMALL(  std::fabs( taiPostLeap - taiPreLeap - ( 1.0 + 2.0E-6 ) ), 1.0E-7 );

        timeScaleConverter->resetTimes< double >( );

        double utcPreLeap = timeScaleConverter->getCurrentTime< double >(
                    tai_scale, utc_scale, taiPreLeap, Eigen::Vector3d::Zero( ) );
        double utcPostLeap = timeScaleConverter->getCurrentTime< double >(
                    tai_scale, utc_scale, taiPostLeap, Eigen::Vector3d::Zero( ) );

        BOOST_CHECK_SMALL( std::fabs( utcPostLeap - utcPreLeap - ( 2.0E-6 ) ), 1.0E-7 );
    }

    // Convert UTC to TAI, and back, in microsecond before and after leap second, and check results in Time precision.
    for( int i = 0; i < leapSecondDays.rows( ); i++ )
    {
        Time utcTimeOfLeapSeconds = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                    leapSecondDays( i, 2 ), leapSecondDays( i, 1 ), leapSecondDays( i, 0 ), 0, 0, 0.0,
                    basic_astrodynamics::JULIAN_DAY_ON_J2000 );

        Time utcInputPreLeap = utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY - 1.0E-6;
        Time utcInputPostLeap = utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY + 1.0E-6;

        Time taiPreLeap = timeScaleConverter->getCurrentTime< Time >(
                    utc_scale, tai_scale, utcInputPreLeap, Eigen::Vector3d::Zero( ) );
        Time taiPostLeap = timeScaleConverter->getCurrentTime< Time >(
                    utc_scale, tai_scale, utcInputPostLeap, Eigen::Vector3d::Zero( ) );

        long double timeDifferenceTai = static_cast< long double >( taiPostLeap - taiPreLeap ) - ( 1.0L + 2.0E-6 );

        BOOST_CHECK_SMALL( std::fabs( timeDifferenceTai ),
                            3600.0L * std::numeric_limits< long double >::epsilon( ) );

        timeScaleConverter->resetTimes< Time >( );

        Time utcPreLeap = timeScaleConverter->getCurrentTime< Time >(
                    tai_scale, utc_scale, taiPreLeap, Eigen::Vector3d::Zero( ) );
        Time utcPostLeap = timeScaleConverter->getCurrentTime< Time >(
                    tai_scale, utc_scale, taiPostLeap, Eigen::Vector3d::Zero( ) );

        long double timeDifferenceUtc = static_cast< long double >( utcPostLeap - utcPreLeap ) - ( 2.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( timeDifferenceUtc ),
                           3600.0L * std::numeric_limits< long double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat

