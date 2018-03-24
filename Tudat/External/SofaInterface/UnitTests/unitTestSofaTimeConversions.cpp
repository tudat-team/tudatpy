/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <Eigen/Core>

#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"
#include "Tudat/Basics/timeType.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::sofa_interface;
using namespace tudat::basic_astrodynamics;

BOOST_AUTO_TEST_SUITE( test_sofa_time_conversions )

//! Function to test Sofa time conversion functions (except TT<->TDB) for templated time scalar type.
template< typename ScalarType >
void testTimeConversions( )
{
    // Define UTC date and time.
    int year = 2006;
    int months = 1;
    int days = 15;
    int hours = 21;
    int minutes = 24;
    ScalarType seconds = static_cast< ScalarType >( 37.5 );

    // Get UTC time.
    ScalarType utcSecondsSinceJ2000 =
            convertCalendarDateToJulianDaysSinceEpoch< ScalarType >(
                year, months, days, hours, minutes, seconds,
                getJulianDayOnJ2000< ScalarType >( ) ) *
            physical_constants::getJulianDay< ScalarType >( );

    // Convert to TAI time.
    ScalarType taiSecondsSinceJ2000 = convertUTCtoTAI( utcSecondsSinceJ2000 );

    // Get number of leap seconds and test
    ScalarType numberOfLeapSeconds = static_cast< ScalarType >( 33 );
    BOOST_CHECK_CLOSE_FRACTION( taiSecondsSinceJ2000 - utcSecondsSinceJ2000, numberOfLeapSeconds,
                                utcSecondsSinceJ2000 * std::numeric_limits< ScalarType >::epsilon( ) );

    // Recover UTC and test
    ScalarType recoveredUtcSecondsSinceJ2000 = convertTAItoUTC( taiSecondsSinceJ2000 );
    BOOST_CHECK_SMALL( recoveredUtcSecondsSinceJ2000 - utcSecondsSinceJ2000,
                       std::numeric_limits< ScalarType >::epsilon( ) );

    // Get TT time and compare against expected result.
    ScalarType ttSecondsSinceJ2000 = convertUTCtoTT( utcSecondsSinceJ2000 );
    ScalarType ttMinusTai = static_cast< ScalarType >( 32.184 );
    BOOST_CHECK_CLOSE_FRACTION( ttSecondsSinceJ2000 - utcSecondsSinceJ2000,
                                numberOfLeapSeconds + ttMinusTai,
                                utcSecondsSinceJ2000 * std::numeric_limits< ScalarType >::epsilon( ) );

    // Recover UTC from TT and compare against original UTC.
    recoveredUtcSecondsSinceJ2000 = convertTTtoUTC( ttSecondsSinceJ2000 );
    BOOST_CHECK_SMALL( recoveredUtcSecondsSinceJ2000 - utcSecondsSinceJ2000,
                       std::numeric_limits< ScalarType >::epsilon( ) );

}

//! Test TDB<->TT conversions using SOFA cookbook data
BOOST_AUTO_TEST_CASE( testSofaTimeConversions )
{

    // Test conversions for double/long double
    testTimeConversions< double >( );
    testTimeConversions< long double >( );


    // Test TDB - TT calculation, using code from SOFA Time Conversions cookbook.
    {
        // Calculate UTC time from date/time
        int year = 2006;
        int months = 1;
        int days = 15;
        int hours = 21;
        int minutes = 24;
        double seconds = 37.5;
        double utcSecondsSinceJ2000 =
                convertCalendarDateToJulianDaysSinceEpoch< double >(
                    year, months, days, hours, minutes, seconds, getJulianDayOnJ2000< double >( ) ) *
                physical_constants::getJulianDay< double >( );

        // Get UT1 time
        double dut = 0.3341;
        double ut1SecondsSinceEpoch = utcSecondsSinceJ2000 + dut;

        // Convert UTC to TT
        double ttSecondsSinceJ2000 = convertUTCtoTT( utcSecondsSinceJ2000 );

        int latnd, latnm, lonwd, lonwm, iy, mo, id, ih, im;
        double slatn, slonw, hm, elon, phi, xyz[3], u, v, sec,
                utc1, utc2, ut11, ut12, ut, tai1, tai2, tt1, tt2,
                tcg1, tcg2;

        // Run code from Sofa cookbook to obtain TDB - TT (and ground station position).
        Eigen::Vector3d referencePoint;
        double dtr;
        {
            /* UTC date and time. */
            iy = 2006;
            mo = 1;
            id = 15;
            ih = 21;
            im = 24;
            sec = 37.5;
            /* Transform into internal format. */
            iauDtf2d ( "UTC", iy, mo, id, ih, im, sec, &utc1, &utc2 );

            /* UTC -> UT1. */
            iauUtcut1 ( utc1, utc2, dut, &ut11, &ut12 );


            /* Extract fraction for TDB-TT calculation, later. */
            ut = fmod ( fmod(ut11,1.0) + fmod(ut12,1.0), 1.0 ) + 0.5;
            /* UTC -> TAI -> TT -> TCG. */
            iauUtctai ( utc1, utc2, &tai1, &tai2 );

            iauTaitt ( tai1, tai2, &tt1, &tt2 );

            iauTttcg ( tt1, tt2, &tcg1, &tcg2 );


            /* Site terrestrial coordinates (WGS84). */
            latnd = 19;
            latnm = 28;
            slatn = 52.5;
            lonwd = 155;
            lonwm = 55;
            slonw = 59.6;
            hm = 0.0;

            /* Transform to geocentric. */
            iauAf2a ( '+', latnd, latnm, slatn, &phi );
            iauAf2a ( '-', lonwd, lonwm, slonw, &elon );
            iauGd2gc ( 1, elon, phi, hm, xyz );
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
            iauDtf2d ( "UTC", iy, mo, id, ih, im, sec, &utc1, &utc2 );


            /* UT1-UTC (s, from IERS). */
            dut = 0.3341;

            /* UTC -> UT1. */
            iauUtcut1 ( utc1, utc2, dut, &ut11, &ut12 );


            /* Extract fraction for TDB-TT calculation, later. */
            ut = fmod ( fmod(ut11,1.0) + fmod(ut12,1.0), 1.0 ) + 0.5;

            /* UTC -> TAI -> TT -> TCG. */
            iauUtctai ( utc1, utc2, &tai1, &tai2 );
            iauTaitt ( tai1, tai2, &tt1, &tt2 );
            iauTttcg ( tt1, tt2, &tcg1, &tcg2 );

            /* TDB-TT (using TT as a substitute for TDB). */
            dtr = iauDtdb ( tt1, tt2, ut, elon, u/1e3, v/1e3 );

            // Define reference point and
            referencePoint = ( Eigen::Vector3d( ) << xyz[ 0 ], xyz[ 1 ], xyz[ 2 ] ).finished( );
        }

        // Calculate UT1 fraction of day
        double utFractionOfDay = std::fmod( ut1SecondsSinceEpoch / physical_constants::JULIAN_DAY + 0.5, 1.0 );

        // Calculate TDB - TT from Sofa and comapre against cookbook result.
        double tdbMinusTt = getTDBminusTT( ttSecondsSinceJ2000, utFractionOfDay, referencePoint );
        BOOST_CHECK_SMALL( tdbMinusTt - dtr, std::numeric_limits< double >::epsilon( ) );

        // Check validity of using TT as subsititute for TDB in conversions
        double dtr2 = iauDtdb ( tt1, tt2 + dtr / physical_constants::JULIAN_DAY, ut, elon, u/1.0e3, v/1.0e3 );
        double dtr3 = iauDtdb ( tt1, tt2 + dtr2 / physical_constants::JULIAN_DAY, ut, elon, u/1.0e3, v/1.0e3 );
        BOOST_CHECK_SMALL( dtr - dtr2, 1.0E-12 );
        BOOST_CHECK_SMALL( dtr2 - dtr3, 1.0E-15 );

        // Check approximate conversion, omitting utc-ut1 correction.
        double tdbMinusTtApproximate = getTDBminusTT( ttSecondsSinceJ2000, referencePoint );
        BOOST_CHECK_SMALL( tdbMinusTtApproximate - tdbMinusTt, 1.0E-10 );    }

}

BOOST_AUTO_TEST_CASE( testLeapSecondIdentification )
{
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

    for( unsigned int i = 0; i < leapSecondDays.rows( ); i++ )
    {
        double utcTimeOfLeapSeconds = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                    leapSecondDays( i, 2 ), leapSecondDays( i, 1 ), leapSecondDays( i, 0 ), 0, 0, 0.0,
                    basic_astrodynamics::JULIAN_DAY_ON_J2000 );
        BOOST_CHECK_EQUAL( sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds - 1.0E-6 ) -
                           sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds + 1.0E-6 ), -1.0 );
        BOOST_CHECK_EQUAL( sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds + 1.0E-6 ) -
                           sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds + 2.0E-6 ), 0.0 );

        double taiPreLeap = tudat::sofa_interface::convertUTCtoTAI< double >(
                    utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY - 1.0E-6 );
        double taiPostLeap = tudat::sofa_interface::convertUTCtoTAI< double >(
                    utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY + 1.0E-6 );
        BOOST_CHECK_SMALL(  std::fabs( taiPostLeap - taiPreLeap - ( 1.0 + 2.0E-6 ) ), 1.0E-7 );

        double utcPreLeap = tudat::sofa_interface::convertTAItoUTC< double >( taiPreLeap );
        double utcPostLeap = tudat::sofa_interface::convertTAItoUTC< double >( taiPostLeap );

        BOOST_CHECK_SMALL( std::fabs( utcPostLeap - utcPreLeap - ( 2.0E-6 ) ), 1.0E-7 );
    }

    for( unsigned int i = 0; i < leapSecondDays.rows( ); i++ )
    {
        Time utcTimeOfLeapSeconds = basic_astrodynamics::convertCalendarDateToJulianDaysSinceEpoch(
                    leapSecondDays( i, 2 ), leapSecondDays( i, 1 ), leapSecondDays( i, 0 ), 0, 0, 0.0,
                    basic_astrodynamics::JULIAN_DAY_ON_J2000 );
        BOOST_CHECK_EQUAL( sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds - 1.0E-6 ) -
                           sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds + 1.0E-6 ), -1.0 );
        BOOST_CHECK_EQUAL( sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds + 1.0E-6 ) -
                           sofa_interface::getDeltaAtFromUtc( utcTimeOfLeapSeconds + 2.0E-6 ), 0.0 );

        Time taiPreLeap = tudat::sofa_interface::convertUTCtoTAI< Time >(
                    utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY - 1.0E-6 );
        Time taiPostLeap = tudat::sofa_interface::convertUTCtoTAI< Time >(
                    utcTimeOfLeapSeconds * physical_constants::JULIAN_DAY + 1.0E-6 );

        long double timeDifferenceTai = static_cast< long double >( taiPostLeap - taiPreLeap ) - ( 1.0L + 2.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( timeDifferenceTai ), 3600.0L * std::numeric_limits< long double >::epsilon( ) );

        Time utcPreLeap = tudat::sofa_interface::convertTAItoUTC< Time >( taiPreLeap );
        Time utcPostLeap = tudat::sofa_interface::convertTAItoUTC< Time >( taiPostLeap );

        long double timeDifferenceUtc = static_cast< long double >( utcPostLeap - utcPreLeap ) - ( 2.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( timeDifferenceUtc ), 3600.0L * std::numeric_limits< long double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




