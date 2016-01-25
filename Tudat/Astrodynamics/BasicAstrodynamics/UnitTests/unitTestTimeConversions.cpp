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
 *      130218    D. Dirkx          File created from personal application.
 *      130301    K. Kumar          Split and refactored unit tests.
 *
 *    References
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;

//! Test the functionality of the time conversion functions.
BOOST_AUTO_TEST_SUITE( test_Time_Conversions )

//! Unit test for Julian day to seconds conversion function.
BOOST_AUTO_TEST_CASE( testJulianDayToSecondsConversions )
{
    // Test conversion from Julian day to seconds since epoch at 0 MJD.
    {
        // Set reference epoch and Julian day for tests.
        const double referenceEpoch = JULIAN_DAY_AT_0_MJD;
        const double julianDay = JULIAN_DAY_AT_0_MJD + 1.0e6 / 86400.0;

        // Set expected seconds since epoch result.
        const double expectedSecondsSinceEpoch = 1.0e6;

        // Compute seconds since epoch given Julian day and reference epoch.
        const double computedSecondsSinceEpoch = convertJulianDayToSecondsSinceEpoch(
                    julianDay, referenceEpoch );

        // Test that computed result matches expected result.
        // Test is run at reduced tolerance, because the final digits of the seconds were lost
        // when converting to Julian day.
        BOOST_CHECK_CLOSE_FRACTION( computedSecondsSinceEpoch, expectedSecondsSinceEpoch,
                                    1.0e-11 );
    }

    // Test conversion from Julian day to seconds since J2000 epoch.
    {
        // Set reference epoch and Julian day for tests.
        const double referenceEpoch = JULIAN_DAY_ON_J2000;
        const double julianDay = JULIAN_DAY_ON_J2000 + 0.5;

        // Set expected seconds since epoch result.
        double expectedSecondsSinceEpoch
                = physical_constants::JULIAN_DAY / 2.0;

        // Compute seconds since epoch given Julian day and reference epoch.
        const double computedSecondsSinceEpoch = convertJulianDayToSecondsSinceEpoch(
                    julianDay, referenceEpoch );

        // Test that computed result matches expected result.
        // Test is run at reduced tolerance, because the final digits of the seconds were lost
        // when converting to Julian day.
        BOOST_CHECK_CLOSE_FRACTION( computedSecondsSinceEpoch, expectedSecondsSinceEpoch,
                                    std::numeric_limits< double >::epsilon( ) );
    }
}

//! Unit test for seconds to Julian day conversion function.
BOOST_AUTO_TEST_CASE( testSecondsSinceEpochToJulianDayConversions )
{
    // Test conversion from seconds since epoch to Julian day.

    // Set reference epoch and seconds since epoch for tests.
    const double referenceEpoch = JULIAN_DAY_AT_0_MJD;
    const double secondsSinceEpoch = 1.0e6;

    // Set expected Julian day result.
    const double expectedJulianDay = JULIAN_DAY_AT_0_MJD + 1.0e6 / 86400.0;

    // Compute Julian day with seconds since reference epoch specified.
    const double computedJulianDay = convertSecondsSinceEpochToJulianDay(
                secondsSinceEpoch, referenceEpoch );

    // Test that computed result matches expected result.
    BOOST_CHECK_CLOSE_FRACTION( computedJulianDay, expectedJulianDay,
                                std::numeric_limits< double >::epsilon( ) );
}

//! Unit test for calendar date to Julian day conversion function.
BOOST_AUTO_TEST_CASE( testConversionCalendarDateToJulianDay )
{
    // Compute the Julian day of the calendar date: January 1st, 2000, at 12h0m0s.
    {
        //Use the function to compute the Julian day.
        const double computedJulianDay = convertCalendarDateToJulianDay ( 2000, 1, 1, 12, 0, 0 );

        //Known Julian day at this calendar date.
        const double expectedJulianDay = basic_astrodynamics::JULIAN_DAY_ON_J2000;

        // Test that computed result matches expected result.
        BOOST_CHECK_CLOSE_FRACTION( computedJulianDay, expectedJulianDay,
                                std::numeric_limits< double >::epsilon( ) );
    }

    //Compute the Julian day of the calendar date: November 17th, 1858, At 0h0m0s.
    {
        //Use the function to compute the Julian day.
        const double computedJulianDay = convertCalendarDateToJulianDay( 1858, 11, 17, 0, 0, 0 );

        //Known Julian day at this calendar date
        const double expectedJulianDay = basic_astrodynamics::JULIAN_DAY_AT_0_MJD;

        // Test that computed result matches expected result.
        BOOST_CHECK_CLOSE_FRACTION( computedJulianDay, expectedJulianDay,
                                std::numeric_limits< double >::epsilon( ) );
    }

    //Test conversion wrapper against boost result.
    {
        const int year = 1749;
        const int month = 3;
        const int day = 30;

        BOOST_CHECK_CLOSE_FRACTION( boost::gregorian::date( year, month, day ).julian_day( ) - 0.5,
                                    convertCalendarDateToJulianDay( year, month, day, 0, 0, 0.0 ),
                                    std::numeric_limits< double >::epsilon( ) );
        const int hour = 4;
        BOOST_CHECK_CLOSE_FRACTION( boost::gregorian::date( year, month, day ).julian_day( ) - 0.5 +
                                    static_cast< double >( hour ) / 24.0,
                                    convertCalendarDateToJulianDay( year, month, day, hour, 0, 0.0 ),
                                    std::numeric_limits< double >::epsilon( ) );

        const int minute = 36;
        BOOST_CHECK_CLOSE_FRACTION( boost::gregorian::date( year, month, day ).julian_day( ) - 0.5 +
                                    static_cast< double >( hour ) / 24.0 + static_cast< double >( minute ) / ( 24.0 * 60.0 ),
                                    convertCalendarDateToJulianDay( year, month, day, hour, minute, 0.0 ),
                                    std::numeric_limits< double >::epsilon( ) );

        const double second = 21.36474854836359;
        BOOST_CHECK_CLOSE_FRACTION( boost::gregorian::date( year, month, day ).julian_day( ) - 0.5 +
                                    static_cast< double >( hour ) / 24.0 + static_cast< double >( minute ) / ( 24.0 * 60.0 ) +
                                    second / ( 24.0 * 60.0 * 60.0 ),
                                    convertCalendarDateToJulianDay( year, month, day, hour, minute, second ),
                                    std::numeric_limits< double >::epsilon( ) );

    }
}

BOOST_AUTO_TEST_CASE( testTimeConversions )
{
    const double testModifiedJulianDay = 54583.87;
    const double testJulianDay = testModifiedJulianDay + JULIAN_DAY_AT_0_MJD;

    BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertModifiedJulianDayToJulianDay( testModifiedJulianDay ),
                                5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay, convertJulianDayToModifiedJulianDay( testJulianDay ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( convertJulianDayToModifiedJulianDay(
                                    convertModifiedJulianDayToJulianDay( testModifiedJulianDay ) ),
                                testModifiedJulianDay, 10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( convertModifiedJulianDayToJulianDay(
                                    convertJulianDayToModifiedJulianDay( testJulianDay ) ),
                                testJulianDay, 5.0 * std::numeric_limits< double >::epsilon( ) );

    double secondsSinceModifedJulianDayZero = testModifiedJulianDay * JULIAN_DAY;

    BOOST_CHECK_CLOSE_FRACTION( secondsSinceModifedJulianDayZero, convertJulianDayToSecondsSinceEpoch(
                                    testJulianDay, JULIAN_DAY_AT_0_MJD ),
                                10.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertSecondsSinceEpochToJulianDay(
                                    secondsSinceModifedJulianDayZero, JULIAN_DAY_AT_0_MJD ),
                                5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay /  JULIAN_YEAR_IN_DAYS, convertSecondsSinceEpochToJulianYearsSinceEpoch(
                                    secondsSinceModifedJulianDayZero ),
                                5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay / (  JULIAN_YEAR_IN_DAYS * 100.0 ),
                                convertSecondsSinceEpochToJulianCenturiesSinceEpoch( secondsSinceModifedJulianDayZero ),
                                5.0 * std::numeric_limits< double >::epsilon( ) );

    int testYear = 2008;
    int testMonth = 4;
    int testDay = 27;
    boost::gregorian::date testCalendarDate( testYear, testMonth, testDay );
    double testFractionOfDay = 0.37;

    boost::gregorian::date calendarDate = convertJulianDayToCalendarDate( testJulianDay );

    BOOST_CHECK_EQUAL( testYear, calendarDate.year( ) );
    BOOST_CHECK_EQUAL( testMonth, calendarDate.month( ) );
    BOOST_CHECK_EQUAL( testDay, calendarDate.day( ) );

    calendarDate = convertJulianDayToCalendarDate( testJulianDay + 0.5 );

    BOOST_CHECK_EQUAL( testYear, calendarDate.year( ) );
    BOOST_CHECK_EQUAL( testMonth, calendarDate.month( ) );
    BOOST_CHECK_EQUAL( testDay + 1, calendarDate.day( ) );

    BOOST_CHECK_CLOSE_FRACTION( calculateJulianDaySinceEpoch(
                testCalendarDate, 0.87, JULIAN_DAY_AT_0_MJD ), testModifiedJulianDay,
                                50.0 * std::numeric_limits< double >::epsilon( ) );

    double secondsSinceJ2000Synchronization = ( TAI_JULIAN_DAY_SINCE_J2000_AT_TIME_SYNCHRONIZATION ) * JULIAN_DAY;

    double testTcg = convertTtToTcg( secondsSinceJ2000Synchronization );
    double testTt = convertTcgToTt( secondsSinceJ2000Synchronization );

    BOOST_CHECK_CLOSE_FRACTION( testTcg, secondsSinceJ2000Synchronization, 5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTt, secondsSinceJ2000Synchronization, 5.0 * std::numeric_limits< double >::epsilon( ) );

    double testTcb = convertTdbToTcb( secondsSinceJ2000Synchronization );
    double testTdb = convertTcbToTdb( secondsSinceJ2000Synchronization );

    BOOST_CHECK_CLOSE_FRACTION( testTdb, secondsSinceJ2000Synchronization + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION, 5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTcb, secondsSinceJ2000Synchronization - TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION, 5.0 * std::numeric_limits< double >::epsilon( ) );

    double testTime = 0.0;

    testTcg = convertTtToTcg( 0.0 );
    testTt = convertTcgToTt( testTcg );

    double expectedTcg = -secondsSinceJ2000Synchronization * LG_TIME_RATE_TERM;
    BOOST_CHECK_CLOSE_FRACTION( testTcg, expectedTcg, 5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( testTt, std::numeric_limits< double >::epsilon( ) );


    testTdb = convertTcbToTdb( 0.0 );
    testTcb = convertTdbToTcb( testTdb );

    double expectedTdb = secondsSinceJ2000Synchronization * LB_TIME_RATE_TERM + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION;

    BOOST_CHECK_CLOSE_FRACTION( testTdb, expectedTdb, 5.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_SMALL( testTcb, std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_CASE( testTimeConversionsLong )
{
    const long double testModifiedJulianDay = static_cast< long double >( 54583.87 );
    const long double testJulianDay = testModifiedJulianDay + JULIAN_DAY_AT_0_MJD_LONG;

    BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertModifiedJulianDayToJulianDay< long double >( testModifiedJulianDay ),
                                5.0 * std::numeric_limits< long double >::epsilon( ));
    BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay, convertJulianDayToModifiedJulianDay< long double >( testJulianDay ),
                                5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( convertJulianDayToModifiedJulianDay< long double >(
                                    convertModifiedJulianDayToJulianDay< long double >( testModifiedJulianDay ) ),
                                testModifiedJulianDay, 5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( convertModifiedJulianDayToJulianDay(
                                    convertJulianDayToModifiedJulianDay( testJulianDay ) ),
                                testJulianDay, 5.0 * std::numeric_limits< long double >::epsilon( ) );

    long double secondsSinceModifedJulianDayZero = testModifiedJulianDay * JULIAN_DAY_LONG;

    BOOST_CHECK_CLOSE_FRACTION( secondsSinceModifedJulianDayZero, convertJulianDayToSecondsSinceEpoch< long double >(
                                    testJulianDay, JULIAN_DAY_AT_0_MJD_LONG ),
                                5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertSecondsSinceEpochToJulianDay< long double >(
                                    secondsSinceModifedJulianDayZero, JULIAN_DAY_AT_0_MJD_LONG ),
                                5.0 * std::numeric_limits< long double >::epsilon( ) );

    long double secondsSinceJ2000Synchronization = ( TAI_JULIAN_DAY_SINCE_J2000_AT_TIME_SYNCHRONIZATION_LONG ) * JULIAN_DAY_LONG;

    long double testTcg = convertTtToTcg< long double >( secondsSinceJ2000Synchronization );
    long double testTt = convertTcgToTt< long double >( secondsSinceJ2000Synchronization );

    BOOST_CHECK_CLOSE_FRACTION( testTcg, secondsSinceJ2000Synchronization, 5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTt, secondsSinceJ2000Synchronization, 5.0 * std::numeric_limits< long double >::epsilon( ) );

    long double testTcb = convertTdbToTcb< long double >( secondsSinceJ2000Synchronization );
    long double testTdb = convertTcbToTdb< long double >( secondsSinceJ2000Synchronization );

    BOOST_CHECK_CLOSE_FRACTION( testTdb, secondsSinceJ2000Synchronization + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION, 5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTcb, secondsSinceJ2000Synchronization - TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION, 5.0 * std::numeric_limits< long double >::epsilon( ) );

    long double testTime = 0.0;

    testTcg = convertTtToTcg< long double >( static_cast< long double >( 0.0 ) );
    testTt = convertTcgToTt< long double >( testTcg );

    long double expectedTcg = -secondsSinceJ2000Synchronization * LG_TIME_RATE_TERM_LONG;
    BOOST_CHECK_CLOSE_FRACTION( testTcg, expectedTcg, 5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTt, 0.0, 5.0 * std::numeric_limits< long double >::epsilon( ) );


    testTdb = convertTcbToTdb< long double >( static_cast< long double >( 0.0 ) );
    testTcb = convertTdbToTcb< long double >( testTdb );

    long double expectedTdb = secondsSinceJ2000Synchronization * LB_TIME_RATE_TERM + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION;

    BOOST_CHECK_CLOSE_FRACTION( testTdb, expectedTdb, 5.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_SMALL( testTcb, 5.0 * std::numeric_limits< long double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
