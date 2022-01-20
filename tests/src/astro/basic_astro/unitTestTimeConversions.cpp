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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/basic_astro/timeConversions.h"

namespace tudat
{
namespace unit_tests
{

using namespace basic_astrodynamics;
using namespace physical_constants;

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

    // Test conversion from Julian day to seconds since epoch at 0 MJD with long doubles.
    {
        // Set reference epoch and Julian day for tests.
        const long double referenceEpoch = JULIAN_DAY_AT_0_MJD_LONG;
        const long double julianDay = JULIAN_DAY_AT_0_MJD_LONG + 1.0e6L / 86400.0L;

        // Set expected seconds since epoch result.
        const long double expectedSecondsSinceEpoch = 1.0e6L;

        // Compute seconds since epoch given Julian day and reference epoch.
        const long double computedSecondsSinceEpoch = convertJulianDayToSecondsSinceEpoch< long double >(
                    julianDay, referenceEpoch );

        // Test that computed result matches expected result.
        // Test is run at reduced tolerance, because the final digits of the seconds were lost
        // when converting to Julian day.
#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
        BOOST_CHECK_CLOSE_FRACTION( computedSecondsSinceEpoch, expectedSecondsSinceEpoch,
                                    1.0e-14 );
#endif
    }

    // Test conversion from Julian day to seconds since J2000 epoch  with long doubles
    {
        // Set reference epoch and Julian day for tests.
        const long double referenceEpoch = JULIAN_DAY_ON_J2000_LONG;
        const long double julianDay = JULIAN_DAY_ON_J2000_LONG + 0.5L;

        // Set expected seconds since epoch result.
        long double expectedSecondsSinceEpoch
                = physical_constants::JULIAN_DAY_LONG / 2.0L;

        // Compute seconds since epoch given Julian day and reference epoch.
        const long double computedSecondsSinceEpoch = convertJulianDayToSecondsSinceEpoch< long double >(
                    julianDay, referenceEpoch );

        // Test that computed result matches expected result.
        // Test is run at reduced tolerance, because the final digits of the seconds were lost
        // when converting to Julian day.
        BOOST_CHECK_CLOSE_FRACTION( computedSecondsSinceEpoch, expectedSecondsSinceEpoch,
                                    std::numeric_limits< long double >::epsilon( ) );
    }
}

//! Unit test for seconds to Julian day conversion function.
BOOST_AUTO_TEST_CASE( testSecondsSinceEpochToJulianDayConversions )
{
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

    {
        // Test conversion from seconds since epoch to Julian day with long doubles.

        // Set reference epoch and seconds since epoch for tests.
        const long double referenceEpoch = JULIAN_DAY_AT_0_MJD_LONG;
        const long double secondsSinceEpoch = 1.0e6L;

        // Set expected Julian day result.
        const long double expectedJulianDay = JULIAN_DAY_AT_0_MJD_LONG + 1.0e6L / 86400.0L;

        // Compute Julian day with seconds since reference epoch specified.
        const long double computedJulianDay = convertSecondsSinceEpochToJulianDay< long double >(
                    secondsSinceEpoch, referenceEpoch );

        // Test that computed result matches expected result.
        BOOST_CHECK_CLOSE_FRACTION( computedJulianDay, expectedJulianDay,
                                    std::numeric_limits< long double >::epsilon( ) );
    }
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
        const double computedJulianDay = convertCalendarDateToJulianDay(
                    1858, 11, 17, 0, 0, 0.0 );

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
                                    static_cast< double >( hour ) / 24.0 +
                                    static_cast< double >( minute ) / ( 24.0 * 60.0 ) +
                                    second / ( 24.0 * 60.0 * 60.0 ),
                                    convertCalendarDateToJulianDay( year, month, day, hour, minute, second ),
                                    std::numeric_limits< double >::epsilon( ) );

    }
}

BOOST_AUTO_TEST_CASE( testTimeConversions )
{
    const double testModifiedJulianDay = 54583.87;
    const double testJulianDay = testModifiedJulianDay + JULIAN_DAY_AT_0_MJD;
    {
        // Test conversions to/from Modified Julian day

        BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertModifiedJulianDayToJulianDay( testModifiedJulianDay ),
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay, convertJulianDayToModifiedJulianDay( testJulianDay ),
                                    10.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( convertJulianDayToModifiedJulianDay(
                                        convertModifiedJulianDayToJulianDay( testModifiedJulianDay ) ),
                                    testModifiedJulianDay, 10.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( convertModifiedJulianDayToJulianDay(
                                        convertJulianDayToModifiedJulianDay( testJulianDay ) ),
                                    testJulianDay, 2.0 * std::numeric_limits< double >::epsilon( ) );

        // Test conversions to/from seconds since epoch
        double secondsSinceModifedJulianDayZero = testModifiedJulianDay * JULIAN_DAY;
        BOOST_CHECK_CLOSE_FRACTION( secondsSinceModifedJulianDayZero, convertJulianDayToSecondsSinceEpoch(
                                        testJulianDay, JULIAN_DAY_AT_0_MJD ),
                                    10.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertSecondsSinceEpochToJulianDay(
                                        secondsSinceModifedJulianDayZero, JULIAN_DAY_AT_0_MJD ),
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay /  JULIAN_YEAR_IN_DAYS,
                                    convertSecondsSinceEpochToJulianYearsSinceEpoch( secondsSinceModifedJulianDayZero ),
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay / (  JULIAN_YEAR_IN_DAYS * 100.0 ),
                                    convertSecondsSinceEpochToJulianCenturiesSinceEpoch( secondsSinceModifedJulianDayZero ),
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
    }

    {
        // Test conversion from Julian day to calendar date
        int testYear = 2008;
        int testMonth = 4;
        int testDay = 27;
        boost::gregorian::date testCalendarDate( testYear, testMonth, testDay );
        double testFractionOfDay = 0.37;

        // Check direct conversion from julian day to calendar date
        boost::gregorian::date calendarDate = convertJulianDayToCalendarDate( testJulianDay );
        BOOST_CHECK_EQUAL( testYear, calendarDate.year( ) );
        BOOST_CHECK_EQUAL( testMonth, calendarDate.month( ) );
        BOOST_CHECK_EQUAL( testDay, calendarDate.day( ) );

        // Check indirect conversion from julian day to calendar date
        BOOST_CHECK_EQUAL( getDaysInMonth( 1, testYear ) +
                           getDaysInMonth( 2, testYear ) +
                           getDaysInMonth( 3, testYear ) +
                           testDay, testCalendarDate.day_of_year( ) );

        calendarDate = convertYearAndDaysInYearToDate(
                    testYear, ( getDaysInMonth( 1, testYear ) +
                                getDaysInMonth( 2, testYear ) +
                                getDaysInMonth( 3, testYear ) +
                                testDay ) - 1 ); // Subtract 1 to go to day 0 for first day of year.
        BOOST_CHECK_EQUAL( testYear, calendarDate.year( ) );
        BOOST_CHECK_EQUAL( testMonth, calendarDate.month( ) );
        BOOST_CHECK_EQUAL( testDay, calendarDate.day( ) );

        // Check if noon reference of Julian day is handled properly
        calendarDate = convertJulianDayToCalendarDate( testJulianDay + 0.5 );
        BOOST_CHECK_EQUAL( testYear, calendarDate.year( ) );
        BOOST_CHECK_EQUAL( testMonth, calendarDate.month( ) );
        BOOST_CHECK_EQUAL( testDay + 1, calendarDate.day( ) );

        BOOST_CHECK_CLOSE_FRACTION( calculateJulianDaySinceEpoch(
                                        testCalendarDate, 0.5 + testFractionOfDay, JULIAN_DAY_AT_0_MJD ),
                                    testModifiedJulianDay,
                                    50.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test relativistic time scale conversions (TCG, TT, TDB, TCB) for double precision.
    {
        double secondsSinceJ2000Synchronization = getTimeOfTaiSynchronizationSinceJ2000< double >( );

        // Test whether TCG and TT are both 0 at synchronization time.
        double testTcg = convertTtToTcg( secondsSinceJ2000Synchronization );
        double testTt = convertTcgToTt( secondsSinceJ2000Synchronization );

        BOOST_CHECK_CLOSE_FRACTION( testTcg, secondsSinceJ2000Synchronization,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testTt, secondsSinceJ2000Synchronization,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        // Test whether TDB and TCB are both at defined relative offset at given time.
        double testTcb = convertTdbToTcb( secondsSinceJ2000Synchronization );
        double testTdb = convertTcbToTdb( secondsSinceJ2000Synchronization );

        BOOST_CHECK_CLOSE_FRACTION( testTdb, secondsSinceJ2000Synchronization + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testTcb, secondsSinceJ2000Synchronization - TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        // Test back and forth TCG<->TT, and expected value of TCG at TT=0..
        double testTime = 0.0;

        testTcg = convertTtToTcg( testTime);
        testTt = convertTcgToTt( testTcg );

        double expectedTcg = -secondsSinceJ2000Synchronization * LG_TIME_RATE_TERM /
                ( 1.0 - physical_constants::LG_TIME_RATE_TERM_LONG );
        BOOST_CHECK_CLOSE_FRACTION( testTcg, expectedTcg, 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( testTt, std::numeric_limits< double >::epsilon( ) );

        // Test back and forth TDB<->TCB, and expected value of TDB at TCB=0.
        testTdb = convertTcbToTdb( testTime );
        testTcb = convertTdbToTcb( testTdb );

        double expectedTdb = secondsSinceJ2000Synchronization * LB_TIME_RATE_TERM + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION;

        BOOST_CHECK_CLOSE_FRACTION( testTdb, expectedTdb, 2.0 * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( testTcb, 2.0 * std::numeric_limits< double >::epsilon( ) );

        // Test back and forth TCG<->TT.
        double secondsSinceModifedJulianDayZero = testModifiedJulianDay * JULIAN_DAY;
        testTime = secondsSinceModifedJulianDayZero;

        testTcg = convertTtToTcg< double >( testTime );
        testTt = convertTcgToTt< double >( testTcg );

        BOOST_CHECK_CLOSE_FRACTION( testTt, testTime, 2.0 * std::numeric_limits< double >::epsilon( ) );

        // Test back and forth TDB<->TCB.
        testTdb = convertTcbToTdb< long double >( testTime );
        testTcb = convertTdbToTcb< long double >( testTdb );

        BOOST_CHECK_CLOSE_FRACTION( testTcb, testTime, 2.0 * std::numeric_limits< double >::epsilon( ) );

        secondsSinceModifedJulianDayZero = 1.0E12;

        // Test rate difference between TCG and TT
        double testTt1 = 0.0;
        double testTt2 = secondsSinceModifedJulianDayZero;
        double testTcg1 = convertTtToTcg< double >( testTt1 );
        double testTcg2 = convertTtToTcg< double >( testTt2 );
        double computedLg = 1.0 - ( testTt2 - testTt1 ) / ( testTcg2 - testTcg1 );
        BOOST_CHECK_SMALL( computedLg - physical_constants::LG_TIME_RATE_TERM,
                                    std::numeric_limits< double >::epsilon( ) );

        testTcg1 = 0.0;
        testTcg2 = secondsSinceModifedJulianDayZero;
        testTt1 = convertTcgToTt< double >( testTcg1 );
        testTt2 = convertTcgToTt< double >( testTcg2 );
        computedLg = 1.0 - ( testTt2 - testTt1 ) / ( testTcg2 - testTcg1 );
        BOOST_CHECK_SMALL( computedLg - physical_constants::LG_TIME_RATE_TERM,
                                    std::numeric_limits< double >::epsilon( ) );


        // Test rate difference between TDB and TCB
        double testTcb1 = 0.0;
        double testTcb2 = secondsSinceModifedJulianDayZero;
        double testTdb1 = convertTcbToTdb< double >( testTcb1 );
        double testTdb2 = convertTcbToTdb< double >( testTcb2 );
        double computedLb = 1.0 - ( testTdb2 - testTdb1 ) / ( testTcb2 - testTcb1 );
        BOOST_CHECK_SMALL( computedLb - physical_constants::LB_TIME_RATE_TERM,
                                    std::numeric_limits< double >::epsilon( ) );

        testTdb1 = 0.0;
        testTdb2 = secondsSinceModifedJulianDayZero;
        testTcb1 = convertTdbToTcb< double >( 0.0 );
        testTcb2 = convertTdbToTcb< double >( secondsSinceModifedJulianDayZero );
        computedLb = 1.0 - ( testTdb2 - testTdb1 ) / ( testTcb2 - testTcb1 );
        BOOST_CHECK_SMALL( computedLb - physical_constants::LB_TIME_RATE_TERM,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Compare TT/TCG conversions against Sofa cookbook (only limited precision).
    {
        double expectedTT = physical_constants::JULIAN_DAY *
                convertCalendarDateToJulianDaysSinceEpoch( 2006, 1, 15, 12, 25, 42.68400, JULIAN_DAY_ON_J2000 );
        double expectedTCG = physical_constants::JULIAN_DAY *
                convertCalendarDateToJulianDaysSinceEpoch( 2006, 1, 15, 12, 25, 43.32269, JULIAN_DAY_ON_J2000 );

        double calculatedTT = convertTcgToTt( expectedTCG );
        double calculatedTCG = convertTtToTcg( expectedTT );

        BOOST_CHECK_SMALL( calculatedTCG - expectedTCG, 5.0E-5 );
        BOOST_CHECK_SMALL( calculatedTT - expectedTT, 5.0E-5 );

    }
}
#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )

// Test relativistic time scale conversions (TCG, TT, TDB, TCB) for long double precision.
BOOST_AUTO_TEST_CASE( testTimeConversionsLong )
{
    // Define test dates (arbitrary).
    const long double testModifiedJulianDay = static_cast< long double >( 54583.87 );
    const long double testJulianDay = testModifiedJulianDay + JULIAN_DAY_AT_0_MJD_LONG;

    // Test whether back and forth conversion between JD and MJD provides correct resulats.
    BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertModifiedJulianDayToJulianDay< long double >( testModifiedJulianDay ),
                                2.0 * std::numeric_limits< long double >::epsilon( ));
    BOOST_CHECK_CLOSE_FRACTION( testModifiedJulianDay, convertJulianDayToModifiedJulianDay< long double >( testJulianDay ),
                                2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( convertJulianDayToModifiedJulianDay< long double >(
                                    convertModifiedJulianDayToJulianDay< long double >( testModifiedJulianDay ) ),
                                testModifiedJulianDay, 2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( convertModifiedJulianDayToJulianDay(
                                    convertJulianDayToModifiedJulianDay( testJulianDay ) ),
                                testJulianDay, 2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test conversion to seconds since Epoch for JD and MJD
    long double secondsSinceModifedJulianDayZero = testModifiedJulianDay * JULIAN_DAY_LONG;

    BOOST_CHECK_CLOSE_FRACTION( secondsSinceModifedJulianDayZero, convertJulianDayToSecondsSinceEpoch< long double >(
                                    testJulianDay, JULIAN_DAY_AT_0_MJD_LONG ),
                                2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testJulianDay, convertSecondsSinceEpochToJulianDay< long double >(
                                    secondsSinceModifedJulianDayZero, JULIAN_DAY_AT_0_MJD_LONG ),
                                2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test whether TCG and TT are both 0 at synchronization time.
    long double secondsSinceJ2000Synchronization = getTimeOfTaiSynchronizationSinceJ2000< long double >( );

    long double testTcg = convertTtToTcg< long double >( secondsSinceJ2000Synchronization );
    long double testTt = convertTcgToTt< long double >( secondsSinceJ2000Synchronization );

    BOOST_CHECK_CLOSE_FRACTION( testTcg, secondsSinceJ2000Synchronization,
                                2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTt, secondsSinceJ2000Synchronization,
                                2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test whether TDB and TCB are both at defined relative offset at given time.
    long double testTcb = convertTdbToTcb< long double >( secondsSinceJ2000Synchronization );
    long double testTdb = convertTcbToTdb< long double >( secondsSinceJ2000Synchronization );

    BOOST_CHECK_CLOSE_FRACTION( testTdb, secondsSinceJ2000Synchronization + TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION,
                                2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTcb, secondsSinceJ2000Synchronization - TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION,
                                2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test back and forth TCG<->TT at t=0, and expected value of TCG at TT=0..
    long double testTime = 0.0L;

    testTcg = convertTtToTcg< long double >( testTime);
    testTt = convertTcgToTt< long double >( testTcg );

    long double expectedTcg = -secondsSinceJ2000Synchronization * LG_TIME_RATE_TERM_LONG /
            ( 1.0L - LG_TIME_RATE_TERM_LONG);
    BOOST_CHECK_CLOSE_FRACTION( testTcg, expectedTcg, 2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( testTt, testTime, 2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test back and forth TDB<->TCB at t=0, and expected value of TDB at TCB=0.
    testTdb = convertTcbToTdb< long double >( testTime );
    testTcb = convertTdbToTcb< long double >( testTdb );

    long double expectedTdb = secondsSinceJ2000Synchronization * LB_TIME_RATE_TERM_LONG +
            TDB_SECONDS_OFFSET_AT_SYNCHRONIZATION_LONG;

    BOOST_CHECK_CLOSE_FRACTION( testTdb, expectedTdb, 2.0 * std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_SMALL( testTcb, 2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test back and forth TCG<->TT.
    testTime = secondsSinceModifedJulianDayZero;

    testTcg = convertTtToTcg< long double >( testTime );
    testTt = convertTcgToTt< long double >( testTcg );

    BOOST_CHECK_CLOSE_FRACTION( testTt, testTime, 2.0 * std::numeric_limits< long double >::epsilon( ) );

    // Test back and forth TDB<->TCB.
    testTdb = convertTcbToTdb< long double >( testTime );
    testTcb = convertTdbToTcb< long double >( testTdb );

    BOOST_CHECK_CLOSE_FRACTION( testTcb, testTime, 2.0 * std::numeric_limits< long double >::epsilon( ) );

    secondsSinceModifedJulianDayZero = 1.0E12;

    // Test rate difference between TCG and TT
    long double testTt1 = 0.0;
    long double testTt2 = secondsSinceModifedJulianDayZero;
    long double testTcg1 = convertTtToTcg< long double >( testTt1 );
    long double testTcg2 = convertTtToTcg< long double >( testTt2 );
    long double computedLg = 1.0L - ( testTt2 - testTt1 ) / ( testTcg2 - testTcg1 );
    BOOST_CHECK_SMALL( computedLg - physical_constants::LG_TIME_RATE_TERM_LONG,
                                std::numeric_limits< long double >::epsilon( ) );

    testTcg1 = 0.0;
    testTcg2 = secondsSinceModifedJulianDayZero;
    testTt1 = convertTcgToTt< long double >( testTcg1 );
    testTt2 = convertTcgToTt< long double >( testTcg2 );
    computedLg = 1.0L - ( testTt2 - testTt1 ) / ( testTcg2 - testTcg1 );
    BOOST_CHECK_SMALL( computedLg - physical_constants::LG_TIME_RATE_TERM_LONG,
                                std::numeric_limits< long double >::epsilon( ) );


    // Test rate difference between TDB and TCB
    long double testTcb1 = 0.0L;
    long double testTcb2 = secondsSinceModifedJulianDayZero;
    long double testTdb1 = convertTcbToTdb< long double >( testTcb1 );
    long double testTdb2 = convertTcbToTdb< long double >( testTcb2 );
    long double computedLb = 1.0L - ( testTdb2 - testTdb1 ) / ( testTcb2 - testTcb1 );
    BOOST_CHECK_SMALL( computedLb - physical_constants::LB_TIME_RATE_TERM_LONG,
                                std::numeric_limits< long double >::epsilon( ) );

    testTdb1 = 0.0;
    testTdb2 = secondsSinceModifedJulianDayZero;
    testTcb1 = convertTdbToTcb< long double >( 0.0 );
    testTcb2 = convertTdbToTcb< long double >( secondsSinceModifedJulianDayZero );
    computedLb = 1.0L - ( testTdb2 - testTdb1 ) / ( testTcb2 - testTcb1 );
    BOOST_CHECK_SMALL( computedLb - physical_constants::LB_TIME_RATE_TERM_LONG,
                                std::numeric_limits< long double >::epsilon( ) );

}
#endif

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
