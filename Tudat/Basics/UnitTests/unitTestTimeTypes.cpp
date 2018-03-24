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

#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Basics/timeType.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_time_type )

using namespace mathematical_constants;

//! Test if Time objects cast to the expected precision
BOOST_AUTO_TEST_CASE( testTimeBasicCasts )
{

    Time testTime( 2, LONG_PI );

    //Test if Time casts to double/long double at the expected precision
    {

        BOOST_CHECK_CLOSE_FRACTION( testTime.getSeconds< long double >( ), 2.0L * TIME_NORMALIZATION_TERM + LONG_PI,
                                    std::numeric_limits< long double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testTime.getSeconds< double >( ), 2.0 * TIME_NORMALIZATION_TERM + PI,
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( static_cast< long double >( testTime ), 2.0L * TIME_NORMALIZATION_TERM + LONG_PI,
                                    std::numeric_limits< long double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( testTime ),  2.0 * TIME_NORMALIZATION_TERM + PI,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test if Time pre-/post-multiplies Eigen vectors at expected level of precision
    {
        Eigen::Vector3d testVector = ( Eigen::Vector3d(  ) << 4.5, 4.5, 4.5 ).finished( );
        Eigen::Matrix< long double, 3, 1 > testVectorLong = ( Eigen::Matrix< long double, 3, 1 >(  ) << 4.5L, 4.5L, 4.5L ).finished( );

        Eigen::Vector3d multipliedTestVector = testTime * testVector;
        Eigen::Matrix< long double, 3, 1 > multipliedTestVectorLong = testTime * testVectorLong;

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( multipliedTestVector( i ), multipliedTestVectorLong( i ),
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( 4.5L * ( 2.0L * TIME_NORMALIZATION_TERM + LONG_PI ), multipliedTestVectorLong( i ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

        multipliedTestVector = testVector * testTime;
        multipliedTestVectorLong = testVectorLong * testTime;

        for( unsigned int i = 0; i < 3; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( multipliedTestVector( i ), multipliedTestVectorLong( i ),
                                        std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( 4.5L * ( 2.0L * TIME_NORMALIZATION_TERM + LONG_PI ), multipliedTestVectorLong( i ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }
    }
}

//! Test if time saves/retrieves entries at the expected level of precision.
BOOST_AUTO_TEST_CASE( testTimeContentsPrecision )
{
    int numberOfDays = 759;
    long double numberOfSeconds = 2.0L * TIME_NORMALIZATION_TERM + LONG_PI;

    Time testTime( numberOfDays, numberOfSeconds );
    BOOST_CHECK_EQUAL( testTime.getFullPeriods( ), numberOfDays + 2 );
    BOOST_CHECK_CLOSE_FRACTION( testTime.getSecondsIntoFullPeriod( ), LONG_PI,
                                TIME_NORMALIZATION_TERM * std::numeric_limits< long double >::epsilon( ) );

    numberOfSeconds = -2.0L * TIME_NORMALIZATION_TERM + LONG_PI;

    Time testTime2 = Time( numberOfDays, numberOfSeconds );
    BOOST_CHECK_EQUAL( testTime2.getFullPeriods( ), numberOfDays - 2 );
    BOOST_CHECK_CLOSE_FRACTION( testTime2.getSecondsIntoFullPeriod( ), LONG_PI,
                                TIME_NORMALIZATION_TERM * std::numeric_limits< long double >::epsilon( ) );

    Time testTime3 = testTime2;
    BOOST_CHECK_EQUAL( testTime2.getSecondsIntoFullPeriod( ), testTime3.getSecondsIntoFullPeriod( ) );
    BOOST_CHECK_EQUAL( testTime2.getFullPeriods( ), testTime3.getFullPeriods( ) );

    numberOfSeconds = -2.0L * TIME_NORMALIZATION_TERM - LONG_PI;

    Time testTime4 = Time( numberOfDays, numberOfSeconds );
    BOOST_CHECK_EQUAL( testTime4.getFullPeriods( ), numberOfDays - 3 );
    BOOST_CHECK_CLOSE_FRACTION( testTime4.getSecondsIntoFullPeriod( ), TIME_NORMALIZATION_TERM - LONG_PI,
                                TIME_NORMALIZATION_TERM * std::numeric_limits< long double >::epsilon( ) );

}

//! Test basic arithmetic operations of time object
BOOST_AUTO_TEST_CASE( testArithmeticOperations )
{
    {
        // Define Time test values
        int numberOfDays1 = 759;
        long double numberOfSeconds1 = 2566.8309405984728595902;
        Time inputTime1( numberOfDays1, numberOfSeconds1 );

        int numberOfDays2 = 2;
        long double numberOfSeconds2 = 1432.48492385475949349;
        Time inputTime2( numberOfDays2, numberOfSeconds2 );

        Time outputTime;

        // Test Time additions
        {
            outputTime = inputTime1 + inputTime2;
            BOOST_CHECK_EQUAL( numberOfDays1 + numberOfDays2 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
            outputTime = inputTime1;
            outputTime += inputTime2;
            BOOST_CHECK_EQUAL( numberOfDays1 + numberOfDays2 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

        // Test addition between doubles and Time
        {
            outputTime = inputTime1 + numberOfSeconds2;
            BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );

            outputTime = inputTime1;
            outputTime += numberOfSeconds2;
            BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );

            outputTime = numberOfSeconds2 + inputTime1;
            BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

        // Test subtractions of Time objects
        {
            outputTime = inputTime2 - inputTime1;
            BOOST_CHECK_EQUAL( numberOfDays2 - numberOfDays1 - 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
            outputTime = inputTime2;
            outputTime -= inputTime1;
            BOOST_CHECK_EQUAL( numberOfDays2 - numberOfDays1 - 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );

            outputTime = inputTime2 - numberOfSeconds1;
            BOOST_CHECK_EQUAL( numberOfDays2 - 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );

            outputTime = numberOfSeconds2 - inputTime1;
            BOOST_CHECK_EQUAL( -numberOfDays1 - 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM,
                                        outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

    }

    // Test division of Time by double/long double values
    {
        // Define Time test values
        Time testTime( 2, LONG_PI );

        Time dividedTime = testTime / 2.0L;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 1 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI / 2.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );

        dividedTime = testTime / 2.0;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 1 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI / 2.0L,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        dividedTime = testTime / 3.0L;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 0 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI / 3.0L + TIME_NORMALIZATION_TERM * 2.0L / 3.0L,
                                    2.0 *std::numeric_limits< long double >::epsilon( ) );
        dividedTime = testTime / 3.0;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 0 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI / 3.0L + TIME_NORMALIZATION_TERM * 2.0L / 3.0L,
                                    2.0 *std::numeric_limits< double >::epsilon( ) );

        testTime = Time( 9, 8.0 * LONG_PI );
        dividedTime = testTime / 2.0L;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 4 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 4.0L + TIME_NORMALIZATION_TERM / 2.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );
        dividedTime = testTime / 2.0;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 4 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 4.0L + TIME_NORMALIZATION_TERM / 2.0L,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        dividedTime = testTime / 3.0L;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 3 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 8.0L / 3.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );
        dividedTime = testTime / 3.0;
        BOOST_CHECK_EQUAL( dividedTime.getFullPeriods( ), 3 );
        BOOST_CHECK_CLOSE_FRACTION( dividedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 8.0L / 3.0L,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );
    }

    // Test multiplication of Time by double/long double values
    {
        Time testTime( 5, 1200.0L + LONG_PI );
        Time multipliedTime = testTime * 2.0L;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 10 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    2400.0L + LONG_PI * 2.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );

        multipliedTime = testTime * 3.0L;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 16 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 3.0L,
                                    200.0 * std::numeric_limits< long double >::epsilon( ) );

        multipliedTime = 2.0L * testTime;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 10 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    2400.0L + LONG_PI * 2.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );

        multipliedTime = 3.0L * testTime;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 16 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 3.0L,
                                    200.0 * std::numeric_limits< long double >::epsilon( ) );

        multipliedTime = testTime * 2.0;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 10 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    2400.0 + LONG_PI * 2.0,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        multipliedTime = testTime * 3.0;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 16 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    PI * 3.0,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        multipliedTime = 2.0 * testTime;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 10 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    2400.0 + LONG_PI * 2.0,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );

        multipliedTime = 3.0 * testTime;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 16 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    PI * 3.0,
                                    2.0 * std::numeric_limits< double >::epsilon( ) );


    }
}

//! Test if the comparison operators defined for the Time object function correctly
BOOST_AUTO_TEST_CASE( testComparisonOperators )
{
    for( unsigned int i = 0; i < 4; i++ )
    {
        // Define test time (same for each test)
        int numberOfDays = 759;
        long double numberOfSeconds = 2.0L * TIME_NORMALIZATION_TERM + LONG_PI;

        Time testTime( numberOfDays, numberOfSeconds );
        double testTimeDouble = testTime.getSeconds< double >( );
        long double testTimeLongDouble = testTime.getSeconds< long double >( );

        // Define test time (different values for each test)
        Time testTime2 = TUDAT_NAN;
        double testTimeDouble2 = TUDAT_NAN;
        long double testTimeLongDouble2 = TUDAT_NAN;

        // Check comparison if only seconds is larger
        if( i == 0 )
        {
            testTime2 = Time( numberOfDays, numberOfSeconds + 1.0L );
            testTimeDouble2 = testTime.getSeconds< double >( ) + 1.0;
            testTimeLongDouble2 = testTime.getSeconds< long double >( ) + 1.0L;
        }
        // Check comparison if only hours is larger
        else if( i == 1 )
        {
            testTime2 = Time( numberOfDays + 3, numberOfSeconds );
            testTimeDouble2 = testTime.getSeconds< double >( ) + 3.0 * TIME_NORMALIZATION_TERM;
            testTimeLongDouble2 = testTime.getSeconds< long double >( ) + 3.0L * TIME_NORMALIZATION_TERM;
        }
        // Check comparison if hours is larger and seconds is larger
        else if( i == 2 )
        {
            testTime2 = Time( numberOfDays + 3, numberOfSeconds + 1.0L );
            testTimeDouble2 = testTime.getSeconds< double >( ) + 3.0 * TIME_NORMALIZATION_TERM + 1.0;
            testTimeLongDouble2 = testTime.getSeconds< long double >( ) + 3.0L * TIME_NORMALIZATION_TERM + 1.0L;
        }
        // Check comparison if hours is larger and seconds is smaller
        else if( i == 3 )
        {
            testTime2 = Time( numberOfDays + 3, numberOfSeconds - 1.0L );
            testTimeDouble2 = testTime.getSeconds< double >( ) + 3.0 * TIME_NORMALIZATION_TERM  - 1.0;
            testTimeLongDouble2 = testTime.getSeconds< long double >( ) + 3.0L * TIME_NORMALIZATION_TERM - 1.0L;
        }

        // Check equals comparison
        BOOST_CHECK( testTime == testTimeDouble );
        BOOST_CHECK( testTime == testTimeLongDouble );
        BOOST_CHECK( testTime == testTime );

        BOOST_CHECK( testTimeDouble  == testTime );
        BOOST_CHECK( testTimeLongDouble  == testTime );
        BOOST_CHECK( testTimeLongDouble  == testTime );

        // Check not-equals comparison
        BOOST_CHECK( testTime2 != testTimeDouble );
        BOOST_CHECK( testTime2 != testTimeLongDouble );
        BOOST_CHECK( testTime2 != testTime );

        BOOST_CHECK( testTimeDouble  != testTime2 );
        BOOST_CHECK( testTimeLongDouble  != testTime2 );

        // Check (strict) greater/less than
        BOOST_CHECK( testTime2 > testTimeDouble );
        BOOST_CHECK( testTime2 > testTimeLongDouble );
        BOOST_CHECK( testTime2 > testTime );


        BOOST_CHECK( testTimeDouble < testTime2 );
        BOOST_CHECK( testTimeLongDouble  < testTime2 );

        BOOST_CHECK( testTime2 >= testTimeDouble );
        BOOST_CHECK( testTime2 >= testTimeLongDouble );
        BOOST_CHECK( testTime2 >= testTime );

        BOOST_CHECK( testTime >= testTime );

        BOOST_CHECK( testTimeDouble <= testTime2 );
        BOOST_CHECK( testTimeLongDouble  <= testTime2 );

        // Check (strict) greater/less than (opposite direction)
        BOOST_CHECK( testTime < testTimeDouble2 );
        BOOST_CHECK( testTime < testTimeLongDouble2 );
        BOOST_CHECK( testTime < testTime2 );

        BOOST_CHECK( testTimeDouble2 > testTime );
        BOOST_CHECK( testTimeLongDouble2  > testTime );

        BOOST_CHECK( testTime <= testTimeDouble2 );
        BOOST_CHECK( testTime <= testTimeLongDouble2 );
        BOOST_CHECK( testTime <= testTime2 );

        BOOST_CHECK( testTimeDouble2 >= testTime );
        BOOST_CHECK( testTimeLongDouble2  >= testTime );

        // Check if comparison picks up small differences
        long double currentFullSeconds = testTime2.getSecondsIntoFullPeriod( );
        BOOST_CHECK( testTime2 != testTime2 + currentFullSeconds * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( testTime2 < testTime2 + currentFullSeconds * std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK( testTime2 != testTime2 - currentFullSeconds * std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( testTime2 > testTime2 - currentFullSeconds * std::numeric_limits< double >::epsilon( ) );

        long double testTimeLongDouble2Rounded = testTimeLongDouble2 *
                ( 1.0L + 2.0L * std::numeric_limits< long double >::epsilon( ) );
        BOOST_CHECK( testTime2 != testTimeLongDouble2Rounded );
        BOOST_CHECK( testTime2 < testTimeLongDouble2Rounded );
    }
}


BOOST_AUTO_TEST_SUITE_END( )


}

}

