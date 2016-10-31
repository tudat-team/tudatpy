#define BOOST_TEST_MAIN

#include <boost/format.hpp>
#include <boost/test/unit_test.hpp>

#include <iostream>
#include <iomanip>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/timeTypes.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_time_type )

using namespace mathematical_constants;

BOOST_AUTO_TEST_CASE( testDoubleLongDoublePrecisions )
{
    {
        Time testTime( 2, LONG_PI );

        BOOST_CHECK_CLOSE_FRACTION( testTime.getSeconds< long double >( ), 2.0L * TIME_NORMALIZATION_TERM + LONG_PI, std::numeric_limits< long double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( testTime.getSeconds< double >( ), 2.0 * TIME_NORMALIZATION_TERM + PI, std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( static_cast< long double >( testTime ), 2.0L * TIME_NORMALIZATION_TERM + LONG_PI, std::numeric_limits< long double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( testTime ),  2.0 * TIME_NORMALIZATION_TERM + PI, std::numeric_limits< double >::epsilon( ) );

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

BOOST_AUTO_TEST_CASE( testDoubleDaySecondCombinations )
{
    int numberOfDays = 759;
    long double numberOfSeconds = 2.0L * TIME_NORMALIZATION_TERM + LONG_PI;

    Time testTime( numberOfDays, numberOfSeconds );
    BOOST_CHECK_EQUAL( testTime.getFullPeriods( ), numberOfDays + 2 );
    BOOST_CHECK_CLOSE_FRACTION( testTime.getSecondsIntoFullPeriod( ), LONG_PI, TIME_NORMALIZATION_TERM * std::numeric_limits< long double >::epsilon( ) );

    numberOfSeconds = -2.0L * TIME_NORMALIZATION_TERM + LONG_PI;

    Time testTime2 = Time( numberOfDays, numberOfSeconds );
    BOOST_CHECK_EQUAL( testTime2.getFullPeriods( ), numberOfDays - 2 );
    BOOST_CHECK_CLOSE_FRACTION( testTime2.getSecondsIntoFullPeriod( ), LONG_PI, TIME_NORMALIZATION_TERM * std::numeric_limits< long double >::epsilon( ) );

    Time testTime3 = testTime2;
    BOOST_CHECK_EQUAL( testTime2.getSecondsIntoFullPeriod( ), testTime3.getSecondsIntoFullPeriod( ) );
    BOOST_CHECK_EQUAL( testTime2.getFullPeriods( ), testTime3.getFullPeriods( ) );

    numberOfSeconds = -2.0L * TIME_NORMALIZATION_TERM - LONG_PI;

    Time testTime4 = Time( numberOfDays, numberOfSeconds );
    BOOST_CHECK_EQUAL( testTime4.getFullPeriods( ), numberOfDays - 3 );
    BOOST_CHECK_CLOSE_FRACTION( testTime4.getSecondsIntoFullPeriod( ), TIME_NORMALIZATION_TERM - LONG_PI, TIME_NORMALIZATION_TERM * std::numeric_limits< long double >::epsilon( ) );

}

BOOST_AUTO_TEST_CASE( testArithmeticOperations )
{
    {
        int numberOfDays1 = 759;
        long double numberOfSeconds1 = 2566.8309405984728595902;
        Time inputTime1( numberOfDays1, numberOfSeconds1 );

        int numberOfDays2 = 2;
        long double numberOfSeconds2 = 1432.48492385475949349;
        Time inputTime2( numberOfDays2, numberOfSeconds2 );

        Time outputTime;
        {
            outputTime = inputTime1 + inputTime2;
            BOOST_CHECK_EQUAL( numberOfDays1 + numberOfDays2 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
            outputTime = inputTime1;
            outputTime += inputTime2;
            BOOST_CHECK_EQUAL( numberOfDays1 + numberOfDays2 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

        {
            outputTime = inputTime1 + numberOfSeconds2;
            BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );

            outputTime = inputTime1;
            outputTime += numberOfSeconds2;
            BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

        {
            outputTime = numberOfSeconds2 + inputTime1;
            BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }

        {
            outputTime = inputTime2 - inputTime1;
            BOOST_CHECK_EQUAL( numberOfDays2 - numberOfDays1 - 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
            outputTime = inputTime2;
            outputTime -= inputTime1;
            BOOST_CHECK_EQUAL( numberOfDays2 - numberOfDays1 - 1, outputTime.getFullPeriods( ) );
            BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                        std::numeric_limits< long double >::epsilon( ) );
        }
        outputTime = inputTime2 - numberOfSeconds1;
        BOOST_CHECK_EQUAL( numberOfDays2 - 1, outputTime.getFullPeriods( ) );
        BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                    std::numeric_limits< long double >::epsilon( ) );

        outputTime = numberOfSeconds2 - inputTime1;
        BOOST_CHECK_EQUAL( -numberOfDays1 - 1, outputTime.getFullPeriods( ) );
        BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                    std::numeric_limits< long double >::epsilon( ) );

    }

    {
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
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );

        multipliedTime = 2.0L * testTime;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 10 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    2400.0L + LONG_PI * 2.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );

        multipliedTime = 3.0L * testTime;
        BOOST_CHECK_EQUAL( multipliedTime.getFullPeriods( ), 16 );
        BOOST_CHECK_CLOSE_FRACTION( multipliedTime.getSecondsIntoFullPeriod( ),
                                    LONG_PI * 3.0L,
                                    2.0 * std::numeric_limits< long double >::epsilon( ) );

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


BOOST_AUTO_TEST_SUITE_END( )


}

}

