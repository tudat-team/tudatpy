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
    Time time( 2, LONG_PI );


    BOOST_CHECK_CLOSE_FRACTION( time.getSeconds< long double >( ), 2.0L * TIME_NORMALIZATION_TERM + LONG_PI, std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( time.getSeconds< double >( ), 2.0 * TIME_NORMALIZATION_TERM + PI, std::numeric_limits< double >::epsilon( ) );

    BOOST_CHECK_CLOSE_FRACTION( static_cast< long double >( time ), 2.0L * TIME_NORMALIZATION_TERM + LONG_PI, std::numeric_limits< long double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION( static_cast< double >( time ),  2.0 * TIME_NORMALIZATION_TERM + PI, std::numeric_limits< double >::epsilon( ) );

    Eigen::Vector3d testVector = ( Eigen::Vector3d(  ) << 4.5, 4.5, 4.5 ).finished( );
    Eigen::Matrix< long double, 3, 1 > testVectorLong = ( Eigen::Matrix< long double, 3, 1 >(  ) << 4.5L, 4.5L, 4.5L ).finished( );

    Eigen::Vector3d multipliedTestVector = time * testVector;
    Eigen::Matrix< long double, 3, 1 > multipliedTestVectorLong = time * testVectorLong;

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( multipliedTestVector( i ), multipliedTestVectorLong( i ),
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( 4.5L * ( 2.0L * TIME_NORMALIZATION_TERM + LONG_PI ), multipliedTestVectorLong( i ),
                                    std::numeric_limits< long double >::epsilon( ) );
    }

    multipliedTestVector = testVector * time;
    multipliedTestVectorLong = testVectorLong * time;

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_CLOSE_FRACTION( multipliedTestVector( i ), multipliedTestVectorLong( i ),
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( 4.5L * ( 2.0L * TIME_NORMALIZATION_TERM + LONG_PI ), multipliedTestVectorLong( i ),
                                    std::numeric_limits< long double >::epsilon( ) );
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

    int numberOfDays1 = 759;
    long double numberOfSeconds1 = 2566.8309405984728595902;
    Time inputTime1( numberOfDays1, numberOfSeconds1 );

    int numberOfDays2 = 2;
    long double numberOfSeconds2 = 1432.48492385475949349;
    Time inputTime2( numberOfDays2, numberOfSeconds2 );

    Time outputTime = inputTime1 + inputTime2;
    BOOST_CHECK_EQUAL( numberOfDays1 + numberOfDays2 + 1, outputTime.getFullPeriods( ) );
    BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                std::numeric_limits< long double >::epsilon( ) );

    outputTime = inputTime1 + numberOfSeconds2;
    BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
    BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                std::numeric_limits< long double >::epsilon( ) );

    outputTime = numberOfSeconds2 + inputTime1;
    BOOST_CHECK_EQUAL( numberOfDays1 + 1, outputTime.getFullPeriods( ) );
    BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds1 + numberOfSeconds2 - TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                std::numeric_limits< long double >::epsilon( ) );

    outputTime = inputTime2 - inputTime1;
    BOOST_CHECK_EQUAL( numberOfDays2 - numberOfDays1 - 1, outputTime.getFullPeriods( ) );
    BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                std::numeric_limits< long double >::epsilon( ) );

    outputTime = inputTime2 - numberOfSeconds1;
    BOOST_CHECK_EQUAL( numberOfDays2 - 1, outputTime.getFullPeriods( ) );
    BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                std::numeric_limits< long double >::epsilon( ) );

    outputTime = numberOfSeconds2 - inputTime1;
    BOOST_CHECK_EQUAL( -numberOfDays1 - 1, outputTime.getFullPeriods( ) );
    BOOST_CHECK_CLOSE_FRACTION( numberOfSeconds2 - numberOfSeconds1 + TIME_NORMALIZATION_TERM, outputTime.getSecondsIntoFullPeriod( ),
                                std::numeric_limits< long double >::epsilon( ) );
}


BOOST_AUTO_TEST_SUITE_END( )


}

}

