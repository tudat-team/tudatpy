/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/basics/timeType.h"
#include "tudat/astro/basic_astro/dateTime.h"

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_date_time )

using namespace mathematical_constants;
using namespace basic_astrodynamics;

std::vector< int > years = { 1899, 1901, 1969, 1970, 1971, 2999, 3000, 3001, 2023, 2373, 1910, 1621, 1900, 2000, 2004 };
// std::vector< int > years = { 1899, 1901, 1969, 1970, 1971, 2999, 3000, 3001, 2023, 2373, 1910, 1621, 1900, 2000, 2004 };
std::vector< std::pair< int, int > > dates = {
    { 5, 17 }, { 1, 1 }, { 8, 31 }, { 12, 17 }, { 12, 31 }, { 2, 29 },
};
std::vector< std::tuple< int, int, long double > > times = { { 8, 34, 30.234567890123456789L },
                                                             { 11, 34, 30.234567890123456789L },
                                                             { 18, 34, 30.234567890123456789L },
                                                             { 23, 34, 30.234567890123456789L },
                                                             { 0, 0, 0.0L },
                                                             { 12, 0, 0.0L },
                                                             { 11, 59, 60.0L - std::numeric_limits< long double >::epsilon( ) * 3600.0L },
                                                             { 23, 59, 60.0L - std::numeric_limits< long double >::epsilon( ) * 3600.0L },
                                                             { 11, 59, std::numeric_limits< long double >::epsilon( ) * 3600.0L },
                                                             { 23, 59, std::numeric_limits< long double >::epsilon( ) * 3600.0L } };

std::string trimTrailingFractionalZeros( const std::string& input )
{
    const std::size_t decimalPointIndex = input.find( "." );
    if( decimalPointIndex == std::string::npos )
    {
        return input;
    }

    std::size_t lastNonZeroIndex = input.find_last_not_of( '0' );
    if( lastNonZeroIndex == std::string::npos || lastNonZeroIndex < decimalPointIndex )
    {
        return input.substr( 0, decimalPointIndex );
    }
    if( input.at( lastNonZeroIndex ) == '.' )
    {
        return input.substr( 0, lastNonZeroIndex );
    }
    return input.substr( 0, lastNonZeroIndex + 1 );
}

BOOST_AUTO_TEST_CASE( testDateTimeConversions )
{
    for( unsigned int i = 0; i < years.size( ); i++ )
    {
        for( unsigned int j = 0; j < dates.size( ); j++ )
        {
            for( unsigned int k = 0; k < times.size( ); k++ )
            {
                bool exceptionCaught = 0;
                try
                {
                    std::cout << "i = " << i << ", j = " << j << ", k = " << k << std::endl;
                    DateTime currentDateTime( years.at( i ),
                                              dates.at( j ).first,
                                              dates.at( j ).second,
                                              std::get< 0 >( times.at( k ) ),
                                              std::get< 1 >( times.at( k ) ),
                                              std::get< 2 >( times.at( k ) ) );
                    Time currentTime = currentDateTime.epoch< Time >( );

                    // Check that hours into current day are calculated correctly
                    int fullPeriodsUntilStartOfCurrentDay = 24 * ( currentTime.getFullPeriods( ) / 24 );
                    if( currentTime.getFullPeriods( ) < 0 && ( currentTime.getFullPeriods( ) % 24 != 0 ) )
                    {
                        fullPeriodsUntilStartOfCurrentDay -= 24;
                    }
                    if( currentDateTime.getHour( ) >= 12 )
                    {
                        BOOST_CHECK_EQUAL( currentTime.getFullPeriods( ) - fullPeriodsUntilStartOfCurrentDay,
                                           currentDateTime.getHour( ) - 12 );
                    }
                    else
                    {
                        BOOST_CHECK_EQUAL( currentTime.getFullPeriods( ) - fullPeriodsUntilStartOfCurrentDay,
                                           currentDateTime.getHour( ) + 12 );
                    }

                    // Check that seconds into current hour are computed correctly
                    long double expectedSecondsIntoFullPeriod = 60.0L * currentDateTime.getMinute( ) + currentDateTime.getSeconds( );
                    BOOST_CHECK_CLOSE_FRACTION( expectedSecondsIntoFullPeriod,
                                                currentTime.getSecondsIntoFullPeriod( ),
                                                std::numeric_limits< long double >::epsilon( ) * 3600.0L );

                    long double secondsSinceMidnight = static_cast< long double >( currentDateTime.getHour( ) ) * 3600.0L +
                            static_cast< long double >( currentDateTime.getMinute( ) ) * 60.0L + currentDateTime.getSeconds( );

                    BOOST_CHECK_SMALL( std::fabs( secondsSinceMidnight - currentTime.secondsSinceMidnight( ) ),
                                       std::numeric_limits< long double >::epsilon( ) );
                    if( currentDateTime.getHour( ) >= 12 )
                    {
                        // TODO: Check why factor 10 multiplication is needed to get it to pass on Windows
                        BOOST_CHECK_SMALL( std::fabs( ( secondsSinceMidnight - currentTime.secondsSinceNoon( ) - 12.0L * 3600.0L ) ),
                                           10.0 * 3600.0 * std::numeric_limits< long double >::epsilon( ) );
                    }
                    else
                    {
                        BOOST_CHECK_SMALL( std::fabs( secondsSinceMidnight + 12.0L * 3600.0L - currentTime.secondsSinceNoon( ) ),
                                           3600.0 * std::numeric_limits< long double >::epsilon( ) );
                    }

                    DateTime reconstructedDateTime = DateTime::fromTime( currentTime );

                    BOOST_CHECK_EQUAL( reconstructedDateTime.getYear( ), currentDateTime.getYear( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getMonth( ), currentDateTime.getMonth( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getDay( ), currentDateTime.getDay( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getHour( ), currentDateTime.getHour( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getMinute( ), currentDateTime.getMinute( ) );
                    BOOST_CHECK_SMALL( std::fabs( reconstructedDateTime.getSeconds( ) - currentDateTime.getSeconds( ) ),
                                       std::numeric_limits< long double >::epsilon( ) * 3600.0L );

                    long double currentJulianDay = julianDayFromTime< long double >( currentTime );
                    Time reconstructedTime = timeFromJulianDay< long double >( currentJulianDay );
                    double timeTolerance = 3.0 * currentJulianDay * 86400.0 * std::numeric_limits< long double >::epsilon( );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( reconstructedTime - currentTime ) ), timeTolerance );

                    long double julianDayFromDateTime = currentDateTime.julianDay< long double >( );
                    double julianDayTolerance = 3.0 * currentJulianDay * std::numeric_limits< long double >::epsilon( );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( julianDayFromDateTime - currentJulianDay ) ), julianDayTolerance );
                    reconstructedDateTime = DateTime::fromJulianDay( currentJulianDay );
                    long double julianDayFromReconstructedDateTime = reconstructedDateTime.julianDay< long double >( );
                    double reconstructedJdTolerance = std::fabs( 3.0 * currentJulianDay * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( julianDayFromReconstructedDateTime - currentJulianDay ) ),
                                       reconstructedJdTolerance );

                    long double currentModifiedJulianDay = modifiedJulianDayFromTime< long double >( currentTime );
                    reconstructedTime = timeFromModifiedJulianDay< long double >( currentModifiedJulianDay );
                    timeTolerance = std::fabs( 3.0 * currentModifiedJulianDay * 86400.0 * std::numeric_limits< long double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( reconstructedTime - currentTime ) ), timeTolerance );
                    reconstructedDateTime = DateTime::fromModifiedJulianDay( currentModifiedJulianDay );
                    long double mjdFromReconstructedDateTime = reconstructedDateTime.modifiedJulianDay< long double >( );
                    double reconstructedMjdTolerance =
                            std::fabs( 3.0 * currentModifiedJulianDay * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( mjdFromReconstructedDateTime - currentModifiedJulianDay ) ),
                                       reconstructedMjdTolerance );

                    long double modifiedJulianDayFromDateTime = currentDateTime.modifiedJulianDay< long double >( );
                    double modifiedJulianDayTolerance = 3.0 * currentJulianDay * std::numeric_limits< long double >::epsilon( );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( modifiedJulianDayFromDateTime - currentModifiedJulianDay ) ),
                                       modifiedJulianDayTolerance );
                }
                catch( std::runtime_error& caughtException )
                {
                    std::cout << "Exception " << caughtException.what( ) << std::endl;
                    exceptionCaught = true;
                }

                if( j == 5 && i < 13 )
                {
                    BOOST_CHECK_EQUAL( exceptionCaught, true );
                }
                else
                {
                    BOOST_CHECK_EQUAL( exceptionCaught, false );
                }
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( testLeapSecondReconstructionFromTime )
{
    DateTime leapSecondDateTime( 2016, 12, 31, 23, 59, 60.0L );
    DateTime reconstructedLeapSecondDateTime = DateTime::fromTime( leapSecondDateTime.epoch< Time >( ) );

    BOOST_CHECK_EQUAL( reconstructedLeapSecondDateTime.getYear( ), leapSecondDateTime.getYear( ) );
    BOOST_CHECK_EQUAL( reconstructedLeapSecondDateTime.getMonth( ), leapSecondDateTime.getMonth( ) );
    BOOST_CHECK_EQUAL( reconstructedLeapSecondDateTime.getDay( ), leapSecondDateTime.getDay( ) );
    BOOST_CHECK_EQUAL( reconstructedLeapSecondDateTime.getHour( ), leapSecondDateTime.getHour( ) );
    BOOST_CHECK_EQUAL( reconstructedLeapSecondDateTime.getMinute( ), leapSecondDateTime.getMinute( ) );
    BOOST_CHECK_SMALL( std::fabs( reconstructedLeapSecondDateTime.getSeconds( ) - leapSecondDateTime.getSeconds( ) ),
                       std::numeric_limits< long double >::epsilon( ) * 3600.0L );

    DateTime fractionalLeapSecondDateTime( 2016, 12, 31, 23, 59, 60.5L );
    DateTime reconstructedFractionalLeapSecondDateTime = DateTime::fromTime( fractionalLeapSecondDateTime.epoch< Time >( ) );

    BOOST_CHECK_EQUAL( reconstructedFractionalLeapSecondDateTime.getYear( ), fractionalLeapSecondDateTime.getYear( ) );
    BOOST_CHECK_EQUAL( reconstructedFractionalLeapSecondDateTime.getMonth( ), fractionalLeapSecondDateTime.getMonth( ) );
    BOOST_CHECK_EQUAL( reconstructedFractionalLeapSecondDateTime.getDay( ), fractionalLeapSecondDateTime.getDay( ) );
    BOOST_CHECK_EQUAL( reconstructedFractionalLeapSecondDateTime.getHour( ), fractionalLeapSecondDateTime.getHour( ) );
    BOOST_CHECK_EQUAL( reconstructedFractionalLeapSecondDateTime.getMinute( ), fractionalLeapSecondDateTime.getMinute( ) );
    BOOST_CHECK_SMALL( std::fabs( reconstructedFractionalLeapSecondDateTime.getSeconds( ) - fractionalLeapSecondDateTime.getSeconds( ) ),
                       std::numeric_limits< long double >::epsilon( ) * 3600.0L );
}

BOOST_AUTO_TEST_CASE( testLeapSecondValidation )
{
    BOOST_CHECK_NO_THROW( DateTime( 2016, 12, 31, 23, 59, 60.0L ) );
    BOOST_CHECK_NO_THROW( DateTime( 2016, 12, 31, 23, 59, 60.5L ) );
    BOOST_CHECK_THROW( DateTime( 2016, 12, 31, 12, 0, 60.0L ), std::runtime_error );
    BOOST_CHECK_THROW( DateTime( 2016, 12, 31, 12, 0, 60.5L ), std::runtime_error );
    BOOST_CHECK_THROW( DateTime( 2017, 1, 1, 23, 59, 60.0L ), std::runtime_error );
    BOOST_CHECK_THROW( DateTime( 2016, 12, 31, 23, 59, 61.0L ), std::runtime_error );
}

BOOST_AUTO_TEST_CASE( testAddSecondsDuringLeapSecond )
{
    const long double epsilon = std::numeric_limits< long double >::epsilon( );
    const long double tolerance = 3600.0L * epsilon;

    DateTime beforeLeapSecond( 2016, 12, 31, 23, 59, 59.0L );
    DateTime leapSecond = beforeLeapSecond.addSecondsToDateTime( 1.0L );

    BOOST_CHECK_EQUAL( leapSecond.getYear( ), 2016 );
    BOOST_CHECK_EQUAL( leapSecond.getMonth( ), 12 );
    BOOST_CHECK_EQUAL( leapSecond.getDay( ), 31 );
    BOOST_CHECK_EQUAL( leapSecond.getHour( ), 23 );
    BOOST_CHECK_EQUAL( leapSecond.getMinute( ), 59 );
    BOOST_CHECK_SMALL( std::fabs( leapSecond.getSeconds( ) - 60.0L ), tolerance );

    DateTime afterLeapSecond = leapSecond.addSecondsToDateTime( 1.0L );
    BOOST_CHECK_EQUAL( afterLeapSecond.getYear( ), 2017 );
    BOOST_CHECK_EQUAL( afterLeapSecond.getMonth( ), 1 );
    BOOST_CHECK_EQUAL( afterLeapSecond.getDay( ), 1 );
    BOOST_CHECK_EQUAL( afterLeapSecond.getHour( ), 0 );
    BOOST_CHECK_EQUAL( afterLeapSecond.getMinute( ), 0 );
    BOOST_CHECK_SMALL( std::fabs( afterLeapSecond.getSeconds( ) - 0.0L ), tolerance );

    DateTime almostOneSecondAfterLeapSecond = leapSecond.addSecondsToDateTime( 1.0L - epsilon );
    BOOST_CHECK_EQUAL( almostOneSecondAfterLeapSecond.getYear( ), 2016 );
    BOOST_CHECK_EQUAL( almostOneSecondAfterLeapSecond.getMonth( ), 12 );
    BOOST_CHECK_EQUAL( almostOneSecondAfterLeapSecond.getDay( ), 31 );
    BOOST_CHECK_EQUAL( almostOneSecondAfterLeapSecond.getHour( ), 23 );
    BOOST_CHECK_EQUAL( almostOneSecondAfterLeapSecond.getMinute( ), 59 );
    BOOST_CHECK_GE( almostOneSecondAfterLeapSecond.getSeconds( ), 60.0L );
    BOOST_CHECK_LT( almostOneSecondAfterLeapSecond.getSeconds( ), 61.0L );

    DateTime midnightNextDay( 2017, 1, 1, 0, 0, 0.0L );
    DateTime epsilonBeforeMidnightNextDay = midnightNextDay.addSecondsToDateTime( -epsilon );
    BOOST_CHECK_EQUAL( epsilonBeforeMidnightNextDay.getYear( ), 2016 );
    BOOST_CHECK_EQUAL( epsilonBeforeMidnightNextDay.getMonth( ), 12 );
    BOOST_CHECK_EQUAL( epsilonBeforeMidnightNextDay.getDay( ), 31 );
    BOOST_CHECK_EQUAL( epsilonBeforeMidnightNextDay.getHour( ), 23 );
    BOOST_CHECK_EQUAL( epsilonBeforeMidnightNextDay.getMinute( ), 59 );
    BOOST_CHECK_GE( epsilonBeforeMidnightNextDay.getSeconds( ), 60.0L );
    BOOST_CHECK_LT( epsilonBeforeMidnightNextDay.getSeconds( ), 61.0L );
}

BOOST_AUTO_TEST_CASE( testDateTimeStringRepresentation )
{
    double timeValue = 0.99999999;
    for( int i = 0; i < 15; i++ )
    {
        DateTime date1 = DateTime( 2008, 4, 1, 4, 6, 59 + timeValue );
        std::string dateString1 = date1.isoString( false, i );

        DateTime date2 = DateTime( 1999, 12, 31, 23, 59, 59 + timeValue );
        std::string dateString2 = date2.isoString( false, i );
        if( i >= 8 )
        {
            std::string testString1 = "2008-04-01 04:06:59.99999999";
            std::string testString2 = "1999-12-31 23:59:59.99999999";

            BOOST_CHECK_EQUAL( testString1, dateString1.substr( 0, testString1.length( ) ) );
            BOOST_CHECK_EQUAL( testString2, dateString2.substr( 0, testString1.length( ) ) );
        }
        else
        {
            std::string testString1 = "2008-04-01 04:07:00";
            std::string testString2 = "2000-01-01 00:00:00";

            BOOST_CHECK_EQUAL( testString1, dateString1.substr( 0, testString1.length( ) ) );
            BOOST_CHECK_EQUAL( testString2, dateString2.substr( 0, testString1.length( ) ) );
        }
    }
}

BOOST_AUTO_TEST_CASE( testIsoInitialization )
{
    std::cout << "Testing ISO initialization" << std::endl;
    std::vector< std::string > testStrings = { "2023-06-20T00:05:23.2817658340294",
                                               "2020-02-29T23:59:59.9999999999998",
                                               "2000-01-01T12:00:00.0000000000000",
                                               "1753-08-09T22:34:10.7295231830836" };

    std::vector< double > julianDays = { 2460116., 2458909., 2451545., 2361551 };
    for( unsigned int i = 0; i < testStrings.size( ); i++ )
    {
        std::cout << "Testing string i = " << i << ": " << testStrings.at( i ) << std::endl;
        DateTime dateTime = DateTime::fromIsoString( testStrings.at( i ) );
        std::string reconstuctedString = dateTime.isoString( true, 17 );

        if( sizeof( long double ) > 8 )
        {
            // Compare canonicalized representations, as ISO formatting may pad trailing zeros.
            BOOST_CHECK_EQUAL( trimTrailingFractionalZeros( testStrings.at( i ) ), trimTrailingFractionalZeros( reconstuctedString ) );
        }
        Time time = timeFromIsoString< Time >( testStrings.at( i ) );
        BOOST_CHECK_SMALL( static_cast< long double >( time - dateTime.epoch< Time >( ) ),
                           ( 3600.0L * std::numeric_limits< long double >::epsilon( ) ) );

        BOOST_CHECK_SMALL( std::fabs( dateTime.julianDay< double >( ) - julianDays.at( i ) ), 0.5 );
        BOOST_CHECK_SMALL( std::fabs( dateTime.julianDay< double >( ) - dateTime.modifiedJulianDay< double >( ) ) - 2400000.5,
                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_CASE( testIsoInitializationLeapSecond )
{
    const std::string leapSecondIsoString = "2016-12-31T23:59:60.5";
    DateTime leapSecondFromIso = DateTime::fromIsoString( leapSecondIsoString );
    DateTime referenceLeapSecond( 2016, 12, 31, 23, 59, 60.5L );

    BOOST_CHECK_EQUAL( leapSecondFromIso.getYear( ), 2016 );
    BOOST_CHECK_EQUAL( leapSecondFromIso.getMonth( ), 12 );
    BOOST_CHECK_EQUAL( leapSecondFromIso.getDay( ), 31 );
    BOOST_CHECK_EQUAL( leapSecondFromIso.getHour( ), 23 );
    BOOST_CHECK_EQUAL( leapSecondFromIso.getMinute( ), 59 );
    BOOST_CHECK_SMALL( std::fabs( leapSecondFromIso.getSeconds( ) - 60.5L ), std::numeric_limits< long double >::epsilon( ) * 3600.0L );
    BOOST_CHECK_SMALL( std::fabs( static_cast< long double >( leapSecondFromIso.epoch< Time >( ) - referenceLeapSecond.epoch< Time >( ) ) ),
                       std::numeric_limits< long double >::epsilon( ) * 3600.0L );
}

BOOST_AUTO_TEST_CASE( testDateTimeDayInYearConversions )
{
    {
        std::cout << "Testing DateTime day in year conversions" << std::endl;
        // Test conversion from Julian day to calendar date
        // same test as in testTimeConversions, but using DateTime class
        int testYear = 2008;
        int testMonth = 4;
        int testDay = 27;
        int testHour = 0;
        int testMinute = 0;
        long double testSeconds = 0.0;
        DateTime testDateTime( testYear, testMonth, testDay, testHour, testMinute, testSeconds );

        // day of year returns 0 for first day of year, thus we need to subtract 1 to get the day of year
        BOOST_CHECK_EQUAL( getDaysInMonth( 1, testYear ) + getDaysInMonth( 2, testYear ) + getDaysInMonth( 3, testYear ) + testDay - 1,
                           testDateTime.dayOfYear( ) );

        DateTime constructedDateTime = DateTime::fromYearAndDaysInYear(
                testYear,
                ( getDaysInMonth( 1, testYear ) + getDaysInMonth( 2, testYear ) + getDaysInMonth( 3, testYear ) + testDay -
                  1 ) );  // Subtract 1 to go to day 0 for first day of year.
        BOOST_CHECK_EQUAL( testYear, constructedDateTime.getYear( ) );
        BOOST_CHECK_EQUAL( testMonth, constructedDateTime.getMonth( ) );
        BOOST_CHECK_EQUAL( testDay, constructedDateTime.getDay( ) );
        BOOST_CHECK_EQUAL( testHour, constructedDateTime.getHour( ) );
        BOOST_CHECK_EQUAL( testMinute, constructedDateTime.getMinute( ) );
        BOOST_CHECK_SMALL( std::fabs( static_cast< double >( testSeconds - constructedDateTime.getSeconds( ) ) ),
                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
