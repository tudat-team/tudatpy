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

BOOST_AUTO_TEST_CASE( testDateTimeConversions )
{
    std::vector< int > years = { 2023, 2373, 1910, 1621, 1900, 2000, 2004 };
    std::vector< std::pair< int, int > > dates = {
        { 5, 17 },
        { 1, 1 },
        { 8, 31 },
        { 12, 17 },
        { 12, 31 },
        { 2, 29 },
    };
    std::vector< std::tuple< int, int, long double > > times = {
        { 8, 34, 30.234567890123456789L },
        { 11, 34, 30.234567890123456789L },
        { 18, 34, 30.234567890123456789L },
        { 23, 34, 30.234567890123456789L },
        { 0, 0, 0.0L },
        { 12, 0, 0.0L },
        { 11, 59, 60.0L - std::numeric_limits< long double >::epsilon( ) * 3600.0L   },
        { 23, 59, 60.0L - std::numeric_limits< long double >::epsilon( ) * 3600.0L   },
        { 11, 59, std::numeric_limits< long double >::epsilon( ) * 3600.0L   },
        { 23, 59, std::numeric_limits< long double >::epsilon( ) * 3600.0L   }
    };

    for( unsigned int i = 0; i < years.size( ); i++ )
    {
        for( unsigned int j = 0; j < dates.size( ); j++ )
        {
            for ( unsigned int k = 0; k < times.size( ); k++ )
            {
                bool exceptionCaught = 0;
                try
                {
                    std::cout << i << " " << j << " " << k << std::endl;
                    DateTime currentDateTime( years.at( i ), dates.at( j ).first, dates.at( j ).second,
                                              std::get<0>( times.at( k )), std::get<1>( times.at( k )),
                                              std::get<2>( times.at( k )));
                    Time currentTime = currentDateTime.epoch< Time >(  );

                    // Check that hours into current day are calculated correctly
                    int fullPeriodsUntilStartOfCurrentDay = 24 * ( currentTime.getFullPeriods( ) / 24 );
                    if ( currentTime.getFullPeriods( ) < 0 && ( currentTime.getFullPeriods( ) % 24 != 0 ))
                    {
                        fullPeriodsUntilStartOfCurrentDay -= 24;
                    }
                    if ( currentDateTime.getHour( ) >= 12 )
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
                    long double expectedSecondsIntoFullPeriod =
                        60.0L * currentDateTime.getMinute( ) + currentDateTime.getSeconds( );
                    BOOST_CHECK_CLOSE_FRACTION( expectedSecondsIntoFullPeriod, currentTime.getSecondsIntoFullPeriod( ),
                                                std::numeric_limits<long double>::epsilon( ) * 3600.0L );

                    long double secondsSinceMidnight =
                        static_cast< long double >( currentDateTime.getHour( ) ) * 3600.0L +
                        static_cast< long double >( currentDateTime.getMinute( ) ) * 60.0L +
                        currentDateTime.getSeconds( );

                    BOOST_CHECK_SMALL( std::fabs( secondsSinceMidnight - currentTime.secondsSinceMidnight( )),
                                       std::numeric_limits<long double>::epsilon( ));
                    if ( currentDateTime.getHour( ) >= 12 )
                    {
                        // TODO: Check why factor 10 multiplication is needed to get it to pass on Windows
                        BOOST_CHECK_SMALL(
                            std::fabs(( secondsSinceMidnight - currentTime.secondsSinceNoon( ) - 12.0L * 3600.0L )),
                            10.0 * 3600.0 * std::numeric_limits<long double>::epsilon( ));
                    }
                    else
                    {
                        BOOST_CHECK_SMALL(
                            std::fabs( secondsSinceMidnight + 12.0L * 3600.0L - currentTime.secondsSinceNoon( )),
                            3600.0 * std::numeric_limits<long double>::epsilon( ));
                    }

                    DateTime reconstructedDateTime = getCalendarDateFromTime( currentTime );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getYear( ), currentDateTime.getYear( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getMonth( ), currentDateTime.getMonth( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getDay( ), currentDateTime.getDay( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getHour( ), currentDateTime.getHour( ) );
                    BOOST_CHECK_EQUAL( reconstructedDateTime.getMinute( ), currentDateTime.getMinute( ) );
                    BOOST_CHECK_SMALL( std::fabs( reconstructedDateTime.getSeconds( ) - currentDateTime.getSeconds( ) ),
                                       std::numeric_limits<long double>::epsilon( ) * 3600.0L );

                    long double currentJulianDay = julianDayFromTime<long double>( currentTime );
                    Time reconstructedTime = timeFromJulianDay< long double >( currentJulianDay );
                    double timeTolerance = 3.0 * currentJulianDay * 86400.0 * std::numeric_limits<long double>::epsilon( );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( reconstructedTime - currentTime )),
                                       timeTolerance );

                    long double julianDayFromDateTime = currentDateTime.julianDay< long double >( );
                    double julianDayTolerance = 3.0 * currentJulianDay * std::numeric_limits<long double>::epsilon( );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( julianDayFromDateTime - currentJulianDay )),
                                       julianDayTolerance );

                    long double currentModifiedJulianDay = modifiedJulianDayFromTime<long double>( currentTime );
                    reconstructedTime = timeFromModifiedJulianDay< long double >( currentModifiedJulianDay );
                    timeTolerance = std::fabs( 3.0 * currentModifiedJulianDay * 86400.0 * std::numeric_limits<long double>::epsilon( ));
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( reconstructedTime - currentTime )),
                                       timeTolerance );

                    long double modifiedJulianDayFromDateTime = currentDateTime.modifiedJulianDay< long double >( );
                    double modifiedJulianDayTolerance = 3.0 * currentJulianDay * std::numeric_limits<long double>::epsilon( );
                    BOOST_CHECK_SMALL( std::fabs( static_cast< double >( modifiedJulianDayFromDateTime - currentModifiedJulianDay )),
                                       modifiedJulianDayTolerance );

                }
                catch( ... )
                {
                    exceptionCaught = true;
                }

                if( j == 5 && i < 5 )
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

BOOST_AUTO_TEST_CASE( testIsoInitialization )
{
    std::vector< std::string > testStrings =
        { "2023-06-20T00:05:23.28176583402943837",
          "2020-02-29T23:59:59.99999999999999998",
          "2000-01-01T12:00:00.00000000000000000",
          "1753-08-09T22:34:10.72952318308363849"  };

    std::vector< double > julianDays =
        { 2460116.,
          2458909.,
          2451545.,
          2361551};
    for( unsigned int i = 0; i < testStrings.size( ); i++ )
    {
        DateTime dateTime = dateTimeFromIsoString( testStrings.at( i ) );
        std::string reconstuctedString = dateTime.isoString( true, 17 );

        // TODO: fix test for long doubles with 64-bit precision
        if (sizeof(long double) > 8)
        {
            BOOST_CHECK_EQUAL( testStrings.at( i ), reconstuctedString );
        }
        Time time = timeFromIsoString<Time>( testStrings.at( i ) );
        BOOST_CHECK_SMALL( static_cast< long double >( time - dateTime.epoch<Time>( )),
                           ( 3600.0L * std::numeric_limits<long double>::epsilon( )));

        BOOST_CHECK_SMALL( std::fabs( dateTime.julianDay< double >( ) - julianDays.at( i ) ), 0.5 );
        BOOST_CHECK_SMALL( std::fabs( dateTime.julianDay< double >( ) - dateTime.modifiedJulianDay< double >( ) ) - 2400000.5,
                           std::numeric_limits< double >::epsilon( ) );

    }

}
BOOST_AUTO_TEST_SUITE_END( )

}

}

