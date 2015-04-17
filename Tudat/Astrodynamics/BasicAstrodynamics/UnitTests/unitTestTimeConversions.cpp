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

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
