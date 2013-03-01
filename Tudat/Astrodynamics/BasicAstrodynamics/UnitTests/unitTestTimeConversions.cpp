/*    Copyright (c) 2010-2013, Delft University of Technology
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

#include <TudatCore/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

#include <Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h>

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
                = basic_astrodynamics::physical_constants::JULIAN_DAY / 2.0;

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

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
