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
 *      110131    E. Iorfida        First creation of the code.
 *      110202    E. Iorfida        Modified structure of the code, unique
 *                                  base class for launch and capture paths.
 *      110208    E. Iorfida        Deleted execute( ) function. Modified getDeltaV into
 *                                  computeDeltaV.
 *      110214    E. Iorfida        Code updated with the modification made in .h/.cpp files
 *                                  about radius of central body.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 *      120416    T. Secretin       Boostified unit test.
 *      120531    P. Musegaas       Complete update of unit test after changing code completely.
 *                                  Changed examples to more accurate ones, verifyable in reader.
 *                                  Added hand calculation test for non-circular parking orbit.
 *      120625    P. Musegaas       Minor changes.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes Astrodynamics II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Astrodynamics/MissionSegments/escapeAndCapture.h"

namespace tudat
{
namespace unit_tests
{

//! Test patched conics implementation.
BOOST_AUTO_TEST_SUITE( test_escape_and_capture )

//! Test delta-V computation for the escape phase.
BOOST_AUTO_TEST_CASE( testDeltaVEscape )
{
    // Set tolerance.
    const double tolerance = 1.0e-15;

    // Expected test result based on Table 18.2 of [Wakker, 2007]. Escape velocity to go to Mars
    // starting from a parking orbit around Earth at 185 km altitude. Result calculated with more
    // precision separately using Equation 18.25 of [Wakker, 2007].
    const double expectedDeltaVEscape = 3614.64460281887;

    // Set required parameters.
    const double gravitationalParameterEarth = 3.986e14;
    const double semiMajorAxis = 6.378e6 + 1.85e5;
    const double eccentricity = 0.0;
    const double excessVelocity = 2944.61246668719;

    // Compute delta-V of escape phase.
    const double deltaVEscape = mission_segments::computeEscapeOrCaptureDeltaV(
                gravitationalParameterEarth, semiMajorAxis, eccentricity, excessVelocity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaVEscape, deltaVEscape, tolerance );
}

//! Test delta-V computation for the capture phase.
BOOST_AUTO_TEST_CASE( testDeltaVCapture )
{
    // Set tolerance.
    const double tolerance = 1.0e-15;

    // Expected test result based on Table 18.2 of [Wakker, 2007]. Capture velocity to arrive at
    // a parking orbit around Mars of 1.1 times the planetary radius starting from Earth. Result
    // calculated with more precision separately using Equation 18.25 of [Wakker, 2007].
    const double expectedDeltaVCapture = 2087.1062716740798;

    // Set required parameters.
    const double gravitationalParameterMars = 4.2830e13;
    const double semiMajorAxis = 1.1 * 3.3895e6;
    const double eccentricity = 0.0;
    const double excessVelocity = 2648.83359973278;

    // Compute delta-V of escape phase.
    const double deltaVCapture = mission_segments::computeEscapeOrCaptureDeltaV(
                gravitationalParameterMars, semiMajorAxis, eccentricity, excessVelocity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaVCapture, deltaVCapture, tolerance );
}

//! Test delta-V computation for an escape from an eccentric parking orbit.
BOOST_AUTO_TEST_CASE( testDeltaVEscapeElliptical )
{
    // Set tolerance.
    const double tolerance = 1.0e-15;

    // Departure from 200 km perigee, 1800 km apogee parking orbit from Earth to Mars transfer.
    // Result calculated to high precision separately using Equation 18.25 of [Wakker, 2007].
    const double expectedDeltaVEscape = 3200.2178506657729;

    // Set required parameters.
    const double gravitationalParameterEarth = 3.986e14;
    const double semiMajorAxis = 6.378e6 + 1.0e6;
    const double eccentricity = 0.108430468961778;
    const double excessVelocity = 2944.61246668719;

    // Compute delta-V of escape phase.
    const double deltaVEscape = mission_segments::computeEscapeOrCaptureDeltaV(
                gravitationalParameterEarth, semiMajorAxis, eccentricity, excessVelocity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaVEscape, deltaVEscape, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
