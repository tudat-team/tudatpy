/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Wakker, K. F. (2007), Lecture Notes astro II (Chapter 18), TU Delft course AE4-874,
 *          Delft University of technology, Delft, The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/mission_segments/escapeAndCapture.h"

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
