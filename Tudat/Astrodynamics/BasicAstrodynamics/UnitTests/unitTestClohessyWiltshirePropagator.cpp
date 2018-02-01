/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications website.
 *        URL http://www.smad.com/vallado, 2001. Last accessed: 25-07-2012.
 *      Wakker, K.F. Astrodynamics I, reader of course AE4874 I. Delft University of Technology,
 *        2007.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/clohessyWiltshirePropagator.h"

namespace tudat
{
namespace unit_tests
{

//! Test suite for Clohessy-Wiltshire propagation.
BOOST_AUTO_TEST_SUITE( test_clohessyWiltshirePropagation )

// Testcase 1: propagation of full state in low-Earth orbit.
// This test benchmarks the propagation of a full state against external data.
BOOST_AUTO_TEST_CASE( test_ClohessyWiltshirePropagation_fullState )
{
    // Set central body gravitational parameter [m^3 s^-2].
    // In this case: central body is Earth.
    const double centralBodyGravitationalParameter1 = 3.986004418e14;

    // Set propagation duration [s].
    // In this case: propagation duration is 30 minutes.
    const double propagationDuration1 = 1800.0;

    // Set reference orbit radius [m].
    // In this case: reference orbit is at 400 km altitude above Earth.
    const double referenceOrbitRadius1 = 6.778137e6;

    // Set initial state [m], [m], [m], [m/s], [m,s], [m/s].
    // In this case: arbitrary non-zero distances and velocities.
    const Eigen::Vector6d initialState1 =
            ( Eigen::Vector6d( ) << 45.0, 37.0, 12.0, 0.08,
              0.03, 0.01 ).finished( );

    // Calculate final state according to Tudat function.
    const Eigen::Vector6d computedFinalState1
            = basic_astrodynamics::propagateClohessyWiltshire(
                initialState1,
                propagationDuration1,
                centralBodyGravitationalParameter1,
                referenceOrbitRadius1 );

    // Set final state according to the MATLAB routine "hillsr" from Vallado [2001].
    const Eigen::Vector6d expectedFinalState1 =
            ( Eigen::Vector6d( ) << 3.806450080201250e2,
              -5.437424675454679e2, 2.509547637285142,
              1.541620605755606e-1, -7.294751390499470e-1,
              -1.662099488431618e-2 ).finished( );

    // Check if computed final state matches expected state.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( computedFinalState1, expectedFinalState1, 1.0e-14 );
}

// Testcase 2: pure harmonic relative motion in low-Earth orbit.
//   According to Wakker [2007] pure harmonic relative motion occurs if two conditions are met:
//   (1) initialRadialVelocity = 0.5 * meanAngularMotion * initialAlongTrackPosition.
//   (2) initialAlongTrackVelocity = -2.0 * meanAngularMotion * initialRadialPosition.
//   Therefore the final state after one orbital revolution must be exactly equal to the initial
//   state. This testcase verifies this.
BOOST_AUTO_TEST_CASE( test_ClohessyWiltshirePropagation_harmonicMotion )
{
    // Set central body gravitational parameter [m^3 s^-2].
    // In this case: central body is Earth.
    const double centralBodyGravitationalParameter2 = 3.986004418e14;

    // Set reference orbit radius [m].
    // In this case: reference orbit is at 400 km altitude above Earth.
    const double referenceOrbitRadius2 = 6.778137e6;

    // Calculate mean angular motion of reference orbit [rad].
    const double meanAngularMotion =
            std::sqrt( centralBodyGravitationalParameter2
                       / ( referenceOrbitRadius2 * referenceOrbitRadius2
                           * referenceOrbitRadius2 ) );

    // Set propagation duration [s].
    // In this case: orbital period of the reference orbit.
    const double propagationDuration2 = 2.0 * mathematical_constants::PI
            / meanAngularMotion;

    // Set initial state [m], [m], [m], [m/s], [m,s], [m/s].
    // In this case: arbitrary values for initial positions and initialCrossTrackVelocity. The
    // initialRadialVelocity and initialAlongTrackVelocity are set to achieve harmonic motion.
    const Eigen::Vector6d initialState2 =
            ( Eigen::Vector6d( )
              << 34.0, 49.0 , 17.0,
              0.5 * meanAngularMotion * 49.0,
              -2.0 * meanAngularMotion * 34.0, 0.04 ).finished( );

    // Calculate final state according to Tudat function.
    const Eigen::Vector6d computedFinalState2
            = basic_astrodynamics::propagateClohessyWiltshire(
                initialState2,
                propagationDuration2,
                centralBodyGravitationalParameter2,
                referenceOrbitRadius2 );

    // Check if computed final state matches the initial state.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( initialState2, computedFinalState2, 1.0e-14 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
