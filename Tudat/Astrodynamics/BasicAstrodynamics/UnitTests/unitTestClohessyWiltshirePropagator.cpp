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
 *      120725    E. Dekens         File created.
 *
 *    References
 *      Vallado, D.A. Fundamentals of Astrodynamics and Applications website.
 *        URL http://www.smad.com/vallado, 2001. Last accessed: 25-07-2012.
 *      Wakker, K.F. Astrodynamics I, reader of course AE4874 I. Delft University of Technology,
 *        2007.
 *
 *    Notes
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

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
    const Eigen::VectorXd initialState1 = ( Eigen::VectorXd( 6 ) << 45.0, 37.0, 12.0, 0.08,
                                                                    0.03, 0.01 ).finished( );

    // Calculate final state according to Tudat function.
    const Eigen::VectorXd computedFinalState1
            = basic_astrodynamics::propagateClohessyWiltshire(
                initialState1,
                propagationDuration1,
                centralBodyGravitationalParameter1,
                referenceOrbitRadius1 );

    // Set final state according to the MATLAB routine "hillsr" from Vallado [2001].
    const Eigen::VectorXd expectedFinalState1 = ( Eigen::VectorXd( 6 ) << 3.806450080201250e2,
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
    const double propagationDuration2 = 2.0 * basic_mathematics::mathematical_constants::PI
            / meanAngularMotion;

    // Set initial state [m], [m], [m], [m/s], [m,s], [m/s].
    // In this case: arbitrary values for initial positions and initialCrossTrackVelocity. The
    // initialRadialVelocity and initialAlongTrackVelocity are set to achieve harmonic motion.
    const Eigen::VectorXd initialState2 = ( Eigen::VectorXd( 6 )
                                            << 34.0, 49.0 , 17.0,
                                            0.5 * meanAngularMotion * 49.0,
                                            -2.0 * meanAngularMotion * 34.0, 0.04 ).finished( );

    // Calculate final state according to Tudat function.
    const Eigen::VectorXd computedFinalState2
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
