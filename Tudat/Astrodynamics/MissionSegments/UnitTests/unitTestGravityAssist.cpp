/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110117    E. Iorfida        First creation of the code.
 *      110202    J. Melman         Some small changes to layout and
 *                                  added some notes and comments.
 *      110203    E. Iorfida        Changed the code respect to the
 *                                  modification inside the main file.
 *      110208    E. Iorfida        Update file with CartesianVelocityElements
 *                                  pointers as input.
 *      110212    J. Melman         Added comments to clarify test case.
 *      110512    K. Kumar          Updated code to not use dynamic memory
 *                                  allocation and new
 *                                  createPredefinedCelestialBody( ) function.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 *      120416    T. Secretin       Boostified unit test.
 *
 *    References
 *
 */

// Temporary notes (move to class/function doxygen):
// Test runs code and verifies result against expected value.
// If the tested code is erroneous, the test function returns a boolean
// true; if the code is correct, the function returns a boolean false.
// A test should be added together with results from the Lambert targeter
// and the Ephemeris class.
// The expected result for the current test should be calculated inside
// this code, not outside by a calculator.
// 

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Bodies/planet.h"
#include "Tudat/Astrodynamics/MissionSegments/gravityAssist.h"

namespace tudat
{
namespace unit_tests
{

//! Test of gravity assist code.
BOOST_AUTO_TEST_SUITE( test_gravity_assist )

//! Test Delta-V computation.
BOOST_AUTO_TEST_CASE( testDeltaV )
{
    // Tolerances.
    // Expected velocity output is defined with an accuracy of 1 m/s.
    const double velocityTolerance = 1.0;

    // In the first test case, the incoming and outgoing inertial
    // velocities are defined such that the hyperbolic excess velocities
    // are equal. In that way, a delta-V is only needed to rotate the
    // velocity vector, which has been calculated by hand.
    // Expected delta-V for a powered swing-by around Mars.
    const double expectedDeltaV = 3.652e3;

    // Define body that is swung by.
    Planet predefinedMars;
    predefinedMars.setPredefinedPlanetSettings( Planet::mars );

    // Define Sun gravitational parameter.
    const double gravitationalParameterSun = 1.32712440018e20;

    // Define planet-Sun distance.
    const double distanceMarsToSun = unit_conversions::convertAstronomicalUnitsToMeters( 1.5 );

    // Define smallest periapsis distance factor.
    const double marsSmallestPeriapsisDistanceFactor = 1.076;

    // Define planet heliocentric velocity vector.
    // The orbit is considered to be circular.
    Eigen::Vector3d marsVelocity;
    marsVelocity.x( ) = 0.0;
    marsVelocity.y( ) = sqrt( gravitationalParameterSun / distanceMarsToSun );
    marsVelocity.z( ) = 0.0;

    // Define pointer to satellite incoming vector.
    using tudat::mathematics::PI;
    CartesianVelocityElements incomingVelocityTest;
    incomingVelocityTest.setCartesianElementXDot( -25.0e3 * std::sin( PI / 6.0 ) );
    incomingVelocityTest.setCartesianElementYDot( 25.0e3 * std::cos( PI / 6.0 ) );
    incomingVelocityTest.setCartesianElementZDot( 0.0 );

    // Define pointer to satellite outgoing vector.
    CartesianVelocityElements outgoingVelocityTest;
    outgoingVelocityTest.setCartesianElementXDot(
                incomingVelocityTest.getCartesianElementXDot( ) );
    outgoingVelocityTest.setCartesianElementYDot(
                2.0 * marsVelocity.y( ) - incomingVelocityTest.getCartesianElementYDot( ) );
    outgoingVelocityTest.setCartesianElementZDot( 0.0 );

    // Create pointers to new Newton Raphson objects.
    NewtonRaphson myNewtonRaphsonGravityAssist;

    // Declare GravityAssist object.
    using astrodynamics::mission_segments::GravityAssist;
    GravityAssist myGravityAssist( predefinedMars.getGravityFieldModel( ),
                                   3398.0e3 * marsSmallestPeriapsisDistanceFactor,
                                   marsVelocity,
                                   &incomingVelocityTest,
                                   &outgoingVelocityTest,
                                   &myNewtonRaphsonGravityAssist );

    // Compute powered gravity-assist implementation.
    const double deltaV = myGravityAssist.computeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, velocityTolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

