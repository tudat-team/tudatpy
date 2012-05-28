/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110207    B. Romgens        File created.
 *      110215    K. Kumar          Minor modifications to layout, comments
 *                                  and variable-naming.
 *      110411    K. Kumar          Added unit test for
 *                                  convertCartesianToSpherical() function.
 *      110701    K. Kumar          Updated failing tests with relative errors.
 *      110708    K. Kumar          Added unit tests for computeSampleMean()
 *                                  and computeSampleVariance() functions.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      111111    K. Kumar          Strange error with convertCylindricalToCartesian function;
 *                                  achieved precision of results is less than machine precision,
 *                                  fixed by using slightly larger precision tolerance.
 *      120127    D. Dirkx          Moved unit to separate file from basic mathematics test; moved
 *                                  tests to separate functions; moved to Tudat Core.
 *      120127    K. Kumar          Transferred unit tests over to Boost unit test framework.
 *      120128    K. Kumar          Changed BOOST_CHECK to BOOST_CHECK_CLOSE_FRACTION and
 *                                  BOOST_CHECK_SMALL for unit test comparisons.
 *      120118    D. Gondelach      Added unit tests for convertCylindricalToCartesian.
 *                                  Removed unit test for old convertCylindricalToCartesian
 *                                  function.
 *      120214    K. Kumar          Branched from old Tudat trunk for new coordinate conversions.
 *      120512    K. Kumar          Boostified unit test.
 *      120514    P. Musegaas       Fixed small error.
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <TudatCore/Basics/testMacros.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_coordinate_conversions )

// Test conversion from cylindrical (r, theta, z) to Cartesian (x, y, z) coordinates.
BOOST_AUTO_TEST_CASE( testCylindricalToCartesianPositionCoordinateConversion )
{
    using tudat::mathematics::PI;
    using tudat::basic_mathematics::coordinate_conversions::xCartesianCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::yCartesianCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::zCartesianCoordinateIndex;

    // Test 1: test conversion of ( 0.0, pi, 1.2 ).
    {
        // Set cylindrical coordinates.
        Eigen::Vector3d cylindricalCoordinates( 0.0, PI, 1.2 );

        // Set expected Cartesian coordinates.
        Eigen::Vector3d expectedCartesianCoordinates( 0.0, 0.0, 1.2 );

        // Convert cylindrical to Cartesian coordinates.
        Eigen::Vector3d computedCartesianCoordinates
                = tudat::basic_mathematics::coordinate_conversions::
                convertCylindricalToCartesian( cylindricalCoordinates );

        // Check if computed Cartesian coordinates match expected values.
        BOOST_CHECK_SMALL( computedCartesianCoordinates( xCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( computedCartesianCoordinates( yCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianCoordinates( zCartesianCoordinateIndex ),
                                    computedCartesianCoordinates( zCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: test conversion of ( 2.3, -pi/2, -3.5 ).
    {
        // Set cylindrical coordinates.
        Eigen::Vector3d cylindricalCoordinates( 2.3, -PI/2, -3.5 );

        // Set expected Cartesian coordinates.
        Eigen::Vector3d expectedCartesianCoordinates( 2.3 * std::cos( -PI / 2.0 ),
                                                      2.3 * std::sin( -PI / 2.0 ),
                                                      -3.5 );

        // Convert cylindrical to Cartesian coordinates.
        Eigen::Vector3d convertedCartesianCoordinates
                = tudat::basic_mathematics::coordinate_conversions::
                convertCylindricalToCartesian( cylindricalCoordinates );

        // Check if converted Cartesian coordinates match expected values.
        BOOST_CHECK_SMALL( convertedCartesianCoordinates ( xCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    expectedCartesianCoordinates.segment( yCartesianCoordinateIndex, 2 ),
                    convertedCartesianCoordinates.segment( yCartesianCoordinateIndex, 2 ),
                    std::numeric_limits< double >::epsilon( ) );
    }
}

// Test conversion from cylindrical (r, theta, z, Vr, Vtheta, Vz) to Cartesian
// (x, y, z, xdot, ydot, zdot) state.
BOOST_AUTO_TEST_CASE( testCylindricalToCartesianPositionAndVelocityCoordinateConversion )
{
    using tudat::mathematics::PI;
    using tudat::basic_mathematics::coordinate_conversions::xCartesianCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::yCartesianCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::zCartesianCoordinateIndex;

    // Test 1: test conversion of (2.1, pi/2.0, 1.2, 5.4, 4.5, -3.9).
    {
        // Set Cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd cylindricalState = ( Eigen::VectorXd( 6 )
                                             << 2.1, PI / 2.0, 1.2, 5.4, 4.5, -3.9 ).finished( );

        // Set expected Cartesian state (x, y, z, xdot, ydot, zdot).
        Eigen::VectorXd expectedCartesianState( 6 );
        expectedCartesianState << ( Eigen::VectorXd( 6 )
                                    << 2.1 * std::cos( PI / 2.0 ),
                                    2.1 * std::sin( PI / 2.0 ),
                                    1.2,
                                    5.4 * std::cos( PI / 2.0 ) - 4.5 * std::sin( PI / 2.0 ),
                                    5.4 * std::sin( PI / 2.0 ) + 4.5 * std::cos( PI / 2.0 ),
                                    -3.9 ).finished( );

        // Convert cylindrical to Cartesian state (x, y, z, xdot, ydot, zdot).
        Eigen::VectorXd convertedCartesianState
                = tudat::basic_mathematics::coordinate_conversions::
                convertCylindricalToCartesian( cylindricalState );

        // Check if converted Cartesian state match expected state.
        BOOST_CHECK_SMALL( convertedCartesianState ( xCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    expectedCartesianState.segment( yCartesianCoordinateIndex, 5 ),
                    convertedCartesianState.segment( yCartesianCoordinateIndex, 5 ),
                    std::numeric_limits< double >::epsilon( ) );

    }

    // Test 2: test conversion of (0.0, 8.2*pi/3.0, -2.5, -5.8, 0.0, 1.7).
    {
        // Set Cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd cylindricalState = ( Eigen::VectorXd( 6 )
                                             << 0.0, 8.2 * PI / 3.0,
                                             -2.5, -5.8, 0.0, 1.7 ).finished( );

        // Set expected Cartesian state (x, y, z, xdot, ydot, zdot).
        Eigen::VectorXd expectedCartesianState( 6 );
        expectedCartesianState << ( Eigen::VectorXd( 6 )
                                    << 0.0,
                                    0.0,
                                    -2.5,
                                    -5.8 * std::cos( 8.2 * PI / 3.0 ),
                                    -5.8 * std::sin( 8.2 * PI / 3.0 ),
                                    1.7 ).finished( );

        // Convert cylindrical to Cartesian state (x, y, z, xdot, ydot, zdot).
        Eigen::VectorXd convertedCartesianState
                = tudat::basic_mathematics::coordinate_conversions::
                convertCylindricalToCartesian( cylindricalState );

        // Check if converted Cartesian state match expected state.
        BOOST_CHECK_SMALL( convertedCartesianState( xCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCartesianState( yCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    expectedCartesianState.segment( zCartesianCoordinateIndex, 4 ),
                    convertedCartesianState.segment( zCartesianCoordinateIndex, 4 ),
                    std::numeric_limits< double >::epsilon( ) );
    }
}

// Test conversion from Cartesian (x, y, z) to cylindrical (r, theta, z) coordinates.
BOOST_AUTO_TEST_CASE( testCartesianToCylindricalPositionCoordinateConversion )
{
    using tudat::mathematics::PI;
    using tudat::basic_mathematics::coordinate_conversions::rCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::thetaCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::zCylindricalCoordinateIndex;

    // Test 1: test conversion of ( 0.0, 0.0, 1.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        Eigen::Vector3d cartesianCoordinates( 0.0, 0.0, 1.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        Eigen::Vector3d expectedCylindricalCoordinates( 0.0, 0.0, 1.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        Eigen::Vector3d convertedCylindricalCoordinates
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianCoordinates );

        // Check if converted cylindrical coordinates match expected values.
        BOOST_CHECK_SMALL( convertedCylindricalCoordinates( rCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCylindricalCoordinates( thetaCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCylindricalCoordinates( zCylindricalCoordinateIndex ),
                                    convertedCylindricalCoordinates( zCylindricalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: test conversion of ( 0.0, 2.0, 1.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        Eigen::Vector3d cartesianCoordinates( 0.0, 2.0, 1.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        Eigen::Vector3d expectedCylindricalCoordinates( 2.0, PI / 2.0, 1.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        Eigen::Vector3d convertedCylindricalCoordinates
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianCoordinates );

        // Check if converted cylindrical coordinates match expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalCoordinates,
                                           convertedCylindricalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: test conversion of ( 0.0, -2.0, -1.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        Eigen::Vector3d cartesianCoordinates( 0.0, -2.0, -1.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        Eigen::Vector3d expectedCylindricalCoordinates( 2.0, 3.0 * PI / 2.0, -1.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        Eigen::Vector3d convertedCylindricalCoordinates
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianCoordinates );

        // Check if converted cylindrical coordinates match expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalCoordinates,
                                           convertedCylindricalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 4: test conversion of ( -5.0, -8.0, 5.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        Eigen::Vector3d cartesianCoordinates( -5.0, -8.0, 5.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        Eigen::Vector3d expectedCylindricalCoordinates( std::sqrt( 25.0 + 64.0 ),
                                                        std::atan2( -8.0,-5.0 ) + 2.0 * PI,
                                                        5.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        Eigen::Vector3d convertedCylindricalCoordinates
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianCoordinates );

        // Check if converted cylindrical coordinates match expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalCoordinates,
                                           convertedCylindricalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

}

// Test conversion from Cartesian (x, y, z, xdot, ydot, zdot) to cylindrical
// (r, theta, z, Vr, Vtheta, Vz) state.
BOOST_AUTO_TEST_CASE( testCartesianToCylindricalPositionAndVelocityCoordinateConversion )
{
    using tudat::mathematics::PI;
    using tudat::basic_mathematics::coordinate_conversions::rCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::thetaCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::zCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::rDotCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::thetaDotCylindricalCoordinateIndex;
    using tudat::basic_mathematics::coordinate_conversions::zDotCylindricalCoordinateIndex;

    // Test 1: test conversion of ( 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).
    {
        // Set Cartesian state (x,y,z,xdot,ydot,zdot).
        Eigen::VectorXd cartesianState = ( Eigen::VectorXd( 6 )
                                           << 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).finished( );

        // Set expected cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd expectedCylindricalState = ( Eigen::VectorXd( 6 )
                                                     << 0.0, 0.0, 1.0,
                                                     std::sqrt( 25.0 + 36.0 ),
                                                     0.0, -9.0 ).finished( );

        // Convert Cartesian to cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd convertedCylindricalState
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianState );

        // Check that converted cylindrical state matches expected state.
        BOOST_CHECK_SMALL( convertedCylindricalState( rCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCylindricalState( thetaCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCylindricalState( zCylindricalCoordinateIndex ),
                                    convertedCylindricalState( zCylindricalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCylindricalState( rDotCylindricalCoordinateIndex ),
                                    convertedCylindricalState( rDotCylindricalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCylindricalState( thetaDotCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCylindricalState( zDotCylindricalCoordinateIndex ),
                                    convertedCylindricalState( zDotCylindricalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: test conversion of ( 2.0, 0.0, -5.0, -4.0, 6.0, -6.0 ).
    {
        // Set Cartesian state (x,y,z,xdot,ydot,zdot).
        Eigen::VectorXd cartesianState = ( Eigen::VectorXd( 6 )
                                           << 2.0, 0.0, -5.0, -4.0, 6.0, -6.0 ).finished( );

        // Set expected cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd expectedCylindricalState = ( Eigen::VectorXd( 6 )
                                                     << 2.0, 0.0, -5.0,
                                                     -4.0, 6.0, -6.0 ).finished( );

        // Convert Cartesian to cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd convertedCylindricalState
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianState );

        // Check that converted cylindrical state matches expected state.
        BOOST_CHECK_CLOSE_FRACTION( expectedCylindricalState( rCylindricalCoordinateIndex ),
                                    convertedCylindricalState( rCylindricalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCylindricalState( thetaCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    expectedCylindricalState.segment( zCylindricalCoordinateIndex, 4 ),
                    convertedCylindricalState.segment( zCylindricalCoordinateIndex, 4 ),
                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: test conversion of ( -7.0, -4.0, 3.0, 5.0, -3.0, 7.0 ).
    {
        // Set Cartesian state (x,y,z,xdot,ydot,zdot).
        Eigen::VectorXd cartesianState = ( Eigen::VectorXd( 6 )
                                           << -7.0, -4.0, 3.0, 5.0, -3.0, 7.0 ).finished( );

        // Set expected cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd expectedCylindricalState = ( Eigen::VectorXd( 6 )
                                                     << std::sqrt( 49.0 + 16.0 ),
                                                     std::atan2( -4.0, -7.0 ) + 2.0 * PI,
                                                     3.0,
                                                     ( -7.0 * 5.0 + ( -4.0 ) * -3.0 )
                                                     / std::sqrt( 49.0 + 16.0 ),
                                                     ( -7.0 * -3.0 - ( -4.0 ) * 5.0 )
                                                     / sqrt( 49.0 + 16.0 ), 7.0 ).finished( );

        // Convert Cartesian to cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        Eigen::VectorXd convertedCylindricalState
                = tudat::basic_mathematics::coordinate_conversions::
                convertCartesianToCylindrical( cartesianState );

        // Check that converted cylindrical state matches expected state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalState, convertedCylindricalState,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
