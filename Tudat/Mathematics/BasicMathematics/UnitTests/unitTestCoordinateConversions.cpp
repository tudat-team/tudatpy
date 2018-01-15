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
 *      Mathworks. gravitysphericalharmonic, implement spherical harmonic representation of
 *        planetary gravity. Help documentation of MATLAB R2012a, 2012.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_coordinate_conversions )

//! Test if spherical-to-Cartesian conversion is working correctly.
BOOST_AUTO_TEST_CASE( testSphericalToCartesianConversion )
{
    using std::sqrt;
    using mathematical_constants::PI;
    using coordinate_conversions::convertSphericalToCartesian;

    // Test 1: test conversion of: ( 0.0, 0.0, 0.0 ).
    {
        Eigen::Vector3d sphericalCoordinates( 0.0, 0.0, 0.0 );

        // Convert spherical coordinates to Cartesian coordinates.
        Eigen::Vector3d cartesianCoordinates = convertSphericalToCartesian( sphericalCoordinates );

        // Expected cartesian coordinates.
        Eigen::Vector3d expectedCartesianCoordinates( 0.0, 0.0, 0.0 );

        // Check if converted Cartesian coordinates are correct.
        TUDAT_CHECK_MATRIX_BASE( cartesianCoordinates, expectedCartesianCoordinates )
                BOOST_CHECK_SMALL( cartesianCoordinates.coeff( row, col ),
                                   std::numeric_limits< double >::min( ) );
    }

    // Test 2: test conversion of: ( 1.0, pi/6, pi/6 ), from Stewart (2003), Exercise 12.7.15.
    {
        Eigen::Vector3d sphericalCoordinates( 1.0, PI / 6.0, PI / 6.0 );

        // Convert spherical coordinates to Cartesian coordinates.
        Eigen::Vector3d cartesianCoordinates = convertSphericalToCartesian( sphericalCoordinates );

        // Expected cartesian coordinates.
        Eigen::Vector3d expectedCartesianCoordinates(
                    sqrt( 3.0 ) / 4.0, 1.0 / 4.0, sqrt( 3.0 ) / 2.0 );

        // Check if converted Cartesian coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                    cartesianCoordinates, expectedCartesianCoordinates, 1.0e-15 );
    }

    // Test 3: test conversion of: ( 2.0, pi/4, pi/3 ), from Stewart (2003), Exercise 12.7.17.
    {
        Eigen::Vector3d sphericalCoordinates( 2.0, PI / 4.0, PI / 3.0 );

        // Convert spherical coordinates to Cartesian coordinates.
        Eigen::Vector3d cartesianCoordinates = convertSphericalToCartesian( sphericalCoordinates );

        // Expected cartesian coordinates.
        Eigen::Vector3d expectedCartesianCoordinates(
                    0.5 * sqrt( 2.0 ), 0.5 * sqrt( 6.0 ), sqrt( 2.0 ) );

        // Check if converted Cartesian coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cartesianCoordinates, expectedCartesianCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 4: test conversion of: ( 2.0, pi/3, pi/4 , from Stewart (2003), Section 12.7,
    //         Example 4.
    {
        Eigen::Vector3d sphericalCoordinates( 2.0, PI / 3.0, PI / 4.0 );

        // Convert spherical coordinates to Cartesian coordinates.
        Eigen::Vector3d cartesianCoordinates = convertSphericalToCartesian( sphericalCoordinates );

        // Expected cartesian coordinates.
        Eigen::Vector3d expectedCartesianCoordinates( sqrt( 1.5 ), sqrt( 1.5 ), 1.0 );

        // Check if converted Cartesian coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( cartesianCoordinates, expectedCartesianCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

//! Test if Cartesian-to-spherical conversion is working correctly.
BOOST_AUTO_TEST_CASE( testCartesianToSphericalConversion )
{
    using std::acos;
    using std::atan2;
    using std::pow;
    using std::sqrt;
    using coordinate_conversions::convertCartesianToSpherical;
    using mathematical_constants::PI;

    // Test 1: Test conversion of: ( 0.0, 0.0, 0.0 ).
    {
        Eigen::Vector3d cartesianCoordinates = Eigen::Vector3d::Zero( );

        // Expected vector in spherical coordinates.
        Eigen::Vector3d expectedSphericalCoordinates = Eigen::Vector3d::Zero( );

        // Result vector in spherical coordinates.
        Eigen::Vector3d sphericalCoordinates = convertCartesianToSpherical( cartesianCoordinates );

        // Check if converted spherical coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalCoordinates, expectedSphericalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: Test conversion of: ( 1.0, sqrt( 3 ), 2 * sqrt(3) ), from Stewart (2003),
    // Exercise 12.7.19.
    {
        Eigen::Vector3d cartesianCoordinates( 1.0, sqrt( 3. ), 2.0 * sqrt( 3. ) );

        // Expected vector in spherical coordinates.
        Eigen::Vector3d expectedSphericalCoordinates( 4.0, PI / 6.0, PI / 3.0 );

        // Result vector in spherical coordinates.
        Eigen::Vector3d sphericalCoordinates = convertCartesianToSpherical( cartesianCoordinates );

        // Check if converted spherical coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalCoordinates, expectedSphericalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: Test conversion of: ( 0.0, -1.0, -1.0 ), from Stewart (2003), Exercise 12.7.21.
    {
        Eigen::Vector3d cartesianCoordinates( 0.0, -1.0, -1.0 );

        // Expected vector in spherical coordinates.
        Eigen::Vector3d expectedSphericalCoordinates( sqrt( 2.0 ), 3.0 * PI / 4.0, -PI / 2.0 );

        // Result vector in spherical coordinates.
        Eigen::Vector3d sphericalCoordinates = Eigen::Vector3d::Zero( );

        // Compute conversions.
        sphericalCoordinates = convertCartesianToSpherical( cartesianCoordinates );

        // Check if converted spherical coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalCoordinates, expectedSphericalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 4: Test conversion of: ( 0.0, 2 sqrt( 3 ), -2.0 ), from Stewart (2003), Section 12.7,
    // Example 5.
    {
        Eigen::Vector3d cartesianCoordinates( 0.0, 2.0 * sqrt( 3. ), -2.0 );

        // Expected vector in spherical coordinates.
        Eigen::Vector3d expectedSphericalCoordinates( 4.0, 2.0 * PI / 3.0, PI / 2.0 );

        // Result vector in spherical coordinates.
        Eigen::Vector3d sphericalCoordinates = convertCartesianToSpherical( cartesianCoordinates );

        // Check if converted spherical coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalCoordinates, expectedSphericalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

// Test conversion from cylindrical (r, theta, z) to Cartesian (x, y, z) coordinates.
BOOST_AUTO_TEST_CASE( testCylindricalToCartesianPositionCoordinateConversion )
{
    using mathematical_constants::PI;
    using coordinate_conversions::xCartesianCoordinateIndex;
    using coordinate_conversions::yCartesianCoordinateIndex;
    using coordinate_conversions::zCartesianCoordinateIndex;

    // Test 1: test conversion of ( 0.0, pi, 1.2 ).
    {
        // Set cylindrical coordinates.
        const Eigen::Vector3d cylindricalCoordinates( 0.0, PI, 1.2 );

        // Set expected Cartesian coordinates.
        const Eigen::Vector3d expectedCartesianCoordinates( 0.0, 0.0, 1.2 );

        // Convert cylindrical to Cartesian coordinates.
        const Eigen::Vector3d computedCartesianCoordinates
                = coordinate_conversions::
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
        const Eigen::Vector3d cylindricalCoordinates( 2.3, -PI/2, -3.5 );

        // Set expected Cartesian coordinates.
        const Eigen::Vector3d expectedCartesianCoordinates( 2.3 * std::cos( -PI / 2.0 ),
                                                            2.3 * std::sin( -PI / 2.0 ),
                                                            -3.5 );

        // Convert cylindrical to Cartesian coordinates.
        const Eigen::Vector3d convertedCartesianCoordinates
                = coordinate_conversions::
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
    using mathematical_constants::PI;
    using coordinate_conversions::xCartesianCoordinateIndex;
    using coordinate_conversions::yCartesianCoordinateIndex;
    using coordinate_conversions::zCartesianCoordinateIndex;

    // Test 1: test conversion of (2.1, pi/2.0, 1.2, 5.4, 4.5, -3.9).
    {
        // Set Cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d cylindricalState =
                ( Eigen::Vector6d( ) << 2.1, PI / 2.0, 1.2, 5.4, 4.5, -3.9 ).finished( );

        // Set expected Cartesian state (x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d expectedCartesianState =
                ( Eigen::Vector6d( )
                  << 2.1 * std::cos( PI / 2.0 ),
                  2.1 * std::sin( PI / 2.0 ),
                  1.2,
                  5.4 * std::cos( PI / 2.0 ) - 4.5 * std::sin( PI / 2.0 ),
                  5.4 * std::sin( PI / 2.0 ) + 4.5 * std::cos( PI / 2.0 ),
                  -3.9 ).finished( );

        // Convert cylindrical to Cartesian state (x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d convertedCartesianState
                = coordinate_conversions::
                convertCylindricalToCartesianState( cylindricalState );

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
        const Eigen::Vector6d cylindricalState = ( Eigen::Vector6d( )
                                                   << 0.0, 8.2 * PI / 3.0,
                                                   -2.5, -5.8, 0.0, 1.7 ).finished( );

        // Set expected Cartesian state (x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d expectedCartesianState = ( Eigen::Vector6d( )
                                                         << 0.0,
                                                         0.0,
                                                         -2.5,
                                                         -5.8 * std::cos( 8.2 * PI / 3.0 ),
                                                         -5.8 * std::sin( 8.2 * PI / 3.0 ),
                                                         1.7 ).finished( );

        // Convert cylindrical to Cartesian state (x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d convertedCartesianState
                = coordinate_conversions::
                convertCylindricalToCartesianState( cylindricalState );

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
    using mathematical_constants::PI;
    using coordinate_conversions::rCylindricalCoordinateIndex;
    using coordinate_conversions::thetaCylindricalCoordinateIndex;
    using coordinate_conversions::zCylindricalCoordinateIndex;

    // Test 1: test conversion of ( 0.0, 0.0, 1.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        const Eigen::Vector3d cartesianCoordinates( 0.0, 0.0, 1.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d expectedCylindricalCoordinates( 0.0, 0.0, 1.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d convertedCylindricalCoordinates
                = coordinate_conversions::
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
        const Eigen::Vector3d cartesianCoordinates( 0.0, 2.0, 1.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d expectedCylindricalCoordinates( 2.0, PI / 2.0, 1.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d convertedCylindricalCoordinates
                = coordinate_conversions::
                convertCartesianToCylindrical( cartesianCoordinates );

        // Check if converted cylindrical coordinates match expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalCoordinates,
                                           convertedCylindricalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: test conversion of ( 0.0, -2.0, -1.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        const Eigen::Vector3d cartesianCoordinates( 0.0, -2.0, -1.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d expectedCylindricalCoordinates( 2.0, 3.0 * PI / 2.0, -1.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d convertedCylindricalCoordinates
                = coordinate_conversions::
                convertCartesianToCylindrical( cartesianCoordinates );

        // Check if converted cylindrical coordinates match expected values.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalCoordinates,
                                           convertedCylindricalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }

    // Test 4: test conversion of ( -5.0, -8.0, 5.0 ).
    {
        // Set Cartesian coordinates (x, y, z).
        const Eigen::Vector3d cartesianCoordinates( -5.0, -8.0, 5.0 );

        // Set expected cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d expectedCylindricalCoordinates( std::sqrt( 25.0 + 64.0 ),
                                                              std::atan2( -8.0,-5.0 ) + 2.0 * PI,
                                                              5.0 );

        // Convert Cartesian to cylindrical coordinates (r, theta, z).
        const Eigen::Vector3d convertedCylindricalCoordinates
                = coordinate_conversions::
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
    using mathematical_constants::PI;
    using coordinate_conversions::rCylindricalCoordinateIndex;
    using coordinate_conversions::thetaCylindricalCoordinateIndex;
    using coordinate_conversions::zCylindricalCoordinateIndex;
    using coordinate_conversions::rDotCylindricalCoordinateIndex;
    using coordinate_conversions::vThetaCylindricalCoordinateIndex;
    using coordinate_conversions::zDotCylindricalCoordinateIndex;

    // Test 1: test conversion of ( 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).
    {
        // Set Cartesian state (x,y,z,xdot,ydot,zdot).
        const Eigen::Vector6d cartesianState = ( Eigen::Vector6d( )
                                                 << 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).finished( );

        // Set expected cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d expectedCylindricalState =
                ( Eigen::Vector6d( )
                                                           << 0.0, 0.0, 1.0,
                                                           std::sqrt( 25.0 + 36.0 ),
                                                           0.0, -9.0 ).finished( );

        // Convert Cartesian to cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d convertedCylindricalState
                = coordinate_conversions::
                convertCartesianToCylindricalState( cartesianState );

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

        BOOST_CHECK_SMALL( convertedCylindricalState( vThetaCylindricalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCylindricalState( zDotCylindricalCoordinateIndex ),
                                    convertedCylindricalState( zDotCylindricalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: test conversion of ( 2.0, 0.0, -5.0, -4.0, 6.0, -6.0 ).
    {
        // Set Cartesian state (x,y,z,xdot,ydot,zdot).
        const Eigen::Vector6d cartesianState = ( Eigen::Vector6d( )
                                                 << 2.0, 0.0, -5.0, -4.0, 6.0, -6.0 ).finished( );

        // Set expected cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d expectedCylindricalState = ( Eigen::Vector6d( )
                                                           << 2.0, 0.0, -5.0,
                                                           -4.0, 6.0, -6.0 ).finished( );

        // Convert Cartesian to cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d convertedCylindricalState
                = coordinate_conversions::
                convertCartesianToCylindricalState( cartesianState );

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
        const Eigen::Vector6d cartesianState = ( Eigen::Vector6d( )
                                                 << -7.0, -4.0, 3.0, 5.0, -3.0, 7.0 ).finished( );

        // Set expected cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d expectedCylindricalState = ( Eigen::Vector6d( )
                                                           << std::sqrt( 49.0 + 16.0 ),
                                                           std::atan2( -4.0, -7.0 ) + 2.0 * PI,
                                                           3.0,
                                                           ( -7.0 * 5.0 + ( -4.0 ) * -3.0 )
                                                           / std::sqrt( 49.0 + 16.0 ),
                                                           ( -7.0 * -3.0 - ( -4.0 ) * 5.0 )
                                                           / sqrt( 49.0 + 16.0 ), 7.0 ).finished( );

        // Convert Cartesian to cylindrical state (r, theta, z, Vr, Vtheta, Vz).
        const Eigen::Vector6d convertedCylindricalState
                = coordinate_conversions::
                convertCartesianToCylindricalState( cartesianState );

        // Check that converted cylindrical state matches expected state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCylindricalState, convertedCylindricalState,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

// Test conversion from spherical to Cartesian gradient.
BOOST_AUTO_TEST_CASE( test_SphericalToCartesianGradientConversion )
{
    // Define an arbitrary spherical gradient.
    const Eigen::Vector3d sphericalGradient( -2.438146967844150, 1.103964186650749e4,
                                             5.932496087870976e1 );

    // Define an arbitrary spherical position vector.
    const Eigen::Vector3d cartesianCoordinates( -7.0e6, 8.5e6, -6.5e6 );

    // Compute Cartesian gradient.
    const Eigen::Vector3d cartesianGradient = coordinate_conversions::
            convertSphericalToCartesianGradient( sphericalGradient, cartesianCoordinates );

    // Define expected Cartesian gradient. These values are obtained from the computations of
    // 'gx', 'gy' and 'gz' in the MATLAB function 'gravitysphericalharmonics'. This function is
    // described by Mathworks [2012].
    const Eigen::Vector3d expectedCartesianGradient( 1.334464111786447, -1.620429182163669,
                                                     1.240151677089496 );

    // Check if the computed Cartesian gradient matches the expected value.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianGradient, cartesianGradient, 1e-15 );
}

// Test derivative of Cartesian gradient, keeping spherical gradient constant.
BOOST_AUTO_TEST_CASE( test_SphericalToCartesianGradientPartialDerivatives )
{
    // Define an arbitrary spherical gradient.
    const Eigen::Vector3d sphericalGradient( -2.438146967844150, 1.103964186650749e4,
                                             5.932496087870976e1 );

    // Define an arbitrary spherical position vector.
    const Eigen::Vector3d cartesianCoordinates( -7.0e6, 8.5e6, -6.5e6 );


    // Define perturbation for numerical computation of partial
    const Eigen::Vector3d cartesianCoordinatePerturbation( 10.0, 10.0, 10.0 );
    Eigen::Vector3d perturbedCartesianCoordinates;

    std::vector< Eigen::Matrix3d > matrixPartials;
    matrixPartials.resize( 3 );

    std::vector< Eigen::Matrix3d > subMatrices;
    subMatrices.resize( 3 );
    // Compute analytical partial (with sub-components)
    coordinate_conversions::getDerivativeOfSphericalToCartesianGradient( sphericalGradient, cartesianCoordinates, subMatrices );

    Eigen::Matrix3d uppPerturbedMatrix;
    for( unsigned int i = 0; i < 3; i++ )
    {
        // Compute partial numerically.
        perturbedCartesianCoordinates = cartesianCoordinates;
        perturbedCartesianCoordinates( i ) += cartesianCoordinatePerturbation( i );
        matrixPartials[ i ] = coordinate_conversions::getSphericalToCartesianGradientMatrix(
                    perturbedCartesianCoordinates );

        perturbedCartesianCoordinates = cartesianCoordinates;
        perturbedCartesianCoordinates( i ) -= cartesianCoordinatePerturbation( i );
        matrixPartials[ i ] -= coordinate_conversions::getSphericalToCartesianGradientMatrix(
                    perturbedCartesianCoordinates );

        matrixPartials[ i ] /= ( 2.0 * cartesianCoordinatePerturbation( i ) );

        // Normalize components for uniform tolerance.
        matrixPartials[ i ].block( 0, 0, 3, 1 ) = matrixPartials[ i ].block( 0, 0, 3, 1 ) / cartesianCoordinates.norm( );
        subMatrices[ i ].block( 0, 0, 3, 1 ) = subMatrices[ i ].block( 0, 0, 3, 1 ) / cartesianCoordinates.norm( );

        for( unsigned k = 0; k < 3; k++ )
        {
            for( unsigned l = 0; l < 3; l++ )
            {
                BOOST_CHECK_SMALL( std::fabs( matrixPartials[ i ]( k, l ) - subMatrices[ i ]( k, l ) ), 1.0E-23 );
            }

        }
    }
}

// Test conversion from Cartesian (x, y, z, xdot, ydot, zdot) to Spherical (radius, azimuth,
// elevation, radial velocity, azimuthal velocity, elevational velocity) state.
BOOST_AUTO_TEST_CASE( testCartesianToSphericalStateConversion )
{
    using mathematical_constants::PI;
    using coordinate_conversions::radiusSphericalCoordinateIndex;
    using coordinate_conversions::azimuthSphericalCoordinateIndex;
    using coordinate_conversions::elevationSphericalCoordinateIndex;
    using coordinate_conversions::radialVelocitySphericalCoordinateIndex;
    using coordinate_conversions::azimuthVelocitySphericalCoordinateIndex;
    using coordinate_conversions::elevationVelocitySphericalCoordinateIndex;
    using coordinate_conversions::convertCartesianToSphericalState;
    using std::sqrt;

    // Test 1: test conversion of ( 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).
    {
        // Set Cartesian state (x, y , z, xdot, ydot, zdot).
        const Eigen::Vector6d cartesianState
                = ( Eigen::Vector6d( ) << 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).finished( );

        // Set expected spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d expectedSphericalState
                = ( Eigen::Vector6d( ) << 1.0, 0.0, PI / 2.0, -9.0, 6.0, -5.0 ).finished( );

        // Convert Cartesian to spherical state.
        const Eigen::Vector6d convertedSphericalState
                = convertCartesianToSphericalState( cartesianState );

        // Check that converted spherical state matches expected spherical state.
        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState( radiusSphericalCoordinateIndex ),
                                    convertedSphericalState( radiusSphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedSphericalState( azimuthSphericalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState( elevationSphericalCoordinateIndex ),
                                    convertedSphericalState( elevationSphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState(
                                        radialVelocitySphericalCoordinateIndex ),
                                    convertedSphericalState(
                                        radialVelocitySphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState(
                                        azimuthVelocitySphericalCoordinateIndex ),
                                    convertedSphericalState(
                                        azimuthVelocitySphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState(
                                        elevationVelocitySphericalCoordinateIndex ),
                                    convertedSphericalState(
                                        elevationVelocitySphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: test conversion of ( 2.0, 0.0, 0.0, -4.0, 6.0, -6.0 ).
    {
        // Set Cartesian state (x, y , z, xdot, ydot, zdot).
        const Eigen::Vector6d cartesianState
                = ( Eigen::Vector6d( ) << 2.0, 0.0, 0.0, -4.0, 6.0, -6.0 ).finished( );

        // Set expected spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d expectedSphericalState
                = ( Eigen::Vector6d( ) << 2.0, 0.0, 0.0, -4.0, 6.0, -6.0 ).finished( );

        // Convert Cartesian to spherical state ( r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d convertedSphericalState
                = convertCartesianToSphericalState( cartesianState );

        // Check that converted spherical state matches expected spherical state.
        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState( radiusSphericalCoordinateIndex ),
                                    convertedSphericalState( radiusSphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedSphericalState( azimuthSphericalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedSphericalState( elevationSphericalCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState(
                                        radialVelocitySphericalCoordinateIndex ),
                                    convertedSphericalState(
                                        radialVelocitySphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState(
                                        azimuthVelocitySphericalCoordinateIndex ),
                                    convertedSphericalState(
                                        azimuthVelocitySphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedSphericalState(
                                        elevationVelocitySphericalCoordinateIndex ),
                                    convertedSphericalState(
                                        elevationVelocitySphericalCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: test conversion of ( -1.0, -1.0, sqrt( 2.0 ), -1.0, 1.0, sqrt( 2.0 ) ).
    {
        // Set Cartesian state ( x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d cartesianState
                = ( Eigen::Vector6d( )
                    << -1.0, -1.0, sqrt( 2.0 ), -1.0, 1.0, sqrt( 2.0 ) ).finished( );

        // Set expected spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d expectedSphericalState
                = ( Eigen::Vector6d( )
                    << 2.0, -PI * 3.0 / 4.0, PI / 4.0, 1.0, -sqrt( 2.0 ), 1.0 ).finished( );

        // Convert Cartesian to spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d convertedSphericalState
                = convertCartesianToSphericalState( cartesianState );

        // Check that converted spherical state matches the expected spherical state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedSphericalState, convertedSphericalState,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

// Test conversion from spherical (radius, azimuth, elevation, radial velocity, azimuthal velocity,
// elevational velocity) to Cartesian state (x, y, z, xdot, ydot, zdot).
BOOST_AUTO_TEST_CASE( testSphericalToCartesianStateConversion )
{
    using mathematical_constants::PI;
    using coordinate_conversions::xCartesianCoordinateIndex;
    using coordinate_conversions::yCartesianCoordinateIndex;
    using coordinate_conversions::zCartesianCoordinateIndex;
    using coordinate_conversions::xDotCartesianCoordinateIndex;
    using coordinate_conversions::yDotCartesianCoordinateIndex;
    using coordinate_conversions::zDotCartesianCoordinateIndex;
    using coordinate_conversions::convertSphericalToCartesianState;
    using std::sqrt;

    // Test 1: test conversion of ( 1.0, 0.0, pi/2, -9.0, 6.0, -5.0 ).
    {
        // Set Spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d sphericalState
                = ( Eigen::Vector6d( ) << 1.0, 0.0, PI / 2.0, -9.0, 6.0, -5.0 ).finished( );

        // Set expected Cartesian state (x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d expectedCartesianState
                = ( Eigen::Vector6d( ) << 0.0, 0.0, 1.0, 5.0, 6.0, -9.0 ).finished( );

        // Convert spherical to Cartesian state.
        const Eigen::Vector6d convertedCartesianState
                = convertSphericalToCartesianState( sphericalState );

        // Check that converted Cartesian state matches expected state.
        BOOST_CHECK_SMALL( convertedCartesianState( xCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCartesianState( yCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( zCartesianCoordinateIndex ),
                                    convertedCartesianState( zCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( xDotCartesianCoordinateIndex ),
                                    convertedCartesianState( xDotCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( yDotCartesianCoordinateIndex ),
                                    convertedCartesianState( yDotCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( zDotCartesianCoordinateIndex ),
                                    convertedCartesianState( zDotCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 2: test conversion of ( 2.0, 0.0, 0.0, -4.0, 6.0, -6.0 ).
    {
        // Set spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d sphericalState
                = ( Eigen::Vector6d( ) << 2.0, 0.0, 0.0, -4.0, 6.0, -6.0 ).finished( );

        // Set expected Cartesian state ( x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d expectedCartesianState
                = ( Eigen::Vector6d( ) << 2.0, 0.0, 0.0, -4.0, 6.0, -6.0 ).finished( );

        // Convert spherical to Cartesian state.
        const Eigen::Vector6d convertedCartesianState
                = convertSphericalToCartesianState( sphericalState );

        // Check that converted Cartesian state matches expected state.
        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( xCartesianCoordinateIndex ),
                                    convertedCartesianState( xCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCartesianState( yCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_SMALL( convertedCartesianState( zCartesianCoordinateIndex ),
                           std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( xDotCartesianCoordinateIndex ),
                                    convertedCartesianState( xDotCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( yDotCartesianCoordinateIndex ),
                                    convertedCartesianState( yDotCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );

        BOOST_CHECK_CLOSE_FRACTION( expectedCartesianState( zDotCartesianCoordinateIndex ),
                                    convertedCartesianState( zDotCartesianCoordinateIndex ),
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 3: test conversion of ( 2.0, -3*pi/4, pi/4, 1.0, -sqrt( 2.0 ), 1.0 ).
    {
        // Set spherical state (r, theta, phi, Vr, Vtheta, Vphi).
        const Eigen::Vector6d sphericalState
                = ( Eigen::Vector6d( )
                    << 2.0, -PI * 3.0 / 4.0, PI / 4.0, 1.0, -sqrt( 2.0 ), 1.0 ).finished( );

        // Set expected Cartesian state (x, y, z, xdot, ydot, zdot).
        const Eigen::Vector6d expectedCartesianState
                = ( Eigen::Vector6d( )
                    << -1.0, -1.0, sqrt( 2.0 ), -1.0, 1.0, sqrt( 2.0 ) ).finished( );

        // Convert spherical to Cartesian state.
        const Eigen::Vector6d convertedCartesianState
                = convertSphericalToCartesianState( sphericalState );

        // Check that converted Cartesian state matches the expected state.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCartesianState, convertedCartesianState,
                                           1.0e-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
