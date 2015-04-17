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
 *      120716    D. Dirkx          Revised and updated unit tests to cover bug in
 *                                  convertCartesianToSpherical().
 *
 *    References
 *      Stewart, J. Calculus: Early Transcendentals, Fourth Edition, Brooks/Cole Publishing
 *          Company, Pacific Grove, CA, USA, 1999.
 *
 *    Notes
 *
 */

#include <cmath>
#include <limits>
#include <iostream>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Basics/testMacros.h"

namespace tudat
{
namespace unit_tests
{

//! Test suite for coordinate conversion functions.
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
        Eigen::VectorXd sphericalCoordinates = convertCartesianToSpherical( cartesianCoordinates );

        // Check if converted spherical coordinates are correct.
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( sphericalCoordinates, expectedSphericalCoordinates,
                                           std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
