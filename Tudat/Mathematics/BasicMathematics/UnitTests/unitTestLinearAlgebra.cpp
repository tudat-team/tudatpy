/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

//! Test suite for unit conversion functions.
BOOST_AUTO_TEST_SUITE( test_unit_conversions )

//! Test if angle between vectors is computed correctly.
BOOST_AUTO_TEST_CASE( testAngleBetweenVectorFunctions )
{
    // Using declarations.
    using std::cos;
    using std::sqrt;
    using linear_algebra::computeAngleBetweenVectors;
    using linear_algebra::computeCosineOfAngleBetweenVectors;

    // Four tests are executed. First, the equality of the caluclated cosineOfAngle and the cosine
    // of the calculated angle is checked. Subsequently, the values of the angle and cosineOfAngle
    // are checked against reference values, which are analytical in the first two cases and
    // taken from Matlab results in the third. The first three tests are written for vectors of length
    // 3. The fourth test is written for a vector of length 5.

    // Test 1: Test values for two equal vectors of length 3.
    {
        Eigen::Vector3d testVector1_ = Eigen::Vector3d( 3.0, 2.1, 4.6 );
        Eigen::Vector3d testVector2_ = Eigen::Vector3d( 3.0, 2.1, 4.6 );

        double angle = computeAngleBetweenVectors( testVector1_, testVector2_ );
        double cosineOfAngle = computeCosineOfAngleBetweenVectors( testVector1_, testVector2_ );

        // Check if computed angle and cosine-of-angle are correct.
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle,
                           std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_GE( cosineOfAngle, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( cosineOfAngle - 1.0, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_LE( angle, std::sqrt( std::numeric_limits< double >::epsilon( ) ) );
    }

    // Test 2: Test values for two equal, but opposite vectors of length 3.
    {
        Eigen::Vector3d testVector1_ = Eigen::Vector3d( 3.0, 2.1, 4.6 );
        Eigen::Vector3d testVector2_ = Eigen::Vector3d( -3.0, -2.1, -4.6 );

        double angle = computeAngleBetweenVectors( testVector1_, testVector2_ );
        double cosineOfAngle = computeCosineOfAngleBetweenVectors( testVector1_, testVector2_ );

        // Check if computed angle and cosine-of-angle are correct.
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle,
                           std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( cosineOfAngle < std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( cosineOfAngle + 1.0, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( angle, mathematical_constants::PI,
                                    std::sqrt( std::numeric_limits< double >::epsilon( ) ) );
    }

    // Test 3: Test values for two vectors of length 3, benchmark values computed using Matlab.
    {
        Eigen::Vector3d testVector1_ = Eigen::Vector3d( 1.0, 2.0, 3.0 );
        Eigen::Vector3d testVector2_ = Eigen::Vector3d( -3.74, 3.7, -4.6 );

        double angle = computeAngleBetweenVectors( testVector1_, testVector2_ );
        double cosineOfAngle = computeCosineOfAngleBetweenVectors( testVector1_, testVector2_ );

        // Check if computed angle and cosine-of-angle are correct.
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle,
                           std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( cosineOfAngle < std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( cosineOfAngle, -0.387790156029810,
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( angle, 1.969029256915446,
                                    std::numeric_limits< double >::epsilon( ) );
    }

    // Test 4: Test values for two vectors of length 5, benchmark values computed using Matlab.
    {
        Eigen::VectorXd testVector1_( 5 );
        testVector1_ << 3.26, 8.66, 1.09, 4.78, 9.92;
        Eigen::VectorXd testVector2_( 5 );
        testVector2_ << 1.05, 0.23, 9.01, 3.25, 7.74;

        double angle = computeAngleBetweenVectors( testVector1_, testVector2_ );
        double cosineOfAngle = computeCosineOfAngleBetweenVectors( testVector1_, testVector2_ );

        // Check if computed angle and cosine-of-angle are correct.
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle,
                           std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( cosineOfAngle > std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( cosineOfAngle, 0.603178944723925,
                                    std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( angle, 0.923315587553074, 1.0e-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
