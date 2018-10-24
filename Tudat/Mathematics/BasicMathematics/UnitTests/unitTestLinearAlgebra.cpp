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

#include <iostream>
#include <cmath>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{

namespace unit_tests
{

Eigen::Matrix3d getTestRotationMatrix(
        const double psiPerturbation,
        const double thetaPerturbation,
        const double phiPerturbation )
{
    return ( Eigen::AngleAxisd( 1.2 + psiPerturbation, Eigen::Vector3d::UnitZ( ) ) *
             Eigen::AngleAxisd( 0.3 + thetaPerturbation, Eigen::Vector3d::UnitX( ) ) *
             Eigen::AngleAxisd( -0.4 + phiPerturbation, Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( );
}

Eigen::Quaterniond perturbQuaternionEntry(
        const Eigen::Quaterniond originalQuaternion,
        const int index,
        const double perturbation )
{
    double perturbedValue, newq0;
    Eigen::Quaterniond newQuaternion;
    if( index == 1 )
    {
        perturbedValue = originalQuaternion.x( ) + perturbation;
        newq0 = std::sqrt( 1.0 - perturbedValue * perturbedValue -
                           originalQuaternion.y( ) * originalQuaternion.y( ) -
                           originalQuaternion.z( ) * originalQuaternion.z( ) );
        newQuaternion = Eigen::Quaterniond( newq0, perturbedValue, originalQuaternion.y( ), originalQuaternion.z( ) );

    }
    else if( index == 2 )
    {
        perturbedValue = originalQuaternion.y( ) + perturbation;
        newq0 = std::sqrt( 1.0 - perturbedValue * perturbedValue -
                           originalQuaternion.x( ) * originalQuaternion.x( ) -
                           originalQuaternion.z( ) * originalQuaternion.z( ) );
        newQuaternion = Eigen::Quaterniond( newq0, originalQuaternion.x( ), perturbedValue, originalQuaternion.z( ) );
    }
    else if( index == 3 )
    {
        perturbedValue = originalQuaternion.z( ) + perturbation;
        newq0 = std::sqrt( 1.0 - perturbedValue * perturbedValue -
                           originalQuaternion.y( ) * originalQuaternion.y( ) -
                           originalQuaternion.x( ) * originalQuaternion.x( ) );
        newQuaternion = Eigen::Quaterniond( newq0, originalQuaternion.x( ), originalQuaternion.y( ), perturbedValue );

    }
    return newQuaternion;
}

Eigen::Matrix3d getTestRotationMatrix(
        const double perturbation,
        const int angleIndex )
{
    if( angleIndex == 0 )
    {
        return getTestRotationMatrix( perturbation, 0.0, 0.0 );
    }
    else if( angleIndex == 1 )
    {
        return getTestRotationMatrix( 0.0, perturbation, 0.0 );
    }
    else if( angleIndex == 2 )
    {
        return getTestRotationMatrix( 0.0, 0.0, perturbation );
    }
    else
    {
        return Eigen::Matrix3d::Zero( );
    }
}

//! Test suite for unit conversion functions.
BOOST_AUTO_TEST_SUITE( test_unit_conversions )

//! Test if quaternion multiplication is computed correctly. The expected results are computed with the
//! built-in MATLAB function quatmultiply.
BOOST_AUTO_TEST_CASE( testQuaternionMultiplication )
{
    // Define tolerance
    double tolerance = 10.0 * std::numeric_limits< double >::epsilon( );

    // Test 1
    {
        // Define quaternions to multiply
        Eigen::Vector4d firstQuaternion;
        firstQuaternion << 1.0, 0.0, 0.0, 0.0;
        Eigen::Vector4d secondQuaternion;
        secondQuaternion << 0.062246819542533, -0.210426307554654, -0.276748204905855, 0.935551459636048;

        // Multiply quaternions
        Eigen::Vector4d computedResultantQuaternion = linear_algebra::quaternionProduct( firstQuaternion,
                                                                                         secondQuaternion );

        // Expected result
        Eigen::Vector4d expectedResultantQuaternion;
        expectedResultantQuaternion << 0.062246819542533, -0.210426307554654, -0.276748204905855, 0.935551459636048;

        // Check if result is correct
        for ( unsigned int i = 0; i < 4; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedResultantQuaternion[ i ], expectedResultantQuaternion[ i ],
                                        tolerance );
        }
    }

    // Test 2
    {
        // Define quaternions to multiply
        Eigen::Vector4d firstQuaternion;
        firstQuaternion << 0.808563085773298, -0.184955787789735, -0.218853569995773, 0.513926266878768;
        Eigen::Vector4d secondQuaternion;
        secondQuaternion << 0.062246819542533, -0.276748204905855, 0.935551459636048, -0.210426307554654;

        // Multiply quaternions
        Eigen::Vector4d computedResultantQuaternion = linear_algebra::quaternionProduct( firstQuaternion,
                                                                                         secondQuaternion );

        // Expected result
        Eigen::Vector4d expectedResultantQuaternion;
        expectedResultantQuaternion << 0.312036681781879, -0.670033212581165, 0.561681701127148, -0.371755658840092;

        // Check if result is correct
        for ( unsigned int i = 0; i < 4; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedResultantQuaternion[ i ], expectedResultantQuaternion[ i ],
                                        tolerance );
        }
    }

    // Test 3
    {
        // Define quaternions to multiply
        Eigen::Vector4d firstQuaternion;
        firstQuaternion << -0.522659045487757, -0.363000140684436, 0.699634559308423, -0.324915225026798;
        Eigen::Vector4d secondQuaternion;
        secondQuaternion << 0.242385025550030, 0.655897240236433, 0.670554074388965, -0.247801418397274;

        // Multiply quaternions
        Eigen::Vector4d computedResultantQuaternion = linear_algebra::quaternionProduct( firstQuaternion,
                                                                                         secondQuaternion );

        // Expected result
        Eigen::Vector4d expectedResultantQuaternion;
        expectedResultantQuaternion << -0.438251193562246, -0.386293632078142, -0.483953161080297, -0.651538532273827;

        // Check if result is correct
        for ( unsigned int i = 0; i < 4; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedResultantQuaternion[ i ], expectedResultantQuaternion[ i ],
                                        tolerance );
        }
    }
}

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

//! Test whether rotation matrix partial derivatives are correctly implemented
BOOST_AUTO_TEST_CASE( testRotationMatrixDerivatives )
{
    using namespace tudat::linear_algebra;
    Eigen::Matrix3d nominalMatrix = getTestRotationMatrix( 0.3, 0.2, 0.6 );
    Eigen::Quaterniond nominalQuaternion = Eigen::Quaterniond( nominalMatrix );
    Eigen::Vector4d nominalQuaternionVector = convertQuaternionToVectorFormat(
                nominalQuaternion );

    std::vector< Eigen::Matrix3d > partialDerivatives;
    partialDerivatives.resize( 4 );
    computePartialDerivativeOfRotationMatrixWrtQuaternion(
                nominalQuaternionVector, partialDerivatives );

    double quaternionPerturbation = 1.0E-7;

    for( unsigned int i = 1; i < 4; i++ )
    {
        Eigen::Quaterniond perturbedQuaternion = perturbQuaternionEntry(
                    nominalQuaternion, i, quaternionPerturbation );
        Eigen::Vector4d perturbedQuaternionVector = convertQuaternionToVectorFormat(
                    perturbedQuaternion );
        Eigen::Matrix3d perturbedRotationMatrix = perturbedQuaternion.toRotationMatrix( );

        Eigen::Matrix3d rotationMatrixChange = perturbedRotationMatrix - nominalMatrix;

        Eigen::Matrix3d computedRotationMatrixChange =
                partialDerivatives[ 0 ] * ( perturbedQuaternionVector[ 0 ] - nominalQuaternionVector[ 0 ] ) +
                partialDerivatives[ i ] * ( perturbedQuaternionVector[ i ] - nominalQuaternionVector[ i ] );

        Eigen::Matrix3d rotationMatrixChangeError = rotationMatrixChange - computedRotationMatrixChange;

        for( unsigned int j = 0; j < 3; j++ )
        {
            for( unsigned int k = 0; k < 3; k++ )
            {
                BOOST_CHECK_SMALL( std::fabs( rotationMatrixChangeError( j, k ) ), 1.0E-6 * quaternionPerturbation );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
