/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/leastSquaresEstimation.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace unit_tests
{

Eigen::Matrix3d getTestRotationMatrix( const double psiPerturbation, const double thetaPerturbation, const double phiPerturbation )
{
    return ( Eigen::AngleAxisd( 1.2 + psiPerturbation, Eigen::Vector3d::UnitZ( ) ) *
             Eigen::AngleAxisd( 0.3 + thetaPerturbation, Eigen::Vector3d::UnitX( ) ) *
             Eigen::AngleAxisd( -0.4 + phiPerturbation, Eigen::Vector3d::UnitZ( ) ) )
            .toRotationMatrix( );
}

Eigen::Quaterniond perturbQuaternionEntry( const Eigen::Quaterniond originalQuaternion, const int index, const double perturbation )
{
    double perturbedValue, newq0;
    Eigen::Quaterniond newQuaternion;
    if( index == 1 )
    {
        perturbedValue = originalQuaternion.x( ) + perturbation;
        newq0 = std::sqrt( 1.0 - perturbedValue * perturbedValue - originalQuaternion.y( ) * originalQuaternion.y( ) -
                           originalQuaternion.z( ) * originalQuaternion.z( ) );
        newQuaternion = Eigen::Quaterniond( newq0, perturbedValue, originalQuaternion.y( ), originalQuaternion.z( ) );
    }
    else if( index == 2 )
    {
        perturbedValue = originalQuaternion.y( ) + perturbation;
        newq0 = std::sqrt( 1.0 - perturbedValue * perturbedValue - originalQuaternion.x( ) * originalQuaternion.x( ) -
                           originalQuaternion.z( ) * originalQuaternion.z( ) );
        newQuaternion = Eigen::Quaterniond( newq0, originalQuaternion.x( ), perturbedValue, originalQuaternion.z( ) );
    }
    else if( index == 3 )
    {
        perturbedValue = originalQuaternion.z( ) + perturbation;
        newq0 = std::sqrt( 1.0 - perturbedValue * perturbedValue - originalQuaternion.y( ) * originalQuaternion.y( ) -
                           originalQuaternion.x( ) * originalQuaternion.x( ) );
        newQuaternion = Eigen::Quaterniond( newq0, originalQuaternion.x( ), originalQuaternion.y( ), perturbedValue );
    }
    return newQuaternion;
}

Eigen::Matrix3d getTestRotationMatrix( const double perturbation, const int angleIndex )
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
        Eigen::Vector4d computedResultantQuaternion = linear_algebra::quaternionProduct( firstQuaternion, secondQuaternion );

        // Expected result
        Eigen::Vector4d expectedResultantQuaternion;
        expectedResultantQuaternion << 0.062246819542533, -0.210426307554654, -0.276748204905855, 0.935551459636048;

        // Check if result is correct
        for( unsigned int i = 0; i < 4; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedResultantQuaternion[ i ], expectedResultantQuaternion[ i ], tolerance );
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
        Eigen::Vector4d computedResultantQuaternion = linear_algebra::quaternionProduct( firstQuaternion, secondQuaternion );

        // Expected result
        Eigen::Vector4d expectedResultantQuaternion;
        expectedResultantQuaternion << 0.312036681781879, -0.670033212581165, 0.561681701127148, -0.371755658840092;

        // Check if result is correct
        for( unsigned int i = 0; i < 4; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedResultantQuaternion[ i ], expectedResultantQuaternion[ i ], tolerance );
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
        Eigen::Vector4d computedResultantQuaternion = linear_algebra::quaternionProduct( firstQuaternion, secondQuaternion );

        // Expected result
        Eigen::Vector4d expectedResultantQuaternion;
        expectedResultantQuaternion << -0.438251193562246, -0.386293632078142, -0.483953161080297, -0.651538532273827;

        // Check if result is correct
        for( unsigned int i = 0; i < 4; i++ )
        {
            BOOST_CHECK_CLOSE_FRACTION( computedResultantQuaternion[ i ], expectedResultantQuaternion[ i ], tolerance );
        }
    }
}

//! Test if angle between vectors is computed correctly.
BOOST_AUTO_TEST_CASE( testAngleBetweenVectorFunctions )
{
    // Using declarations.
    using linear_algebra::computeAngleBetweenVectors;
    using linear_algebra::computeCosineOfAngleBetweenVectors;
    using std::cos;
    using std::sqrt;

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
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle, std::numeric_limits< double >::epsilon( ) );
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
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( cosineOfAngle < std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_SMALL( cosineOfAngle + 1.0, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( angle, mathematical_constants::PI, std::sqrt( std::numeric_limits< double >::epsilon( ) ) );
    }

    // Test 3: Test values for two vectors of length 3, benchmark values computed using Matlab.
    {
        Eigen::Vector3d testVector1_ = Eigen::Vector3d( 1.0, 2.0, 3.0 );
        Eigen::Vector3d testVector2_ = Eigen::Vector3d( -3.74, 3.7, -4.6 );

        double angle = computeAngleBetweenVectors( testVector1_, testVector2_ );
        double cosineOfAngle = computeCosineOfAngleBetweenVectors( testVector1_, testVector2_ );

        // Check if computed angle and cosine-of-angle are correct.
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( cosineOfAngle < std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( cosineOfAngle, -0.387790156029810, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( angle, 1.969029256915446, std::numeric_limits< double >::epsilon( ) );
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
        BOOST_CHECK_SMALL( cos( angle ) - cosineOfAngle, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK( cosineOfAngle > std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( cosineOfAngle, 0.603178944723925, std::numeric_limits< double >::epsilon( ) );
        BOOST_CHECK_CLOSE_FRACTION( angle, 0.923315587553074, 1.0e-15 );
    }
}

//! Test whether rotation matrix partial derivatives are correctly implemented
BOOST_AUTO_TEST_CASE( testRotationMatrixDerivatives )
{
    using namespace tudat::linear_algebra;
    Eigen::Matrix3d nominalMatrix = getTestRotationMatrix( 0.3, 0.2, 0.6 );
    Eigen::Quaterniond nominalQuaternion = Eigen::Quaterniond( nominalMatrix );
    Eigen::Vector4d nominalQuaternionVector = convertQuaternionToVectorFormat( nominalQuaternion );

    std::vector< Eigen::Matrix3d > partialDerivatives;
    partialDerivatives.resize( 4 );
    computePartialDerivativeOfRotationMatrixWrtQuaternion( nominalQuaternionVector, partialDerivatives );

    double quaternionPerturbation = 1.0E-7;

    for( unsigned int i = 1; i < 4; i++ )
    {
        Eigen::Quaterniond perturbedQuaternion = perturbQuaternionEntry( nominalQuaternion, i, quaternionPerturbation );
        Eigen::Vector4d perturbedQuaternionVector = convertQuaternionToVectorFormat( perturbedQuaternion );
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

// Mixed observation units must not cause SVD to drop parameters with smaller
// weighted column norms. Exercise both weight representations and augmented systems.
BOOST_AUTO_TEST_CASE( test_LeastSquaresMixedObservationScales )
{
    using namespace tudat::linear_algebra;
    Eigen::MatrixXd design( 6, 3 );
    design << 1, 0, 0, 2, 0, 0, 0, 1, 1, 0, 2, -1, 0, -1, 2, 0, 1, 3;
    const Eigen::Vector3d expected( 3.0, -2.0, 4.0 );
    const Eigen::VectorXd residuals = design * expected;
    Eigen::VectorXd weights( 6 );
    weights << 1.0E-12, 2.0E-12, 1.0E12, 2.0E12, 3.0E12, 4.0E12;

    for( bool correlated : { false, true } )
    {
        Eigen::SparseMatrix< double > weightMatrix( 6, 6 );
        for( int i = 0; i < 6; ++i )
        {
            weightMatrix.insert( i, i ) = weights( i );
        }
        if( correlated )
        {
            weightMatrix.insert( 0, 1 ) = weightMatrix.insert( 1, 0 ) = 0.2E-12;
            weightMatrix.insert( 2, 3 ) = weightMatrix.insert( 3, 2 ) = 0.2E12;
        }
        for( bool withPrior : { false, true } )
        {
            Eigen::MatrixXd prior = Eigen::MatrixXd::Zero( 3, 3 );
            Eigen::MatrixXd softNormal = Eigen::MatrixXd::Zero( 3, 3 );
            if( withPrior )
            {
                prior.diagonal( ) << 5.0E-12, 6.0E12, 7.0E12;
                softNormal.diagonal( ) << 2.0E-12, 3.0E12, 4.0E12;
            }
            for( bool constrained : { false, true } )
            {
                Eigen::MatrixXd constraint( constrained ? 1 : 0, 3 );
                if( constrained )
                {
                    constraint << 1.0, 1.0E-12, 0.0;
                }
                const Eigen::VectorXd constraintRhs = constraint * expected;
                const Eigen::MatrixXd noConsider( 0, 0 );
                const Eigen::VectorXd noDeviation( 0 );
                const auto result = correlated ? performLeastSquaresAdjustmentFromDesignMatrix( design,
                                                                                                residuals,
                                                                                                weightMatrix,
                                                                                                prior,
                                                                                                1.0E8,
                                                                                                constraint,
                                                                                                constraintRhs,
                                                                                                noConsider,
                                                                                                noDeviation,
                                                                                                softNormal,
                                                                                                softNormal * expected,
                                                                                                -expected )
                                               : performLeastSquaresAdjustmentFromDesignMatrix( design,
                                                                                                residuals,
                                                                                                weights,
                                                                                                prior,
                                                                                                1.0E8,
                                                                                                constraint,
                                                                                                constraintRhs,
                                                                                                noConsider,
                                                                                                noDeviation,
                                                                                                softNormal,
                                                                                                softNormal * expected,
                                                                                                -expected );
                BOOST_CHECK_SMALL( ( result.first.head( 3 ) - expected ).norm( ), 1.0E-11 );
                // The returned information matrix retains the caller's coordinates.
                const Eigen::MatrixXd information = design.transpose( ) * weightMatrix * design + prior + softNormal;
                BOOST_CHECK_SMALL( ( result.second.topLeftCorner( 3, 3 ) - information ).norm( ), 1.0E-20 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
