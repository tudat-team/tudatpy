/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Mathworks. Legendre - Associated Legendre functions. Help documentation of MATLAB R2012a,
 *        2012a.
 *      Mathworks. Gravitysphericalharmonic - Implement spherical harmonic representation of
 *        planetary gravity. Help documentation of Aerospace Toolbox of MATLAB R2012a, 2012b.
 *      Eberly, D. Spherical Harmonics. Help documentation of Geometric Tools, 2008. Available at
 *        URL http://www.geometrictools.com/Documentation/SphericalHarmonics.pdf. Last access:
 *        09-09-2012.
 *      Heiskanen, W.A., Moritz, H. Physical geodesy. Freeman, 1967.
 *
 */

#define BOOST_TEST_MAIN

#include <array>
#include <iostream>
#include <iomanip>
#include <memory>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/included/unit_test.hpp>
#include <boost/math/special_functions/legendre.hpp>
#include <boost/math/special_functions/factorials.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/basics/testMacros.h"

#include "tudat/math/basic/cayleyKleinParameters.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/sphericalHarmonicTransformations.h"
#include "tudat/math/basic/wignerDMatrices.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_WignerDMatrices )

void updateWignerCacheFromQuaternionVector( WignerDMatricesCache& wignerDMatrixCache, const Eigen::Vector4d& quaternionVector )
{
    const std::complex< double > cayleyKleinA( quaternionVector( 0 ), -quaternionVector( 3 ) );
    const std::complex< double > cayleyKleinB( quaternionVector( 2 ), -quaternionVector( 1 ) );
    wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );
}

//! Compute values of Wigner D-matrix components at edge of blocks, Varschalovich et al., p. 115, Section 4.17, Eq. (8)
void getWignerDValueAtBoundary( const int degree,
                                const int orderToEvaluate,
                                const double angleAlpha,
                                const double angleBeta,
                                const double angleGamma,
                                std::complex< double >& valueOrderEqualsDegree,
                                std::complex< double >& valueOrderEqualsMinusDegree,
                                std::complex< double >& valueOrderPrimeEqualsDegree,
                                std::complex< double >& valueOrderPrimeEqualsMinusDegree )
{
    double factorialTerm = std::sqrt( boost::math::factorial< double >( 2 * degree ) /
                                      ( boost::math::factorial< double >( degree + orderToEvaluate ) *
                                        boost::math::factorial< double >( degree - orderToEvaluate ) ) );

    valueOrderEqualsDegree = factorialTerm * std::pow( std::cos( angleBeta / 2.0 ), degree + orderToEvaluate ) *
            std::pow( -std::sin( angleBeta / 2.0 ), degree - orderToEvaluate ) *
            std::exp( -mathematical_constants::COMPLEX_I * ( orderToEvaluate * angleAlpha + degree * angleGamma ) );

    valueOrderEqualsMinusDegree = factorialTerm * std::pow( std::cos( angleBeta / 2.0 ), degree - orderToEvaluate ) *
            std::pow( std::sin( angleBeta / 2.0 ), degree + orderToEvaluate ) *
            std::exp( -mathematical_constants::COMPLEX_I * ( orderToEvaluate * angleAlpha - degree * angleGamma ) );

    valueOrderPrimeEqualsDegree = factorialTerm * std::pow( std::cos( angleBeta / 2.0 ), degree + orderToEvaluate ) *
            std::pow( std::sin( angleBeta / 2.0 ), degree - orderToEvaluate ) *
            std::exp( -mathematical_constants::COMPLEX_I * ( degree * angleAlpha + orderToEvaluate * angleGamma ) );

    valueOrderPrimeEqualsMinusDegree = factorialTerm * std::pow( std::cos( angleBeta / 2.0 ), degree - orderToEvaluate ) *
            std::pow( -std::sin( angleBeta / 2.0 ), degree + orderToEvaluate ) *
            std::exp( -mathematical_constants::COMPLEX_I * ( -degree * angleAlpha + orderToEvaluate * angleGamma ) );
}

// Varshalovich, p. 112, Eq. (1) and (2)
BOOST_AUTO_TEST_CASE( test_Wigner_D_Matrices )
{
    // Create Wigner D-matrix calculation object
    int maximumDegree = 32;
    WignerDMatricesCache wignerDMatrixCache( maximumDegree );

    // Define angles
    double angleAlpha = 0.0, angleBeta = 0.0, angleGamma = 0.0;
    std::complex< double > cayleyKleinA, cayleyKleinB;

    // Update Wigner D-matrices for zero angles
    convert323EulerAnglesToCayleyKleinParameters( -angleAlpha, -angleBeta, -angleGamma, cayleyKleinA, cayleyKleinB );
    wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );

    // Check that Wigner D-matrices are real unit matrices: Varschalovich et al., p. 112, Section 4.16, Eq. (1)
    for( int i = 0; i <= maximumDegree; i++ )
    {
        Eigen::MatrixXcd dMatrixError = wignerDMatrixCache.getWignerDMatrix( i ) - Eigen::MatrixXcd::Identity( 2 * i + 1, 2 * i + 1 );
        for( int j = 0; j < 2 * i + 1; j++ )
        {
            for( int k = 0; k < 2 * i + 1; k++ )
            {
                BOOST_CHECK_SMALL( dMatrixError( j, k ).real( ), 100.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( dMatrixError( j, k ).imag( ), 100.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }
    }

    // Reset Angles: no intermediate rotation about y-axis
    angleAlpha = 0.543;
    angleBeta = 0.0;
    angleGamma = -1.073762;

    convert323EulerAnglesToCayleyKleinParameters( -angleAlpha, -angleBeta, -angleGamma, cayleyKleinA, cayleyKleinB );
    wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );

    // Check values accorging to Varschalovich et al., p. 112, Section 4.16, Eq. (2)
    for( int i = 0; i <= maximumDegree; i++ )
    {
        Eigen::MatrixXcd currentDMatrix = wignerDMatrixCache.getWignerDMatrix( i );
        for( int j = 0; j < 2 * i + 1; j++ )
        {
            int m = j - i;
            for( int k = 0; k < 2 * i + 1; k++ )
            {
                // Off-diagonal values are zero
                if( j != k )
                {
                    BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).real( ) ), 100.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).imag( ) ), 100.0 * std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    // Compute expected value on diagonal
                    std::complex< double > expectedDMatrixValue = std::exp(
                            -static_cast< double >( m ) * tudat::mathematical_constants::COMPLEX_I * ( angleAlpha + angleGamma ) );
                    BOOST_CHECK_CLOSE_FRACTION( currentDMatrix( j, k ).real( ),
                                                expectedDMatrixValue.real( ),
                                                100.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_CLOSE_FRACTION( currentDMatrix( j, k ).imag( ),
                                                expectedDMatrixValue.imag( ),
                                                100.0 * std::numeric_limits< double >::epsilon( ) );
                }
            }
        }
    }

    // Check Wigner D-Matrix values with intermediate y-rotation, both with and without z-rotation
    for( int useZRotations = 0; useZRotations < 2; useZRotations++ )
    {
        // Set angles.
        angleBeta = 0.543;
        if( useZRotations == 0 )
        {
            angleAlpha = 0.0;
            angleGamma = 0.0;
        }
        else
        {
            angleAlpha = 0.513483;
            angleGamma = -1.073762;
        }
        double cosBeta = std::cos( angleBeta );
        double sinBeta = std::sin( angleBeta );

        // Update Wigner D-Matrices
        convert323EulerAnglesToCayleyKleinParameters( -angleAlpha, -angleBeta, -angleGamma, cayleyKleinA, cayleyKleinB );
        wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );

        for( int i = 0; i <= maximumDegree; i++ )
        {
            Eigen::MatrixXcd currentDMatrix = wignerDMatrixCache.getWignerDMatrix( i );

            // For small orders, and no z-rotations, compute explicit formulations from Varschalovich et al.,
            // Tables 4.4 and 4.6, p. 119
            Eigen::MatrixXd testMatrix = Eigen::MatrixXd( 2 * i + 1, 2 * i + 1 );
            if( useZRotations == 0 )
            {
                if( i == 1 )
                {
                    testMatrix << ( 1.0 + cosBeta ) / 2.0, sinBeta / std::sqrt( 2.0 ), ( 1.0 - cosBeta ) / 2.0, -sinBeta / std::sqrt( 2.0 ),
                            cosBeta, sinBeta / std::sqrt( 2.0 ), ( 1.0 - cosBeta ) / 2.0, -sinBeta / std::sqrt( 2.0 ),
                            ( 1.0 + cosBeta ) / 2.0;
                }
                else if( i == 2 )
                {
                    testMatrix << std::pow( 1.0 + cosBeta, 2 ) / 4.0, sinBeta * ( 1.0 + cosBeta ) / 2.0,
                            0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, sinBeta * ( 1.0 - cosBeta ) / 2.0,
                            std::pow( 1.0 - cosBeta, 2 ) / 4.0, -sinBeta * ( 1.0 + cosBeta ) / 2.0,
                            ( 2.0 * cosBeta * cosBeta + cosBeta - 1.0 ) / 2.0, std::sqrt( 1.5 ) * sinBeta * cosBeta,
                            -( 2.0 * cosBeta * cosBeta - cosBeta - 1.0 ) / 2.0, sinBeta * ( 1.0 - cosBeta ) / 2.0,
                            0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, -std::sqrt( 1.5 ) * sinBeta * cosBeta,
                            ( 3.0 * cosBeta * cosBeta - 1.0 ) / 2.0, std::sqrt( 1.5 ) * sinBeta * cosBeta,
                            0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, -sinBeta * ( 1.0 - cosBeta ) / 2.0,
                            -( 2.0 * cosBeta * cosBeta - cosBeta - 1.0 ) / 2.0, -std::sqrt( 1.5 ) * sinBeta * cosBeta,
                            ( 2.0 * cosBeta * cosBeta + cosBeta - 1.0 ) / 2.0, sinBeta * ( 1.0 + cosBeta ) / 2.0,
                            std::pow( 1.0 - cosBeta, 2 ) / 4.0, -sinBeta * ( 1.0 - cosBeta ) / 2.0,
                            0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, -sinBeta * ( 1.0 + cosBeta ) / 2.0,
                            std::pow( 1.0 + cosBeta, 2 ) / 4.0;
                }
            }

            for( int j = 0; j < 2 * i + 1; j++ )
            {
                int m = j - i;

                for( int k = 0; k < 2 * i + 1; k++ )
                {
                    // Check if diagonal values from Varschalovich et al., p. 114, Section 4.17, Eq (2)
                    if( j == i && k == i )
                    {
                        double testValue = boost::math::legendre_p< double >( i, std::cos( angleBeta ) );
                        BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).real( ) - testValue ),
                                           10000.0 * std::numeric_limits< double >::epsilon( ) );
                        BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).imag( ) ),
                                           10000.0 * std::numeric_limits< double >::epsilon( ) );
                    }

                    if( useZRotations == 0 )
                    {
                        // Check against explicit matrices
                        if( i == 1 || i == 2 )
                        {
                            BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).real( ) - testMatrix( j, k ) ),
                                               10.0 * std::numeric_limits< double >::epsilon( ) );
                        }

                        // Check that imaginary parts are zero for no z-rotations
                        BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).imag( ) ), 10.0 * std::numeric_limits< double >::epsilon( ) );
                    }
                }
            }

            // Check components at edge of blocks, from Varschalovich et al., p. 115, Section 4.17, Eq. (8)
            for( int j = -i; j <= i; j++ )
            {
                std::complex< double > valueOrderEqualsDegree;
                std::complex< double > valueOrderEqualsMinusDegree;
                std::complex< double > valueOrderPrimeEqualsDegree;
                std::complex< double > valueOrderPrimeEqualsMinusDegree;

                getWignerDValueAtBoundary( i,
                                           j,
                                           angleAlpha,
                                           angleBeta,
                                           angleGamma,
                                           valueOrderEqualsDegree,
                                           valueOrderEqualsMinusDegree,
                                           valueOrderPrimeEqualsDegree,
                                           valueOrderPrimeEqualsMinusDegree );

                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( 0, j + i ) - valueOrderEqualsMinusDegree ).real( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( 0, j + i ) - valueOrderEqualsMinusDegree ).imag( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( 2 * i, j + i ) - valueOrderEqualsDegree ).real( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( 2 * i, j + i ) - valueOrderEqualsDegree ).imag( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( j + i, 2 * i ) - valueOrderPrimeEqualsDegree ).real( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( j + i, 2 * i ) - valueOrderPrimeEqualsDegree ).imag( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( j + i, 0 ) - valueOrderPrimeEqualsMinusDegree ).real( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( std::fabs( ( currentDMatrix( j + i, 0 ) - valueOrderPrimeEqualsMinusDegree ).imag( ) ),
                                   10.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( test_Wigner_D_Matrix_Quaternion_Derivatives )
{
    const int maximumDegree = 8;
    const double perturbation = 1.0E-7;

    Eigen::Quaterniond rotation = Eigen::AngleAxisd( 0.37, Eigen::Vector3d( 1.0, 2.0, -0.5 ).normalized( ) ) *
            Eigen::AngleAxisd( -0.28, Eigen::Vector3d( -0.2, 0.7, 1.0 ).normalized( ) );
    rotation.normalize( );

    const Eigen::Vector4d quaternionVector = linear_algebra::convertQuaternionToVectorFormat( rotation );
    std::shared_ptr< WignerDMatricesCache > wignerDMatrixCache = std::make_shared< WignerDMatricesCache >( maximumDegree );
    updateWignerCacheFromQuaternionVector( *wignerDMatrixCache, quaternionVector );

    std::array< std::vector< Eigen::MatrixXcd >, 4 > analyticalDerivatives;
    computeDerivativeOfWignerDMatricesWrtQuaternion( rotation, wignerDMatrixCache, analyticalDerivatives );

    double maximumAbsoluteError = 0.0;
    double maximumRelativeError = 0.0;
    int maximumErrorQuaternionIndex = -1;
    int maximumErrorDegree = -1;

    for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
    {
        Eigen::Vector4d perturbedQuaternionVector = quaternionVector;
        perturbedQuaternionVector( quaternionIndex ) += perturbation;
        WignerDMatricesCache upperWignerDMatrixCache( maximumDegree );
        updateWignerCacheFromQuaternionVector( upperWignerDMatrixCache, perturbedQuaternionVector );

        perturbedQuaternionVector = quaternionVector;
        perturbedQuaternionVector( quaternionIndex ) -= perturbation;
        WignerDMatricesCache lowerWignerDMatrixCache( maximumDegree );
        updateWignerCacheFromQuaternionVector( lowerWignerDMatrixCache, perturbedQuaternionVector );

        for( int degree = 0; degree <= maximumDegree; degree++ )
        {
            const Eigen::MatrixXcd numericalDerivative =
                    ( upperWignerDMatrixCache.getWignerDMatrix( degree ) - lowerWignerDMatrixCache.getWignerDMatrix( degree ) ) /
                    ( 2.0 * perturbation );
            const double absoluteError = ( analyticalDerivatives.at( quaternionIndex ).at( degree ) - numericalDerivative ).norm( );
            const double relativeError = absoluteError / std::max( 1.0, numericalDerivative.norm( ) );

            if( absoluteError > maximumAbsoluteError )
            {
                maximumAbsoluteError = absoluteError;
                maximumRelativeError = relativeError;
                maximumErrorQuaternionIndex = quaternionIndex;
                maximumErrorDegree = degree;
            }
        }
    }

    BOOST_CHECK_MESSAGE( maximumAbsoluteError < 1.0E-7,
                         "Maximum absolute Wigner-D derivative error is " << maximumAbsoluteError << " for q" << maximumErrorQuaternionIndex
                                                                          << " at degree " << maximumErrorDegree );
    BOOST_CHECK_MESSAGE( maximumRelativeError < 1.0E-8,
                         "Maximum relative Wigner-D derivative error is " << maximumRelativeError << " for q" << maximumErrorQuaternionIndex
                                                                          << " at degree " << maximumErrorDegree );
}

BOOST_AUTO_TEST_CASE( test_Wigner_D_TransformedCoefficient_Quaternion_Derivatives )
{
    const int maximumDegree = 5;
    const double perturbation = 1.0E-7;

    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumDegree + 1 );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd::Zero( maximumDegree + 1, maximumDegree + 1 );
    for( int degree = 0; degree <= maximumDegree; degree++ )
    {
        for( int order = 0; order <= degree; order++ )
        {
            cosineCoefficients( degree, order ) = 0.03 * ( degree + 1 ) - 0.02 * order;
            if( order > 0 )
            {
                sineCoefficients( degree, order ) = -0.04 * degree + 0.015 * ( order + 1 );
            }
        }
    }

    Eigen::Quaterniond rotation = Eigen::AngleAxisd( -0.41, Eigen::Vector3d( 0.3, -1.0, 0.8 ).normalized( ) ) *
            Eigen::AngleAxisd( 0.22, Eigen::Vector3d( 1.0, 0.4, -0.2 ).normalized( ) );
    rotation.normalize( );
    const Eigen::Vector4d quaternionVector = linear_algebra::convertQuaternionToVectorFormat( rotation );

    std::shared_ptr< WignerDMatricesCache > wignerDMatrixCache = std::make_shared< WignerDMatricesCache >( maximumDegree );
    updateWignerCacheFromQuaternionVector( *wignerDMatrixCache, quaternionVector );

    std::array< std::vector< Eigen::MatrixXcd >, 4 > analyticalWignerDerivatives;
    computeDerivativeOfWignerDMatricesWrtQuaternion( rotation, wignerDMatrixCache, analyticalWignerDerivatives );

    double maximumCosineError = 0.0;
    double maximumSineError = 0.0;
    for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
    {
        Eigen::MatrixXd analyticalCosinePartial;
        Eigen::MatrixXd analyticalSinePartial;
        transformSphericalHarmonicCoefficientsWithWignerD( cosineCoefficients,
                                                           sineCoefficients,
                                                           analyticalWignerDerivatives.at( quaternionIndex ),
                                                           analyticalCosinePartial,
                                                           analyticalSinePartial );

        Eigen::Vector4d perturbedQuaternionVector = quaternionVector;
        perturbedQuaternionVector( quaternionIndex ) += perturbation;
        WignerDMatricesCache upperWignerDMatrixCache( maximumDegree );
        updateWignerCacheFromQuaternionVector( upperWignerDMatrixCache, perturbedQuaternionVector );
        Eigen::MatrixXd upperCosineCoefficients;
        Eigen::MatrixXd upperSineCoefficients;
        transformSphericalHarmonicCoefficientsWithWignerD( cosineCoefficients,
                                                           sineCoefficients,
                                                           upperWignerDMatrixCache.getWignerDMatrices( ),
                                                           upperCosineCoefficients,
                                                           upperSineCoefficients );

        perturbedQuaternionVector = quaternionVector;
        perturbedQuaternionVector( quaternionIndex ) -= perturbation;
        WignerDMatricesCache lowerWignerDMatrixCache( maximumDegree );
        updateWignerCacheFromQuaternionVector( lowerWignerDMatrixCache, perturbedQuaternionVector );
        Eigen::MatrixXd lowerCosineCoefficients;
        Eigen::MatrixXd lowerSineCoefficients;
        transformSphericalHarmonicCoefficientsWithWignerD( cosineCoefficients,
                                                           sineCoefficients,
                                                           lowerWignerDMatrixCache.getWignerDMatrices( ),
                                                           lowerCosineCoefficients,
                                                           lowerSineCoefficients );

        const Eigen::MatrixXd numericalCosinePartial = ( upperCosineCoefficients - lowerCosineCoefficients ) / ( 2.0 * perturbation );
        const Eigen::MatrixXd numericalSinePartial = ( upperSineCoefficients - lowerSineCoefficients ) / ( 2.0 * perturbation );
        maximumCosineError = std::max( maximumCosineError, ( analyticalCosinePartial - numericalCosinePartial ).cwiseAbs( ).maxCoeff( ) );
        maximumSineError = std::max( maximumSineError, ( analyticalSinePartial - numericalSinePartial ).cwiseAbs( ).maxCoeff( ) );
    }

    BOOST_CHECK_MESSAGE( maximumCosineError < 2.0E-8, "Maximum transformed cosine coefficient derivative error is " << maximumCosineError );
    BOOST_CHECK_MESSAGE( maximumSineError < 2.0E-8, "Maximum transformed sine coefficient derivative error is " << maximumSineError );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
