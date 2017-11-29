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

#include <iostream>
#include <iomanip>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>
#include <boost/math/special_functions/legendre.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/BasicMathematics/cayleyKleinParameters.h"
#include "Tudat/Mathematics/BasicMathematics/wignerDMatrices.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_WignerDMatrices )

// Varshalovich, p. 112, Eq. (1) and (2)
BOOST_AUTO_TEST_CASE( test_Wigner_D_Matrices )
{
    int maximumDegree = 2;
    WignerDMatricesCache wignerDMatrixCache( maximumDegree );
    double angleAlpha = 0.0, angleBeta = 0.0, angleGamma = 0.0;
    std::complex< double > cayleyKleinA, cayleyKleinB;

    convert323EulerAnglesToCayleyKleinParameters( -angleAlpha, -angleBeta, -angleGamma, cayleyKleinA, cayleyKleinB );
    wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );

    for( unsigned int i = 0; i <= maximumDegree; i++ )
    {
        Eigen::MatrixXcd dMatrixError = wignerDMatrixCache.getWignerDMatrix( i ) -
                Eigen::MatrixXcd::Identity( 2 * i + 1, 2 * i + 1 );
        for( unsigned int j = 0; j < 2 * i + 1; j++ )
        {
            for( unsigned int k = 0; k < 2 * i + 1; k++ )
            {
                BOOST_CHECK_SMALL( dMatrixError( j, k ).real( ), 100.0 * std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_SMALL( dMatrixError( j, k ).imag( ), 100.0 * std::numeric_limits< double >::epsilon( ) );

            }
        }
    }

    angleAlpha = 0.543;
    angleBeta = 0.0;
    angleGamma = -1.073762;

    convert323EulerAnglesToCayleyKleinParameters( -angleAlpha, -angleBeta, -angleGamma, cayleyKleinA, cayleyKleinB );
    wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );

    for( unsigned int i = 0; i <= maximumDegree; i++ )
    {
        Eigen::MatrixXcd currentDMatrix = wignerDMatrixCache.getWignerDMatrix( i );
        for( unsigned int j = 0; j < 2 * i + 1; j++ )
        {
            int m = j - i;

            for( unsigned int k = 0; k < 2 * i + 1; k++ )
            {
                std::complex< double > expectedDMatrixValue;
                if( j != k )
                {
                    BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).real( ) ), 100.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).imag( ) ), 100.0 * std::numeric_limits< double >::epsilon( ) );
                }
                else
                {
                    expectedDMatrixValue = std::exp( static_cast< double >( m ) * tudat::mathematical_constants::COMPLEX_I * (
                                                         angleAlpha + angleGamma ) );
                    //std::cout<<i<<" "<<j<<" "<<k<<" "<<currentDMatrix( j, k )<<" "<<expectedDMatrixValue<<std::endl;

                    BOOST_CHECK_CLOSE_FRACTION( currentDMatrix( j, k ).real( ), expectedDMatrixValue.real( ),
                                                100.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_CLOSE_FRACTION( currentDMatrix( j, k ).imag( ), expectedDMatrixValue.imag( ),
                                                100.0 * std::numeric_limits< double >::epsilon( ) );
                }
            }
        }

        //std::cout<<std::endl;
    }

    angleAlpha = 0.0;//0.543;
    angleBeta = 0.513483;
    angleGamma = 0.0;//-1.073762;

    double cosBeta = std::cos( angleBeta );
    double sinBeta = std::sin( angleBeta );

    convert323EulerAnglesToCayleyKleinParameters( -angleAlpha, -angleBeta, -angleGamma, cayleyKleinA, cayleyKleinB );
    wignerDMatrixCache.updateMatrices( cayleyKleinA, cayleyKleinB );

    for( unsigned int i = 0; i <= maximumDegree; i++ )
    {
        Eigen::MatrixXcd currentDMatrix = wignerDMatrixCache.getWignerDMatrix( i );
        std::cout<<std::setprecision( 3 )<<wignerDMatrixCache.getWignerDMatrix( i )<<std::endl<<std::endl;

        Eigen::MatrixXd testMatrix = Eigen::MatrixXd( 2 * i + 1, 2 * i + 1);
        if( i == 1 )
        {
            testMatrix << ( 1.0 + cosBeta ) / 2.0, sinBeta / std::sqrt( 2.0 ), ( 1.0 - cosBeta ) / 2.0,
                    -sinBeta / std::sqrt( 2.0 ), cosBeta, sinBeta / std::sqrt( 2.0 ),
                    ( 1.0 - cosBeta ) / 2.0, -sinBeta / std::sqrt( 2.0 ), ( 1.0 + cosBeta ) / 2.0;
        }
        else if( i == 2 )
        {
            testMatrix << std::pow( 1.0 + cosBeta, 2 ) / 4.0, sinBeta * ( 1.0 + cosBeta )/ 2.0, 0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, sinBeta * ( 1.0 - cosBeta ) / 2.0, std::pow( 1.0 - cosBeta, 2 ) / 4.0,
                    -sinBeta * ( 1.0 + cosBeta )/ 2.0, ( 2.0 * cosBeta * cosBeta + cosBeta - 1.0 ) / 2.0, std::sqrt( 1.5 ) * sinBeta * cosBeta, -( 2.0 * cosBeta * cosBeta - cosBeta - 1.0 ) / 2.0, sinBeta * ( 1.0 - cosBeta )/ 2.0,
                    0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, -std::sqrt( 1.5 ) * sinBeta * cosBeta, ( 3.0 * cosBeta * cosBeta - 1.0 ) / 2.0, std::sqrt( 1.5 ) * sinBeta * cosBeta, 0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta,
            -sinBeta * ( 1.0 - cosBeta )/ 2.0, -( 2.0 * cosBeta * cosBeta - cosBeta - 1.0 ) / 2.0, -std::sqrt( 1.5 ) * sinBeta * cosBeta, ( 2.0 * cosBeta * cosBeta + cosBeta - 1.0 ) / 2.0, sinBeta * ( 1.0 + cosBeta )/ 2.0,
                    std::pow( 1.0 - cosBeta, 2 ) / 4.0, - sinBeta * ( 1.0 - cosBeta )/ 2.0, 0.5 * std::sqrt( 1.5 ) * sinBeta * sinBeta, -sinBeta * ( 1.0 + cosBeta ) / 2.0, std::pow( 1.0 + cosBeta, 2 ) / 4.0;

        }

        std::cout<<std::setprecision( 3 )<<testMatrix<<std::endl<<std::endl;


        for( unsigned int j = 0; j < 2 * i + 1; j++ )
        {
            int m = j - i;

            for( unsigned int k = 0; k < 2 * i + 1; k++ )
            {
                if( j == i && k == i )
                {
                    double testValue =  boost::math::legendre_p< double >( i, std::cos( angleBeta ) );
                    std::cout<<"Test: "<<testValue<<std::endl<<std::endl;;
                    BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).real( ) - testValue ),
                                       10000.0 * std::numeric_limits< double >::epsilon( ) );
                    BOOST_CHECK_SMALL( std::fabs( currentDMatrix( j, k ).imag( ) ),
                                       10000.0 * std::numeric_limits< double >::epsilon( ) );

                }
            }
        }
    }


}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
