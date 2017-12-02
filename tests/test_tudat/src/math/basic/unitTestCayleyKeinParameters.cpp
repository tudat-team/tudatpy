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

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/testMacros.h"

#include "Tudat/Mathematics/BasicMathematics/cayleyKleinParameters.h"

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_CayleyKleinParameters )

BOOST_AUTO_TEST_CASE( test_CayleyKleinConversions )
{
    double angleAlpha, angleBeta, angleGamma;

    {
        angleAlpha = 0.21;
        angleBeta = -0.321;
        angleGamma = 0.87;

        for( int angleTest = 0; angleTest < 2; angleTest++ )
        {
            std::complex< double > cayleyKleinA, cayleyKleinB;
            Eigen::Matrix3d directRotationMatrix;

            // Test 3-2-3 rotation convention
            if( angleTest == 0 )
            {
                directRotationMatrix =
                        Eigen::AngleAxisd( -angleGamma, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( ) *
                        Eigen::AngleAxisd( -angleBeta, Eigen::Vector3d::UnitY( ) ).toRotationMatrix( ) *
                        Eigen::AngleAxisd( -angleAlpha, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );
                basic_mathematics::convert323EulerAnglesToCayleyKleinParameters(
                            angleAlpha, angleBeta, angleGamma, cayleyKleinA, cayleyKleinB );
            }
            // Test 3-1-3 rotation convention
            else
            {
                basic_mathematics::convert313EulerAnglesToCayleyKleinParameters(
                            angleAlpha, angleBeta, angleGamma, cayleyKleinA, cayleyKleinB );
                directRotationMatrix =
                        Eigen::AngleAxisd( -angleGamma, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( ) *
                        Eigen::AngleAxisd( -angleBeta, Eigen::Vector3d::UnitX( ) ).toRotationMatrix( ) *
                        Eigen::AngleAxisd( -angleAlpha, Eigen::Vector3d::UnitZ( ) ).toRotationMatrix( );
            }
            std::complex< double > cayleyKleinC = -std::conj( cayleyKleinB );
            std::complex< double > cayleyKleinD = std::conj( cayleyKleinA );


            Eigen::Quaterniond directQuaternion = Eigen::Quaterniond( directRotationMatrix );
            std::complex< double > cayleyKleinATest, cayleyKleinBTest;
            basic_mathematics::convertQuaterionToCayleyKleinParameters(
                        directQuaternion, cayleyKleinATest, cayleyKleinBTest );

            BOOST_CHECK_CLOSE_FRACTION( cayleyKleinA.real( ), cayleyKleinATest.real( ), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( cayleyKleinA.imag( ), cayleyKleinATest.imag( ), 10.0 * std::numeric_limits< double >::epsilon( ) );

            BOOST_CHECK_CLOSE_FRACTION( cayleyKleinB.real( ), cayleyKleinBTest.real( ), 10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( cayleyKleinB.imag( ), cayleyKleinBTest.imag( ), 10.0 * std::numeric_limits< double >::epsilon( ) );



            BOOST_CHECK_CLOSE_FRACTION( directQuaternion.w( ), ( 0.5 * ( cayleyKleinA + cayleyKleinD ) ).real( ),
                                        10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( directQuaternion.x( ),
                                        ( 0.5 * mathematical_constants::COMPLEX_I * ( cayleyKleinB + cayleyKleinC ) ).real( ),
                                        10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( directQuaternion.y( ), ( 0.5 * ( cayleyKleinB - cayleyKleinC ) ).real( ),
                                        10.0 * std::numeric_limits< double >::epsilon( ) );
            BOOST_CHECK_CLOSE_FRACTION( directQuaternion.z( ),
                                        ( 0.5 * mathematical_constants::COMPLEX_I * ( cayleyKleinA - cayleyKleinD ) ).real( ),
                                        10.0 * std::numeric_limits< double >::epsilon( ) );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat
