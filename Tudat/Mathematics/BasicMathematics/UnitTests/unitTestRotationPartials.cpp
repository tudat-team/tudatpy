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
 *
 *    References
 *
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <iostream>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/rotationRepresentations.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_mathematics;

BOOST_AUTO_TEST_SUITE( test_rotation_partials )

Eigen::Matrix3d getTestRotationMatrix(
        const double perturbation )
{   
    return ( Eigen::AngleAxisd( -( -1.2 + perturbation ), Eigen::Vector3d::UnitZ( ) ) *
             Eigen::AngleAxisd( -( 0.3 + perturbation ), Eigen::Vector3d::UnitX( ) ) *
             Eigen::AngleAxisd( -( -0.4 + perturbation ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( );
}

Eigen::Matrix3d getTestRotationMatrix(
        const double psiPerturbation,
        const double thetaPerturbation,
        const double phiPerturbation )
{
    return ( Eigen::AngleAxisd( -( -1.2 + psiPerturbation ), Eigen::Vector3d::UnitZ( ) ) *
             Eigen::AngleAxisd( -( 0.3 + thetaPerturbation ), Eigen::Vector3d::UnitX( ) ) *
             Eigen::AngleAxisd( -( -0.4 + phiPerturbation ), Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( );
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

BOOST_AUTO_TEST_CASE( testEulerAngles )
{
    Eigen::Quaterniond quaternion  = Eigen::Quaterniond( getTestRotationMatrix( 0.0 ) );
    Eigen::Vector3d eulerAngles = get313EulerAnglesFromQuaternion( quaternion );


    BOOST_CHECK_SMALL( std::fabs( eulerAngles( 0 ) + 0.4 ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( eulerAngles( 1 ) - 0.3 ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( eulerAngles( 2 ) + 1.2 ), 1.0E-15 );


    Eigen::Quaterniond recomputedQuaternion  = getQuaternionFrom313EulerAngles( eulerAngles );

    Eigen::Quaterniond manualQuaternion = Eigen::Quaterniond(
                Eigen::AngleAxisd( -eulerAngles( 2 ), Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( -eulerAngles( 1 ), Eigen::Vector3d::UnitX( ) ) *
                Eigen::AngleAxisd( -eulerAngles( 0 ), Eigen::Vector3d::UnitZ( ) ) );
    Eigen::Matrix3d rotationMatrix = manualQuaternion.toRotationMatrix( );

    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.w( ) - quaternion.w( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.x( ) - quaternion.x( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.y( ) - quaternion.y( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.z( ) - quaternion.z( ) ), 1.0E-15 );

    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.w( ) - recomputedQuaternion.w( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.x( ) - recomputedQuaternion.x( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.y( ) - recomputedQuaternion.y( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.z( ) - recomputedQuaternion.z( ) ), 1.0E-15 );
}

BOOST_AUTO_TEST_CASE( testEulerAnglePartials )
{
    Eigen::Quaterniond quaternion  = Eigen::Quaterniond( getTestRotationMatrix( 0.0 ) );

    Eigen::Matrix< double, 3, 4 > eulerAnglePartial;
    eulerAnglePartial = calculateEulerAngle313WrtQuaternionPartial( quaternion );


    Eigen::Matrix< double, 3, 4 > eulerAnglePartial2;
    eulerAnglePartial2 = calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
                get313EulerAnglesFromQuaternion( quaternion ) );

    Eigen::Vector4d quaternionInVectorForm = linear_algebra::convertQuaternionToVectorFormat( quaternion );

    Eigen::Matrix< double, 3, 4 > manualAnglePartial = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix< double, 3, 4 > manualQ0PartialContribution = Eigen::Matrix< double, 3, 4 >::Zero( );


    double perturbation = 1.0E-6;
    Eigen::Quaterniond upperturbedQuaternion, downperturbedQuaternion;
    for( int i = 1; i < 4; i++ )
    {
        upperturbedQuaternion = perturbQuaternionEntry( quaternion, i, perturbation );
        downperturbedQuaternion = perturbQuaternionEntry( quaternion, i, -perturbation );

        manualAnglePartial.block( 0, i, 3, 1 ) =
                ( get313EulerAnglesFromQuaternion( upperturbedQuaternion ) -
                  get313EulerAnglesFromQuaternion( downperturbedQuaternion ) ) / ( 2.0 * perturbation );
        manualQ0PartialContribution.block( 0, i, 3, 1 ) = eulerAnglePartial.block( 0, 0, 3, 1 ) *
                -quaternionInVectorForm( i ) / quaternionInVectorForm( 0 );
    }

//    std::cout<<eulerAnglePartial<<std::endl<<std::endl<<
//               manualQ0PartialContribution<<std::endl<<std::endl<<
//               manualAnglePartial<<std::endl<<std::endl<<
//               eulerAnglePartial <<std::endl<<std::endl<< manualQ0PartialContribution - manualAnglePartial<<std::endl<<std::endl<<
//               eulerAnglePartial + ( manualAnglePartial - manualQ0PartialContribution )<<std::endl<<std::endl<<
//               eulerAnglePartial - ( manualAnglePartial - manualQ0PartialContribution )<<std::endl;


    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 1; j < 4; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( eulerAnglePartial( i, j )  -eulerAnglePartial2( i, j ) ), 1.0E-15 );
            BOOST_CHECK_SMALL( std::fabs( eulerAnglePartial( i, j ) -
                                          ( manualAnglePartial( i, j ) - manualQ0PartialContribution( i, j ) ) ), 1.0E-9 );
        }
        BOOST_CHECK_SMALL( std::fabs( eulerAnglePartial( i, 0 )  -eulerAnglePartial2( i, 0 ) ), 1.0E-15 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




