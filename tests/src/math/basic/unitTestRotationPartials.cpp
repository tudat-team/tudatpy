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

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>
#include <iostream>

#include <boost/test/unit_test.hpp>


#include <Eigen/Core>

#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/rotationRepresentations.h"

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

    BOOST_CHECK_SMALL( std::fabs( eulerAngles( 0 ) + 1.2 ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( eulerAngles( 1 ) - 0.3 ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( eulerAngles( 2 ) + 0.4 ), 1.0E-15 );


    Eigen::Quaterniond recomputedQuaternion  = getQuaternionFrom313EulerAngles( eulerAngles );

    Eigen::Quaterniond manualQuaternion = Eigen::Quaterniond(
                Eigen::AngleAxisd( -eulerAngles( 0 ), Eigen::Vector3d::UnitZ( ) ) *
                Eigen::AngleAxisd( -eulerAngles( 1 ), Eigen::Vector3d::UnitX( ) ) *
                Eigen::AngleAxisd( -eulerAngles( 2 ), Eigen::Vector3d::UnitZ( ) ) );
//    Eigen::Matrix3d rotationMatrix = manualQuaternion.toRotationMatrix( );

    Eigen::Vector3d newEulerAngles = get313EulerAnglesFromQuaternion( recomputedQuaternion );
//    std::cout<<eulerAngles.transpose( )<<std::endl;
//    std::cout<<newEulerAngles.transpose( )<<std::endl;

//    std::cout<<quaternion.w( )<<" "<<recomputedQuaternion.w( )<<" "<<manualQuaternion.w( )<<std::endl;
//    std::cout<<quaternion.x( )<<" "<<recomputedQuaternion.x( )<<" "<<manualQuaternion.x( )<<std::endl;
//    std::cout<<quaternion.y( )<<" "<<recomputedQuaternion.y( )<<" "<<manualQuaternion.y( )<<std::endl;
//    std::cout<<quaternion.z( )<<" "<<recomputedQuaternion.z( )<<" "<<manualQuaternion.z( )<<std::endl;

    BOOST_CHECK_SMALL( std::fabs( newEulerAngles.x( ) - eulerAngles.x( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( newEulerAngles.y( ) - eulerAngles.y( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( newEulerAngles.z( ) - eulerAngles.z( ) ), 1.0E-15 );

    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.w( ) - quaternion.w( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.x( ) - quaternion.x( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.y( ) - quaternion.y( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( recomputedQuaternion.z( ) - quaternion.z( ) ), 1.0E-15 );

    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.w( ) - recomputedQuaternion.w( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.x( ) - recomputedQuaternion.x( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.y( ) - recomputedQuaternion.y( ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( manualQuaternion.z( ) - recomputedQuaternion.z( ) ), 1.0E-15 );
}


BOOST_AUTO_TEST_CASE( testEulerAngleRetrieval )
{
    {
        const double angleX = 2.1;
        const double angleZ = -0.3;
        const double angleY = 1.45;

        Eigen::Matrix3d rotationMatrix =
                ( Eigen::AngleAxisd( -angleX, Eigen::Vector3d::UnitX( ) ) *
                  Eigen::AngleAxisd( -angleZ, Eigen::Vector3d::UnitZ( ) ) *
                  Eigen::AngleAxisd( -angleY, Eigen::Vector3d::UnitY( ) ) ).toRotationMatrix( );

        Eigen::Vector3d eulerAngles = get132EulerAnglesFromRotationMatrix(
                    rotationMatrix );

        BOOST_CHECK_CLOSE_FRACTION( angleX, eulerAngles( 0 ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( angleZ, eulerAngles( 1 ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( angleY, eulerAngles( 2 ), 1.0E-15 );
    }

    {
        const double angleZ2 = -2.1;
        const double angleX = 0.3;
        const double angleZ1 = 1.45;

        Eigen::Matrix3d rotationMatrix =
                ( Eigen::AngleAxisd( -angleZ2, Eigen::Vector3d::UnitZ( ) ) *
                  Eigen::AngleAxisd( -angleX, Eigen::Vector3d::UnitX( ) ) *
                  Eigen::AngleAxisd( -angleZ1, Eigen::Vector3d::UnitZ( ) ) ).toRotationMatrix( );
        Eigen::Matrix3d manualMatrix =
                ( Eigen::Matrix3d( ) <<
                  std::cos( angleZ2 ), std::sin( angleZ2 ), 0.0,
                  -std::sin( angleZ2 ), std::cos( angleZ2 ), 0.0,
                  0.0, 0.0, 1.0 ).finished( ) *
                ( Eigen::Matrix3d( ) <<
                  1.0, 0.0, 0.0,
                  0.0, std::cos( angleX ), std::sin( angleX ),
                  0.0, -std::sin( angleX ), std::cos( angleX ) ).finished( ) *
                ( Eigen::Matrix3d( ) <<
                  std::cos( angleZ1 ), std::sin( angleZ1 ), 0.0,
                  -std::sin( angleZ1 ), std::cos( angleZ1 ), 0.0,
                  0.0, 0.0, 1.0 ).finished( );
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rotationMatrix, manualMatrix, 1.0E-15 );

        Eigen::Quaterniond rotationQuaternion =
                Eigen::Quaterniond( manualMatrix );
        double quaternionEntryW = 0.5 * std::sqrt( 1.0 + manualMatrix( 0, 0 ) + manualMatrix( 1, 1 ) + manualMatrix( 2, 2 ) );
        Eigen::Quaterniond manualQuaternion = Eigen::Quaterniond(
                    quaternionEntryW,
                    0.5 * ( manualMatrix( 1, 2 ) - manualMatrix( 2, 1 ) ) / ( 2.0 * quaternionEntryW ),
                    0.5 * ( manualMatrix( 2, 0 ) - manualMatrix( 0, 2 ) ) / ( 2.0 *quaternionEntryW ),
                    0.5 * ( manualMatrix( 0, 1 ) - manualMatrix( 1, 0 ) ) / ( 2.0 *quaternionEntryW ) );
        BOOST_CHECK_CLOSE_FRACTION( manualQuaternion.w( ), rotationQuaternion.w( ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( manualQuaternion.x( ), -rotationQuaternion.x( ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( manualQuaternion.y( ), -rotationQuaternion.y( ), 1.0E-15 );
        BOOST_CHECK_CLOSE_FRACTION( manualQuaternion.z( ), -rotationQuaternion.z( ), 1.0E-15 );


        Eigen::Vector3d eulerAngles = get313EulerAnglesFromRotationMatrix( rotationMatrix );

        BOOST_CHECK_CLOSE_FRACTION( angleZ2, eulerAngles( 0 ), 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( angleX, eulerAngles( 1 ), 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( angleZ1, eulerAngles( 2 ), 1.0E-14 );

        eulerAngles = get313EulerAnglesFromQuaternion( rotationQuaternion );

        BOOST_CHECK_CLOSE_FRACTION( angleZ2, eulerAngles( 0 ), 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( angleX, eulerAngles( 1 ), 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( angleZ1, eulerAngles( 2 ), 1.0E-14 );

    }
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

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 1; j < 4; j++ )
        {
            BOOST_CHECK_LT( std::fabs( eulerAnglePartial( i, j )  - eulerAnglePartial2( i, j ) ), 5.0E-14 );
            BOOST_CHECK_LT( std::fabs( eulerAnglePartial( i, j ) -
                                          ( manualAnglePartial( i, j ) - manualQ0PartialContribution( i, j ) ) ), 1.0E-9 );
        }
        BOOST_CHECK_LT( std::fabs( eulerAnglePartial( i, 0 )  - eulerAnglePartial2( i, 0 ) ), 1.0E-14 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat




