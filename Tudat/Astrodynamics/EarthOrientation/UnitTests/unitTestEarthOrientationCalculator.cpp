/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>


#include "Tudat/Basics/testMacros.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/Astrodynamics/EarthOrientation/UnitTests/sofaEarthOrientationCookbookExamples.cpp"

namespace tudat
{
namespace unit_tests
{

using namespace earth_orientation;

BOOST_AUTO_TEST_SUITE( test_earth_orientation )

// Compare to SOFA Cookbook
BOOST_AUTO_TEST_CASE( testEarthOrientationRotationSetupAgainstSofa )
{
    std::cout<<"A: "<<std::endl<<getSofaEarthOrientationExamples( 2 ) - getSofaEarthOrientationExamples( 1 )<<std::endl<<std::endl;
    std::cout<<"B: "<<std::endl<<getSofaEarthOrientationExamples( 4 ) - getSofaEarthOrientationExamples( 3 )<<std::endl;


    double terrestrialTimeDaysSinceMjd0 = 54195.50075444445;
    double terrestrialTimeSecondsSinceJ2000 = ( terrestrialTimeDaysSinceMjd0-
                                                ( basic_astrodynamics::JULIAN_DAY_ON_J2000 -
                                                  basic_astrodynamics::JULIAN_DAY_AT_0_MJD ) ) * physical_constants::JULIAN_DAY;


    std::pair< Eigen::Vector2d, double > positionOfCipInGcrs = sofa_interface::getPositionOfCipInGcrs(
                terrestrialTimeSecondsSinceJ2000,
                basic_astrodynamics::JULIAN_DAY_ON_J2000,
                iau_2006 );

    double arcSecondToRadian = 4.848136811095359935899141E-6;

    double uncorrectedX = 0.0007122638811749685;
    double uncorrectedY = 4.438634561981791e-05;
    double s = -0.002200475 * arcSecondToRadian;

    BOOST_CHECK_SMALL( std::fabs( uncorrectedX - positionOfCipInGcrs.first( 0 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY - positionOfCipInGcrs.first( 1 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( s - positionOfCipInGcrs.second ), 1.0E-15 );

    //  CIP offsets wrt IAU 2006/2000A (mas->radians).
    double DX06 =  0.1750 * arcSecondToRadian / 1000.0;
    double DY06 = -0.2259 * arcSecondToRadian / 1000.0;

    //double X = uncorrectedX + DX06;
    //double Y = uncorrectedY + DY06;

    double X = 0.0007122647295252042;
    double Y = 4.43852488746973e-05;

    Eigen::Matrix3d cirsToGcrsRotation = calculateRotationFromCirsToGcrs( X, Y, s ).toRotationMatrix( );

    Eigen::Matrix3d expectedCirsToGcrsRotation;
    expectedCirsToGcrsRotation<<0.999999746339445, -5.138822464778592e-09, -0.0007122647299987151,
            -2.647522615722986e-08, 0.9999999990149748, -4.438524127611243e-05,
            0.0007122647295252042, 4.43852488746973e-05, 0.9999997453544199;


    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedCirsToGcrsRotation.transpose( )( i, j ) -
                                          cirsToGcrsRotation( i, j ) ), 1.0E-15 );
        }
    }

    std::cout<<"Testing: "<<( expectedCirsToGcrsRotation.transpose( )-cirsToGcrsRotation )<<std::endl<<std::endl;

    double era = 13.318492966097 * mathematical_constants::PI / 180.0;
    Eigen::Matrix3d tirsToCirsRotation = Eigen::Matrix3d ( calculateRotationFromTirsToCirs( era ) );

    double xPole = 0.034928200 * arcSecondToRadian;
    double yPole = 0.483316300 * arcSecondToRadian;

    Eigen::Matrix3d itrsToTirsRotation = calculateRotationFromItrsToTirs( xPole, yPole, getApproximateTioLocator(
                                                                              terrestrialTimeSecondsSinceJ2000 ) ).toRotationMatrix( );

    Eigen::Matrix3d itrsToGcrsRotation = cirsToGcrsRotation * tirsToCirsRotation * itrsToTirsRotation;

    Eigen::Matrix3d itrsToGcrsRotationDirect = Eigen::Matrix3d(
                calculateRotationFromItrsToGcrs( X, Y, s, era, xPole, yPole, getApproximateTioLocator(
                                                     terrestrialTimeSecondsSinceJ2000 ) ) );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( itrsToGcrsRotationDirect, itrsToGcrsRotation, 1.0E-10 );

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( itrsToGcrsRotationDirect( i, j ) -
                               itrsToGcrsRotation( i, j ) ), 1.0E-14 );
        }
    }

    Eigen::Matrix3d expectedGcrsToTirs;

    expectedGcrsToTirs<<+0.973104317573127,+0.230363826247709,-0.000703332818416,
            -0.230363798804181,0.973104570735574,0.000120888551078,
            +0.000712264729525,+0.000044385248875,+0.999999745354420;

    std::cout<<( expectedGcrsToTirs.transpose( ) - cirsToGcrsRotation * tirsToCirsRotation )<<std::endl<<std::endl;

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedGcrsToTirs.transpose( )( i, j ) -
                               ( cirsToGcrsRotation * tirsToCirsRotation )( i, j ) ), 1.0E-14 );
        }
    }

    Eigen::Matrix3d expectedGcrsToItrs;
    expectedGcrsToItrs<< 0.973104317697536,0.230363826239128,-0.000703163481769,
            -0.230363800456036,0.973104570632801,0.000118545368117,
            0.000711560162594,0.000046626402444,0.999999745754024;

    std::cout<<"Testing total: "<< ( expectedGcrsToItrs.transpose( ) - itrsToGcrsRotation )<<std::endl;

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( expectedGcrsToItrs.transpose( )( i, j ) -
                               itrsToGcrsRotation( i, j ), 1.0E-14 );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat



