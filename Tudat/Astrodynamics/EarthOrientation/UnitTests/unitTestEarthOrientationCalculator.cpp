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

//! Test whether Tudat interfaces and SOFA (orientation cookbook case 3) provide same rotation matrices for same input.
BOOST_AUTO_TEST_CASE( testEarthOrientationRotationSetupAgainstSofa )
{
    double arcSecondToRadian = 4.848136811095359935899141E-6;

    // Define test time.
    double terrestrialTimeDaysSinceMjd0 = 54195.50075444445;
    double terrestrialTimeSecondsSinceJ2000 = ( terrestrialTimeDaysSinceMjd0-
                                                ( basic_astrodynamics::JULIAN_DAY_ON_J2000 -
                                                  basic_astrodynamics::JULIAN_DAY_AT_0_MJD ) ) * physical_constants::JULIAN_DAY;

    // Define precession/nutation parameters
    double X = 0.0007122647295989105;
    double Y = 4.438525042571229e-05;
    double s = -0.002200475 * arcSecondToRadian;

    // Compute CIRS->GCRS rotation matrix in Tudat, and compare against Sofa example.
    Eigen::Matrix3d cirsToGcrsRotation = calculateRotationFromCirsToGcrs( X, Y, s ).toRotationMatrix( );
    Eigen::Matrix3d expectedCirsToGcrsRotation;
    expectedCirsToGcrsRotation<<     0.999999746339445, -5.138822464778592e-09, -0.0007122647300724212,
            -2.647522726051399e-08,     0.9999999990149748, -4.438524282712702e-05,
            0.0007122647295989105,  4.438525042571229e-05,     0.9999997453544198;

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedCirsToGcrsRotation.transpose( )( i, j ) -
                                          cirsToGcrsRotation( i, j ) ), 1.0E-15 );
        }
    }

    // Define Earth rotation angle
    double era = 0.2324515536620879;

    // Compute TIRS->CIRS rotation matrix in Tudat
    Eigen::Matrix3d tirsToCirsRotation = Eigen::Matrix3d ( calculateRotationFromTirsToCirs( era ) );

    // Define polar motion values
    double xPole = 0.034928200 * arcSecondToRadian;
    double yPole = 0.483316300 * arcSecondToRadian;

    // Compute ITRS->TIRS rotation matrix in Tudat
    Eigen::Matrix3d itrsToTirsRotation = calculateRotationFromItrsToTirs(
                xPole, yPole, getApproximateTioLocator(
                    terrestrialTimeSecondsSinceJ2000 ) ).toRotationMatrix( );

    // Compute ITRS -> GCRS by concatenation and directly, and compare against each other.
    Eigen::Matrix3d itrsToGcrsRotation = cirsToGcrsRotation * tirsToCirsRotation * itrsToTirsRotation;
    Eigen::Matrix3d itrsToGcrsRotationDirect = Eigen::Matrix3d(
                calculateRotationFromItrsToGcrs( X, Y, s, era, xPole, yPole, getApproximateTioLocator(
                                                     terrestrialTimeSecondsSinceJ2000 ) ) );
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( itrsToGcrsRotationDirect( i, j ) -
                                          itrsToGcrsRotation( i, j ) ), 1.0E-14 );
        }
    }

    // Compare Tudat GCRS -> TIRS against Sofa
    Eigen::Matrix3d expectedGcrsToTirs;
    expectedGcrsToTirs<<0.9731043175731277,     0.2303638262477064, -0.0007033328188453794,
            -0.2303637988041795,     0.9731045707355742,  0.0001208885495858678,
            0.0007122647295989105,  4.438525042571229e-05,     0.9999997453544198;
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( expectedGcrsToTirs.transpose( )( i, j ) -
                                          ( cirsToGcrsRotation * tirsToCirsRotation )( i, j ) ), 1.0E-15 );
        }
    }

    // Compare Tudat GCRS -> ITRS against Sofa
    Eigen::Matrix3d expectedGcrsToItrs;
    expectedGcrsToItrs<<    0.973104317697536,     0.2303638262391256, -0.0007031634821983242,
            -0.2303638004560344,     0.9731045706328012,   0.000118545366624876,
            0.000711560162667892,  4.662640399540082e-05,     0.9999997457540244;
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( expectedGcrsToItrs.transpose( )( i, j ) -
                               itrsToGcrsRotation( i, j ), 1.0E-15 );
        }
    }
}

//! Test whether SOFA functions and Tudat functions provide same X, Y, ERA and s as a function of time.
BOOST_AUTO_TEST_CASE( testEarthOrientationAngleFunctionsAgainstSofa )
{
    double arcSecondToRadian = 4.848136811095359935899141E-6;

    // Define input time.
    double terrestrialTimeDaysSinceMjd0 = 54195.50075444445;
    double terrestrialTimeSecondsSinceJ2000 = ( terrestrialTimeDaysSinceMjd0-
                                                ( basic_astrodynamics::JULIAN_DAY_ON_J2000 -
                                                  basic_astrodynamics::JULIAN_DAY_AT_0_MJD ) ) * physical_constants::JULIAN_DAY;

    // Set SOFA values
    double uncorrectedX = 0.0007122638811749685;
    double uncorrectedY = 4.438634561981791e-05;
    double s = -0.002200475 * arcSecondToRadian;
    double era = 0.2324515536620879;

    // Use Tudat functions to compute precession/nutation parameters
    std::pair< Eigen::Vector2d, double > positionOfCipInGcrs = sofa_interface::getPositionOfCipInGcrs(
                terrestrialTimeSecondsSinceJ2000, basic_astrodynamics::JULIAN_DAY_ON_J2000, iau_2006 );

    // Compare SOFA values against Tudat function output
    BOOST_CHECK_SMALL( std::fabs( uncorrectedX - positionOfCipInGcrs.first( 0 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY - positionOfCipInGcrs.first( 1 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( s - positionOfCipInGcrs.second ), 1.0E-15 );

    // Compute Earth rotation angles from EarthOrientationAnglesCalculator object and compare against SOFA
    boost::shared_ptr< tudat::earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            tudat::earth_orientation::createStandardEarthOrientationCalculator( );
    std::pair< Eigen::Vector5d, double > rotationAngles = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs< double >(
                terrestrialTimeSecondsSinceJ2000, basic_astrodynamics::tt_scale );

    // Compare SOFA values against EarthOrientationAnglesCalculator output
    BOOST_CHECK_SMALL( std::fabs( uncorrectedX - positionOfCipInGcrs.first( 0 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY - positionOfCipInGcrs.first( 1 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( s - positionOfCipInGcrs.second ), 1.0E-15 );
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat



