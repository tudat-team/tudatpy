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

#include <limits>
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
    expectedCirsToGcrsRotation <<     0.999999746339445, -5.138822464778592e-09, -0.0007122647300724212,
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

    //getSofaEarthOrientationExamples( 3 );

    // Define Earth rotation angle
    double era = 0.2324515536620879;

    double ut1JulianDay = 2454195.5;
    double ut1FractionOfDay = 0.499999165813831;

    double ut1 = ( ut1JulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 + ut1FractionOfDay ) * physical_constants::JULIAN_DAY;
    Time ut1Long = tudat::Time( ( ut1JulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 24,
                                ut1FractionOfDay * physical_constants::JULIAN_DAY_LONG );

    // Compute TIRS->CIRS rotation matrix in Tudat
    Eigen::Matrix3d tirsToCirsRotation = Eigen::Matrix3d ( calculateRotationFromTirsToCirs( era ) );

    // Define polar motion values
    double xPole = 0.034928200 * arcSecondToRadian;
    double yPole = 0.483316300 * arcSecondToRadian;

    // Compute ITRS->TIRS rotation matrix in Tudat
    Eigen::Matrix3d itrsToTirsRotation = calculateRotationFromItrsToTirs(
                xPole, yPole, getApproximateTioLocator(
                    terrestrialTimeSecondsSinceJ2000 ) ).toRotationMatrix( );

    // Compute ITRS -> GCRS by concatenation and directly, and compare against each other (accurate UT1 representation).
    Eigen::Matrix3d itrsToGcrsRotation = cirsToGcrsRotation * tirsToCirsRotation * itrsToTirsRotation;
    Eigen::Matrix3d itrsToGcrsRotationDirect = Eigen::Matrix3d(
                calculateRotationFromItrsToGcrs( X, Y, s, ut1Long, xPole, yPole, getApproximateTioLocator(
                                                     terrestrialTimeSecondsSinceJ2000 ) ) );
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( itrsToGcrsRotationDirect( i, j ) -
                                          itrsToGcrsRotation( i, j ) ), 1.0E-14 );
        }
    }

    // Compute ITRS -> GCRS by directly, and compare against concatenated computation (inaccurate UT1 representation).
    Eigen::Matrix3d itrsToGcrsRotationDirectInaccurate = Eigen::Matrix3d(
                calculateRotationFromItrsToGcrs( X, Y, s, ut1, xPole, yPole, getApproximateTioLocator(
                                                     terrestrialTimeSecondsSinceJ2000 ) ) );

    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( itrsToGcrsRotationDirectInaccurate( i, j ) -
                                          itrsToGcrsRotation( i, j ) ), 1.0E-12 );
        }
    }

    // Compare Tudat GCRS -> TIRS against Sofa
    Eigen::Matrix3d expectedGcrsToTirs;
    expectedGcrsToTirs << 0.9731043175731277,     0.2303638262477064, -0.0007033328188453794,
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
    expectedGcrsToItrs <<    0.973104317697536,     0.2303638262391256, -0.0007031634821983242,
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
    double terrestrialTimeFullDaysSinceMjd0 = 54195;
    double terrestrialTimeDayFractionsSinceMjd0 = 0.50075444445;

    double terrestrialTimeDaysSinceMjd0 = terrestrialTimeFullDaysSinceMjd0 + terrestrialTimeDayFractionsSinceMjd0;
    double terrestrialTimeSecondsSinceJ2000Inaccruate = ( terrestrialTimeDaysSinceMjd0 -
                                                ( basic_astrodynamics::JULIAN_DAY_ON_J2000 -
                                                  basic_astrodynamics::JULIAN_DAY_AT_0_MJD ) ) * physical_constants::JULIAN_DAY;

    double ut1JulianDay = 2454195.5;
    double ut1FractionOfDay = 0.499999165813831;

    Time ut1TimeSecondsSinceJ2000 = tudat::Time( ( ut1JulianDay - basic_astrodynamics::JULIAN_DAY_ON_J2000 ) * 24,
                                                        ut1FractionOfDay * physical_constants::JULIAN_DAY_LONG );

    // Set SOFA values
    double sofaXValue = 0.0007122647295989105;
    double sofaYValue = 4.438525042571229e-05;

    double sofaXCorrection = 8.484239419416879e-10;
    double sofaYCorrection = -1.095194105626442e-09;

    double uncorrectedX = sofaXValue - sofaXCorrection;
    double uncorrectedY = sofaYValue - sofaYCorrection;

    double s = -0.002200475 * arcSecondToRadian;
    double era = 0.2324515536620879;

    // Use Tudat functions to compute precession/nutation parameters
    std::pair< Eigen::Vector2d, double > positionOfCipInGcrs = sofa_interface::getPositionOfCipInGcrs(
                terrestrialTimeSecondsSinceJ2000Inaccruate, basic_astrodynamics::JULIAN_DAY_ON_J2000,
                basic_astrodynamics::iau_2006 );

    // Compare SOFA values against Tudat function output
    BOOST_CHECK_SMALL( std::fabs( uncorrectedX - positionOfCipInGcrs.first( 0 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY - positionOfCipInGcrs.first( 1 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( s - positionOfCipInGcrs.second ), 1.0E-15 );

    Eigen::Vector2d tudatXYCorrection = createStandardEarthOrientationCalculator( )->getPrecessionNutationCalculator( )->getDailyCorrectionInterpolator( )->interpolate(
                terrestrialTimeSecondsSinceJ2000Inaccruate );

    // Compute Earth rotation angles from EarthOrientationAnglesCalculator object and compare against SOFA
    std::shared_ptr< tudat::earth_orientation::EarthOrientationAnglesCalculator > earthOrientationCalculator =
            tudat::earth_orientation::createStandardEarthOrientationCalculator( );

    // Compare SOFA values against EarthOrientationAnglesCalculator  with double input
    std::pair< Eigen::Vector5d, double > rotationAnglesInaccurate = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs< double >(
                terrestrialTimeSecondsSinceJ2000Inaccruate, basic_astrodynamics::tt_scale );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedX + tudatXYCorrection.x( ) - rotationAnglesInaccurate.first( 0 ) ), 1.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY + tudatXYCorrection.y( ) - rotationAnglesInaccurate.first( 1 ) ), 1.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( s - rotationAnglesInaccurate.first( 2 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( era - tudat::sofa_interface::calculateEarthRotationAngleTemplated< double >( rotationAnglesInaccurate.second ) ), 1.0E-9 );

    rotationAnglesInaccurate = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs< double >(
                    ut1TimeSecondsSinceJ2000.getSeconds< double >( ), basic_astrodynamics::ut1_scale );
    double eraDouble = tudat::sofa_interface::calculateEarthRotationAngleTemplated< double >( rotationAnglesInaccurate.second );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedX + tudatXYCorrection.x( ) - rotationAnglesInaccurate.first( 0 ) ), 1.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY + tudatXYCorrection.y( ) - rotationAnglesInaccurate.first( 1 ) ), 1.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( era - eraDouble ), 1.0E-12 );

    // Compare SOFA values against EarthOrientationAnglesCalculator output with Time input
    std::pair< Eigen::Vector5d, Time > rotationAngles = earthOrientationCalculator->getRotationAnglesFromItrsToGcrs< Time >(
                ut1TimeSecondsSinceJ2000, basic_astrodynamics::ut1_scale );
    double eraTime = tudat::sofa_interface::calculateEarthRotationAngleTemplated< Time >( rotationAngles.second );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedX + tudatXYCorrection.x( ) - rotationAngles.first( 0 ) ), 1.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( uncorrectedY + tudatXYCorrection.y( )- rotationAngles.first( 1 ) ), 1.0E-13 );
    BOOST_CHECK_SMALL( std::fabs( s - rotationAngles.first( 2 ) ), 1.0E-15 );
    BOOST_CHECK_SMALL( std::fabs( era - eraTime ), 1.0E-12 );

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat



