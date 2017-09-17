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

using namespace tudat::earth_orientation;
using mathematical_constants::PI;

BOOST_AUTO_TEST_SUITE( test_earth_orientation )

// Compare to SOFA Cookbook
BOOST_AUTO_TEST_CASE( testEarthOrientationRotationSetupAgainstSofa )
{
    std::cout<<"A: "<<std::endl<<getSofaEarthOrientationExamples( 2 ) - getSofaEarthOrientationExamples( 1 )<<std::endl<<std::endl;
    std::cout<<"B: "<<std::endl<<getSofaEarthOrientationExamples( 4 ) - getSofaEarthOrientationExamples( 3 )<<std::endl;
    double terrestrialTimeDaysSinceMjd0 = 54195.500754444445192;
    double terrestrialTimeSecondsSinceJ2000 = ( 54195.500754444445192 -
                                                ( basic_astrodynamics::JULIAN_DAY_ON_J2000 -
                                                  basic_astrodynamics::JULIAN_DAY_AT_0_MJD ) ) * physical_constants::JULIAN_DAY;


    std::pair< Eigen::Vector2d, double > positionOfCipInGcrs = sofa_interface::getPositionOfCipInGcrs(
                terrestrialTimeSecondsSinceJ2000,
                basic_astrodynamics::JULIAN_DAY_ON_J2000,
                iau_2006 );

    double arcSecondToRadian = 4.848136811095359935899141E-6;

    double X = 0.000712264729525;
    double Y = 0.000044385248875;
    double s = -0.002200475 * arcSecondToRadian;

    Eigen::Matrix3d cirsToGcrsRotation = calculateRotationFromCirsToGcrs( X, Y, s ).toRotationMatrix( );

    Eigen::Matrix3d expectedCirsToGcrsRotation;
    expectedCirsToGcrsRotation<<0.999999746339445,-0.000000005138822,-0.000712264729999,
            -0.000000026475226,0.999999999014975,-0.000044385241276,
            0.000712264729525,0.000044385248875,0.999999745354420;

    {
        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedCirsToGcrsRotation.transpose( ), cirsToGcrsRotation, 1.0E-7 );

        TUDAT_CHECK_MATRIX_BASE( expectedCirsToGcrsRotation.transpose( ), cirsToGcrsRotation )
                BOOST_CHECK_SMALL( expectedCirsToGcrsRotation.transpose( ).coeff(row, col) -
                                   cirsToGcrsRotation.coeff(row, col), 1.0E-15 );
    }

    std::cout<<"Testing: "<<( expectedCirsToGcrsRotation.transpose( )-cirsToGcrsRotation )<<std::endl<<std::endl;

    double era = 13.318492966097 * PI / 180.0;
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

    Eigen::Matrix3d expectedGcrsToTirs;

    expectedGcrsToTirs<<+0.973104317573127,+0.230363826247709,-0.000703332818416,
            -0.230363798804181,0.973104570735574,0.000120888551078,
            +0.000712264729525,+0.000044385248875,+0.999999745354420;

    std::cout<<( expectedGcrsToTirs.transpose( ) - cirsToGcrsRotation * tirsToCirsRotation )<<std::endl<<std::endl;

    Eigen::Matrix3d expectedGcrsToItrs;
    expectedGcrsToItrs<< 0.973104317697536,0.230363826239128,-0.000703163481769,
            -0.230363800456036,0.973104570632801,0.000118545368117,
            0.000711560162594,0.000046626402444,0.999999745754024;

    std::cout<<"Testing total: "<< ( expectedGcrsToItrs.transpose( ) - itrsToGcrsRotation )<<std::endl;

    {
        TUDAT_CHECK_MATRIX_BASE( expectedGcrsToItrs.transpose( ), itrsToGcrsRotation )
                BOOST_CHECK_SMALL( expectedGcrsToItrs.transpose( ).coeff(row, col) -
                                   itrsToGcrsRotation.coeff(row, col), 1.0E-14 );
    }
}

BOOST_AUTO_TEST_CASE( testEarthOrientationRotationSetupAgainstSpice )
{
    double testTime = 1.0E1;
    tudat::spice_interface::loadStandardSpiceKernels( );
    tudat::spice_interface::loadSpiceKernelInTudat( tudat::input_output::getSpiceKernelPath( ) + "earth_latest_high_prec.bpc" );
    tudat::spice_interface::loadSpiceKernelInTudat( tudat::input_output::getSpiceKernelPath( ) + "earth_fixed.tf" );

    std::cout<<tudat::spice_interface::computeRotationQuaternionBetweenFrames(
                 "J2000", "EARTH_FIXED", testTime ).toRotationMatrix( )<<std::endl;
}


//BOOST_AUTO_TEST_CASE( testDefaultEarthOrientationCalculator )
//{
//    boost::shared_ptr< EarthOrientationAnglesCalculator > earthOrientationCalculator =
//            createStandardEarthOrientationCalculator( );
//    boost::gregorian::date currentDate( 2002, 04, 05 );
//    double currentFractionOfDay = 0.2;
//    double secondsSinceEpoch =
//            basic_astrodynamics::calculateJulianDaySinceEpoch( currentDate, currentFractionOfDay ) * physical_constants::JULIAN_DAY;
//    std::cout<<earthOrientationCalculator->getRotationAnglesFromItrsToGcrs( secondsSinceEpoch )<<std::endl<<std::endl;
//    std::cout<<earthOrientationCalculator->getRotationAnglesFromItrsToGcrs( secondsSinceEpoch ) * 3600.0 * 180.0 / (
//                   mathematical_constants::PI ) <<std::endl<<std::endl;

//}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat



