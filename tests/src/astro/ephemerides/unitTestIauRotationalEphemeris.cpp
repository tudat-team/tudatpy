/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ephemerides/iauRotationModel.h"
#include "tudat/interface/spice/spiceInterface.h"

namespace tudat
{
namespace unit_tests
{

using namespace ephemerides;
using namespace spice_interface;

BOOST_AUTO_TEST_SUITE( test_iau_rotation_model )


BOOST_AUTO_TEST_CASE( test_IauUranusRotationModel )
{
    loadSpiceKernelInTudat( paths::getSpiceKernelPath( ) + "/pck00010.tpc" );

    std::string baseFrameOrientation = "J2000";
    std::string targetFrameOrientation = "IAU_Uranus";



    double degreeToRadian = unit_conversions::convertDegreesToRadians( 1.0 );
    double nominalMeridian = 203.81 * degreeToRadian;
    Eigen::Vector2d nominalPole = ( Eigen::Vector2d( ) << 257.311, -15.175 ).finished( ) * degreeToRadian;
    double rotationRate = -501.1600928 * degreeToRadian / physical_constants::JULIAN_DAY;
    Eigen::Vector2d polePrecession = Eigen::Vector2d::Zero( );

    std::map< double, std::pair< double, double > > meridianPeriodicTerms;
    std::map< double, std::pair< Eigen::Vector2d, double > > polePeriodicTerms;

    // Create rotation model
    std::shared_ptr< tudat::ephemerides::IauRotationModel > iauRotationModel =
            std::make_shared< ephemerides::IauRotationModel >(
                baseFrameOrientation, targetFrameOrientation,
                nominalMeridian, nominalPole, rotationRate, polePrecession,
                meridianPeriodicTerms, polePeriodicTerms );

    double testTime = 1.0E9;
    Eigen::Matrix3d spiceRotation =
            computeRotationMatrixBetweenFrames( "J2000", "IAU_Uranus", testTime );
    Eigen::Matrix3d tudatRotation = iauRotationModel->getRotationMatrixToTargetFrame( testTime );

    Eigen::Matrix3d spiceRotationDerivative =
            computeRotationMatrixDerivativeBetweenFrames( "J2000", "IAU_Uranus", testTime );
    Eigen::Matrix3d tudatRotationDerivative = iauRotationModel->getDerivativeOfRotationToTargetFrame( testTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceRotation, tudatRotation, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceRotationDerivative, tudatRotationDerivative, 1.0E-10 );

}

BOOST_AUTO_TEST_CASE( test_IauJupiterRotationModel )
{
    loadSpiceKernelInTudat( paths::getSpiceKernelPath( ) + "/pck00010.tpc" );

    std::string baseFrameOrientation = "J2000";
    std::string targetFrameOrientation = "IAU_Jupiter";

    double degreeToRadian = unit_conversions::convertDegreesToRadians( 1.0 );
    double nominalMeridian = 284.95 * degreeToRadian;
    Eigen::Vector2d nominalPole = ( Eigen::Vector2d( ) << 268.056595, 64.495303 ).finished( ) * degreeToRadian;
    double rotationRate = 870.5360000 * degreeToRadian / physical_constants::JULIAN_DAY;
    Eigen::Vector2d polePrecession = ( Eigen::Vector2d( ) << -0.006499, 0.002413 ).finished( ) * degreeToRadian / physical_constants::JULIAN_CENTURY;


    std::map< double, std::pair< double, double > > meridianPeriodicTerms;
    std::map< double, std::pair< Eigen::Vector2d, double > > polePeriodicTerms;
    polePeriodicTerms[ 4850.4046 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.000117, 0.000050 ).finished( ), degreeToRadian * 99.360714 );
    polePeriodicTerms[ 1191.9605 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.000938, 0.000404 ).finished( ), degreeToRadian * 175.895369 );
    polePeriodicTerms[ 262.5475 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.001432, 0.000617 ).finished( ), degreeToRadian * 300.323162 );
    polePeriodicTerms[ 6070.2476 * degreeToRadian / physical_constants::JULIAN_CENTURY ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.000030, -0.000013 ).finished( ), degreeToRadian * 114.012305 );
    polePeriodicTerms[ 64.3000 * degreeToRadian / physical_constants::JULIAN_CENTURY  ] =
            std::make_pair( degreeToRadian * ( Eigen::Vector2d( ) << 0.002150, 0.000926 ).finished( ), degreeToRadian * 49.511251 );


    // Create rotation model
    std::shared_ptr< tudat::ephemerides::IauRotationModel > iauRotationModel =
            std::make_shared< ephemerides::IauRotationModel >(
                baseFrameOrientation, targetFrameOrientation,
                nominalMeridian, nominalPole, rotationRate, polePrecession,
                meridianPeriodicTerms, polePeriodicTerms );

    double testTime = 1.0E9;
    Eigen::Matrix3d spiceRotation =
            computeRotationMatrixBetweenFrames( "J2000", "IAU_Jupiter", testTime );
    Eigen::Matrix3d tudatRotation = iauRotationModel->getRotationMatrixToTargetFrame( testTime );

    Eigen::Matrix3d spiceRotationDerivative =
            computeRotationMatrixDerivativeBetweenFrames( "J2000", "IAU_Jupiter", testTime );
    Eigen::Matrix3d tudatRotationDerivative = iauRotationModel->getDerivativeOfRotationToTargetFrame( testTime );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceRotation, tudatRotation, 1.0E-10 );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( spiceRotationDerivative, tudatRotationDerivative, 1.0E-10 );

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
