/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <memory>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readPsfFile.h"

namespace tudat
{
namespace unit_tests
{

namespace
{

void checkClose( const double actual, const double expected, const double tolerance = 1.0E-12 )
{
    BOOST_CHECK_SMALL( std::fabs( actual - expected ), tolerance );
}

void checkVector2Close( const Eigen::Vector2d& actual, const double expectedFirst, const double expectedSecond )
{
    checkClose( actual( 0 ), expectedFirst );
    checkClose( actual( 1 ), expectedSecond );
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_psf_file_reader )

BOOST_AUTO_TEST_CASE( testSinglePsfFileReader )
{
    const std::string file = "/home/dominic/Downloads/psf_vgr2_neptune.txt";
    const input_output::psf::RawPsfFileContents psfFile = input_output::psf::readPsfFile( file );

    BOOST_CHECK_EQUAL( psfFile.spacecraftId_, "VGR2" );
    BOOST_CHECK_EQUAL( psfFile.numberOfCameras_, 2 );
    BOOST_CHECK_EQUAL( psfFile.equinox_, 2000 );
    BOOST_CHECK_EQUAL( psfFile.psfId_, "NEPTUNE" );
    BOOST_CHECK_EQUAL( psfFile.psfProgram_, "ODAP" );
    BOOST_CHECK_EQUAL( psfFile.psfGenerationTimeUtcString_, "2010 MAY 18 14:52:03" );
    BOOST_REQUIRE_EQUAL( psfFile.psfComments_.size( ), 3 );
    BOOST_CHECK_EQUAL( psfFile.psfComments_.at( 0 ), "Voyager PSF used in creating the Neptunian" );

    BOOST_REQUIRE_EQUAL( psfFile.cameraModels_.size( ), 2 );
    BOOST_REQUIRE( psfFile.cameraModels_.count( "A" ) == 1 );
    BOOST_REQUIRE( psfFile.cameraModels_.count( "B" ) == 1 );

    const input_output::psf::CameraModel& cameraA = psfFile.cameraModels_.at( "A" );
    BOOST_CHECK_EQUAL( cameraA.cameraId_, "A" );
    checkClose( cameraA.focalLengthMm_, 1503.4900 );
    checkVector2Close( cameraA.principalPoint_, 398.030, 399.390 );
    checkClose( cameraA.fieldOfViewBounds_( 0 ), 0.000 );
    checkClose( cameraA.fieldOfViewBounds_( 1 ), 800.000 );
    checkClose( cameraA.fieldOfViewBounds_( 2 ), 0.000 );
    checkClose( cameraA.fieldOfViewBounds_( 3 ), 800.000 );
    checkClose( cameraA.kMatrix_( 0, 0 ), 7.1808000E+01 );
    checkClose( cameraA.kMatrix_( 1, 0 ), -1.0938000E+00 );
    checkClose( cameraA.kMatrix_( 0, 1 ), 7.6998000E-01 );
    checkClose( cameraA.kMatrix_( 1, 1 ), 7.2951000E+01 );
    checkClose( cameraA.kMatrix_( 0, 2 ), -9.3650000E-03 );
    checkClose( cameraA.kMatrix_( 1, 2 ), 6.5570000E-03 );

    const input_output::psf::CameraModel& cameraB = psfFile.cameraModels_.at( "B" );
    BOOST_CHECK_EQUAL( cameraB.cameraId_, "B" );
    checkClose( cameraB.focalLengthMm_, 200.7260 );
    checkVector2Close( cameraB.principalPoint_, 402.631, 404.101 );
    checkClose( cameraB.mountingOffsetsDegrees_( 0 ), -3.0845000E-02 );
    checkClose( cameraB.mountingOffsetsDegrees_( 1 ), -6.7920000E-03 );
    checkClose( cameraB.mountingOffsetsDegrees_( 2 ), 1.7111900E-01 );

    BOOST_REQUIRE_EQUAL( psfFile.images_.size( ), 845 );
    const input_output::psf::RawPsfFileImageContents& firstImage = psfFile.images_.front( );
    BOOST_CHECK_EQUAL( firstImage.pictureName_, "02890B+34" );
    BOOST_CHECK_EQUAL( firstImage.pictureNumber_, 1 );
    BOOST_CHECK_EQUAL( firstImage.endOfExposureTimeUtcString_, "1988 NOV 14 20:09:38.128" );
    BOOST_CHECK_EQUAL( firstImage.cameraId_, "A" );
    checkClose( firstImage.exposureTimeSeconds_, 1.9200 );
    BOOST_CHECK_EQUAL( firstImage.pictureDeletionFlag_, 0 );
    checkClose( firstImage.rightAscensionDegrees_, 298.030911 );
    checkClose( firstImage.declinationDegrees_, -18.310783 );
    checkClose( firstImage.twistDegrees_, -92.712239 );

    BOOST_REQUIRE_EQUAL( firstImage.measurements_.size( ), 3 );
    const std::shared_ptr< input_output::psf::RawPsfMeasurement > neptuneMeasurement = firstImage.measurements_.at( 0 );
    BOOST_REQUIRE( neptuneMeasurement != nullptr );
    BOOST_CHECK_EQUAL( static_cast< int >( neptuneMeasurement->opticalImageType_ ),
                       static_cast< int >( input_output::psf::OpticalImageType::planet ) );
    BOOST_CHECK_EQUAL( neptuneMeasurement->imageName_, "NEPTUNE" );
    BOOST_CHECK_EQUAL( neptuneMeasurement->imageId_, 8 );
    BOOST_CHECK_EQUAL( neptuneMeasurement->useFlag_, 0 );
    checkVector2Close( neptuneMeasurement->observedPixelLine_, 493.159, 176.480 );
    checkVector2Close( neptuneMeasurement->localCorrection_, 1.791, -2.065 );
    checkVector2Close( neptuneMeasurement->sigmaPixelLine_, 2.000, 2.000 );
    checkVector2Close( neptuneMeasurement->getEffectivePixelLine( ), 491.368, 178.545 );

    const std::shared_ptr< input_output::psf::RawPsfMeasurement > tritonMeasurement = firstImage.measurements_.at( 1 );
    BOOST_REQUIRE( tritonMeasurement != nullptr );
    BOOST_CHECK_EQUAL( static_cast< int >( tritonMeasurement->opticalImageType_ ),
                       static_cast< int >( input_output::psf::OpticalImageType::satellite ) );
    BOOST_CHECK_EQUAL( tritonMeasurement->imageName_, "TRITON" );
    BOOST_CHECK_EQUAL( tritonMeasurement->imageId_, 801 );
    checkVector2Close( tritonMeasurement->observedPixelLine_, 400.493, 182.554 );
    checkVector2Close( tritonMeasurement->localCorrection_, 1.400, -2.313 );
    checkVector2Close( tritonMeasurement->sigmaPixelLine_, 0.500, 0.500 );

    const std::shared_ptr< input_output::psf::RawPsfStarMeasurement > starMeasurement =
            std::dynamic_pointer_cast< input_output::psf::RawPsfStarMeasurement >( firstImage.measurements_.at( 2 ) );
    BOOST_REQUIRE( starMeasurement != nullptr );
    BOOST_CHECK_EQUAL( static_cast< int >( starMeasurement->opticalImageType_ ),
                       static_cast< int >( input_output::psf::OpticalImageType::star ) );
    BOOST_CHECK_EQUAL( starMeasurement->imageName_, "170603" );
    BOOST_CHECK_EQUAL( starMeasurement->imageId_, 170603 );
    checkVector2Close( starMeasurement->observedPixelLine_, 192.038, 454.355 );
    checkVector2Close( starMeasurement->localCorrection_, -1.471, -0.428 );
    checkVector2Close( starMeasurement->sigmaPixelLine_, 0.250, 0.250 );
    checkClose( starMeasurement->rightAscensionDegrees_, 298.1430650 );
    checkClose( starMeasurement->declinationDegrees_, -18.3433350 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
