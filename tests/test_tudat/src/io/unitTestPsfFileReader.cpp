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

#include <cctype>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readPsfFile.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/estimation_setup/processPsfFile.h"

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

std::string normalizePsfPictureName( const std::string& rawPictureName )
{
    std::string normalizedName;
    for( const char character : rawPictureName )
    {
        if( !std::isspace( static_cast< unsigned char >( character ) ) )
        {
            normalizedName.push_back( static_cast< char >( std::toupper( static_cast< unsigned char >( character ) ) ) );
        }
    }

    const std::string::size_type plusPosition = normalizedName.find( '+' );
    if( plusPosition == std::string::npos )
    {
        return normalizedName;
    }

    std::string prefix = normalizedName.substr( 0, plusPosition );
    while( prefix.size( ) > 1 && prefix.front( ) == '0' )
    {
        prefix.erase( prefix.begin( ) );
    }

    std::string suffix = normalizedName.substr( plusPosition + 1 );
    while( suffix.size( ) < 2 )
    {
        suffix.insert( suffix.begin( ), '0' );
    }

    return prefix + "+" + suffix;
}

const input_output::psf::RawPsfFileImageContents& getPsfImage( const input_output::psf::RawPsfFileContents& psfFile,
                                                               const std::string& normalizedPictureName )
{
    for( const input_output::psf::RawPsfFileImageContents& image : psfFile.images_ )
    {
        if( normalizePsfPictureName( image.pictureName_ ) == normalizedPictureName )
        {
            return image;
        }
    }
    throw std::runtime_error( "Could not find PSF image '" + normalizedPictureName + "'." );
}

std::shared_ptr< input_output::psf::RawPsfMeasurement > getPsfMeasurement( const input_output::psf::RawPsfFileImageContents& image,
                                                                           const std::string& targetName )
{
    const std::string normalizedTargetName = observation_models::convertStringToUpperCase( targetName );
    for( const std::shared_ptr< input_output::psf::RawPsfMeasurement >& measurement : image.measurements_ )
    {
        if( measurement != nullptr && observation_models::convertStringToUpperCase( measurement->imageName_ ) == normalizedTargetName )
        {
            return measurement;
        }
    }
    throw std::runtime_error( "Could not find PSF measurement for target '" + targetName + "'." );
}

Eigen::Vector2d convertUnitVectorToRightAscensionAndDeclinationDegrees( const Eigen::Vector3d& unitVector )
{
    double rightAscension = std::atan2( unitVector.y( ), unitVector.x( ) );
    if( rightAscension < 0.0 )
    {
        rightAscension += 2.0 * mathematical_constants::PI;
    }

    return ( Eigen::Vector2d( ) << rightAscension, std::asin( unitVector.z( ) ) ).finished( ) * 180.0 / mathematical_constants::PI;
}

double getSmallestAngularDifferenceDegrees( const double firstAngle, const double secondAngle )
{
    double difference = firstAngle - secondAngle;
    while( difference > 180.0 )
    {
        difference -= 360.0;
    }
    while( difference < -180.0 )
    {
        difference += 360.0;
    }
    return difference;
}

Eigen::Vector2d getRightAscensionAndDeclinationFromPsfPixelLine( const input_output::psf::RawPsfFileContents& psfFile,
                                                                 const input_output::psf::RawPsfFileImageContents& image,
                                                                 const input_output::psf::RawPsfMeasurement& measurement )
{
    const std::shared_ptr< system_models::PsfCameraProjectionModel > cameraProjectionModel =
            input_output::psf::createPsfCameraProjectionModel( psfFile.cameraProperties_.at( image.cameraId_ ) );
    const Eigen::Vector3d cameraFrameUnitVector = cameraProjectionModel->pixelLineToUnitVector( measurement.getEffectivePixelLine( ) );
    const Eigen::Vector3d inertialUnitVector =
            observation_models::getPsfPictureRotationFromInertialToCameraFrame( image ).inverse( ) * cameraFrameUnitVector;
    const double observationTime = observation_models::getPsfPictureObservationTime< double >( image, true );
    // Jacobson's 1991 tabulated astrometry is in the B1950/FK4 direction frame, while PSF pointing is stored as J2000.
    const Eigen::Vector3d b1950UnitVector =
            spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "B1950", observationTime ) * inertialUnitVector;
    return convertUnitVectorToRightAscensionAndDeclinationDegrees( b1950UnitVector.normalized( ) );
}

Eigen::Vector2d getRightAscensionAndDeclinationFromSpiceRelativePosition( const Eigen::Vector3d& relativePosition )
{
    return convertUnitVectorToRightAscensionAndDeclinationDegrees( relativePosition.normalized( ) );
}

struct JacobsonAstrometricReference {
    std::string pictureName_;
    std::string targetName_;
    double julianEphemerisDate_;
    double rightAscensionDegrees_;
    double declinationDegrees_;
    double lightTimeSeconds_;
    Eigen::Vector3d neptuneBarycentricVoyagerPositionKilometers_;
};

}  // namespace

BOOST_AUTO_TEST_SUITE( test_psf_file_reader )

BOOST_AUTO_TEST_CASE( testSinglePsfFileReader )
{
    const std::string file = paths::getTudatTestDataPath( ) + "/psf/psf_vgr2_neptune.txt";
    const input_output::psf::RawPsfFileContents psfFile = input_output::psf::readPsfFile( file );

    BOOST_CHECK_EQUAL( psfFile.spacecraftId_, "VGR2" );
    BOOST_CHECK_EQUAL( psfFile.numberOfCameras_, 2 );
    BOOST_CHECK_EQUAL( psfFile.equinox_, 2000 );
    BOOST_CHECK_EQUAL( psfFile.psfId_, "NEPTUNE" );
    BOOST_CHECK_EQUAL( psfFile.psfProgram_, "ODAP" );
    BOOST_CHECK_EQUAL( psfFile.psfGenerationTimeUtcString_, "2010 MAY 18 14:52:03" );
    BOOST_REQUIRE_EQUAL( psfFile.psfComments_.size( ), 3 );
    BOOST_CHECK_EQUAL( psfFile.psfComments_.at( 0 ), "Voyager PSF used in creating the Neptunian" );

    BOOST_REQUIRE_EQUAL( psfFile.cameraProperties_.size( ), 2 );
    BOOST_REQUIRE( psfFile.cameraProperties_.count( "A" ) == 1 );
    BOOST_REQUIRE( psfFile.cameraProperties_.count( "B" ) == 1 );

    const input_output::psf::RawPsfCameraProperties& cameraA = psfFile.cameraProperties_.at( "A" );
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

    const input_output::psf::RawPsfCameraProperties& cameraB = psfFile.cameraProperties_.at( "B" );
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

BOOST_AUTO_TEST_CASE( testPsfFileReaderJacobson1991Consistency )
{
    const std::string testDataPath = paths::getTudatTestDataPath( );
    const input_output::psf::RawPsfFileContents psfFile = input_output::psf::readPsfFile( testDataPath + "/psf/psf_vgr2_neptune.txt" );

    spice_interface::loadStandardSpiceKernels( std::vector< std::string >( ) );
    spice_interface::loadSpiceKernelInTudat( testDataPath + "/spice/vgr2_nep097.bsp" );

    const std::vector< JacobsonAstrometricReference > astrometricReferences = {
        { "2890B+34",
          "TRITON",
          2447480.34066378,
          297.31160,
          -18.32821,
          1368.242,
          Eigen::Vector3d( -178334616.8, 346108595.1, 128939829.6 ) },
        { "4994B+58",
          "TRITON",
          2447550.48720449,
          297.32997,
          -18.31099,
          1029.258,
          Eigen::Vector3d( -134196154.2, 260416600.6, 97025687.3 ) },
        { "5912B+21", "TRITON", 2447581.06662806, 297.30440, -18.37665, 882.443, Eigen::Vector3d( -114954294.2, 223070441.2, 83117675.7 ) },
        { "7773B+53", "TRITON", 2447643.11764157, 297.24355, -18.24112, 581.563, Eigen::Vector3d( -75908145.9, 147301365.2, 54901517.8 ) },
        { "8585B+59", "TRITON", 2447670.18760590, 297.18558, -18.43828, 452.283, Eigen::Vector3d( -58873770.6, 114250315.9, 42592968.9 ) },
        { "9180B+55", "TRITON", 2447690.01869003, 297.20556, -18.19753, 355.223, Eigen::Vector3d( -46394034.5, 90037807.3, 33576066.8 ) }
    };

    const double maximumAstrometricDifference = 2.0E-3;
    for( const JacobsonAstrometricReference& reference : astrometricReferences )
    {
        const input_output::psf::RawPsfFileImageContents& image = getPsfImage( psfFile, reference.pictureName_ );
        const std::shared_ptr< input_output::psf::RawPsfMeasurement > measurement = getPsfMeasurement( image, reference.targetName_ );

        const Eigen::Vector2d pixelLineRightAscensionAndDeclination =
                getRightAscensionAndDeclinationFromPsfPixelLine( psfFile, image, *measurement );
        BOOST_CHECK_SMALL(
                getSmallestAngularDifferenceDegrees( pixelLineRightAscensionAndDeclination.x( ), reference.rightAscensionDegrees_ ),
                maximumAstrometricDifference );
        BOOST_CHECK_SMALL( pixelLineRightAscensionAndDeclination.y( ) - reference.declinationDegrees_, maximumAstrometricDifference );

        const double psfObservationJulianEphemerisDate = basic_astrodynamics::JULIAN_DAY_ON_J2000 +
                observation_models::getPsfPictureObservationTime< double >( image, true ) / physical_constants::JULIAN_DAY;
        BOOST_CHECK_SMALL( psfObservationJulianEphemerisDate - reference.julianEphemerisDate_, 1.0E-6 );
    }

    const double maximumPaperStateReconstructionDifference = 1.0E-3;
    for( const JacobsonAstrometricReference& reference : astrometricReferences )
    {
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );
        const Eigen::Vector3d tritonPositionAtEmission = spice_interface::getBodyCartesianPositionAtEpoch(
                "801", "8", "B1950", "NONE", receptionTime - reference.lightTimeSeconds_ );
        const Eigen::Vector3d paperNeptuneBarycentricVoyagerPosition = 1000.0 * reference.neptuneBarycentricVoyagerPositionKilometers_;
        const Eigen::Vector3d targetPositionFromPaperVoyagerState = tritonPositionAtEmission - paperNeptuneBarycentricVoyagerPosition;
        const Eigen::Vector2d stateReconstructedRightAscensionAndDeclination =
                getRightAscensionAndDeclinationFromSpiceRelativePosition( targetPositionFromPaperVoyagerState );

        BOOST_CHECK_SMALL( getSmallestAngularDifferenceDegrees( stateReconstructedRightAscensionAndDeclination.x( ),
                                                                reference.rightAscensionDegrees_ ),
                           maximumPaperStateReconstructionDifference );
        BOOST_CHECK_SMALL( stateReconstructedRightAscensionAndDeclination.y( ) - reference.declinationDegrees_,
                           maximumPaperStateReconstructionDifference );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
