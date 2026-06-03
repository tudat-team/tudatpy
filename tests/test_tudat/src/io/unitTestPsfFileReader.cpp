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

#include <algorithm>
#include <cctype>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceEphemeris.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/io/readPsfFile.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
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
                                                                 const Eigen::Vector2d& pixelLine )
{
    const std::shared_ptr< system_models::PsfCameraProjectionModel > cameraProjectionModel =
            input_output::psf::createPsfCameraProjectionModel( psfFile.cameraProperties_.at( image.cameraId_ ) );
    const Eigen::Vector3d cameraFrameUnitVector = cameraProjectionModel->pixelLineToUnitVector( pixelLine );
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

Eigen::Vector3d computeApparentDirectionWithStellarAberration( const Eigen::Vector3d& actualDirection,
                                                               const Eigen::Vector3d& observerVelocityDividedBySpeedOfLight )
{
    const double directionDotVelocity = actualDirection.dot( observerVelocityDividedBySpeedOfLight );
    const double velocitySquared = observerVelocityDividedBySpeedOfLight.squaredNorm( );
    const double scaleFactor = -directionDotVelocity + std::sqrt( 1.0 - velocitySquared + directionDotVelocity * directionDotVelocity );
    return ( scaleFactor * actualDirection + observerVelocityDividedBySpeedOfLight ).normalized( );
}

Eigen::Vector2d getRightAscensionAndDeclinationFromAberrationCorrectedPsfPixelLine( const input_output::psf::RawPsfFileContents& psfFile,
                                                                                    const input_output::psf::RawPsfFileImageContents& image,
                                                                                    const Eigen::Vector2d& pixelLine,
                                                                                    const Eigen::Vector3d& observerBarycentricVelocity )
{
    const std::shared_ptr< system_models::PsfCameraProjectionModel > cameraProjectionModel =
            input_output::psf::createPsfCameraProjectionModel( psfFile.cameraProperties_.at( image.cameraId_ ) );
    const Eigen::Vector3d cameraFrameUnitVector = cameraProjectionModel->pixelLineToUnitVector( pixelLine );
    const Eigen::Vector3d apparentInertialUnitVector =
            observation_models::getPsfPictureRotationFromInertialToCameraFrame( image ).inverse( ) * cameraFrameUnitVector;
    const Eigen::Vector3d astrometricInertialUnitVector = computeApparentDirectionWithStellarAberration(
            apparentInertialUnitVector.normalized( ), -observerBarycentricVelocity / physical_constants::SPEED_OF_LIGHT );
    const double observationTime = observation_models::getPsfPictureObservationTime< double >( image, true );
    const Eigen::Vector3d b1950UnitVector =
            spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "B1950", observationTime ) * astrometricInertialUnitVector;
    return convertUnitVectorToRightAscensionAndDeclinationDegrees( b1950UnitVector.normalized( ) );
}

Eigen::Vector2d getRightAscensionAndDeclinationFromObservationModelPixelLine( const input_output::psf::RawPsfFileContents& psfFile,
                                                                              const input_output::psf::RawPsfFileImageContents& image,
                                                                              const Eigen::Vector2d& modelPixelLine )
{
    return getRightAscensionAndDeclinationFromPsfPixelLine( psfFile, image, modelPixelLine );
}

Eigen::Vector2d getObservationFromSetAtTime(
        const std::shared_ptr< observation_models::SingleObservationSet< double, double > >& observationSet,
        const double observationTime )
{
    for( unsigned int i = 0; i < observationSet->getNumberOfObservables( ); ++i )
    {
        if( std::fabs( observationSet->getObservationTime( i ) - observationTime ) < 1.0E-6 )
        {
            return observationSet->getObservation( i );
        }
    }

    throw std::runtime_error( "Could not find PSF observation at requested time." );
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

Eigen::Vector3d estimatePaperVoyagerVelocityKilometersPerSecond( const std::vector< JacobsonAstrometricReference >& astrometricReferences,
                                                                 const unsigned int referenceIndex )
{
    if( astrometricReferences.size( ) < 2 )
    {
        throw std::runtime_error( "Cannot estimate paper Voyager velocity from fewer than two reference states." );
    }

    unsigned int lowerIndex = referenceIndex;
    unsigned int upperIndex = referenceIndex;
    if( referenceIndex == 0 )
    {
        upperIndex = 1;
    }
    else if( referenceIndex == astrometricReferences.size( ) - 1 )
    {
        lowerIndex = referenceIndex - 1;
    }
    else
    {
        lowerIndex = referenceIndex - 1;
        upperIndex = referenceIndex + 1;
    }

    const double lowerTime =
            spice_interface::convertJulianDateToEphemerisTime( astrometricReferences.at( lowerIndex ).julianEphemerisDate_ );
    const double upperTime =
            spice_interface::convertJulianDateToEphemerisTime( astrometricReferences.at( upperIndex ).julianEphemerisDate_ );
    return ( astrometricReferences.at( upperIndex ).neptuneBarycentricVoyagerPositionKilometers_ -
             astrometricReferences.at( lowerIndex ).neptuneBarycentricVoyagerPositionKilometers_ ) /
            ( upperTime - lowerTime );
}

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

    const double maximumAstrometricDifference = 1.0 / 3600.0;
    double maximumPsfObservationJulianDateDifference = 0.0;
    for( unsigned int referenceIndex = 0; referenceIndex < astrometricReferences.size( ); ++referenceIndex )
    {
        const JacobsonAstrometricReference& reference = astrometricReferences.at( referenceIndex );
        const input_output::psf::RawPsfFileImageContents& image = getPsfImage( psfFile, reference.pictureName_ );
        const std::shared_ptr< input_output::psf::RawPsfMeasurement > measurement = getPsfMeasurement( image, reference.targetName_ );
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );
        const Eigen::Quaterniond rotationFromB1950ToJ2000 =
                spice_interface::computeRotationQuaternionBetweenFrames( "B1950", "J2000", receptionTime );
        const Eigen::Vector3d observerBarycentricVelocity =
                spice_interface::getBodyCartesianStateAtEpoch( "8", "SSB", "J2000", "NONE", receptionTime ).segment( 3, 3 ) +
                rotationFromB1950ToJ2000 *
                        ( 1000.0 * estimatePaperVoyagerVelocityKilometersPerSecond( astrometricReferences, referenceIndex ) );

        const Eigen::Vector2d pixelLineRightAscensionAndDeclination = getRightAscensionAndDeclinationFromAberrationCorrectedPsfPixelLine(
                psfFile, image, measurement->getEffectivePixelLine( ), observerBarycentricVelocity );
        BOOST_CHECK_SMALL(
                getSmallestAngularDifferenceDegrees( pixelLineRightAscensionAndDeclination.x( ), reference.rightAscensionDegrees_ ),
                maximumAstrometricDifference );
        BOOST_CHECK_SMALL( pixelLineRightAscensionAndDeclination.y( ) - reference.declinationDegrees_, maximumAstrometricDifference );

        const double psfObservationJulianEphemerisDate = basic_astrodynamics::JULIAN_DAY_ON_J2000 +
                observation_models::getPsfPictureObservationTime< double >( image, true ) / physical_constants::JULIAN_DAY;
        maximumPsfObservationJulianDateDifference =
                std::max( maximumPsfObservationJulianDateDifference,
                          std::fabs( psfObservationJulianEphemerisDate - reference.julianEphemerisDate_ ) );
        BOOST_CHECK_SMALL( psfObservationJulianEphemerisDate - reference.julianEphemerisDate_, 5.0E-9 );
    }

    std::cout << "\nMaximum PSF observation JED difference for selected rows: " << std::scientific << std::setprecision( 16 )
              << maximumPsfObservationJulianDateDifference << " days\n";

    std::cout << "\nManual Jacobson paper-state residuals for Triton:\n"
              << std::setw( 12 ) << "Picture" << std::setw( 14 ) << "dRA [arcsec]" << std::setw( 16 ) << "dDec [arcsec]" << std::setw( 15 )
              << "Sep [arcsec]" << "\n";
    for( const JacobsonAstrometricReference& reference : astrometricReferences )
    {
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );
        const double emissionTime = receptionTime - reference.lightTimeSeconds_;
        const Eigen::Vector3d tritonPositionAtEmission =
                spice_interface::getBodyCartesianPositionAtEpoch( "801", "8", "B1950", "NONE", emissionTime );
        const Eigen::Vector3d paperNeptuneBarycentricVoyagerPosition = 1000.0 * reference.neptuneBarycentricVoyagerPositionKilometers_;
        const Eigen::Vector3d targetPositionFromPaperVoyagerState = tritonPositionAtEmission - paperNeptuneBarycentricVoyagerPosition;
        const Eigen::Vector2d stateReconstructedRightAscensionAndDeclination =
                getRightAscensionAndDeclinationFromSpiceRelativePosition( targetPositionFromPaperVoyagerState );
        const double rightAscensionDifference = getSmallestAngularDifferenceDegrees( stateReconstructedRightAscensionAndDeclination.x( ),
                                                                                     reference.rightAscensionDegrees_ );
        const double declinationDifference = stateReconstructedRightAscensionAndDeclination.y( ) - reference.declinationDegrees_;
        const double angularDifference = std::sqrt(
                std::pow( std::cos( reference.declinationDegrees_ * mathematical_constants::PI / 180.0 ) * rightAscensionDifference, 2.0 ) +
                std::pow( declinationDifference, 2.0 ) );
        std::cout << std::setw( 12 ) << reference.pictureName_ << std::setw( 14 ) << std::fixed << std::setprecision( 3 )
                  << 3600.0 * rightAscensionDifference << std::setw( 16 ) << 3600.0 * declinationDifference << std::setw( 15 )
                  << 3600.0 * angularDifference << "\n";

        BOOST_CHECK_SMALL( rightAscensionDifference, maximumAstrometricDifference );
        BOOST_CHECK_SMALL( declinationDifference, maximumAstrometricDifference );
    }

    observation_models::PsfFileObservationConversionSettings conversionSettings(
            "VGR2", std::map< std::string, std::string >( { { "TRITON", "TRITON" } } ) );
    conversionSettings.useRawImageNameAsBodyNameIfUnmapped_ = false;
    const std::shared_ptr< observation_models::ObservationCollection< double, double > > observationCollection =
            observation_models::createPsfFileObservationCollection< double, double >( psfFile, conversionSettings );

    simulation_setup::SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "VGR2", false );
    bodies.createEmptyBody< double, double >( "TRITON", false );
    const std::shared_ptr< ephemerides::ConstantEphemeris > voyagerEphemeris =
            std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    const std::shared_ptr< ephemerides::ConstantEphemeris > tritonEphemeris =
            std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" );
    bodies.at( "VGR2" )->setEphemeris( voyagerEphemeris );
    bodies.at( "TRITON" )->setEphemeris( tritonEphemeris );
    bodies.at( "VGR2" )->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
            Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "VGR2_Fixed" ) );
    bodies.processBodyFrameDefinitions< double, double >( );
    observation_models::addPsfCamerasToBodies( psfFile, bodies, conversionSettings );

    std::map< std::string, std::shared_ptr< observation_models::ObservationModel< 2, double, double > > > observationModelsByCamera;
    const double maximumModelPixelLineDifference = 1.0;
    const double maximumObservationModelToManualDifference = 1.0E-8;
    std::cout << "\nFull observation model residuals for Triton:\n"
              << std::setw( 12 ) << "Picture" << std::setw( 18 ) << "dRA model-man" << std::setw( 18 ) << "dDec model-man"
              << std::setw( 18 ) << "Sep model-man" << std::setw( 18 ) << "dRA model-data" << std::setw( 18 ) << "dDec model-data"
              << std::setw( 18 ) << "Sep model-data" << std::setw( 13 ) << "dPixel" << std::setw( 13 ) << "dLine" << std::setw( 13 )
              << "dPL norm" << "\n";
    for( unsigned int referenceIndex = 0; referenceIndex < astrometricReferences.size( ); ++referenceIndex )
    {
        const JacobsonAstrometricReference& reference = astrometricReferences.at( referenceIndex );
        const input_output::psf::RawPsfFileImageContents& image = getPsfImage( psfFile, reference.pictureName_ );
        const std::shared_ptr< input_output::psf::RawPsfMeasurement > measurement = getPsfMeasurement( image, reference.targetName_ );
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );

        observation_models::LinkEnds linkEnds;
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId( reference.targetName_ );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId( "VGR2", image.cameraId_ );
        const observation_models::LinkDefinition linkDefinition( linkEnds );
        const std::vector< std::shared_ptr< observation_models::SingleObservationSet< double, double > > > observationSets =
                observationCollection->getSingleLinkAndTypeObservationSets( observation_models::pixel_coordinates, linkDefinition );
        BOOST_REQUIRE_EQUAL( observationSets.size( ), 1 );

        const double observationTime = observation_models::getPsfPictureObservationTime< double >( image, true );
        const Eigen::Vector2d observationSetPixelLine = getObservationFromSetAtTime( observationSets.at( 0 ), observationTime );
        BOOST_CHECK_SMALL( ( observationSetPixelLine - measurement->getEffectivePixelLine( ) ).norm( ), 1.0E-12 );

        if( observationModelsByCamera.count( image.cameraId_ ) == 0 )
        {
            observationModelsByCamera[ image.cameraId_ ] =
                    observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                            observation_models::pixelCoordinatesSettings(
                                    linkDefinition,
                                    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > >( ),
                                    nullptr,
                                    std::make_shared< observation_models::LightTimeConvergenceCriteria >( ),
                                    true ),
                            bodies );
        }

        const double manualEmissionTime = receptionTime - reference.lightTimeSeconds_;
        const Eigen::Vector3d paperNeptuneBarycentricVoyagerPosition = 1000.0 * reference.neptuneBarycentricVoyagerPositionKilometers_;
        const Eigen::Vector3d manualTritonPositionAtEmission =
                spice_interface::getBodyCartesianPositionAtEpoch( "801", "8", "B1950", "NONE", manualEmissionTime );
        const Eigen::Vector3d manualRelativePosition = manualTritonPositionAtEmission - paperNeptuneBarycentricVoyagerPosition;

        const Eigen::Quaterniond rotationFromB1950ToJ2000 =
                spice_interface::computeRotationQuaternionBetweenFrames( "B1950", "J2000", receptionTime );
        const Eigen::Vector6d neptuneStateAtReception =
                spice_interface::getBodyCartesianStateAtEpoch( "8", "SSB", "J2000", "NONE", receptionTime );
        Eigen::Vector6d paperVoyagerState = Eigen::Vector6d::Zero( );
        paperVoyagerState.segment( 0, 3 ) =
                neptuneStateAtReception.segment( 0, 3 ) + rotationFromB1950ToJ2000 * paperNeptuneBarycentricVoyagerPosition;
        paperVoyagerState.segment( 3, 3 ) = neptuneStateAtReception.segment( 3, 3 ) +
                rotationFromB1950ToJ2000 *
                        ( 1000.0 * estimatePaperVoyagerVelocityKilometersPerSecond( astrometricReferences, referenceIndex ) );
        Eigen::Vector6d tritonStateAtEmission = Eigen::Vector6d::Zero( );
        tritonStateAtEmission.segment( 0, 3 ) =
                neptuneStateAtReception.segment( 0, 3 ) + rotationFromB1950ToJ2000 * manualTritonPositionAtEmission;
        voyagerEphemeris->updateConstantState( paperVoyagerState );
        tritonEphemeris->updateConstantState( tritonStateAtEmission );
        bodies.at( "VGR2" )->updateConstantEphemerisDependentMemberQuantities( );
        bodies.at( "VGR2" )->recomputeStateOnNextCall( );
        bodies.at( "TRITON" )->updateConstantEphemerisDependentMemberQuantities( );
        bodies.at( "TRITON" )->recomputeStateOnNextCall( );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d modelPixelLine = observationModelsByCamera.at( image.cameraId_ )
                                                       ->computeIdealObservationsWithLinkEndData(
                                                               observationTime, observation_models::receiver, linkEndTimes, linkEndStates );
        BOOST_REQUIRE_EQUAL( linkEndStates.size( ), 2 );
        const Eigen::Vector3d modelRelativePositionJ2000 = linkEndStates.at( 0 ).segment( 0, 3 ) - linkEndStates.at( 1 ).segment( 0, 3 );
        const Eigen::Vector3d manualRelativePositionJ2000 = rotationFromB1950ToJ2000 * manualRelativePosition;
        BOOST_CHECK_SMALL( ( modelRelativePositionJ2000.normalized( ) - manualRelativePositionJ2000.normalized( ) ).norm( ), 1.0E-14 );

        const Eigen::Vector2d directManualModelPixelLine =
                bodies.at( "VGR2" )
                        ->getVehicleSystems( )
                        ->getCamera( image.cameraId_ )
                        ->calculateObservableFromInertial( computeApparentDirectionWithStellarAberration(
                                                                   manualRelativePositionJ2000.normalized( ),
                                                                   paperVoyagerState.segment( 3, 3 ) / physical_constants::SPEED_OF_LIGHT ),
                                                           observationTime );
        BOOST_CHECK_SMALL( ( modelPixelLine - directManualModelPixelLine ).norm( ), 1.0E-9 );
        const Eigen::Vector2d geometricModelRightAscensionAndDeclination = convertUnitVectorToRightAscensionAndDeclinationDegrees(
                ( spice_interface::computeRotationQuaternionBetweenFrames( "J2000", "B1950", observationTime ) *
                  modelRelativePositionJ2000 )
                        .normalized( ) );
        const Eigen::Vector2d apparentModelRightAscensionAndDeclination =
                getRightAscensionAndDeclinationFromObservationModelPixelLine( psfFile, image, modelPixelLine );
        BOOST_CHECK_SMALL(
                getSmallestAngularDifferenceDegrees( geometricModelRightAscensionAndDeclination.x( ), reference.rightAscensionDegrees_ ),
                maximumAstrometricDifference );
        BOOST_CHECK_SMALL( geometricModelRightAscensionAndDeclination.y( ) - reference.declinationDegrees_, maximumAstrometricDifference );

        const Eigen::Vector2d manualRightAscensionAndDeclination =
                getRightAscensionAndDeclinationFromSpiceRelativePosition( manualRelativePosition );
        const double modelToManualRightAscensionDifference = getSmallestAngularDifferenceDegrees(
                geometricModelRightAscensionAndDeclination.x( ), manualRightAscensionAndDeclination.x( ) );
        const double modelToManualDeclinationDifference =
                geometricModelRightAscensionAndDeclination.y( ) - manualRightAscensionAndDeclination.y( );
        const double modelToManualAngularDifference =
                std::sqrt( std::pow( std::cos( reference.declinationDegrees_ * mathematical_constants::PI / 180.0 ) *
                                             modelToManualRightAscensionDifference,
                                     2.0 ) +
                           std::pow( modelToManualDeclinationDifference, 2.0 ) );
        BOOST_CHECK_SMALL( modelToManualRightAscensionDifference, maximumObservationModelToManualDifference );
        BOOST_CHECK_SMALL( modelToManualDeclinationDifference, maximumObservationModelToManualDifference );
        BOOST_CHECK_SMALL( modelToManualAngularDifference, maximumObservationModelToManualDifference );

        const Eigen::Vector2d observationSetRightAscensionAndDeclination =
                getRightAscensionAndDeclinationFromPsfPixelLine( psfFile, image, observationSetPixelLine );
        const double modelToObservationSetRightAscensionDifference = getSmallestAngularDifferenceDegrees(
                apparentModelRightAscensionAndDeclination.x( ), observationSetRightAscensionAndDeclination.x( ) );
        const double modelToObservationSetDeclinationDifference =
                apparentModelRightAscensionAndDeclination.y( ) - observationSetRightAscensionAndDeclination.y( );
        const double modelToObservationSetAngularDifference =
                std::sqrt( std::pow( std::cos( reference.declinationDegrees_ * mathematical_constants::PI / 180.0 ) *
                                             modelToObservationSetRightAscensionDifference,
                                     2.0 ) +
                           std::pow( modelToObservationSetDeclinationDifference, 2.0 ) );
        const Eigen::Vector2d modelPixelLineDifference = modelPixelLine - observationSetPixelLine;
        std::cout << std::setw( 12 ) << reference.pictureName_ << std::setw( 14 ) << std::fixed << std::setprecision( 3 )
                  << 3600.0 * modelToManualRightAscensionDifference << std::setw( 18 ) << 3600.0 * modelToManualDeclinationDifference
                  << std::setw( 18 ) << 3600.0 * modelToManualAngularDifference << std::setw( 18 )
                  << 3600.0 * modelToObservationSetRightAscensionDifference << std::setw( 18 )
                  << 3600.0 * modelToObservationSetDeclinationDifference << std::setw( 18 )
                  << 3600.0 * modelToObservationSetAngularDifference << std::setw( 13 ) << modelPixelLineDifference.x( ) << std::setw( 13 )
                  << modelPixelLineDifference.y( ) << std::setw( 13 ) << modelPixelLineDifference.norm( ) << "\n";
        BOOST_CHECK_SMALL( modelPixelLineDifference.norm( ), maximumModelPixelLineDifference );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
