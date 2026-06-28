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
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
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

Eigen::Vector2d getRightAscensionAndDeclinationFromSpiceRelativePosition( const Eigen::Vector3d& relativePosition )
{
    return convertUnitVectorToRightAscensionAndDeclinationDegrees( relativePosition.normalized( ) );
}

Eigen::Vector3d computeApparentDirectionWithStellarAberration( const Eigen::Vector3d& actualDirection,
                                                               const Eigen::Vector3d& observerVelocityDividedBySpeedOfLight )
{
    // Apply the same first-order astrometric-to-apparent direction convention used by the PSF image reductions.
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

Eigen::Vector2d getObservationFromDatasetAtTime(
        const std::shared_ptr< observation_models::ObservationDataset< double, double > >& observationDataset,
        const observation_models::LinkDefinition& linkDefinition,
        const double observationTime )
{
    const std::vector< observation_models::ObservationId > matchingObservationIds = observationDataset->getObservationIdsMatchingCondition(
            observation_models::ObservationCondition< double, double >::observableType( observation_models::pixel_coordinates ) &&
            observation_models::ObservationCondition< double, double >::linkDefinition( linkDefinition ) );

    for( const observation_models::ObservationId observationId : matchingObservationIds )
    {
        if( std::fabs( observationDataset->getObservationTime( observationId ) - observationTime ) < 1.0E-6 )
        {
            return observationDataset->getObservationValue( observationId );
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

    // Jacobson provides discrete Voyager positions, but aberration requires velocity. A local finite difference is sufficient
    // for the sub-arcsecond consistency check because it only feeds the stellar-aberration correction.
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

BOOST_AUTO_TEST_SUITE( test_psf_pixel_coordinates_observation_model )

BOOST_AUTO_TEST_CASE( testPsfFileReaderJacobson1991Consistency )
{
    const std::string testDataPath = paths::getTudatTestDataPath( );
    const input_output::psf::RawPsfFileContents psfFile = input_output::psf::readPsfFile( testDataPath + "/psf/psf_vgr2_neptune.txt" );

    spice_interface::loadStandardSpiceKernels( std::vector< std::string >( ) );
    spice_interface::loadSpiceKernelInTudat( testDataPath + "/spice/vgr2_nep097.bsp" );

    // Reference rows transcribed from Jacobson (1991): image id, published B1950/FK4 astrometry, light time,
    // and Voyager position with respect to the Neptune barycenter in the same frame as the paper table.
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
    for( unsigned int referenceIndex = 0; referenceIndex < astrometricReferences.size( ); ++referenceIndex )
    {
        // First validate the file-level interpretation without involving Tudat bodies or observation models.
        const JacobsonAstrometricReference& reference = astrometricReferences.at( referenceIndex );
        const input_output::psf::RawPsfFileImageContents& image = getPsfImage( psfFile, reference.pictureName_ );
        const std::shared_ptr< input_output::psf::RawPsfMeasurement > measurement = getPsfMeasurement( image, reference.targetName_ );
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );

        // The paper-state velocity estimate is in the paper frame, while Tudat's aberration correction uses J2000 states.
        const Eigen::Quaterniond rotationFromB1950ToJ2000 =
                spice_interface::computeRotationQuaternionBetweenFrames( "B1950", "J2000", receptionTime );
        const Eigen::Vector3d observerBarycentricVelocity =
                spice_interface::getBodyCartesianStateAtEpoch( "8", "SSB", "J2000", "NONE", receptionTime ).segment( 3, 3 ) +
                rotationFromB1950ToJ2000 *
                        ( 1000.0 * estimatePaperVoyagerVelocityKilometersPerSecond( astrometricReferences, referenceIndex ) );

        // Direct check: PSF pixel/line + PSF camera calibration + picture pointing should reproduce Jacobson's
        // published astrometry after undoing the stellar-aberration convention used in the image reduction.
        const Eigen::Vector2d pixelLineRightAscensionAndDeclination = getRightAscensionAndDeclinationFromAberrationCorrectedPsfPixelLine(
                psfFile, image, measurement->getEffectivePixelLine( ), observerBarycentricVelocity );
        BOOST_CHECK_SMALL(
                getSmallestAngularDifferenceDegrees( pixelLineRightAscensionAndDeclination.x( ), reference.rightAscensionDegrees_ ),
                maximumAstrometricDifference );
        BOOST_CHECK_SMALL( pixelLineRightAscensionAndDeclination.y( ) - reference.declinationDegrees_, maximumAstrometricDifference );

        // This checks the UTC parsing and internal conversion to ephemeris time used by PSF observation construction.
        const double psfObservationJulianEphemerisDate = basic_astrodynamics::JULIAN_DAY_ON_J2000 +
                observation_models::getPsfPictureObservationTime< double >( image, true ) / physical_constants::JULIAN_DAY;
        BOOST_CHECK_SMALL( psfObservationJulianEphemerisDate - reference.julianEphemerisDate_, 5.0E-9 );
    }

    std::cout << "\nManual Jacobson paper-state residuals for Triton:\n"
              << std::setw( 12 ) << "Picture" << std::setw( 14 ) << "dRA [arcsec]" << std::setw( 16 ) << "dDec [arcsec]" << std::setw( 15 )
              << "Sep [arcsec]" << "\n";
    for( const JacobsonAstrometricReference& reference : astrometricReferences )
    {
        // Independent manual reconstruction: combine Triton's SPICE state at Jacobson's emission epoch with the
        // paper's Voyager position. This intentionally does not use the Tudat pixel-coordinate observation model.
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );
        const double emissionTime = receptionTime - reference.lightTimeSeconds_;

        // Triton's position is sampled at photon emission, but Voyager's tabulated position is at reception.
        const Eigen::Vector3d tritonPositionAtEmission =
                spice_interface::getBodyCartesianPositionAtEpoch( "801", "8", "B1950", "NONE", emissionTime );
        const Eigen::Vector3d paperNeptuneBarycentricVoyagerPosition = 1000.0 * reference.neptuneBarycentricVoyagerPositionKilometers_;
        const Eigen::Vector3d targetPositionFromPaperVoyagerState = tritonPositionAtEmission - paperNeptuneBarycentricVoyagerPosition;

        // Jacobson's published RA/DEC is the geometric astrometric direction in the paper frame.
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

    // Exercise the PSF-to-observation conversion separately from the model calculation below.
    const std::shared_ptr< observation_models::ObservationDataset< double, double > > observationDataset =
            observation_models::createPsfFileObservationDataset< double, double >( psfFile, conversionSettings );

    // The model test below forces only Voyager's state from the paper table. Triton uses a normal Neptune-barycentric
    // SPICE ephemeris, so the model computes the target state and light-time geometry through the standard path.
    simulation_setup::SystemOfBodies bodies( "8", "J2000" );
    bodies.createEmptyBody< double, double >( "8", false );
    bodies.createEmptyBody< double, double >( "VGR2", false );
    bodies.createEmptyBody< double, double >( "TRITON", false );
    const std::shared_ptr< ephemerides::ConstantEphemeris > voyagerEphemeris =
            std::make_shared< ephemerides::ConstantEphemeris >( Eigen::Vector6d::Zero( ), "8", "J2000" );
    bodies.at( "8" )->setEphemeris( std::make_shared< ephemerides::SpiceEphemeris >( "8", "SSB", false, false, false, "J2000" ) );
    bodies.at( "VGR2" )->setEphemeris( voyagerEphemeris );
    bodies.at( "TRITON" )->setEphemeris( std::make_shared< ephemerides::SpiceEphemeris >( "801", "8", false, false, false, "J2000" ) );

    // Picture pointing from the PSF file is camera-specific; the spacecraft body frame can therefore remain identity.
    bodies.at( "VGR2" )->setRotationalEphemeris( std::make_shared< ephemerides::ConstantRotationalEphemeris >(
            Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "VGR2_Fixed" ) );
    bodies.processBodyFrameDefinitions< double, double >( );
    observation_models::addPsfCamerasToBodies( psfFile, bodies, conversionSettings );

    std::map< std::string, std::shared_ptr< observation_models::ObservationModel< 2, double, double > > > observationModelsByCamera;
    const double maximumModelPixelLineDifference = 1.0;
    std::cout << "\nFull observation model residuals for Triton:\n"
              << std::setw( 12 ) << "Picture" << std::setw( 13 ) << "dPixel" << std::setw( 13 ) << "dLine" << std::setw( 13 ) << "dPL norm"
              << "\n";
    for( unsigned int referenceIndex = 0; referenceIndex < astrometricReferences.size( ); ++referenceIndex )
    {
        const JacobsonAstrometricReference& reference = astrometricReferences.at( referenceIndex );
        const input_output::psf::RawPsfFileImageContents& image = getPsfImage( psfFile, reference.pictureName_ );
        const std::shared_ptr< input_output::psf::RawPsfMeasurement > measurement = getPsfMeasurement( image, reference.targetName_ );
        const double receptionTime = spice_interface::convertJulianDateToEphemerisTime( reference.julianEphemerisDate_ );

        // The receiver link end includes the camera id, matching the body/camera created from the PSF camera block.
        observation_models::LinkEnds linkEnds;
        linkEnds[ observation_models::transmitter ] = observation_models::LinkEndId( reference.targetName_ );
        linkEnds[ observation_models::receiver ] = observation_models::LinkEndId( "VGR2", image.cameraId_ );
        const observation_models::LinkDefinition linkDefinition( linkEnds );
        const double observationTime = observation_models::getPsfPictureObservationTime< double >( image, true );
        const Eigen::Vector2d observationSetPixelLine =
                getObservationFromDatasetAtTime( observationDataset, linkDefinition, observationTime );

        // The observation set should store the PSF effective pixel/line exactly; dynamics only enter when the model is evaluated.
        BOOST_CHECK_SMALL( ( observationSetPixelLine - measurement->getEffectivePixelLine( ) ).norm( ), 1.0E-12 );

        // Observation models are camera-specific, so reuse them when multiple selected pictures use the same camera.
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

        const Eigen::Vector3d paperNeptuneBarycentricVoyagerPosition = 1000.0 * reference.neptuneBarycentricVoyagerPositionKilometers_;

        // Convert the paper geometry into J2000 so the observation model uses the same Neptune-barycentric
        // Voyager position as the paper, while Triton's state is provided by its regular SPICE ephemeris.
        const Eigen::Quaterniond rotationFromB1950ToJ2000 =
                spice_interface::computeRotationQuaternionBetweenFrames( "B1950", "J2000", receptionTime );
        const Eigen::Vector6d neptuneStateAtReception =
                spice_interface::getBodyCartesianStateAtEpoch( "8", "SSB", "J2000", "NONE", receptionTime );
        Eigen::Vector6d paperVoyagerState = Eigen::Vector6d::Zero( );

        // Positions must share the Neptune-barycentric origin used by Jacobson's table. The velocity remains barycentric
        // because the stellar-aberration correction depends on the observer's inertial velocity.
        paperVoyagerState.segment( 0, 3 ) = rotationFromB1950ToJ2000 * paperNeptuneBarycentricVoyagerPosition;
        paperVoyagerState.segment( 3, 3 ) = neptuneStateAtReception.segment( 3, 3 ) +
                rotationFromB1950ToJ2000 *
                        ( 1000.0 * estimatePaperVoyagerVelocityKilometersPerSecond( astrometricReferences, referenceIndex ) );

        // Voyager's constant ephemeris is reused for all rows; force cached body states to be recomputed after each update.
        voyagerEphemeris->updateConstantState( paperVoyagerState );
        bodies.at( "VGR2" )->updateConstantEphemerisDependentMemberQuantities( );
        bodies.at( "VGR2" )->recomputeStateOnNextCall( );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;
        const Eigen::Vector2d modelPixelLine = observationModelsByCamera.at( image.cameraId_ )
                                                       ->computeIdealObservationsWithLinkEndData(
                                                               observationTime, observation_models::receiver, linkEndTimes, linkEndStates );
        BOOST_REQUIRE_EQUAL( linkEndStates.size( ), 2 );
        const Eigen::Vector3d modelRelativePositionJ2000 = linkEndStates.at( 0 ).segment( 0, 3 ) - linkEndStates.at( 1 ).segment( 0, 3 );

        // Compare the full observation-model path against a direct camera call using the model's own link-end geometry.
        const Eigen::Vector2d directCameraModelPixelLine =
                bodies.at( "VGR2" )
                        ->getVehicleSystems( )
                        ->getCamera( image.cameraId_ )
                        ->calculateObservableFromInertial( computeApparentDirectionWithStellarAberration(
                                                                   modelRelativePositionJ2000.normalized( ),
                                                                   paperVoyagerState.segment( 3, 3 ) / physical_constants::SPEED_OF_LIGHT ),
                                                           observationTime );
        BOOST_CHECK_SMALL( ( modelPixelLine - directCameraModelPixelLine ).norm( ), 1.0E-9 );
        const Eigen::Vector2d modelPixelLineDifference = modelPixelLine - observationSetPixelLine;

        // Keep the printed table focused on the observable tested by this block.
        std::cout << std::setw( 12 ) << reference.pictureName_ << std::setw( 14 ) << std::fixed << std::setprecision( 3 )
                  << modelPixelLineDifference.x( ) << std::setw( 13 ) << modelPixelLineDifference.y( ) << std::setw( 13 )
                  << modelPixelLineDifference.norm( ) << "\n";
        BOOST_CHECK_SMALL( modelPixelLineDifference.norm( ), maximumModelPixelLineDifference );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
