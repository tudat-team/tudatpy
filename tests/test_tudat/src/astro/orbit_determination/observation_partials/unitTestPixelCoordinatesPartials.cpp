/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/cameraPointingCorrection.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialRotationalState.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/simulation/estimation_setup/createEstimatableParametersFactory.h"
#include "tudat/support/sumLmkTestHelpers.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/support/numericalObservationPartial.h"
#include "tudat/simulation/environment_setup/createCameras.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::gravitation;
using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;
using namespace tudat::observation_partials;
using namespace tudat::estimatable_parameters;

namespace
{

SystemOfBodies createSyntheticSumLmkBodies( )
{
    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );

    bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
    spacecraftState( 2 ) = -10000.0;
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );

    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );

    bodies.processBodyFrameDefinitions< double, double >( );
    return bodies;
}

//! Minimal in-memory SUM image with synthetic (+Z-boresight) geometry: spacecraft behind a
//! landmark at the origin, identity camera axes, so projection succeeds. Used to exercise the
//! conversion plumbing (SIGMA_PTG validation, etc.) independent of the SUM/LMK file parser.
input_output::sum_lmk::SumImageData makeSyntheticSumImage( const std::string& imageId, const Eigen::Vector3d& pointingSigma )
{
    input_output::sum_lmk::SumImageData image;
    image.imageId_ = imageId;
    image.utcEpochString_ = "2015 JUN 05 07:24:42.053";
    image.imageSize_ = Eigen::Vector2i( 1024, 1024 );
    image.focalLengthMm_ = 100.0;
    image.opticalCenter_ = Eigen::Vector2d( 512.0, 512.0 );
    // SCOBJ = spacecraft-to-object: with the spacecraft at z = -10000 and target at the origin,
    // the spacecraft-to-object vector is +10000 along z.
    image.spacecraftObjectVector_ = Eigen::Vector3d( 0.0, 0.0, 10000.0 );
    image.cameraAxes_ = Eigen::Matrix3d::Identity( );
    image.kMatrix_ << 10.0, 0.0, 0.0, 0.0, 10.0, 0.0;
    image.pointingSigma_ = pointingSigma;
    input_output::sum_lmk::SumLandmarkObservation observation;
    observation.landmarkId_ = "LMK0001";
    observation.pixelCoordinates_ = Eigen::Vector2d( 512.0, 512.0 );
    image.landmarkObservations_.push_back( observation );
    return image;
}

std::map< std::string, input_output::sum_lmk::LmkLandmarkData > makeSyntheticLandmarks( )
{
    input_output::sum_lmk::LmkLandmarkData landmark;
    landmark.landmarkId_ = "LMK0001";
    landmark.bodyFixedPosition_ = Eigen::Vector3d::Zero( );
    return { { "LMK0001", landmark } };
}

}  // namespace

BOOST_AUTO_TEST_SUITE( test_pixel_coordinates_partials )

//! Test partial derivatives of pixel coordinate observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testPixelCoordinatesPartials )
{
    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );
    parameterPerturbationMultipliers( 2 ) = 10.0;

    // Test partials with constant ephemerides (allows test of position partials).
    {
        const double stateEvaluationTime = 1.1E7;

        // Create environment.
        SystemOfBodies bodies = setupEnvironment( groundStations, 1.0E7, 1.2E7, stateEvaluationTime, true );

        // Create a camera on Earth with non-identity mounting angles, pointed along the body-fixed Earth-to-Mars line of sight.
        const std::string cameraName = "Camera";
        const Eigen::Vector3d cameraBodyFixedPosition = ( Eigen::Vector3d( ) << 1.7E6, -6.2E6, 1.3E5 ).finished( );
        const Eigen::Vector3d marsStationPositionInInertialFrame =
                bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ).segment( 0, 3 ) +
                bodies.at( "Mars" )->getRotationalEphemeris( )->getRotationToTargetFrame( stateEvaluationTime ).inverse( ) *
                        bodies.at( "Mars" )
                                ->getGroundStation( groundStations[ 1 ].second )
                                ->getNominalStationState( )
                                ->getNominalCartesianPosition( );
        const Eigen::Vector3d relativePositionInInertialFrame = marsStationPositionInInertialFrame -
                bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ).segment( 0, 3 );
        const Eigen::Vector3d relativePositionInBodyFixedFrame =
                bodies.at( "Earth" )->getRotationalEphemeris( )->getRotationToTargetFrame( stateEvaluationTime ) *
                        relativePositionInInertialFrame -
                cameraBodyFixedPosition;
        const Eigen::Vector3d cameraBoresightDirection = relativePositionInBodyFixedFrame.normalized( );
        const Eigen::Vector3d cameraBoresightEulerAngles =
                ( Eigen::Vector3d( ) << std::atan2( cameraBoresightDirection.y( ), cameraBoresightDirection.x( ) ),
                  std::asin( cameraBoresightDirection.z( ) ),
                  0.3 )
                        .finished( );

        createCamera( bodies.at( "Earth" ),
                      cameraName,
                      cameraBoresightEulerAngles,
                      std::make_pair( 1200.0, 900.0 ),
                      std::make_pair( 320.0, 240.0 ),
                      cameraBodyFixedPosition );

        // Set link ends for observation model.
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = groundStations[ 1 ];
        linkEnds[ receiver ] = std::make_pair( std::string( "Earth" ), cameraName );

        // Generate pixel coordinates model.
        std::vector< std::string > perturbingBodies;
        perturbingBodies.push_back( "Earth" );
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings = {
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies )
        };
        std::shared_ptr< ObservationModel< 2 > > pixelCoordinatesModel =
                observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                        pixelCoordinatesSettings( linkEnds, lightTimeCorrectionSettings ), bodies );

        // Create parameter objects.
        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > initialStateParameters;
        initialStateParameters.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
                "Earth",
                bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ),
                "SSB",
                "ECLIPJ2000" ) );
        initialStateParameters.push_back( std::make_shared< InitialTranslationalStateParameter< double > >(
                "Mars",
                bodies.at( "Mars" )->getStateInBaseFrameFromEphemeris< double, double >( stateEvaluationTime ),
                "SSB",
                "ECLIPJ2000" ) );

        std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                std::make_shared< EstimatableParameterSet< double > >(
                        std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ),
                        initialStateParameters );

        testObservationPartials( pixelCoordinatesModel,
                                 bodies,
                                 fullEstimatableParameterSet,
                                 linkEnds,
                                 pixel_coordinates,
                                 1.0E-5,
                                 true,
                                 false,
                                 1.0,
                                 parameterPerturbationMultipliers,
                                 nullptr,
                                 stateEvaluationTime );
    }
}

BOOST_AUTO_TEST_CASE( testSumLmkConversionAndCameraPointingPartials )
{
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/sample.sum" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks =
            input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", testDataPath + "/LMK0002.lmk" } );

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    LoadedPixelTestData conversionResult =
            loadPixelTestData( sumImages, landmarks, bodies );

    BOOST_REQUIRE( conversionResult.observationCollection_ != nullptr );
    const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( "IMG001" );
    BOOST_CHECK_EQUAL( cameraName, "Camera_IMG001" );
    BOOST_CHECK( bodies.at( "Spacecraft" )->getVehicleSystems( )->getCameraMap( ).count( cameraName ) == 1 );
    BOOST_CHECK( bodies.at( "Target" )->getGroundStationMap( ).count( "LMK0001" ) == 1 );
    BOOST_CHECK( bodies.at( "Target" )->getGroundStationMap( ).count( "LMK0002" ) == 1 );
    BOOST_REQUIRE_EQUAL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.size( ), 1 );
    BOOST_CHECK_EQUAL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).first.first, camera_pointing_correction );
    BOOST_CHECK_SMALL( conversionResult.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).second( 0 ) - 1.0E8, 1.0E-6 );

    LinkEnds linkEndsLmk1;
    linkEndsLmk1[ transmitter ] = LinkEndId( "Target", "LMK0001" );
    linkEndsLmk1[ receiver ] = LinkEndId( "Spacecraft", cameraName );
    LinkEnds linkEndsLmk2;
    linkEndsLmk2[ transmitter ] = LinkEndId( "Target", "LMK0002" );
    linkEndsLmk2[ receiver ] = LinkEndId( "Spacecraft", cameraName );

    const std::vector< std::shared_ptr< SingleObservationSet< double, double > > > lmk2ObservationSets =
            conversionResult.observationCollection_->getSingleLinkAndTypeObservationSets( pixel_coordinates,
                                                                                          LinkDefinition( linkEndsLmk2 ) );
    BOOST_REQUIRE_EQUAL( lmk2ObservationSets.size( ), 1 );
    BOOST_CHECK_SMALL( lmk2ObservationSets.at( 0 )->getObservation( 0 )( 0 ) - 612.0, 1.0E-12 );
    BOOST_CHECK_SMALL( lmk2ObservationSets.at( 0 )->getObservation( 0 )( 1 ) - 512.0, 1.0E-12 );
    BOOST_CHECK( linkEndsLmk1.at( receiver ) == linkEndsLmk2.at( receiver ) );

    std::shared_ptr< ObservationModel< 2, double, double > > pixelCoordinatesModel =
            ObservationModelCreator< 2, double, double >::createObservationModel(
                    pixelCoordinatesSettings( LinkDefinition( linkEndsLmk2 ) ), bodies );

    const double observationTime = lmk2ObservationSets.at( 0 )->getObservationTime( 0 );
    std::vector< Eigen::Vector6d > linkEndStates;
    std::vector< double > linkEndTimes;
    Eigen::Vector2d computedObservation =
            pixelCoordinatesModel->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimes, linkEndStates, nullptr );
    BOOST_CHECK_SMALL( computedObservation( 0 ) - 612.0, 1.0E-10 );
    BOOST_CHECK_SMALL( computedObservation( 1 ) - 512.0, 1.0E-10 );

    std::shared_ptr< CameraPointingCorrection > pointingParameter =
            std::make_shared< CameraPointingCorrection >( bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraName );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( { pointingParameter } ) );

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 2 > > >, std::shared_ptr< PositionPartialScaling > >
            partialsAndScaling = ObservationPartialCreator< 2, double, double >::createObservationPartials(
                    pixelCoordinatesModel, bodies, parameterSet );
    BOOST_REQUIRE_EQUAL( partialsAndScaling.first.size( ), 1 );
    std::shared_ptr< ObservationPartial< 2 > > pointingPartial = partialsAndScaling.first.begin( )->second;
    BOOST_REQUIRE( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( pointingPartial ) != nullptr );

    auto computePixelObservation = [ & ]( const Eigen::Vector3d& correction ) {
        pointingParameter->setParameterValue( correction );
        std::vector< Eigen::Vector6d > currentStates;
        std::vector< double > currentTimes;
        return pixelCoordinatesModel->computeObservationsWithLinkEndData( observationTime, receiver, currentTimes, currentStates, nullptr );
    };

    auto computeAnalyticalPointingPartial = [ & ]( const Eigen::Vector3d& correction ) {
        pointingParameter->setParameterValue( correction );
        std::vector< Eigen::Vector6d > currentStates;
        std::vector< double > currentTimes;
        const Eigen::Vector2d currentObservation = pixelCoordinatesModel->computeObservationsWithLinkEndData(
                observationTime, receiver, currentTimes, currentStates, nullptr );
        partialsAndScaling.second->update( currentStates, currentTimes, receiver, currentObservation );
        return pointingPartial->calculatePartial( currentStates, currentTimes, receiver, nullptr, currentObservation ).at( 0 ).first;
    };

    auto computeNumericalPointingPartial = [ & ]( const Eigen::Vector3d& correction ) {
        const double perturbation = 1.0E-7;
        Eigen::Matrix< double, 2, 3 > numericalPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
        for( int i = 0; i < 3; ++i )
        {
            Eigen::Vector3d perturbationVector = Eigen::Vector3d::Zero( );
            perturbationVector( i ) = perturbation;
            numericalPartial.col( i ) = ( computePixelObservation( correction + perturbationVector ) -
                                          computePixelObservation( correction - perturbationVector ) ) /
                    ( 2.0 * perturbation );
        }
        pointingParameter->setParameterValue( correction );
        return numericalPartial;
    };

    const Eigen::Vector3d zeroCorrection = Eigen::Vector3d::Zero( );
    const Eigen::Matrix< double, 2, 3 > zeroAnalyticalPartial = computeAnalyticalPointingPartial( zeroCorrection );
    const Eigen::Matrix< double, 2, 3 > zeroNumericalPartial = computeNumericalPointingPartial( zeroCorrection );
    BOOST_CHECK_SMALL( ( zeroAnalyticalPartial - zeroNumericalPartial ).norm( ), 1.0E-5 );

    const Eigen::Vector3d nonZeroCorrection = ( Eigen::Vector3d( ) << 2.0E-4, -1.0E-4, 3.0E-4 ).finished( );
    const Eigen::Matrix< double, 2, 3 > nonZeroAnalyticalPartial = computeAnalyticalPointingPartial( nonZeroCorrection );
    const Eigen::Matrix< double, 2, 3 > nonZeroNumericalPartial = computeNumericalPointingPartial( nonZeroCorrection );
    BOOST_CHECK_SMALL( ( nonZeroAnalyticalPartial - nonZeroNumericalPartial ).norm( ), 1.0E-5 );

    const Eigen::Vector2d shiftedObservation = computePixelObservation( nonZeroCorrection );
    BOOST_CHECK( ( shiftedObservation - computedObservation ).norm( ) > 1.0E-3 );
    pointingParameter->setParameterValue( Eigen::Vector3d::Zero( ) );
}

//! UTC epoch string -> TDB seconds since J2000 must match an independent SPICE str2et result.
BOOST_AUTO_TEST_CASE( testSumUtcToTdbMatchesSpice )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const input_output::sum_lmk::SumImageData image = input_output::sum_lmk::readSumFile( testDataPath + "/W48230079013.SUM" );

    const double tdbSecondsSinceJ2000 = sumImageTdbTimeForTest( image );
    // str2et interprets the calendar string as UTC and returns ephemeris time (TDB) seconds past J2000.
    const double expected = spice_interface::convertDateStringToEphemerisTime( "2015-04-14T16:26:27.974" );
    BOOST_CHECK_SMALL( tdbSecondsSinceJ2000 - expected, 1.0E-3 );
}

//! Image-ID sanitization into Camera_<id>: safe characters preserved, others replaced, empty throws.
BOOST_AUTO_TEST_CASE( testSumCameraNameSanitization )
{
    BOOST_CHECK_EQUAL( input_output::sum_lmk::getSumImageCameraName( "W48230079013" ), "Camera_W48230079013" );
    BOOST_CHECK_EQUAL( input_output::sum_lmk::getSumImageCameraName( "  N123 " ), "Camera_N123" );
    BOOST_CHECK_EQUAL( input_output::sum_lmk::getSumImageCameraName( "IMG/01:a" ), "Camera_IMG_01_a" );
    BOOST_CHECK_THROW( input_output::sum_lmk::getSumImageCameraName( "   " ), std::runtime_error );
}

//! SIGMA_PTG -> inverse a-priori covariance: positivity and all-or-nothing validation.
BOOST_AUTO_TEST_CASE( testSumPointingSigmaValidation )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = makeSyntheticLandmarks( );

    // Negative sigma component must throw.
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage(
                "IMGNEG", ( Eigen::Vector3d( ) << -1.0E-4, -1.0E-4, -1.0E-4 ).finished( ) ) };
        BOOST_CHECK_THROW( ( loadPixelTestData( images, landmarks, bodies ) ),
                           std::runtime_error );
    }
    // Partially-specified sigma (one finite, two NaN) must throw.
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage(
                "IMGPART", ( Eigen::Vector3d( ) << 1.0E-4, TUDAT_NAN, TUDAT_NAN ).finished( ) ) };
        BOOST_CHECK_THROW( ( loadPixelTestData( images, landmarks, bodies ) ),
                           std::runtime_error );
    }
    // Valid uniform sigma yields one inverse-variance entry of 1/sigma^2 per component.
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage( "IMGOK",
                                                                                                   Eigen::Vector3d::Constant( 1.0E-4 ) ) };
        LoadedPixelTestData result =
                loadPixelTestData( images, landmarks, bodies );
        BOOST_REQUIRE_EQUAL( result.inverseAprioriCovarianceDiagonalEntries_.size( ), 1 );
        BOOST_CHECK_CLOSE( result.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).second( 0 ), 1.0E8, 1.0E-6 );
        BOOST_CHECK_CLOSE( result.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).second( 2 ), 1.0E8, 1.0E-6 );
    }
}

//! SPC SUM's four DISTORTION values are not sufficient to reconstruct its separate Owen
//! distortion model. The converter must reject a non-zero row before it mutates the body system,
//! rather than silently projecting with a zero-distortion camera.
BOOST_AUTO_TEST_CASE( testSumNonZeroDistortionIsRejected )
{
    spice_interface::loadStandardSpiceKernels( );
    input_output::sum_lmk::SumImageData image = makeSyntheticSumImage( "IMGDIST", Eigen::Vector3d::Constant( 1.0E-4 ) );
    image.distortionCoefficients_( 2 ) = 1.0E-8;

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    try
    {
        loadPixelTestData( { image }, makeSyntheticLandmarks( ), bodies );
        BOOST_FAIL( "Expected non-zero SUM DISTORTION validation to throw." );
    }
    catch( const std::runtime_error& exception )
    {
        const std::string message = exception.what( );
        BOOST_CHECK( message.find( "DISTORTION" ) != std::string::npos );
        BOOST_CHECK( message.find( "INIT_LITHOS" ) != std::string::npos );
    }

    BOOST_CHECK( bodies.at( "Target" )->getGroundStationMap( ).empty( ) );
    BOOST_CHECK( bodies.at( "Spacecraft" )->getVehicleSystems( )->getCameraMap( ).empty( ) );
}

//! Programmatic SUM inputs bypass the file reader, so calibration and attitude must be revalidated
//! before the converter adds landmarks or cameras to the caller's body system.
BOOST_AUTO_TEST_CASE( testSumCalibrationAndAttitudeAreValidatedBeforeEnvironmentMutation )
{
    spice_interface::loadStandardSpiceKernels( );
    const auto checkRejectedImage = [ ]( const input_output::sum_lmk::SumImageData& image,
                                                             const std::string& expectedDiagnostic ) {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        bool threw = false;
        try
        {
            loadPixelTestData( { image }, makeSyntheticLandmarks( ), bodies );
        }
        catch( const std::runtime_error& exception )
        {
            threw = true;
            BOOST_CHECK( std::string( exception.what( ) ).find( expectedDiagnostic ) != std::string::npos );
        }
        BOOST_CHECK( threw );
        BOOST_CHECK( bodies.at( "Target" )->getGroundStationMap( ).empty( ) );
        BOOST_CHECK( bodies.at( "Spacecraft" )->getVehicleSystems( )->getCameraMap( ).empty( ) );
    };

    input_output::sum_lmk::SumImageData invalidFocalLength = makeSyntheticSumImage( "IMGBADFOCAL", Eigen::Vector3d::Constant( 1.0E-4 ) );
    invalidFocalLength.focalLengthMm_ = 0.0;
    checkRejectedImage( invalidFocalLength, "MMFL" );

    input_output::sum_lmk::SumImageData singularKMatrix = makeSyntheticSumImage( "IMGBADKMAT", Eigen::Vector3d::Constant( 1.0E-4 ) );
    singularKMatrix.kMatrix_.row( 1 ) = 2.0 * singularKMatrix.kMatrix_.row( 0 );
    checkRejectedImage( singularKMatrix, "K-MATRIX" );

    input_output::sum_lmk::SumImageData nonOrthonormalAttitude = makeSyntheticSumImage( "IMGBADATT", Eigen::Vector3d::Constant( 1.0E-4 ) );
    nonOrthonormalAttitude.cameraAxes_( 0, 0 ) = 2.0;
    checkRejectedImage( nonOrthonormalAttitude, "CX/CY/CZ" );
}

//! Missing LMK definitions must be reported in one aggregated diagnostic.
BOOST_AUTO_TEST_CASE( testSumLmkMissingLandmarkValidation )
{
    spice_interface::loadStandardSpiceKernels( );
    const Eigen::Vector3d pointingSigma = Eigen::Vector3d::Constant( 1.0E-4 );
    input_output::sum_lmk::SumImageData firstImage = makeSyntheticSumImage( "IMGMISS1", pointingSigma );
    input_output::sum_lmk::SumImageData secondImage = makeSyntheticSumImage( "IMGMISS2", pointingSigma );
    firstImage.landmarkObservations_.at( 0 ).landmarkId_ = "LMK_MISSING_A";
    secondImage.landmarkObservations_.at( 0 ).landmarkId_ = "LMK_MISSING_B";

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    try
    {
        loadPixelTestData( { firstImage, secondImage },
                                                             std::map< std::string, input_output::sum_lmk::LmkLandmarkData >( ),
                                                             bodies );
        BOOST_FAIL( "Expected missing LMK validation to throw." );
    }
    catch( const std::runtime_error& exception )
    {
        const std::string message = exception.what( );
        BOOST_CHECK( message.find( "LMK_MISSING_A" ) != std::string::npos );
        BOOST_CHECK( message.find( "LMK_MISSING_B" ) != std::string::npos );
        BOOST_CHECK( message.find( "IMGMISS1" ) != std::string::npos );
        BOOST_CHECK( message.find( "IMGMISS2" ) != std::string::npos );
    }
}

//! The optional missing-LMK skip mode drops only unavailable landmark observations and keeps valid ones.
BOOST_AUTO_TEST_CASE( testSumLmkMissingLandmarkSkipMode )
{
    spice_interface::loadStandardSpiceKernels( );
    input_output::sum_lmk::SumImageData image = makeSyntheticSumImage( "IMGMIXED", Eigen::Vector3d::Constant( 1.0E-4 ) );
    input_output::sum_lmk::SumLandmarkObservation missingObservation;
    missingObservation.landmarkId_ = "LMK_MISSING";
    missingObservation.pixelCoordinates_ = Eigen::Vector2d( 600.0, 600.0 );
    image.landmarkObservations_.push_back( missingObservation );

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    LoadedPixelTestData result =
            loadPixelTestData( { image }, makeSyntheticLandmarks( ), bodies, "Target", "Spacecraft", true );

    BOOST_REQUIRE_EQUAL( result.observationCollection_->getConcatenatedObservations( ).size( ), 2 );
    BOOST_REQUIRE_EQUAL( result.observationModelSettings_.size( ), 1 );
    BOOST_CHECK_EQUAL( result.observationModelSettings_.at( 0 )->linkEnds_.linkEnds_.at( transmitter ).stationName_, "LMK0001" );
}

//! The conversion helper must require a rotational ephemeris on the receiver body and must not
//! attach one implicitly.
BOOST_AUTO_TEST_CASE( testSumReceiverEphemerisGuard )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = makeSyntheticLandmarks( );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage( "IMGEPH",
                                                                                               Eigen::Vector3d::Constant( 1.0E-4 ) ) };

    auto makeBodies = []( const bool attachReceiverRotation ) {
        SystemOfBodies bodies( "SSB", "J2000" );
        bodies.createEmptyBody< double, double >( "Target", false );
        bodies.createEmptyBody< double, double >( "Spacecraft", false );
        bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
        Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
        spacecraftState( 2 ) = -10000.0;
        bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );
        bodies.at( "Target" )
                ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                        Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Target_Fixed" ) );
        if( attachReceiverRotation )
        {
            bodies.at( "Spacecraft" )
                    ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                            Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
        }
        return bodies;
    };

    // Loading metadata does not invent a receiver rotation. The observation model requires one.
    {
        SystemOfBodies bodies = makeBodies( false );
        const auto loaded = loadPixelTestData( images, landmarks, bodies );
        BOOST_CHECK( bodies.at( "Spacecraft" )->getRotationalEphemeris( ) == nullptr );
        BOOST_CHECK_THROW( ( createObservationSimulators< double, double >( loaded.observationModelSettings_, bodies ) ), std::runtime_error );
    }
    // Caller-supplied placeholder rotation -> succeeds.
    {
        SystemOfBodies bodies = makeBodies( true );
        BOOST_CHECK_NO_THROW( ( loadPixelTestData( images, landmarks, bodies ) ) );
    }
}

//! Generic time-bounds filtering works on pixel observation collections and model settings can be
//! regenerated from the filtered collection, matching the Doppler workflow.
BOOST_AUTO_TEST_CASE( testSumLmkTimeBoundsFilteringWorkflow )
{
    spice_interface::loadStandardSpiceKernels( );
    input_output::sum_lmk::SumImageData inWindowImage = makeSyntheticSumImage( "IMGIN", Eigen::Vector3d::Constant( 1.0E-4 ) );
    input_output::sum_lmk::SumImageData outOfWindowImage = makeSyntheticSumImage( "IMGOUT", Eigen::Vector3d::Constant( 1.0E-4 ) );
    outOfWindowImage.utcEpochString_ = "2015 JUN 05 08:24:42.053";

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    LoadedPixelTestData conversionResult = loadPixelTestData(
            { inWindowImage, outOfWindowImage }, makeSyntheticLandmarks( ), bodies );

    BOOST_REQUIRE_EQUAL( conversionResult.observationCollection_->getConcatenatedTimeVector( ).size( ), 4 );
    BOOST_REQUIRE_EQUAL( conversionResult.observationModelSettings_.size( ), 2 );

    const double inWindowTime = sumImageTdbTimeForTest( inWindowImage );
    const std::shared_ptr< ObservationFilterBase > timeBoundsFilter =
            observationFilter( time_bounds_filtering, inWindowTime - 1.0, inWindowTime + 1.0, true, true );
    conversionResult.observationCollection_->filterObservations( timeBoundsFilter );
    conversionResult.observationCollection_->removeEmptySingleObservationSets( );

    const std::vector< double > filteredTimes = conversionResult.observationCollection_->getConcatenatedTimeVector( );
    BOOST_REQUIRE_EQUAL( filteredTimes.size( ), 2 );
    for( const double filteredTime : filteredTimes )
    {
        BOOST_CHECK_SMALL( filteredTime - inWindowTime, 1.0E-9 );
    }

    const std::vector< std::shared_ptr< ObservationModelSettings > > filteredModelSettings =
            pixelModelSettingsForTest( conversionResult.observationCollection_ );
    BOOST_REQUIRE_EQUAL( filteredModelSettings.size( ), 1 );
    BOOST_CHECK_EQUAL( filteredModelSettings.at( 0 )->linkEnds_.linkEnds_.at( receiver ).stationName_, "Camera_IMGIN" );
}

//! End-to-end regression on real SPC SUM/LMK data: the observation model must reproduce the
//! measured landmark pixels to within the SPC solution residual.
//!
//! SCOBJ is the spacecraft-to-object vector in the target body-fixed frame, so the spacecraft
//! body-fixed position is -SCOBJ and the landmark position relative to the spacecraft is
//! VLM + SCOBJ. With the target at the origin under identity rotation, the camera axes CX/CY/CZ
//! (rows of R_targetBodyFixed->camera) place the landmark on the +Z boresight and the projection
//! through MMFL/CTR/K-MATRIX matches the measured pixels (~1 px residual). A sign/convention
//! regression (e.g. VLM - SCOBJ) would put the landmark at negative Z and fail by tens of pixels.
BOOST_AUTO_TEST_CASE( testRealSumLmkReprojection )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/W48230079013.SUM" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = input_output::sum_lmk::readLmkFiles(
            { testDataPath + "/FA0115.LMK", testDataPath + "/FC0098.LMK", testDataPath + "/GD0036.LMK" } );
    BOOST_REQUIRE_EQUAL( sumImages.size( ), 1 );
    const input_output::sum_lmk::SumImageData& image = sumImages.at( 0 );

    // Body-fixed = inertial here: Target at the origin under identity rotation, spacecraft at -SCOBJ.
    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );
    bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
    spacecraftState.segment( 0, 3 ) = -image.spacecraftObjectVector_;
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    // Real geometry now projects to positive camera-frame Z.
    for( const auto& entry : landmarks )
    {
        const Eigen::Vector3d relativeCameraFrame =
                image.cameraAxes_ * ( entry.second.bodyFixedPosition_ - spacecraftState.segment( 0, 3 ) );
        BOOST_CHECK_MESSAGE( relativeCameraFrame.z( ) > 0.0, "Expected positive camera-frame Z for real SPC landmark " + entry.first );
    }

    // Conversion succeeds on real geometry.
    LoadedPixelTestData conversionResult =
            loadPixelTestData( sumImages, landmarks, bodies );
    const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( "W48230079013" );

    // Each modelled pixel matches the measured SUM pixel within the SPC solution residual.
    for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Target", observation.landmarkId_ );
        linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
        std::shared_ptr< ObservationModel< 2, double, double > > model =
                ObservationModelCreator< 2, double, double >::createObservationModel(
                        pixelCoordinatesSettings( LinkDefinition( linkEnds ) ), bodies );
        std::vector< Eigen::Vector6d > states;
        std::vector< double > times;
        const Eigen::Vector2d modelled = model->computeObservationsWithLinkEndData( 0.0, receiver, times, states, nullptr );
        BOOST_CHECK_MESSAGE( ( modelled - observation.pixelCoordinates_ ).norm( ) < 3.0,
                             "Reprojection residual too large for landmark " + observation.landmarkId_ );
    }
}

//! The camera_pointing_correction partial must be created only for the matching receiver camera.
BOOST_AUTO_TEST_CASE( testCameraPointingPartialRoutingGuard )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/sample.sum" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks =
            input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", testDataPath + "/LMK0002.lmk" } );

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    LoadedPixelTestData conversionResult =
            loadPixelTestData( sumImages, landmarks, bodies );
    const std::string cameraName = input_output::sum_lmk::getSumImageCameraName( "IMG001" );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Target", "LMK0001" );
    linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
    std::shared_ptr< ObservationModel< 2, double, double > > pixelCoordinatesModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( pixelCoordinatesSettings( LinkDefinition( linkEnds ) ),
                                                                                  bodies );

    // Add a second, unrelated camera so a non-matching pointing parameter can be requested.
    createCamera( bodies.at( "Spacecraft" ),
                  "Camera_OTHER",
                  Eigen::Vector3d::Zero( ),
                  std::make_pair( 1200.0, 900.0 ),
                  std::make_pair( 320.0, 240.0 ),
                  Eigen::Vector3d::Zero( ) );

    // Matching parameter -> exactly one pointing partial.
    {
        std::shared_ptr< CameraPointingCorrection > parameter =
                std::make_shared< CameraPointingCorrection >( bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraName );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet = std::make_shared< EstimatableParameterSet< double > >(
                std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( { parameter } ) );
        auto partials =
                ObservationPartialCreator< 2, double, double >::createObservationPartials( pixelCoordinatesModel, bodies, parameterSet );
        BOOST_CHECK_EQUAL( partials.first.size( ), 1 );
    }
    // Non-matching camera -> no pointing partial.
    {
        std::shared_ptr< CameraPointingCorrection > parameter = std::make_shared< CameraPointingCorrection >(
                bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", "Camera_OTHER" );
        std::shared_ptr< EstimatableParameterSet< double > > parameterSet = std::make_shared< EstimatableParameterSet< double > >(
                std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( { parameter } ) );
        auto partials =
                ObservationPartialCreator< 2, double, double >::createObservationPartials( pixelCoordinatesModel, bodies, parameterSet );
        BOOST_CHECK_EQUAL( partials.first.size( ), 0 );
    }
}

//! Regression test for the pointing correction being applied by the Camera itself rather than by the
//! picture-specific rotation function that the SUM/LMK conversion installs. A camera created through the
//! ordinary environment-setup path has no such function and reaches the camera frame via the host body's
//! rotational ephemeris instead; the estimated correction must affect that path identically, otherwise the
//! pointing partial would be non-zero while the observable is insensitive to the parameter.
BOOST_AUTO_TEST_CASE( testPointingCorrectionAppliedForBodyFixedCamera )
{
    spice_interface::loadStandardSpiceKernels( );

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );

    // Landmark off the camera boresight, so that all three rotation degrees of freedom (including the one
    // about the boresight) move the pixel and the partial has no structurally vanishing column.
    createGroundStation( bodies.at( "Target" ),
                         "LMK0001",
                         ( Eigen::Vector3d( ) << 400.0, -250.0, 100.0 ).finished( ),
                         coordinate_conversions::cartesian_position,
                         std::vector< std::shared_ptr< GroundStationMotionSettings > >( ) );

    // Ordinary camera creation: boresight Euler angles chosen such that the body-fixed-to-camera rotation is
    // the identity, so the camera looks along +Z body-fixed, towards the target.
    const std::string cameraName = "Camera_BodyFixed";
    createCamera( bodies.at( "Spacecraft" ),
                  cameraName,
                  ( Eigen::Vector3d( ) << -mathematical_constants::PI / 2.0, mathematical_constants::PI / 2.0, 0.0 ).finished( ),
                  std::make_pair( 1200.0, 900.0 ),
                  std::make_pair( 320.0, 240.0 ),
                  Eigen::Vector3d::Zero( ) );

    std::shared_ptr< system_models::Camera > camera = bodies.at( "Spacecraft" )->getVehicleSystems( )->getCamera( cameraName );
    BOOST_REQUIRE( !camera->hasRotationFromInertialToCameraFrameFunction( ) );
    BOOST_CHECK_SMALL(
            ( camera->getNominalRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( ) - Eigen::Matrix3d::Identity( ) ).norm( ),
            1.0E-12 );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Target", "LMK0001" );
    linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
    std::shared_ptr< ObservationModel< 2, double, double > > pixelCoordinatesModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( pixelCoordinatesSettings( LinkDefinition( linkEnds ) ),
                                                                                  bodies );

    std::shared_ptr< CameraPointingCorrection > pointingParameter =
            std::make_shared< CameraPointingCorrection >( bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraName );
    std::shared_ptr< EstimatableParameterSet< double > > parameterSet = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( { pointingParameter } ) );

    auto partialsAndScaling =
            ObservationPartialCreator< 2, double, double >::createObservationPartials( pixelCoordinatesModel, bodies, parameterSet );
    BOOST_REQUIRE_EQUAL( partialsAndScaling.first.size( ), 1 );
    std::shared_ptr< ObservationPartial< 2 > > pointingPartial = partialsAndScaling.first.begin( )->second;

    const double observationTime = 0.0;
    auto computePixelObservation = [ & ]( const Eigen::Vector3d& correction ) {
        pointingParameter->setParameterValue( correction );
        std::vector< Eigen::Vector6d > currentStates;
        std::vector< double > currentTimes;
        return pixelCoordinatesModel->computeObservationsWithLinkEndData( observationTime, receiver, currentTimes, currentStates, nullptr );
    };

    auto computeAnalyticalPointingPartial = [ & ]( const Eigen::Vector3d& correction ) {
        pointingParameter->setParameterValue( correction );
        std::vector< Eigen::Vector6d > currentStates;
        std::vector< double > currentTimes;
        const Eigen::Vector2d currentObservation = pixelCoordinatesModel->computeObservationsWithLinkEndData(
                observationTime, receiver, currentTimes, currentStates, nullptr );
        partialsAndScaling.second->update( currentStates, currentTimes, receiver, currentObservation );
        return pointingPartial->calculatePartial( currentStates, currentTimes, receiver, nullptr, currentObservation ).at( 0 ).first;
    };

    auto computeNumericalPointingPartial = [ & ]( const Eigen::Vector3d& correction ) {
        const double perturbation = 1.0E-7;
        Eigen::Matrix< double, 2, 3 > numericalPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
        for( int i = 0; i < 3; ++i )
        {
            Eigen::Vector3d perturbationVector = Eigen::Vector3d::Zero( );
            perturbationVector( i ) = perturbation;
            numericalPartial.col( i ) = ( computePixelObservation( correction + perturbationVector ) -
                                          computePixelObservation( correction - perturbationVector ) ) /
                    ( 2.0 * perturbation );
        }
        pointingParameter->setParameterValue( correction );
        return numericalPartial;
    };

    const Eigen::Vector2d nominalObservation = computePixelObservation( Eigen::Vector3d::Zero( ) );

    // The observable must actually respond to the correction on this code path.
    const Eigen::Vector3d correction = ( Eigen::Vector3d( ) << 2.0E-4, -1.0E-4, 3.0E-4 ).finished( );
    BOOST_CHECK( ( computePixelObservation( correction ) - nominalObservation ).norm( ) > 1.0E-2 );

    // ... and the analytical partial must agree with the observable's actual sensitivity, both at zero and at
    // a non-zero correction (where the left Jacobian of SO(3) differs from the identity).
    for( const Eigen::Vector3d& evaluationPoint : { Eigen::Vector3d( Eigen::Vector3d::Zero( ) ), correction } )
    {
        const Eigen::Matrix< double, 2, 3 > analyticalPartial = computeAnalyticalPointingPartial( evaluationPoint );
        const Eigen::Matrix< double, 2, 3 > numericalPartial = computeNumericalPointingPartial( evaluationPoint );
        BOOST_CHECK_SMALL( ( analyticalPartial - numericalPartial ).norm( ), 1.0E-4 );
        BOOST_CHECK( analyticalPartial.norm( ) > 1.0E-2 );
    }

    pointingParameter->setParameterValue( Eigen::Vector3d::Zero( ) );

    // Requesting the parameter on a body without vehicle systems must fail cleanly rather than segfault.
    bodies.createEmptyBody< double, double >( "NoSystemsBody", false );
    const std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterSettingsWithoutSystems = { cameraPointingCorrection(
            "NoSystemsBody", cameraName ) };
    BOOST_CHECK_THROW( createParametersToEstimate( parameterSettingsWithoutSystems, bodies ), std::runtime_error );
}

//! The conversion result carries one pixel-coordinate observation model setting per (image, landmark).
BOOST_AUTO_TEST_CASE( testSumLmkObservationModelSettings )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/sample.sum" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks =
            input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", testDataPath + "/LMK0002.lmk" } );

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    LoadedPixelTestData conversionResult =
            loadPixelTestData( sumImages, landmarks, bodies );

    // Two (image, landmark) link ends -> two settings, all pixel_coordinates.
    BOOST_REQUIRE_EQUAL( conversionResult.observationModelSettings_.size( ), 2 );
    for( const auto& settings : conversionResult.observationModelSettings_ )
    {
        BOOST_CHECK_EQUAL( settings->observableType_, pixel_coordinates );
    }
    // The standalone helper reproduces the same count from the collection alone.
    const std::vector< std::shared_ptr< ObservationModelSettings > > standaloneSettings =
            pixelModelSettingsForTest( conversionResult.observationCollection_ );
    BOOST_CHECK_EQUAL( standaloneSettings.size( ), 2 );
}

//! Zero-residual round trip: the synthetic fixture's stored pixels equal the modelled pixels, so the
//! O-C residuals are ~0. Exercises the full simulator/residual plumbing.
BOOST_AUTO_TEST_CASE( testSumLmkResidualsZeroForSelfConsistentData )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/sample.sum" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks =
            input_output::sum_lmk::readLmkFiles( { testDataPath + "/LMK0001.lmk", testDataPath + "/LMK0002.lmk" } );

    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    LoadedPixelTestData conversionResult =
            loadPixelTestData( sumImages, landmarks, bodies );

    const auto simulators = createObservationSimulators< double, double >( conversionResult.observationModelSettings_, bodies );
    computeResidualsAndDependentVariables( conversionResult.observationCollection_, simulators, bodies );
    const Eigen::VectorXd residuals = conversionResult.observationCollection_->getConcatenatedResiduals( );
    BOOST_REQUIRE_EQUAL( residuals.size( ), 4 );  // 2 landmarks x 2 components
    BOOST_CHECK_SMALL( residuals.norm( ), 1.0E-6 );

    const auto camera = bodies.at( "Spacecraft" )->getVehicleSystems( )->getCamera( "Camera_IMG001" );
    const Eigen::Matrix3d nominal = camera->getNominalRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( );
    camera->setPointingCorrection( Eigen::Vector3d( 1.0E-3, 0.0, 0.0 ) );
    BOOST_CHECK_SMALL( ( camera->getPointingCorrection( ) - Eigen::Vector3d( 1.0E-3, 0.0, 0.0 ) ).norm( ), 1.0E-15 );
    BOOST_CHECK_SMALL( ( camera->getNominalRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( ) - nominal ).norm( ), 1.0E-15 );
    BOOST_CHECK_GT( ( camera->getRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( ) - nominal ).norm( ), 1.0E-4 );
    computeResidualsAndDependentVariables( conversionResult.observationCollection_, simulators, bodies );
    BOOST_CHECK_GT( ( conversionResult.observationCollection_->getConcatenatedResiduals( ) - residuals ).norm( ), 1.0E-3 );
    camera->setPointingCorrection( Eigen::Vector3d::Zero( ) );
    computeResidualsAndDependentVariables( conversionResult.observationCollection_, simulators, bodies );
    BOOST_CHECK_SMALL( ( conversionResult.observationCollection_->getConcatenatedResiduals( ) - residuals ).norm( ), 1.0E-10 );
}

//! The same direct partial handles pointing with and without stellar aberration, and receiver velocity
//! affects the pixels only when aberration is enabled.
BOOST_AUTO_TEST_CASE( testPixelPointingAndReceiverVelocityWithAberration )
{
    spice_interface::loadStandardSpiceKernels( );
    for( const bool aberration : { false, true } )
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        Eigen::Vector6d receiverState;
        receiverState << 200.0, -300.0, -10000.0, 30000.0, -20000.0, 10000.0;
        bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >(
                [ &receiverState ]( ) { return receiverState; }, "SSB", "J2000" ) );
        const auto loaded = loadPixelTestData(
                { makeSyntheticSumImage( "IMGABERRATION", Eigen::Vector3d::Constant( 1.0E-4 ) ) }, makeSyntheticLandmarks( ), bodies );
        const auto link = loaded.observationCollection_->getLinkDefinitionsForSingleObservable( pixel_coordinates ).front( );
        const auto model = ObservationModelCreator< 2, double, double >::createObservationModel(
                pixelCoordinatesSettings( link, {}, nullptr, std::make_shared< LightTimeConvergenceCriteria >( ), aberration ), bodies );
        const auto parameter = std::make_shared< CameraPointingCorrection >(
                bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", "Camera_IMGABERRATION" );
        const auto parameters = std::make_shared< EstimatableParameterSet< double > >(
                std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >{ parameter } );
        const auto partials = ObservationPartialCreator< 2, double, double >::createObservationPartials( model, bodies, parameters );
        BOOST_REQUIRE_EQUAL( partials.first.size( ), 1 );
        BOOST_REQUIRE( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( partials.first.begin( )->second ) != nullptr );
        const double epoch = loaded.observationCollection_->getConcatenatedTimeVector( ).front( );
        const std::function< Eigen::VectorXd( const double ) > observable = [ & ]( const double time ) {
            // The receiver ephemeris changes between numerical evaluations at the same epoch.
            bodies.at( "Spacecraft" )->recomputeStateOnNextCall( );
            return Eigen::VectorXd( model->computeObservations( time, receiver ) );
        };
        for( const Eigen::Vector3d& correction : std::vector< Eigen::Vector3d >{
                     Eigen::Vector3d::Zero( ), Eigen::Vector3d( 1.0E-10, -2.0E-10, 3.0E-10 ),
                     Eigen::Vector3d( 1.0E-6, -2.0E-6, 3.0E-6 ),
                     Eigen::Vector3d( 0.02, -0.01, 0.03 ) } )
        {
            parameter->setParameterValue( correction );
            std::vector< double > times;
            std::vector< Eigen::Vector6d > states;
            const Eigen::Vector2d observation = model->computeObservationsWithLinkEndData( epoch, receiver, times, states );
            partials.second->update( states, times, receiver, observation );
            const auto analytical = partials.first.begin( )->second->calculatePartial( states, times, receiver, nullptr, observation );
            BOOST_REQUIRE_EQUAL( analytical.size( ), 1 );
            BOOST_CHECK_EQUAL( analytical.front( ).second, epoch );
            const Eigen::MatrixXd numerical = calculateNumericalObservationParameterPartial(
                    parameter, Eigen::Vector3d::Constant( 1.0E-7 ), observable, epoch );
            BOOST_CHECK_SMALL( ( analytical.front( ).first - numerical ).norm( ), 1.0E-5 );
        }
        parameter->setParameterValue( Eigen::Vector3d::Zero( ) );
        std::vector< double > times;
        std::vector< Eigen::Vector6d > states;
        const Eigen::Vector2d observation = model->computeObservationsWithLinkEndData( epoch, receiver, times, states );
        partials.second->update( states, times, receiver, observation );
        const auto scaling = std::dynamic_pointer_cast< PixelCoordinatesScaling >( partials.second );
        for( int component = 0; component < 6; ++component )
        {
            const double original = receiverState( component );
            const double step = component < 3 ? 0.01 : 10.0;
            receiverState( component ) = original + step;
            const Eigen::Vector2d up = observable( epoch );
            receiverState( component ) = original - step;
            const Eigen::Vector2d down = observable( epoch );
            receiverState( component ) = original;
            const Eigen::Vector2d numerical = ( up - down ) / ( 2.0 * step );
            const Eigen::Vector2d analytical = component < 3 ?
                    Eigen::Vector2d( scaling->getPositionScalingFactor( receiver ).col( component ) ) :
                    Eigen::Vector2d( scaling->getVelocityScalingFactor( receiver ).col( component - 3 ) );
            BOOST_CHECK_SMALL( ( analytical - numerical ).norm( ), component < 3 ? 1.0E-8 : 1.0E-11 );
        }
    }
}

//! A target rotation changes both the landmark position and the camera's reference frame.
BOOST_AUTO_TEST_CASE( testPixelReferenceBodyRotationPartial )
{
    spice_interface::loadStandardSpiceKernels( );
    for( const bool aberration : { false, true } )
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const auto image = makeSyntheticSumImage( "IMGROTATION", Eigen::Vector3d::Constant( 1.0E-4 ) );
        const double epoch = sumImageTdbTimeForTest( image );
        const auto rotation = std::make_shared< SimpleRotationalEphemeris >(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitX( ) ) ), 1.0E-4, epoch - 1000.0,
                "J2000", "Target_Fixed" );
        bodies.at( "Target" )->setRotationalEphemeris( rotation );
        Eigen::Vector6d receiverState;
        receiverState << 1500.0, -700.0, -10000.0, 30000.0, -20000.0, 10000.0;
        bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( receiverState, "SSB", "J2000" ) );
        auto landmarks = makeSyntheticLandmarks( );
        landmarks.at( "LMK0001" ).bodyFixedPosition_ = Eigen::Vector3d( 400.0, -250.0, 100.0 );
        const auto loaded = loadPixelTestData( { image }, landmarks, bodies );
        const auto model = ObservationModelCreator< 2, double, double >::createObservationModel(
                pixelCoordinatesSettings( loaded.observationCollection_->getLinkDefinitionsForSingleObservable( pixel_coordinates ).front( ),
                                          {}, nullptr, std::make_shared< LightTimeConvergenceCriteria >( ), aberration ), bodies );
        const auto parameter = std::make_shared< RotationRate >( rotation, "Target" );
        const auto parameters = std::make_shared< EstimatableParameterSet< double > >(
                std::vector< std::shared_ptr< EstimatableParameter< double > > >{ parameter },
                std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ) );
        const auto partials = ObservationPartialCreator< 2, double, double >::createObservationPartials( model, bodies, parameters );
        BOOST_REQUIRE_EQUAL( partials.first.size( ), 1 );
        std::vector< double > times;
        std::vector< Eigen::Vector6d > states;
        const Eigen::Vector2d observation = model->computeObservationsWithLinkEndData( epoch, receiver, times, states );
        partials.second->update( states, times, receiver, observation );
        const auto contributions = partials.first.begin( )->second->calculatePartial( states, times, receiver, nullptr, observation );
        BOOST_REQUIRE_EQUAL( contributions.size( ), 2 );
        BOOST_CHECK_EQUAL( contributions.back( ).second, epoch );
        const Eigen::Vector2d analytical = contributions.at( 0 ).first + contributions.at( 1 ).first;
        const std::function< Eigen::VectorXd( const double ) > observable = [ & ]( const double time ) {
            return Eigen::VectorXd( model->computeObservations( time, receiver ) );
        };
        const Eigen::VectorXd numerical = calculateNumericalObservationParameterPartial( parameter, 1.0E-8, observable, epoch );
        BOOST_CHECK_SMALL( ( analytical - numerical ).norm( ) / numerical.norm( ), 1.0E-7 );
    }
}

BOOST_AUTO_TEST_CASE( testPixelReferenceBodyRotationalStatePartial )
{
    spice_interface::loadStandardSpiceKernels( );
    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
    rotationalState.head< 4 >( ) = linear_algebra::convertQuaternionToVectorFormat(
            Eigen::Quaterniond( Eigen::AngleAxisd( 0.15, Eigen::Vector3d( 1.0, 2.0, -1.0 ).normalized( ) ) ) );
    const auto rotation = std::make_shared< ConstantRotationalEphemeris >( rotationalState, "J2000", "Target_Fixed" );
    bodies.at( "Target" )->setRotationalEphemeris( rotation );
    auto landmarks = makeSyntheticLandmarks( );
    landmarks.at( "LMK0001" ).bodyFixedPosition_ = Eigen::Vector3d( 400.0, -250.0, 100.0 );
    const auto loaded = loadPixelTestData(
            { makeSyntheticSumImage( "IMGSTATE", Eigen::Vector3d::Constant( 1.0E-4 ) ) }, landmarks, bodies );
    const auto model = ObservationModelCreator< 2, double, double >::createObservationModel( loaded.observationModelSettings_.front( ), bodies );
    const auto parameter = std::make_shared< InitialRotationalStateParameter< double > >(
            "Target", rotationalState, []( ) { return Eigen::Matrix3d::Identity( ); }, "SSB", "J2000" );
    const auto parameters = std::make_shared< EstimatableParameterSet< double > >(
            std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ),
            std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >{ parameter } );
    const auto partials = ObservationPartialCreator< 2, double, double >::createObservationPartials( model, bodies, parameters );
    BOOST_REQUIRE_EQUAL( partials.first.size( ), 1 );
    const double epoch = loaded.observationCollection_->getConcatenatedTimeVector( ).front( );
    std::vector< double > times;
    std::vector< Eigen::Vector6d > states;
    const Eigen::Vector2d observation = model->computeObservationsWithLinkEndData( epoch, receiver, times, states );
    partials.second->update( states, times, receiver, observation );
    const auto contributions = partials.first.begin( )->second->calculatePartial( states, times, receiver, nullptr, observation );
    BOOST_REQUIRE_EQUAL( contributions.size( ), 2 );
    for( const auto& contribution : contributions )
    {
        BOOST_REQUIRE_EQUAL( contribution.first.cols( ), 7 );
    }
    const Eigen::Matrix< double, 2, 7 > analytical = contributions.at( 0 ).first + contributions.at( 1 ).first;
    BOOST_CHECK_SMALL( analytical.rightCols< 3 >( ).norm( ), 1.0E-15 );
    const Eigen::Matrix4d normalizationPartial =
            Eigen::Matrix4d::Identity( ) - rotationalState.head< 4 >( ) * rotationalState.head< 4 >( ).transpose( );
    const Eigen::Matrix< double, 2, 4 > normalizedAnalytical = analytical.leftCols< 4 >( ) * normalizationPartial;
    const double step = 1.0E-7;
    for( int component = 0; component < 4; ++component )
    {
        Eigen::Vector7d perturbed = rotationalState;
        perturbed( component ) += step;
        rotation->updateConstantState( perturbed );
        const Eigen::Vector2d up = model->computeObservations( epoch, receiver );
        perturbed( component ) -= 2.0 * step;
        rotation->updateConstantState( perturbed );
        const Eigen::Vector2d down = model->computeObservations( epoch, receiver );
        BOOST_CHECK_SMALL( ( normalizedAnalytical.col( component ) - ( up - down ) / ( 2.0 * step ) ).norm( ), 1.0E-5 );
    }
    rotation->updateConstantState( rotationalState );
}

//! Real-data O-C residuals against a reference trajectory (spacecraft at -SCOBJ): the SPICE-trajectory
//! analogue. Residuals must be at the SPC solution level (~1 px), not the tens-of-px a convention
//! regression would produce.
BOOST_AUTO_TEST_CASE( testSumLmkResidualsRealData )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/W48230079013.SUM" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = input_output::sum_lmk::readLmkFiles(
            { testDataPath + "/FA0115.LMK", testDataPath + "/FC0098.LMK", testDataPath + "/GD0036.LMK" } );
    const input_output::sum_lmk::SumImageData& image = sumImages.at( 0 );

    SystemOfBodies bodies( "SSB", "J2000" );
    bodies.createEmptyBody< double, double >( "Target", false );
    bodies.createEmptyBody< double, double >( "Spacecraft", false );
    bodies.at( "Target" )->setEphemeris( std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ), "SSB", "J2000" ) );
    Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
    spacecraftState.segment( 0, 3 ) = -image.spacecraftObjectVector_;
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Target_Fixed" ) );
    bodies.at( "Spacecraft" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( Eigen::Matrix3d::Identity( ) ), "J2000", "Spacecraft_Fixed" ) );
    bodies.processBodyFrameDefinitions< double, double >( );

    LoadedPixelTestData conversionResult =
            loadPixelTestData( sumImages, landmarks, bodies );

    const auto simulators = createObservationSimulators< double, double >( conversionResult.observationModelSettings_, bodies );
    computeResidualsAndDependentVariables( conversionResult.observationCollection_, simulators, bodies );
    const Eigen::VectorXd residuals = conversionResult.observationCollection_->getConcatenatedResiduals( );
    BOOST_REQUIRE_EQUAL( residuals.size( ), 6 );  // 3 landmarks x 2 components
    for( int i = 0; i < 3; ++i )
    {
        BOOST_CHECK_MESSAGE( residuals.segment( 2 * i, 2 ).norm( ) < 3.0,
                             "Real-data pixel residual too large for observation index " + std::to_string( i ) );
    }

    // Per-set residual statistics are populated and finite.
    Eigen::VectorXd startTimes, durations, meanValues, rmsValues;
    getResidualStatistics( conversionResult.observationCollection_, startTimes, durations, meanValues, rmsValues );
    BOOST_REQUIRE_EQUAL( rmsValues.size( ), 3 );
    for( int i = 0; i < rmsValues.size( ); ++i )
    {
        BOOST_CHECK( std::isfinite( rmsValues( i ) ) && rmsValues( i ) < 3.0 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
