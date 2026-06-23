/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/io/basicInputOutput.h"
#include "tudat/interface/spice/spiceInterface.h"

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/cameraPointingCorrection.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/constantRotationRate.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/simulation/estimation_setup/processSumLmkFiles.h"
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
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( sumImages, landmarks, bodies, conversionSettings );

    BOOST_REQUIRE( conversionResult.observationCollection_ != nullptr );
    BOOST_REQUIRE( conversionResult.imageIdToCameraName_.count( "IMG001" ) == 1 );
    const std::string cameraName = conversionResult.imageIdToCameraName_.at( "IMG001" );
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

    const double tdbSecondsSinceJ2000 = observation_models::detail::convertSumUtcStringToSecondsSinceJ2000< double >( image );
    // str2et interprets the calendar string as UTC and returns ephemeris time (TDB) seconds past J2000.
    const double expected = spice_interface::convertDateStringToEphemerisTime( "2015-04-14T16:26:27.974" );
    BOOST_CHECK_SMALL( tdbSecondsSinceJ2000 - expected, 1.0E-3 );
}

//! Image-ID sanitization into Camera_<id>: safe characters preserved, others replaced, empty throws.
BOOST_AUTO_TEST_CASE( testSumCameraNameSanitization )
{
    namespace detail = observation_models::detail;
    BOOST_CHECK_EQUAL( detail::sanitizeSumImageIdToCameraName( "W48230079013" ), "Camera_W48230079013" );
    BOOST_CHECK_EQUAL( detail::sanitizeSumImageIdToCameraName( "  N123 " ), "Camera_N123" );
    BOOST_CHECK_EQUAL( detail::sanitizeSumImageIdToCameraName( "IMG/01:a" ), "Camera_IMG_01_a" );
    BOOST_CHECK_THROW( detail::sanitizeSumImageIdToCameraName( "   " ), std::runtime_error );
}

//! SIGMA_PTG -> inverse a-priori covariance: positivity and all-or-nothing validation.
BOOST_AUTO_TEST_CASE( testSumPointingSigmaValidation )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = makeSyntheticLandmarks( );
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );

    // Negative sigma component must throw.
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage(
                "IMGNEG", ( Eigen::Vector3d( ) << -1.0E-4, -1.0E-4, -1.0E-4 ).finished( ) ) };
        BOOST_CHECK_THROW( ( createSumLmkObservationCollection< double, double >( images, landmarks, bodies, conversionSettings ) ),
                           std::runtime_error );
    }
    // Partially-specified sigma (one finite, two NaN) must throw.
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage(
                "IMGPART", ( Eigen::Vector3d( ) << 1.0E-4, TUDAT_NAN, TUDAT_NAN ).finished( ) ) };
        BOOST_CHECK_THROW( ( createSumLmkObservationCollection< double, double >( images, landmarks, bodies, conversionSettings ) ),
                           std::runtime_error );
    }
    // Valid uniform sigma yields one inverse-variance entry of 1/sigma^2 per component.
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage( "IMGOK",
                                                                                                   Eigen::Vector3d::Constant( 1.0E-4 ) ) };
        SumLmkObservationConversionResult< double, double > result =
                createSumLmkObservationCollection< double, double >( images, landmarks, bodies, conversionSettings );
        BOOST_REQUIRE_EQUAL( result.inverseAprioriCovarianceDiagonalEntries_.size( ), 1 );
        BOOST_CHECK_CLOSE( result.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).second( 0 ), 1.0E8, 1.0E-6 );
        BOOST_CHECK_CLOSE( result.inverseAprioriCovarianceDiagonalEntries_.at( 0 ).second( 2 ), 1.0E8, 1.0E-6 );
    }
}

//! The conversion helper must require a rotational ephemeris on the receiver body and must not
//! attach one implicitly.
BOOST_AUTO_TEST_CASE( testSumReceiverEphemerisGuard )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = makeSyntheticLandmarks( );
    const std::vector< input_output::sum_lmk::SumImageData > images = { makeSyntheticSumImage( "IMGEPH",
                                                                                               Eigen::Vector3d::Constant( 1.0E-4 ) ) };
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );

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

    // No rotational ephemeris on the receiver -> throws (no implicit attitude attached).
    {
        SystemOfBodies bodies = makeBodies( false );
        BOOST_CHECK_THROW( ( createSumLmkObservationCollection< double, double >( images, landmarks, bodies, conversionSettings ) ),
                           std::runtime_error );
    }
    // Caller-supplied placeholder rotation -> succeeds.
    {
        SystemOfBodies bodies = makeBodies( true );
        BOOST_CHECK_NO_THROW( ( createSumLmkObservationCollection< double, double >( images, landmarks, bodies, conversionSettings ) ) );
    }
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
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( sumImages, landmarks, bodies, conversionSettings );
    const std::string cameraName = conversionResult.imageIdToCameraName_.at( "W48230079013" );

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
    SumLmkObservationConversionSettings conversionSettings( "Target", "Spacecraft" );
    SumLmkObservationConversionResult< double, double > conversionResult =
            createSumLmkObservationCollection< double, double >( sumImages, landmarks, bodies, conversionSettings );
    const std::string cameraName = conversionResult.imageIdToCameraName_.at( "IMG001" );

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

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
