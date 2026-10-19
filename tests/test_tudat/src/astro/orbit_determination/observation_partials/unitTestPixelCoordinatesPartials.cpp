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

//! Minimal in-memory SUM image with synthetic (+Z-boresight) geometry: spacecraft behind a landmark at the origin, and
//! identity camera axes.
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

//! End-to-end check of the SUM/LMK processing on real SPC data (Rosetta/OSIRIS image W48230079013), with the target rotated
//! w.r.t. the inertial frame. Each step is compared to an independent reference: the landmark positions to the LMK VLM records,
//! the observation epoch to SPICE, the camera orientation to CX/CY/CZ combined with the target rotation, and the modelled pixels
//! to the SPC projection evaluated by hand from the SUM/LMK values (and to the measured pixels, within the SPC solution residual).
BOOST_AUTO_TEST_CASE( testSumLmkProcessingRealData )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/W48230079013.SUM" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = input_output::sum_lmk::readLmkFiles(
            { testDataPath + "/FA0115.LMK", testDataPath + "/FC0098.LMK", testDataPath + "/GD0036.LMK" } );
    const input_output::sum_lmk::SumImageData& image = sumImages.at( 0 );
    const std::string cameraName = "Camera_W48230079013";

    // Target at the origin with a constant, non-trivial rotation; spacecraft at -SCOBJ in the target body-fixed frame.
    SystemOfBodies bodies = createSyntheticSumLmkBodies( );
    const Eigen::Matrix3d rotationToInertialFrame =
            Eigen::AngleAxisd( 0.7, Eigen::Vector3d( 1.0, -2.0, 0.5 ).normalized( ) ).toRotationMatrix( );
    bodies.at( "Target" )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    Eigen::Quaterniond( rotationToInertialFrame ), "J2000", "Target_Fixed" ) );
    Eigen::Vector6d spacecraftState = Eigen::Vector6d::Zero( );
    spacecraftState.segment( 0, 3 ) = -rotationToInertialFrame * image.spacecraftObjectVector_;
    bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( spacecraftState, "SSB", "J2000" ) );

    LoadedPixelTestData loaded = loadPixelTestData( sumImages, landmarks, bodies );

    // Landmarks are ground stations at the LMK positions.
    for( const auto& landmark : landmarks )
    {
        const Eigen::Vector3d stationPosition =
                bodies.at( "Target" )->getGroundStation( landmark.first )->getNominalStationState( )->getNominalCartesianPosition( );
        BOOST_CHECK_SMALL( ( stationPosition - landmark.second.bodyFixedPosition_ ).norm( ), 1.0E-9 );
    }

    // The observation epoch is the SUM UTC epoch, converted to TDB.
    const double epoch = loaded.observationCollection_->getConcatenatedTimeVector( ).front( );
    BOOST_CHECK_SMALL( epoch - spice_interface::convertDateStringToEphemerisTime( "2015-04-14T16:26:27.974" ), 1.0E-3 );

    // The camera is oriented by CX/CY/CZ w.r.t. the rotating target body-fixed frame.
    const std::shared_ptr< system_models::Camera > camera = bodies.at( "Spacecraft" )->getVehicleSystems( )->getCamera( cameraName );
    BOOST_CHECK_SMALL( ( camera->getRotationFromInertialToCameraFrame( epoch ).toRotationMatrix( ) -
                         image.cameraAxes_ * rotationToInertialFrame.transpose( ) )
                               .norm( ),
                       1.0E-6 );

    BOOST_REQUIRE_EQUAL( image.landmarkObservations_.size( ), 3 );
    for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Target", observation.landmarkId_ );
        linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
        const std::shared_ptr< ObservationModel< 2, double, double > > model =
                ObservationModelCreator< 2, double, double >::createObservationModel(
                        pixelCoordinatesSettings( LinkDefinition( linkEnds ) ), bodies );
        const Eigen::Vector2d modelledPixel = model->computeObservations( epoch, receiver );

        // SPC projection: landmark relative to the spacecraft at VLM + SCOBJ, rotated by CX/CY/CZ, projected with MMFL/CTR/K-MATRIX.
        const Eigen::Vector3d cameraFramePosition =
                image.cameraAxes_ * ( landmarks.at( observation.landmarkId_ ).bodyFixedPosition_ + image.spacecraftObjectVector_ );
        const double x = image.focalLengthMm_ * cameraFramePosition.x( ) / cameraFramePosition.z( );
        const double y = image.focalLengthMm_ * cameraFramePosition.y( ) / cameraFramePosition.z( );
        const Eigen::Vector2d expectedPixel = image.opticalCenter_ + image.kMatrix_ * Eigen::Vector3d( x, y, x * y );

        BOOST_CHECK_SMALL( ( modelledPixel - expectedPixel ).norm( ), 1.0E-3 );
        BOOST_CHECK_SMALL( ( modelledPixel - observation.pixelCoordinates_ ).norm( ), 3.0 );
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
    LoadedPixelTestData result = loadPixelTestData( { image }, makeSyntheticLandmarks( ), bodies, "Target", "Spacecraft", true );

    BOOST_REQUIRE_EQUAL( result.observationCollection_->getConcatenatedObservations( ).size( ), 2 );
    BOOST_REQUIRE_EQUAL( result.observationModelSettings_.size( ), 1 );
    BOOST_CHECK_EQUAL( result.observationModelSettings_.at( 0 )->linkEnds_.linkEnds_.at( transmitter ).stationName_, "LMK0001" );
}

//! The pointing correction must also be applied to a camera fixed to its host body (rather than one with picture-specific
//! pointing, as created from SUM data), and its partial must match a numerical partial on that path.
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
    BOOST_CHECK_SMALL( ( camera->getRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( ) - Eigen::Matrix3d::Identity( ) ).norm( ),
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
        bodies.at( "Spacecraft" )
                ->setEphemeris( std::make_shared< ConstantEphemeris >( [ &receiverState ]( ) { return receiverState; }, "SSB", "J2000" ) );
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
        for( const Eigen::Vector3d& correction : std::vector< Eigen::Vector3d >{ Eigen::Vector3d::Zero( ),
                                                                                 Eigen::Vector3d( 1.0E-10, -2.0E-10, 3.0E-10 ),
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
            const Eigen::MatrixXd numerical =
                    calculateNumericalObservationParameterPartial( parameter, Eigen::Vector3d::Constant( 1.0E-7 ), observable, epoch );
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
            const Eigen::Vector2d analytical = component < 3
                    ? Eigen::Vector2d( scaling->getPositionScalingFactor( receiver ).col( component ) )
                    : Eigen::Vector2d( scaling->getVelocityScalingFactor( receiver ).col( component - 3 ) );
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
        const auto rotation =
                std::make_shared< SimpleRotationalEphemeris >( Eigen::Quaterniond( Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitX( ) ) ),
                                                               1.0E-4,
                                                               epoch - 1000.0,
                                                               "J2000",
                                                               "Target_Fixed" );
        bodies.at( "Target" )->setRotationalEphemeris( rotation );
        Eigen::Vector6d receiverState;
        receiverState << 1500.0, -700.0, -10000.0, 30000.0, -20000.0, 10000.0;
        bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( receiverState, "SSB", "J2000" ) );
        auto landmarks = makeSyntheticLandmarks( );
        landmarks.at( "LMK0001" ).bodyFixedPosition_ = Eigen::Vector3d( 400.0, -250.0, 100.0 );
        const auto loaded = loadPixelTestData( { image }, landmarks, bodies );
        const auto model = ObservationModelCreator< 2, double, double >::createObservationModel(
                pixelCoordinatesSettings(
                        loaded.observationCollection_->getLinkDefinitionsForSingleObservable( pixel_coordinates ).front( ),
                        {},
                        nullptr,
                        std::make_shared< LightTimeConvergenceCriteria >( ),
                        aberration ),
                bodies );
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
    const auto loaded =
            loadPixelTestData( { makeSyntheticSumImage( "IMGSTATE", Eigen::Vector3d::Constant( 1.0E-4 ) ) }, landmarks, bodies );
    const auto model =
            ObservationModelCreator< 2, double, double >::createObservationModel( loaded.observationModelSettings_.front( ), bodies );
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

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
