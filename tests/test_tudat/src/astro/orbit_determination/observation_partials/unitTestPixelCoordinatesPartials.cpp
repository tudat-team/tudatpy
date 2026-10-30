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

std::map< std::string, input_output::sum_lmk::LmkLandmarkData > makeSyntheticLandmarks( )
{
    input_output::sum_lmk::LmkLandmarkData landmark;
    landmark.landmarkId_ = "LMK0001";
    landmark.bodyFixedPosition_ = Eigen::Vector3d::Zero( );
    return { { "LMK0001", landmark } };
}

}  // namespace

//! Set the target rotation, and create the pixel model of the image and its summed analytical partial w.r.t. the given parameters.
std::pair< std::shared_ptr< ObservationModel< 2, double, double > >, Eigen::MatrixXd > createModelAndRotationPartial(
        SystemOfBodies& bodies,
        const input_output::sum_lmk::SumImageData& image,
        const std::map< std::string, input_output::sum_lmk::LmkLandmarkData >& landmarks,
        const std::shared_ptr< RotationalEphemeris >& rotation,
        const std::shared_ptr< EstimatableParameterSet< double > >& parameters,
        const bool aberration )
{
    const double epoch = sumImageTdbTimeForTest( image );
    bodies.at( "Target" )->setRotationalEphemeris( rotation );
    const auto loaded = loadPixelTestData( { image }, landmarks, bodies );
    const auto model = ObservationModelCreator< 2, double, double >::createObservationModel(
            pixelCoordinatesSettings( loaded.observationCollection_->getLinkDefinitionsForSingleObservable( pixel_coordinates ).front( ),
                                      {},
                                      nullptr,
                                      std::make_shared< LightTimeConvergenceCriteria >( ),
                                      aberration ),
            bodies );
    const auto partials = ObservationPartialCreator< 2, double, double >::createObservationPartials( model, bodies, parameters );
    BOOST_REQUIRE_EQUAL( partials.first.size( ), 1 );
    std::vector< double > times;
    std::vector< Eigen::Vector6d > states;
    const Eigen::Vector2d observation = model->computeObservationsWithLinkEndData( epoch, receiver, times, states );
    partials.second->update( states, times, receiver, observation );
    const auto contributions = partials.first.begin( )->second->calculatePartial( states, times, receiver, nullptr, observation );
    BOOST_REQUIRE_EQUAL( contributions.size( ), 2 );
    BOOST_CHECK_EQUAL( contributions.back( ).second, epoch );
    return std::make_pair( model, Eigen::MatrixXd( contributions.at( 0 ).first + contributions.at( 1 ).first ) );
}

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

    // Test partials with constant ephemerides (allows test of position partials), without and with stellar aberration (which
    // adds a dependency on the receiver velocity).
    for( const bool aberration : { false, true } )
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
                        pixelCoordinatesSettings( linkEnds,
                                                  lightTimeCorrectionSettings,
                                                  nullptr,
                                                  std::make_shared< LightTimeConvergenceCriteria >( ),
                                                  aberration ),
                        bodies );

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
//! A 2x2-binned copy of the image, taken by the same camera, checks that each image is projected with its own calibration.
BOOST_AUTO_TEST_CASE( testSumLmkProcessingRealData )
{
    spice_interface::loadStandardSpiceKernels( );
    const std::string testDataPath = paths::getTudatTestDataPath( ) + "/sum_lmk";
    const std::vector< input_output::sum_lmk::SumImageData > sumImages =
            input_output::sum_lmk::readSumFiles( { testDataPath + "/W48230079013.SUM" } );
    const std::map< std::string, input_output::sum_lmk::LmkLandmarkData > landmarks = input_output::sum_lmk::readLmkFiles(
            { testDataPath + "/FA0115.LMK", testDataPath + "/FC0098.LMK", testDataPath + "/GD0036.LMK" } );
    const input_output::sum_lmk::SumImageData& image = sumImages.at( 0 );
    const std::string cameraName = sumLmkTestCameraName;

    // Binning halves the pixel scale, and moves the optical centre to the centre of the binned pixel grid.
    input_output::sum_lmk::SumImageData binnedImage = image;
    binnedImage.imageId_ = "W48230079013_BINNED";
    binnedImage.utcEpochString_ = "2015 APR 14 16:36:27.974";
    binnedImage.kMatrix_ = 0.5 * image.kMatrix_;
    binnedImage.opticalCenter_ = ( image.opticalCenter_ - Eigen::Vector2d::Constant( 0.5 ) ) / 2.0;
    binnedImage.imageSize_ = image.imageSize_ / 2;

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

    LoadedPixelTestData loaded = loadPixelTestData( { image, binnedImage }, landmarks, bodies );

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
    for( const input_output::sum_lmk::SumImageData& currentImage : { image, binnedImage } )
    {
        for( const input_output::sum_lmk::SumLandmarkObservation& observation : currentImage.landmarkObservations_ )
        {
            LinkEnds linkEnds;
            linkEnds[ transmitter ] = LinkEndId( "Target", observation.landmarkId_ );
            linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
            const std::shared_ptr< ObservationModel< 2, double, double > > model =
                    ObservationModelCreator< 2, double, double >::createObservationModel(
                            pixelCoordinatesSettings( LinkDefinition( linkEnds ) ), bodies );
            const Eigen::Vector2d modelledPixel = model->computeObservations( sumImageTdbTimeForTest( currentImage ), receiver );

            // SPC projection: landmark relative to the spacecraft at VLM + SCOBJ, rotated by CX/CY/CZ, projected with
            // MMFL/CTR/K-MATRIX.
            const Eigen::Vector3d cameraFramePosition = currentImage.cameraAxes_ *
                    ( landmarks.at( observation.landmarkId_ ).bodyFixedPosition_ + currentImage.spacecraftObjectVector_ );
            const double x = currentImage.focalLengthMm_ * cameraFramePosition.x( ) / cameraFramePosition.z( );
            const double y = currentImage.focalLengthMm_ * cameraFramePosition.y( ) / cameraFramePosition.z( );
            const Eigen::Vector2d expectedPixel = currentImage.opticalCenter_ + currentImage.kMatrix_ * Eigen::Vector3d( x, y, x * y );
            BOOST_CHECK_SMALL( ( modelledPixel - expectedPixel ).norm( ), 1.0E-3 );
        }
    }

    // The modelled pixels of the full-resolution image agree with the measured ones, within the SPC solution residual.
    for( const input_output::sum_lmk::SumLandmarkObservation& observation : image.landmarkObservations_ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = LinkEndId( "Target", observation.landmarkId_ );
        linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
        const Eigen::Vector2d modelledPixel = ObservationModelCreator< 2, double, double >::createObservationModel(
                                                      pixelCoordinatesSettings( LinkDefinition( linkEnds ) ), bodies )
                                                      ->computeObservations( epoch, receiver );
        BOOST_CHECK_SMALL( ( modelledPixel - observation.pixelCoordinates_ ).norm( ), 3.0 );
    }
}

//! The pointing partial, with and without stellar aberration, at zero and non-zero corrections, of a camera with picture-specific
//! pointing (as created from SUM data) and of a camera fixed to its host body. The SUM camera takes two images with different
//! calibrations, and an observation of each image must depend only on the correction of that image. The first image also observes
//! a landmark without LMK data, which the optional skip mode must drop, keeping the valid observations. Receiver state partials
//! with aberration are tested in testPixelCoordinatesPartials.
BOOST_AUTO_TEST_CASE( testPixelPointingPartials )
{
    spice_interface::loadStandardSpiceKernels( );
    for( const bool aberration : { false, true } )
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        Eigen::Vector6d receiverState;
        receiverState << 200.0, -300.0, -10000.0, 30000.0, -20000.0, 10000.0;
        bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( receiverState, "SSB", "J2000" ) );
        input_output::sum_lmk::SumImageData firstImage = makeSyntheticSumImage( "IMGFIRST" );
        input_output::sum_lmk::SumLandmarkObservation missingObservation;
        missingObservation.landmarkId_ = "LMK_MISSING";
        missingObservation.pixelCoordinates_ = Eigen::Vector2d( 600.0, 600.0 );
        firstImage.landmarkObservations_.push_back( missingObservation );
        input_output::sum_lmk::SumImageData secondImage = makeSyntheticSumImage( "IMGSECOND" );
        secondImage.utcEpochString_ = "2015 JUN 05 07:34:42.053";
        secondImage.kMatrix_ *= 0.5;
        const auto loaded =
                loadPixelTestData( { firstImage, secondImage }, makeSyntheticLandmarks( ), bodies, "Target", "Spacecraft", true );

        // Only the observations of the landmark with LMK data are kept, one per image.
        BOOST_REQUIRE_EQUAL( loaded.observationCollection_->getConcatenatedObservations( ).size( ), 4 );
        BOOST_REQUIRE_EQUAL( loaded.observationModelSettings_.size( ), 1 );
        BOOST_CHECK_EQUAL( loaded.observationModelSettings_.at( 0 )->linkEnds_.linkEnds_.at( transmitter ).stationName_, "LMK0001" );

        // The camera finds the picture of an observation from its epoch.
        const std::vector< double > epochs = { sumImageTdbTimeForTest( firstImage ), sumImageTdbTimeForTest( secondImage ) };
        const std::shared_ptr< system_models::Camera > camera =
                bodies.at( "Spacecraft" )->getVehicleSystems( )->getCamera( sumLmkTestCameraName );
        BOOST_CHECK( camera->getPictureIds( ) == ( std::vector< std::string >{ "IMGFIRST", "IMGSECOND" } ) );
        BOOST_CHECK_EQUAL( camera->getPictureIndex( epochs.at( 0 ) ), 0 );
        BOOST_CHECK_EQUAL( camera->getPictureIndex( epochs.at( 1 ) ), 1 );

        // Ordinary camera creation: boresight Euler angles chosen such that the body-fixed-to-camera rotation is the identity, so
        // the camera looks along +Z body-fixed, towards the target.
        const std::string bodyFixedCameraName = "Camera_BodyFixed";
        createCamera( bodies.at( "Spacecraft" ),
                      bodyFixedCameraName,
                      ( Eigen::Vector3d( ) << -mathematical_constants::PI / 2.0, mathematical_constants::PI / 2.0, 0.0 ).finished( ),
                      std::make_pair( 1200.0, 900.0 ),
                      std::make_pair( 320.0, 240.0 ),
                      Eigen::Vector3d::Zero( ) );
        const std::shared_ptr< system_models::Camera > bodyFixedCamera =
                bodies.at( "Spacecraft" )->getVehicleSystems( )->getCamera( bodyFixedCameraName );
        BOOST_REQUIRE( !bodyFixedCamera->hasRotationFromInertialToCameraFrameFunction( ) );
        BOOST_CHECK_SMALL(
                ( bodyFixedCamera->getRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( ) - Eigen::Matrix3d::Identity( ) ).norm( ),
                1.0E-12 );

        for( const std::string& cameraName : { sumLmkTestCameraName, bodyFixedCameraName } )
        {
            LinkEnds linkEnds;
            linkEnds[ transmitter ] = LinkEndId( "Target", "LMK0001" );
            linkEnds[ receiver ] = LinkEndId( "Spacecraft", cameraName );
            const auto model = ObservationModelCreator< 2, double, double >::createObservationModel(
                    pixelCoordinatesSettings(
                            LinkDefinition( linkEnds ), {}, nullptr, std::make_shared< LightTimeConvergenceCriteria >( ), aberration ),
                    bodies );
            const auto parameter = std::make_shared< CameraPointingCorrection >(
                    bodies.at( "Spacecraft" )->getVehicleSystems( ), "Spacecraft", cameraName );
            const auto parameters = std::make_shared< EstimatableParameterSet< double > >(
                    std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                    std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >{ parameter } );
            const auto partials = ObservationPartialCreator< 2, double, double >::createObservationPartials( model, bodies, parameters );
            BOOST_REQUIRE_EQUAL( partials.first.size( ), 1 );
            BOOST_REQUIRE( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( partials.first.begin( )->second ) != nullptr );

            // Two pictures for the SUM camera, a single correction for the body-fixed camera.
            const int numberOfCorrections = parameter->getParameterSize( ) / 3;
            BOOST_REQUIRE_EQUAL( numberOfCorrections, ( cameraName == bodyFixedCameraName ) ? 1 : 2 );
            const std::function< Eigen::VectorXd( const double ) > observable = [ & ]( const double time ) {
                return Eigen::VectorXd( model->computeObservations( time, receiver ) );
            };
            // Observations of both images, which have different calibrations: each depends only on the correction of its own image.
            for( unsigned int imageIndex = 0; imageIndex < epochs.size( ); ++imageIndex )
            {
                const double epoch = epochs.at( imageIndex );
                const int pictureIndex = ( numberOfCorrections == 1 ) ? 0 : imageIndex;
                for( const Eigen::Vector3d& correction : std::vector< Eigen::Vector3d >{ Eigen::Vector3d::Zero( ),
                                                                                         Eigen::Vector3d( 1.0E-10, -2.0E-10, 3.0E-10 ),
                                                                                         Eigen::Vector3d( 1.0E-6, -2.0E-6, 3.0E-6 ),
                                                                                         Eigen::Vector3d( 0.02, -0.01, 0.03 ) } )
                {
                    Eigen::VectorXd parameterValue( 3 * numberOfCorrections );
                    for( int i = 0; i < numberOfCorrections; ++i )
                    {
                        parameterValue.segment( 3 * i, 3 ) = ( 1.0 - 3.0 * i ) * correction;
                    }
                    parameter->setParameterValue( parameterValue );
                    std::vector< double > times;
                    std::vector< Eigen::Vector6d > states;
                    const Eigen::Vector2d observation = model->computeObservationsWithLinkEndData( epoch, receiver, times, states );
                    partials.second->update( states, times, receiver, observation );
                    const auto analytical =
                            partials.first.begin( )->second->calculatePartial( states, times, receiver, nullptr, observation );
                    BOOST_REQUIRE_EQUAL( analytical.size( ), 1 );
                    BOOST_CHECK_EQUAL( analytical.front( ).second, epoch );
                    const Eigen::MatrixXd numerical = calculateNumericalObservationParameterPartial(
                            parameter, Eigen::VectorXd::Constant( 3 * numberOfCorrections, 1.0E-7 ), observable, epoch );
                    BOOST_CHECK_SMALL( ( analytical.front( ).first - numerical ).norm( ), 1.0E-5 );
                    Eigen::MatrixXd otherPicturesPartial = analytical.front( ).first;
                    BOOST_CHECK( otherPicturesPartial.middleCols( 3 * pictureIndex, 3 ).norm( ) > 1.0E-2 );
                    otherPicturesPartial.middleCols( 3 * pictureIndex, 3 ).setZero( );
                    BOOST_CHECK( otherPicturesPartial.isZero( 0.0 ) );
                }
            }
            parameter->setParameterValue( Eigen::VectorXd::Zero( 3 * numberOfCorrections ) );
        }
    }
}

//! A target rotation changes both the landmark position and the camera's reference frame: the sum of both contributions is
//! tested for a rotation model parameter (with and without stellar aberration), and for the target's rotational state.
BOOST_AUTO_TEST_CASE( testPixelReferenceBodyRotationPartials )
{
    spice_interface::loadStandardSpiceKernels( );
    const auto image = makeSyntheticSumImage( "IMGROTATION" );
    const double epoch = sumImageTdbTimeForTest( image );
    auto landmarks = makeSyntheticLandmarks( );
    landmarks.at( "LMK0001" ).bodyFixedPosition_ = Eigen::Vector3d( 400.0, -250.0, 100.0 );

    // Rotation rate, with a moving receiver so that aberration is relevant
    for( const bool aberration : { false, true } )
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        Eigen::Vector6d receiverState;
        receiverState << 1500.0, -700.0, -10000.0, 30000.0, -20000.0, 10000.0;
        bodies.at( "Spacecraft" )->setEphemeris( std::make_shared< ConstantEphemeris >( receiverState, "SSB", "J2000" ) );
        const auto rotation =
                std::make_shared< SimpleRotationalEphemeris >( Eigen::Quaterniond( Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitX( ) ) ),
                                                               1.0E-4,
                                                               epoch - 1000.0,
                                                               "J2000",
                                                               "Target_Fixed" );
        const auto parameter = std::make_shared< RotationRate >( rotation, "Target" );
        const auto modelAndPartial =
                createModelAndRotationPartial( bodies,
                                               image,
                                               landmarks,
                                               rotation,
                                               std::make_shared< EstimatableParameterSet< double > >(
                                                       std::vector< std::shared_ptr< EstimatableParameter< double > > >{ parameter },
                                                       std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ) ),
                                               aberration );
        const std::function< Eigen::VectorXd( const double ) > observable = [ & ]( const double time ) {
            return Eigen::VectorXd( modelAndPartial.first->computeObservations( time, receiver ) );
        };
        const Eigen::VectorXd numerical = calculateNumericalObservationParameterPartial( parameter, 1.0E-8, observable, epoch );
        BOOST_CHECK_SMALL( ( modelAndPartial.second - numerical ).norm( ) / numerical.norm( ), 1.0E-7 );
    }

    // Initial rotational state: only the quaternion entries contribute, compared to central differences in the normalized quaternion
    {
        SystemOfBodies bodies = createSyntheticSumLmkBodies( );
        Eigen::Vector7d rotationalState = Eigen::Vector7d::Zero( );
        rotationalState.head< 4 >( ) = linear_algebra::convertQuaternionToVectorFormat(
                Eigen::Quaterniond( Eigen::AngleAxisd( 0.15, Eigen::Vector3d( 1.0, 2.0, -1.0 ).normalized( ) ) ) );
        const auto rotation = std::make_shared< ConstantRotationalEphemeris >( rotationalState, "J2000", "Target_Fixed" );
        const auto parameter = std::make_shared< InitialRotationalStateParameter< double > >(
                "Target", rotationalState, []( ) { return Eigen::Matrix3d::Identity( ); }, "SSB", "J2000" );
        const auto modelAndPartial = createModelAndRotationPartial(
                bodies,
                image,
                landmarks,
                rotation,
                std::make_shared< EstimatableParameterSet< double > >(
                        std::vector< std::shared_ptr< EstimatableParameter< double > > >( ),
                        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ),
                        std::vector< std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >{ parameter } ),
                false );
        BOOST_REQUIRE_EQUAL( modelAndPartial.second.cols( ), 7 );
        const Eigen::Matrix< double, 2, 7 > analytical = modelAndPartial.second;
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
            const Eigen::Vector2d up = modelAndPartial.first->computeObservations( epoch, receiver );
            perturbed( component ) -= 2.0 * step;
            rotation->updateConstantState( perturbed );
            const Eigen::Vector2d down = modelAndPartial.first->computeObservations( epoch, receiver );
            BOOST_CHECK_SMALL( ( normalizedAnalytical.col( component ) - ( up - down ) / ( 2.0 * step ) ).norm( ), 1.0E-5 );
        }
        rotation->updateConstantState( rotationalState );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
