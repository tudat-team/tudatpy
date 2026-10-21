/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <memory>
#include <string>
#include <vector>

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/ephemerides/constantRotationalEphemeris.h"
#include "tudat/math/basic/rotationRepresentations.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createCameras.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::ephemerides;
using namespace tudat::observation_models;
using namespace tudat::physical_constants;
using namespace tudat::simulation_setup;
using namespace tudat::system_models;

BOOST_AUTO_TEST_SUITE( test_pixel_coordinates_observation_model )

Eigen::Quaterniond getTestRotationFromInertialToBodyFixed( )
{
    return Eigen::Quaterniond( Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitY( ) ) *
                               Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitX( ) ) );
}

BOOST_AUTO_TEST_CASE( testPixelCoordinatesObservationModelFactory )
{
    const std::string targetName = "Target";
    const std::string observerName = "Observer";
    const std::string cameraName = "Camera";
    const std::string globalFrame = "ECLIPJ2000";
    const std::string observerFixedFrame = "Observer_fixed";

    Eigen::Vector6d observerState = Eigen::Vector6d::Zero( );
    observerState.segment( 0, 3 ) = ( Eigen::Vector3d( ) << 1000.0, -2000.0, 3000.0 ).finished( );

    const Eigen::Vector3d cameraBodyFixedPosition = ( Eigen::Vector3d( ) << 10.0, -20.0, 30.0 ).finished( );
    const Eigen::Vector3d relativeTargetPositionInCameraFrame = ( Eigen::Vector3d( ) << 700.0, -300.0, 4.0E6 ).finished( );
    const Eigen::Vector3d cameraBoresightEulerAngles = ( Eigen::Vector3d( ) << 0.2, -0.1, 0.3 ).finished( );
    std::shared_ptr< CameraSettings > cameraSettings = std::make_shared< CameraSettings >( cameraName,
                                                                                           cameraBoresightEulerAngles,
                                                                                           std::make_pair( 1200.0, 900.0 ),
                                                                                           std::make_pair( 320.0, 240.0 ),
                                                                                           cameraBodyFixedPosition );
    const Eigen::Quaterniond rotationFromBodyFixedToCameraFrame =
            basic_mathematics::getQuaternionFrom313EulerAngles( cameraSettings->getCamera313EulerAngles( ) );
    const Eigen::Quaterniond rotationFromInertialToBodyFixed = getTestRotationFromInertialToBodyFixed( );
    const Eigen::Vector3d cameraInertialPositionOffset = rotationFromInertialToBodyFixed.inverse( ) * cameraBodyFixedPosition;
    const Eigen::Vector3d relativeTargetPositionInInertialFrame = rotationFromInertialToBodyFixed.inverse( ) *
            rotationFromBodyFixedToCameraFrame.inverse( ) * relativeTargetPositionInCameraFrame;

    Eigen::Vector6d targetState = Eigen::Vector6d::Zero( );
    targetState.segment( 0, 3 ) = observerState.segment( 0, 3 ) + cameraInertialPositionOffset + relativeTargetPositionInInertialFrame;

    SystemOfBodies bodies( "SSB", globalFrame );
    bodies.createEmptyBody( targetName, false );
    bodies.createEmptyBody( observerName, false );
    bodies.at( targetName )->setEphemeris( std::make_shared< ConstantEphemeris >( targetState, "SSB", globalFrame ) );
    bodies.at( observerName )->setEphemeris( std::make_shared< ConstantEphemeris >( observerState, "SSB", globalFrame ) );
    bodies.processBodyFrameDefinitions( );
    bodies.at( observerName )
            ->setRotationalEphemeris( std::make_shared< ConstantRotationalEphemeris >(
                    rotationFromInertialToBodyFixed.inverse( ), globalFrame, observerFixedFrame ) );

    createCamera( bodies.at( observerName ), cameraSettings );
    std::shared_ptr< Camera > camera = bodies.at( observerName )->getVehicleSystems( )->getCamera( cameraName );
    BOOST_CHECK_SMALL( ( camera->getRotationFromBodyFixedToCameraFrame( ).toRotationMatrix( ) -
                         rotationFromBodyFixedToCameraFrame.toRotationMatrix( ) )
                               .norm( ),
                       10.0 * std::numeric_limits< double >::epsilon( ) );

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( targetName, "" );
    linkEnds[ receiver ] = LinkEndId( observerName, cameraName );

    const Eigen::Vector2d absoluteBias = ( Eigen::Vector2d( ) << 0.25, -0.5 ).finished( );
    std::shared_ptr< ObservationModelSettings > idealObservableSettings = pixelCoordinatesSettings( linkEnds );
    std::shared_ptr< ObservationModelSettings > biasedObservableSettings =
            pixelCoordinatesSettings( linkEnds,
                                      std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
                                      std::make_shared< ConstantObservationBiasSettings >( absoluteBias, true ) );

    std::shared_ptr< ObservationModel< 2, double, double > > idealObservationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( idealObservableSettings, bodies );
    std::shared_ptr< ObservationModel< 2, double, double > > biasedObservationModel =
            ObservationModelCreator< 2, double, double >::createObservationModel( biasedObservableSettings, bodies );

    const double receiverObservationTime = 12345.0;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;
    Eigen::Vector2d idealObservation = idealObservationModel->computeIdealObservationsWithLinkEndData(
            receiverObservationTime, receiver, linkEndTimes, linkEndStates );

    const Eigen::Vector2d expectedObservation =
            ( Eigen::Vector2d( ) << cameraSettings->getOpticalCenter( ).first +
                              cameraSettings->getFocalLengths( ).first * relativeTargetPositionInCameraFrame.x( ) /
                                      relativeTargetPositionInCameraFrame.z( ),
              cameraSettings->getOpticalCenter( ).second +
                      cameraSettings->getFocalLengths( ).second * relativeTargetPositionInCameraFrame.y( ) /
                              relativeTargetPositionInCameraFrame.z( ) )
                    .finished( );

    BOOST_CHECK_EQUAL( linkEndTimes.size( ), 2 );
    BOOST_CHECK_EQUAL( linkEndStates.size( ), 2 );
    const double comparisonTolerance = 1.0E-12;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( idealObservation, expectedObservation, comparisonTolerance );

    const double expectedLightTime = relativeTargetPositionInCameraFrame.norm( ) / SPEED_OF_LIGHT;
    const double timeTolerance = 100.0 * std::numeric_limits< double >::epsilon( ) * receiverObservationTime;
    const double stateTolerance = 1.0E-8;
    BOOST_CHECK_SMALL( linkEndTimes.at( 0 ) - ( receiverObservationTime - expectedLightTime ), timeTolerance );
    BOOST_CHECK_SMALL( linkEndTimes.at( 1 ) - receiverObservationTime, timeTolerance );
    BOOST_CHECK_SMALL( ( linkEndStates.at( 0 ) - targetState ).norm( ), stateTolerance );

    Eigen::Vector6d expectedReceiverState = observerState;
    expectedReceiverState.segment( 0, 3 ) += cameraInertialPositionOffset;
    BOOST_CHECK_SMALL( ( linkEndStates.at( 1 ) - expectedReceiverState ).norm( ), stateTolerance );

    Eigen::Vector2d biasedObservation = biasedObservationModel->computeObservations( receiverObservationTime, receiver );
    Eigen::Vector2d expectedBiasedObservation = expectedObservation + absoluteBias;
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( biasedObservation, expectedBiasedObservation, comparisonTolerance );

    BOOST_CHECK_THROW( idealObservationModel->computeIdealObservationsWithLinkEndData(
                               receiverObservationTime, transmitter, linkEndTimes, linkEndStates ),
                       std::runtime_error );

    BOOST_CHECK_THROW( camera->calculateObservableFromBodyFixed( rotationFromBodyFixedToCameraFrame.inverse( ) *
                                                                 ( Eigen::Vector3d( ) << 1.0, 1.0, -1.0 ).finished( ) ),
                       std::runtime_error );
    BOOST_CHECK_THROW( camera->calculateObservableFromBodyFixed(
                               rotationFromBodyFixedToCameraFrame.inverse( ) *
                               ( Eigen::Vector3d( ) << 1.0, 1.0, 0.5 * std::numeric_limits< double >::epsilon( ) ).finished( ) ),
                       std::runtime_error );

    LinkEnds linkEndsWithoutCamera;
    linkEndsWithoutCamera[ transmitter ] = LinkEndId( targetName, "" );
    linkEndsWithoutCamera[ receiver ] = LinkEndId( observerName, "" );
    BOOST_CHECK_THROW( pixelCoordinatesSettings( linkEndsWithoutCamera ), std::runtime_error );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
