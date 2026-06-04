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

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include <boost/test/unit_test.hpp>
#include <Eigen/Geometry>

#include "tudat/basics/testMacros.h"
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/observation_models/azimuthElevationObservationModel.h"
#include "tudat/astro/orbit_determination/observation_partials/azimuthElevationPartial.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createObservationPartials.h"
#include "tudat/support/observationPartialTestFunctions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::estimatable_parameters;
using namespace tudat::observation_models;
using namespace tudat::observation_partials;
using namespace tudat::simulation_setup;
using namespace tudat::spice_interface;

BOOST_AUTO_TEST_SUITE( test_azimuth_elevation_partials )

Eigen::Vector2d calculateAzimuthElevationFromAngularPosition( const Eigen::Vector2d& angularPosition,
                                                              const Eigen::Matrix3d& rotationFromInertialToTopocentricFrame,
                                                              const bool invertLineOfSight )
{
    const double rightAscension = angularPosition.x( );
    const double declination = angularPosition.y( );
    const double lineOfSightSign = invertLineOfSight ? -1.0 : 1.0;
    Eigen::Vector3d inertialLineOfSight = ( Eigen::Vector3d( ) << std::cos( declination ) * std::cos( rightAscension ),
                                            std::cos( declination ) * std::sin( rightAscension ),
                                            std::sin( declination ) )
                                                  .finished( );
    Eigen::Vector3d topocentricLineOfSight = lineOfSightSign * rotationFromInertialToTopocentricFrame * inertialLineOfSight;

    return ( Eigen::Vector2d( ) << std::atan2( topocentricLineOfSight.x( ), topocentricLineOfSight.y( ) ),
             std::atan2( topocentricLineOfSight.z( ), topocentricLineOfSight.segment( 0, 2 ).norm( ) ) )
            .finished( );
}

Eigen::Vector2d calculateAzimuthElevationFromLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< ground_stations::PointingAnglesCalculator >& pointingAnglesCalculator,
        const LinkEndType stationLinkEndType )
{
    const int stationIndex = ( stationLinkEndType == transmitter ) ? 0 : 1;
    const int targetIndex = ( stationLinkEndType == transmitter ) ? 1 : 0;
    std::pair< double, double > elevationAzimuth = pointingAnglesCalculator->calculatePointingAngles(
            ( linkEndStates.at( targetIndex ) - linkEndStates.at( stationIndex ) ).segment( 0, 3 ), linkEndTimes.at( stationIndex ) );
    return ( Eigen::Vector2d( ) << elevationAzimuth.second, elevationAzimuth.first ).finished( );
}

BOOST_AUTO_TEST_CASE( testAzimuthElevationWrtAngularPositionPartial )
{
    Eigen::Matrix3d rotationFromInertialToTopocentricFrame =
            ( Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitZ( ) ) * Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitY( ) ) *
              Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitX( ) ) )
                    .toRotationMatrix( );
    Eigen::Vector2d angularPosition = ( Eigen::Vector2d( ) << 0.7, 0.3 ).finished( );

    for( int invertLineOfSight = 0; invertLineOfSight < 2; invertLineOfSight++ )
    {
        Eigen::Matrix2d analyticalPartial = calculatePartialOfAzimuthElevationWrtAngularPosition(
                angularPosition, rotationFromInertialToTopocentricFrame, static_cast< bool >( invertLineOfSight ) );
        Eigen::Matrix2d numericalPartial = Eigen::Matrix2d::Zero( );

        const double perturbation = 1.0E-7;
        for( int i = 0; i < 2; i++ )
        {
            Eigen::Vector2d perturbedAngularPosition = angularPosition;
            perturbedAngularPosition( i ) += perturbation;
            Eigen::Vector2d upPerturbedObservation = calculateAzimuthElevationFromAngularPosition(
                    perturbedAngularPosition, rotationFromInertialToTopocentricFrame, static_cast< bool >( invertLineOfSight ) );

            perturbedAngularPosition = angularPosition;
            perturbedAngularPosition( i ) -= perturbation;
            Eigen::Vector2d downPerturbedObservation = calculateAzimuthElevationFromAngularPosition(
                    perturbedAngularPosition, rotationFromInertialToTopocentricFrame, static_cast< bool >( invertLineOfSight ) );

            numericalPartial.col( i ) = ( upPerturbedObservation - downPerturbedObservation ) / ( 2.0 * perturbation );
        }

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalPartial, numericalPartial, 1.0E-7 );
    }
}

BOOST_AUTO_TEST_CASE( testAzimuthElevationScaling )
{
    Eigen::Quaterniond rotationFromInertialToBodyFixedFrame( Eigen::Quaterniond::Identity( ) );
    Eigen::Quaterniond rotationFromBodyFixedToTopocentricFrame( Eigen::AngleAxisd( 0.4, Eigen::Vector3d::UnitZ( ) ) *
                                                                Eigen::AngleAxisd( -0.2, Eigen::Vector3d::UnitY( ) ) *
                                                                Eigen::AngleAxisd( 0.1, Eigen::Vector3d::UnitX( ) ) );

    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            std::make_shared< ground_stations::PointingAnglesCalculator >(
                    [ = ]( const double ) { return rotationFromInertialToBodyFixedFrame; },
                    [ = ]( const double ) { return rotationFromBodyFixedToTopocentricFrame; } );

    std::vector< Eigen::Vector6d > linkEndStates( 2, Eigen::Vector6d::Zero( ) );
    linkEndStates.at( 0 ).segment( 0, 3 ) = ( Eigen::Vector3d( ) << -2.0E7, 1.1E8, 4.0E7 ).finished( );
    linkEndStates.at( 1 ).segment( 0, 3 ) = ( Eigen::Vector3d( ) << 8.0E6, -3.0E7, 1.5E7 ).finished( );
    std::vector< double > linkEndTimes = { 10.0, 20.0 };

    const double positionPerturbation = 10.0;
    for( int stationLinkEndIndex = 0; stationLinkEndIndex < 2; stationLinkEndIndex++ )
    {
        LinkEndType stationLinkEndType = ( stationLinkEndIndex == 0 ) ? receiver : transmitter;
        AzimuthElevationScaling scaling( pointingAnglesCalculator, stationLinkEndType );
        scaling.update( linkEndStates, linkEndTimes, receiver, Eigen::Vector2d::Zero( ) );

        for( LinkEndType linkEndType : { transmitter, receiver } )
        {
            const int linkEndIndex = ( linkEndType == transmitter ) ? 0 : 1;
            Eigen::Matrix< double, 2, 3 > numericalScaling = Eigen::Matrix< double, 2, 3 >::Zero( );

            for( int i = 0; i < 3; i++ )
            {
                std::vector< Eigen::Vector6d > perturbedLinkEndStates = linkEndStates;
                perturbedLinkEndStates.at( linkEndIndex )( i ) += positionPerturbation;
                Eigen::Vector2d upPerturbedObservation = calculateAzimuthElevationFromLinkEndStates(
                        perturbedLinkEndStates, linkEndTimes, pointingAnglesCalculator, stationLinkEndType );

                perturbedLinkEndStates = linkEndStates;
                perturbedLinkEndStates.at( linkEndIndex )( i ) -= positionPerturbation;
                Eigen::Vector2d downPerturbedObservation = calculateAzimuthElevationFromLinkEndStates(
                        perturbedLinkEndStates, linkEndTimes, pointingAnglesCalculator, stationLinkEndType );

                numericalScaling.col( i ) = ( upPerturbedObservation - downPerturbedObservation ) / ( 2.0 * positionPerturbation );
            }

            Eigen::Matrix< double, 2, 3 > analyticalScaling = scaling.getFixedTimePositionScalingFactor( linkEndType );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( analyticalScaling, numericalScaling, 1.0E-7 );
        }
    }
}

BOOST_AUTO_TEST_CASE( testAzimuthElevationPartials )
{
    Eigen::VectorXd parameterPerturbationMultipliers = Eigen::VectorXd::Constant( 4, 1.0 );
    parameterPerturbationMultipliers( 2 ) = 10.0;

    // Define and create ground stations.
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.resize( 2 );
    groundStations[ 0 ] = std::make_pair( "Earth", "Graz" );
    groundStations[ 1 ] = std::make_pair( "Mars", "MSL" );

    // Test partials with constant ephemerides (allows test of position partials)
    {
        for( bool useOblateEarthShape : { false, true } )
        {
            BOOST_TEST_CONTEXT( "Oblate Earth shape: " << useOblateEarthShape )
            {
                // Create environment
                SystemOfBodies bodies =
                        setupEnvironment( groundStations, 1.0E7, 1.2E7, 1.1E7, true, 1.0, false, false, useOblateEarthShape );

                // Set link ends for observation model
                LinkEnds linkEnds;
                linkEnds[ transmitter ] = groundStations[ 1 ];
                linkEnds[ receiver ] = groundStations[ 0 ];

                // Generate azimuth/elevation model
                std::vector< std::string > perturbingBodies;
                perturbingBodies.push_back( "Earth" );
                std::shared_ptr< ObservationModel< 2 > > azimuthElevationModel =
                        observation_models::ObservationModelCreator< 2, double, double >::createObservationModel(
                                observation_models::azimuthElevationSettings(
                                        linkEnds,
                                        { std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( perturbingBodies ) } ),
                                bodies );

                // Create parameter objects.
                std::shared_ptr< EstimatableParameterSet< double > > fullEstimatableParameterSet =
                        createEstimatableParameters( bodies, 1.1E7 );

                for( double positionPerturbationMultiplier : { 1.0, 1.0E-1, 1.0E-2, 1.0E-3 } )
                {
                    BOOST_TEST_CONTEXT( "Position perturbation multiplier: " << positionPerturbationMultiplier )
                    {
                        testObservationPartials( azimuthElevationModel,
                                                 bodies,
                                                 fullEstimatableParameterSet,
                                                 linkEnds,
                                                 azimuth_elevation_angle,
                                                 1.0E-4,
                                                 true,
                                                 false,
                                                 positionPerturbationMultiplier,
                                                 parameterPerturbationMultipliers );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests
}  // namespace tudat
