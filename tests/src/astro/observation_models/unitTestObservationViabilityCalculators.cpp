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

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>


#include "tudat/basics/testMacros.h"

#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::ground_stations;
using namespace tudat::unit_conversions;

BOOST_AUTO_TEST_SUITE( test_observation_viability_calculators )

BOOST_AUTO_TEST_CASE( testSeparateObservationViabilityCalculators )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create default Earth
    std::vector< std::string > bodyNames = { "Earth" };
    BodyListSettings bodySettings = getDefaultBodySettings(
                bodyNames, "Earth", "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define ground station
    double stationRadius = spice_interface::getAverageRadius( "Earth" );
    double stationLatitude = convertDegreesToRadians( 20.0 );
    double stationLongitude = convertDegreesToRadians( 50.0 );
    createGroundStation( bodies.at( "Earth" ), "Station",
                         ( Eigen::Vector3d( ) << stationRadius, stationLatitude, stationLongitude ).finished( ),
                         coordinate_conversions::spherical_position );

    // Get inertial ground station state function
    std::function< Eigen::Vector6d( const double ) > groundStationStateFunction =
            getLinkEndCompleteEphemerisFunction( std::make_pair< std::string, std::string >( "Earth", "Station" ), bodies );

    // Get Earth-fixed ground station position
    Eigen::Vector3d earthFixedGroundStationState =
            bodies.at( "Earth" )->getGroundStation( "Station" )->getNominalStationState( )->getCartesianStateInTime( 0.0, "Earth" ).segment( 0, 3 );

    // Define limiting elevation angle for test
    double testAngle = 20.0 * mathematical_constants::PI / 180.0;

    // Retrieve pointing angle calculator
    std::shared_ptr< PointingAnglesCalculator > pointingAngleCalculator =
            bodies.at( "Earth" )->getGroundStation( "Station" )->getPointingAnglesCalculator( );

    // Define test body for occultation/avoidance tests
    Eigen::Vector6d testBodyState = Eigen::Vector6d::Zero( );
    testBodyState( 0 ) = 8000.0E3;
    double testBodyRadius = 1000.0E3;

    // Declare test variables
    Eigen::Vector3d manualGroundStationInertialPosition;
    Eigen::Vector3d vectorToTarget;
    Eigen::Vector3d vectorToTestBody;
    Eigen::Vector6d targetState = Eigen::Vector6d::Zero( );
//    Eigen::Vector6d groundStationState = Eigen::Vector6d::Zero( );
//    Eigen::Vector3d localVerticalVector = Eigen::Vector3d::UnitZ( );

    // Run different test cases:
    // 0: one-way uplink
    // 1: one-way downlink
    for( unsigned int i = 0; i < 2; i++ )
    {
        int targetIndex = -1;
        int groundStationIndex = -1;
        std::vector< std::pair< int, int > > linkEndIndices;
        if( i == 0 )
        {
            targetIndex = 1;
            groundStationIndex = 0;
            linkEndIndices = { { 0, 1 } };
        }
        else if( i == 1 )
        {
            targetIndex = 0;
            groundStationIndex = 1;
            linkEndIndices = { { 1, 0 } };
        }
        std::shared_ptr< ObservationViabilityCalculator > minimumElevationCalculator =
                std::make_shared< MinimumElevationAngleCalculator >(
                    linkEndIndices, testAngle, pointingAngleCalculator );
        std::shared_ptr< ObservationViabilityCalculator > bodyAvoidanceCalculator =
                std::make_shared< BodyAvoidanceAngleCalculator >(
                    linkEndIndices, testAngle, [=](const double){return testBodyState;}, "TestBody" );
        std::shared_ptr< ObservationViabilityCalculator > bodyOccultationCalculator =
                std::make_shared< OccultationCalculator >(
                    linkEndIndices, [=](const double){return testBodyState;}, testBodyRadius );


        std::vector< double > linkEndTimes;
        linkEndTimes.resize( 2 );

        std::vector< Eigen::Vector6d > linkEndStates;
        linkEndStates.resize( 2 );

        // Test viability settings for ~100 orbits
        for( int j = 0; j < 20000; j++ )
        {
            // Set spacraft state in circular orbit (velocity not set)
            targetState( 0 ) = std::sin( ( static_cast< double >( j ) / 200.0 ) * 2.0 * mathematical_constants::PI ) * 10000.0E3;
            targetState( 2 ) = std::cos( ( static_cast< double >( j ) / 200.0 ) * 2.0 * mathematical_constants::PI ) * 10000.0E3;
            // Set link end times (exagerated light time for testing purposes)
            linkEndTimes[ 0 ] = 2.0 * 3600.0 * ( static_cast< double >( j ) / 200.0 );
            linkEndTimes[ 1 ] = 2.0 * 3600.0 * ( static_cast< double >( j ) / 200.0 ) + 100.0;

            linkEndStates[ groundStationIndex ] = groundStationStateFunction( linkEndTimes[ groundStationIndex ] );
            linkEndStates[ targetIndex ] = targetState;

            // Manually calculate elevation angle
            manualGroundStationInertialPosition = spice_interface::computeRotationMatrixBetweenFrames(
                        "IAU_Earth", "J2000", linkEndTimes[ groundStationIndex ] ) * earthFixedGroundStationState;
            vectorToTarget = ( targetState.segment( 0, 3 ) - manualGroundStationInertialPosition );
            vectorToTestBody = ( testBodyState.segment( 0, 3 ) - manualGroundStationInertialPosition );

            // Test elevation angle
            {
                double manualElevationAngle = mathematical_constants::PI / 2.0 -
                        std::acos( manualGroundStationInertialPosition.normalized( ).dot( vectorToTarget.normalized( ) ) );

                // Compute viability and check against manual calculation
                bool isObservationViable = minimumElevationCalculator->isObservationViable(
                            linkEndStates, linkEndTimes );
                if( manualElevationAngle > testAngle )
                {
                    BOOST_CHECK_EQUAL( isObservationViable, 1 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( isObservationViable, 0 );
                }
            }

            // Test avoidance angle
            {
                double manualAvoidanceAngle =
                        std::acos( ( vectorToTestBody.normalized( ) ).dot( vectorToTarget.normalized( ) ) );

                // Compute viability and check against manual calculation
                bool isObservationViable = bodyAvoidanceCalculator->isObservationViable(
                            linkEndStates, linkEndTimes );
                if( manualAvoidanceAngle > testAngle )
                {
                    BOOST_CHECK_EQUAL( isObservationViable, 1 );
                }
                else
                {
                    BOOST_CHECK_EQUAL( isObservationViable, 0 );
                }
            }

            // Test occultation
            {

                bool manualIsObservationViable = true;
                if( computeOccultation(
                            manualGroundStationInertialPosition,
                            targetState.segment( 0, 3 ),
                            testBodyState.segment( 0, 3 ),
                            testBodyRadius ) )
                {
                    manualIsObservationViable = false;
                }

                // Compute viability and check against manual calculation
                bool isObservationViable = bodyOccultationCalculator->isObservationViable(
                            linkEndStates, linkEndTimes );
                BOOST_CHECK_EQUAL( isObservationViable, manualIsObservationViable );
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( testStationAngleCalculations )
{

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Create default Earth
    std::vector< std::string > bodyNames = { "Earth", "Moon" };
    BodyListSettings bodySettings = getDefaultBodySettings(
                bodyNames, "Earth", "J2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Define ground station
    double stationRadius = spice_interface::getAverageRadius( "Earth" );
    double stationLatitude = convertDegreesToRadians( 20.0 );
    double stationLongitude = convertDegreesToRadians( 50.0 );
    createGroundStation( bodies.at( "Earth" ), "Station",
                         ( Eigen::Vector3d( ) << stationRadius, stationLatitude, stationLongitude ).finished( ),
                         coordinate_conversions::spherical_position );

    // Get inertial ground station state function
    std::function< Eigen::Vector6d( const double ) > groundStationStateFunction =
            getLinkEndCompleteEphemerisFunction( std::make_pair< std::string, std::string >( "Earth", "Station" ), bodies );

    // Get Earth-fixed ground station position
//    Eigen::Vector3d earthFixedGroundStationState =
//            bodies.at( "Earth" )->getGroundStation( "Station" )->getNominalStationState( )->getCartesianStateInTime( 0.0 ).segment( 0, 3 );

    // Define limiting elevation angle for test
    double testAngle = 20.0 * mathematical_constants::PI / 180.0;

    // Retrieve pointing angle calculator
    std::shared_ptr< PointingAnglesCalculator > pointingAngleCalculator =
            bodies.at( "Earth" )->getGroundStation( "Station" )->getPointingAnglesCalculator( );

    std::vector< double > times;
    for( int i = 0; i < 1000; i++ )
    {
        times.push_back( static_cast< double >( i ) * 3600.0 );
    }

    std::map< double, Eigen::VectorXd > targetAnglesAndRange = getTargetAnglesAndRange(
            bodies, std::make_pair< std::string, std::string >( "Earth", "Station" ),
            "Moon", times, true );
    std::map< double, Eigen::VectorXd > targetAnglesAndRange2 = getTargetAnglesAndRange(
            bodies, std::make_pair< std::string, std::string >( "Earth", "Station" ),
            "Moon", times, false );
    for( auto it : targetAnglesAndRange )
    {
        double time = it.first;
        Eigen::Vector6d stateOfMoon = bodies.at( "Moon" )->getStateInBaseFrameFromEphemeris( time );
        Eigen::Vector6d relativeState =  stateOfMoon - groundStationStateFunction( time );
        std::pair< double, double > angles = pointingAngleCalculator->calculatePointingAngles( relativeState.segment( 0, 3 ), time );
        double elevationAngleDifference = angles.first - ( it.second( 0 ) + targetAnglesAndRange2.at( time )( 0 ) ) / 2;
        double azimuthAngleDifference = angles.second - ( it.second( 1 ) + targetAnglesAndRange2.at( time )( 1 ) ) / 2;
        double rangeDifference = relativeState.segment( 0, 3 ).norm( ) - ( it.second( 2 ) + targetAnglesAndRange2.at( time )( 2 ) ) / 2;
        BOOST_CHECK_SMALL( elevationAngleDifference, 1.0E-7 );
        BOOST_CHECK_SMALL( azimuthAngleDifference, 1.0E-7 );
        BOOST_CHECK_SMALL( rangeDifference, 1.0E-3 );

    }
}


//! Test function to manually compute elevation angle(s) of observation link at given body
std::vector< double > getBodyLinkElevationAngles(
        const LinkEnds linkEnds,
        const ObservableType observableType,
        const std::string referenceBody,
        const std::vector< Eigen::Vector6d > linkEndStates,
        const std::vector< double > linkEndTimes,
        const SystemOfBodies& bodies )
{
    std::shared_ptr< ground_stations::PointingAnglesCalculator > currentPointingAnglesCalculator;
    std::vector< double > elevationAngles;
    switch( observableType )
    {
    case one_way_range:
    case one_way_doppler:
    case angular_position:
    {
        if( linkEnds.at( transmitter ).bodyName_ == referenceBody )
        {
            currentPointingAnglesCalculator = bodies.at( referenceBody )->getGroundStation(
                        linkEnds.at( transmitter ).stationName_ )->getPointingAnglesCalculator( );
            elevationAngles.push_back( currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
                                           ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ), linkEndTimes.at( 0 ) ) );

        }
        else if( linkEnds.at( receiver ).bodyName_ == referenceBody )
        {
            currentPointingAnglesCalculator = bodies.at( referenceBody )->getGroundStation(
                        linkEnds.at( receiver ).stationName_ )->getPointingAnglesCalculator( );
            elevationAngles.push_back( currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
                                           ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment( 0, 3 ), linkEndTimes.at( 1 ) ) );
        }
        break;
    }
//    case one_way_differenced_range:
//    {
//        if( linkEnds.at( transmitter ).bodyName_ == referenceBody )
//        {
//            currentPointingAnglesCalculator = bodies.at( referenceBody )->getGroundStation(
//                        linkEnds.at( transmitter ).stationName_ )->getPointingAnglesCalculator( );
//            elevationAngles.push_back( currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
//                                           ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ), linkEndTimes.at( 0 ) ) );
//            elevationAngles.push_back( currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
//                                           ( linkEndStates.at( 3 ) - linkEndStates.at( 2 ) ).segment( 0, 3 ), linkEndTimes.at( 2 ) ) );
//
//        }
//        else if( linkEnds.at( receiver ).bodyName_ == referenceBody )
//        {
//            currentPointingAnglesCalculator = bodies.at( referenceBody )->getGroundStation(
//                        linkEnds.at( receiver ).stationName_ )->getPointingAnglesCalculator( );
//            elevationAngles.push_back( currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
//                                           ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment( 0, 3 ), linkEndTimes.at( 1 ) ) );
//            elevationAngles.push_back( currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
//                                           ( linkEndStates.at( 2 ) - linkEndStates.at( 3 ) ).segment( 0, 3 ), linkEndTimes.at( 3 ) ) );
//        }
//        break;
//    }
    case n_way_range:
    {
        int linkEndIndex = 0;
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
             linkEndIterator++ )
        {
            if( linkEndIterator->second.bodyName_ == referenceBody )
            {
                currentPointingAnglesCalculator = bodies.at( referenceBody )->getGroundStation(
                            linkEndIterator->second.stationName_ )->getPointingAnglesCalculator( );
                if( linkEndIndex != 0 )
                {
                    elevationAngles.push_back(
                                currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
                                    ( linkEndStates.at( 2 * ( linkEndIndex - 1 ) ) -
                                      linkEndStates.at( 2 * ( linkEndIndex - 1 ) + 1 ) ).segment( 0, 3 ),
                                    linkEndTimes.at( 2 * ( linkEndIndex - 1 ) + 1 ) ) );
                }

                if( linkEndIndex != static_cast< int >( linkEnds.size( ) ) - 1 )
                {
                    elevationAngles.push_back(
                                currentPointingAnglesCalculator->calculateElevationAngleFromInertialVector(
                                    ( linkEndStates.at( 2 * linkEndIndex + 1) -
                                      linkEndStates.at( 2 * linkEndIndex ) ).segment( 0, 3 ),
                                    linkEndTimes.at( 2 * linkEndIndex ) ) );
                }

            }
            linkEndIndex++;
        }
        break;
    }
    default:
        throw std::runtime_error( "Error when testing elevation angle viability, observable not recognized" );

    }

    return elevationAngles;
}

//! Test function to manually compute body avoidance angle(s) of observation link at given body
std::vector< double > getBodyCosineAvoidanceAngles(
        const LinkEnds linkEnds,
        const ObservableType observableType,
        const std::string referenceBody,
        const std::string bodyToAvoid,
        const std::vector< Eigen::Vector6d > linkEndStates,
        const std::vector< double > linkEndTimes,
        const SystemOfBodies& bodies )
{
    std::vector< double > cosineAvoidanceAngles;
    switch( observableType )
    {
    case one_way_range:
    case one_way_doppler:
    case angular_position:
    {
        if( linkEnds.at( transmitter ).bodyName_ == referenceBody )
        {
            cosineAvoidanceAngles.push_back(
                        linear_algebra::computeCosineOfAngleBetweenVectors(
                            ( ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ) ), (
                                bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
                                    ( linkEndTimes.at( 1 ) + linkEndTimes.at( 0 ) ) / 2.0 ).segment( 0, 3 ) -
                                linkEndStates.at( 0 ).segment( 0, 3 ) ) ) );
        }
        else if( linkEnds.at( receiver ).bodyName_ == referenceBody )
        {
            cosineAvoidanceAngles.push_back(
                        linear_algebra::computeCosineOfAngleBetweenVectors(
                            ( ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment( 0, 3 ) ), (
                                bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
                                    ( linkEndTimes.at( 1 ) + linkEndTimes.at( 0 ) ) / 2.0 ).segment( 0, 3 ) -
                                linkEndStates.at( 1 ).segment( 0, 3 ) ) ) );
        }
        break;
    }
//    case one_way_differenced_range:
//    {
//        if( linkEnds.at( transmitter ).bodyName_ == referenceBody )
//        {
//            cosineAvoidanceAngles.push_back(
//                        linear_algebra::computeCosineOfAngleBetweenVectors(
//                            ( ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ) ), (
//                                bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
//                                    ( linkEndTimes.at( 1 ) + linkEndTimes.at( 0 ) ) / 2.0 ).segment( 0, 3 ) -
//                                linkEndStates.at( 0 ).segment( 0, 3 ) ) ) );
//            cosineAvoidanceAngles.push_back(
//                        linear_algebra::computeCosineOfAngleBetweenVectors(
//                            ( ( linkEndStates.at( 3 ) - linkEndStates.at( 2 ) ).segment( 0, 3 ) ), (
//                                bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
//                                    ( linkEndTimes.at( 3 ) + linkEndTimes.at( 2 ) ) / 2.0 ).segment( 0, 3 ) -
//                                linkEndStates.at( 2 ).segment( 0, 3 ) ) ) );
//        }
//        else if( linkEnds.at( receiver ).bodyName_ == referenceBody )
//        {
//            cosineAvoidanceAngles.push_back(
//                        linear_algebra::computeCosineOfAngleBetweenVectors(
//                            ( ( linkEndStates.at( 0 ) - linkEndStates.at( 1 ) ).segment( 0, 3 ) ), (
//                                bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
//                                    ( linkEndTimes.at( 1 ) + linkEndTimes.at( 0 ) ) / 2.0 ).segment( 0, 3 ) -
//                                linkEndStates.at( 1 ).segment( 0, 3 ) ) ) );
//            cosineAvoidanceAngles.push_back(
//                        linear_algebra::computeCosineOfAngleBetweenVectors(
//                            ( ( linkEndStates.at( 2 ) - linkEndStates.at( 3 ) ).segment( 0, 3 ) ), (
//                                bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
//                                    ( linkEndTimes.at( 3 ) + linkEndTimes.at( 2 ) ) / 2.0 ).segment( 0, 3 ) -
//                                linkEndStates.at( 3 ).segment( 0, 3 ) ) ) );
//        }
//        break;
//    }
    case n_way_range:
    {
        int linkEndIndex = 0;
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
             linkEndIterator++ )
        {
            if( linkEndIterator->second.bodyName_ == referenceBody )
            {
                if( linkEndIndex != 0 )
                {
                    cosineAvoidanceAngles.push_back(
                                linear_algebra::computeCosineOfAngleBetweenVectors(
                                    ( ( linkEndStates.at( 2 * ( linkEndIndex - 1 ) ) -
                                        linkEndStates.at( 2 * ( linkEndIndex - 1 ) + 1 ) ).segment( 0, 3 ) ), (
                                        bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
                                            ( linkEndTimes.at( 2 * ( linkEndIndex - 1 ) + 1  ) +
                                              linkEndTimes.at( 2 * ( linkEndIndex - 1 ) ) ) / 2.0 ).segment( 0, 3 ) -
                                        linkEndStates.at( 2 * ( linkEndIndex - 1 ) + 1 ).segment( 0, 3 ) ) ) );
                }

                if( linkEndIndex != static_cast< int >( linkEnds.size( ) ) - 1 )
                {
                    cosineAvoidanceAngles.push_back(
                                linear_algebra::computeCosineOfAngleBetweenVectors(
                                    ( linkEndStates.at( 2 * linkEndIndex + 1 ) -
                                      linkEndStates.at( 2 * linkEndIndex ) ).segment( 0, 3 ), (
                                        bodies.at( bodyToAvoid )->getStateInBaseFrameFromEphemeris< double, double >(
                                            ( linkEndTimes.at( 2 * linkEndIndex + 1  ) +
                                              linkEndTimes.at( 2 * linkEndIndex  ) ) / 2.0 ).segment( 0, 3 ) -
                                        linkEndStates.at( 2 * linkEndIndex ).segment( 0, 3 ) ) ) );
                }

            }
            linkEndIndex++;
        }
        break;
    }
    default:
        throw std::runtime_error( "Error when testing avoidance angle viability, observable not recognized" );

    }

    return cosineAvoidanceAngles;
}

//! Test function to manually compute body distance between body center of mass and observation link vector
std::vector< double > getDistanceBetweenLineOfSightVectorAndPoint(
        const std::string bodyToAnalyze,
        const std::vector< Eigen::Vector6d > linkEndStates,
        const std::vector< double > linkEndTimes,
        const SystemOfBodies& bodies )
{
    std::vector< double > pointDistances;
    for( unsigned int i = 0; i < linkEndStates.size( ); i += 2 )
    {
        Eigen::Vector3d linkEnd1Position = linkEndStates.at( i ).segment( 0, 3 );
        Eigen::Vector3d linkEnd2Position = linkEndStates.at( i + 1 ).segment( 0, 3 );

        Eigen::Vector3d pointPosition = bodies.at( bodyToAnalyze )->getStateInBaseFrameFromEphemeris< double, double >(
                    ( linkEndTimes.at( i ) + linkEndTimes.at( i ) ) / 2.0 ).segment( 0, 3 );

        pointDistances.push_back( ( ( linkEnd2Position - pointPosition ).cross( linkEnd1Position - pointPosition ) ).norm( ) /
                                  ( linkEnd2Position - linkEnd1Position ).norm( ) );

    }
    return pointDistances;
}

//! Test if observation viability calculators correctly constrain observations
/*!
 *  Test if observation viability calculators correctly constrain observations. Test is performed for a list of observables,
 *  and a list of link ends per observable. The elevation angle, occultation and avoidance angle conditions are all checked.
 */
BOOST_AUTO_TEST_CASE( testObservationViabilityCalculators )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define environment settings
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames );

    // Set simplified rotation models for Earth/Mars (spice rotation retrieval is slow)
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Earth",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Earth", 0.0 ),
                0.0, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY ) );
    bodySettings.at( "Mars" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Mars",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Mars", 0.0 ),
                0.0, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) );

    // Set unrealistically large radius of Moon to make iss occultation mode significant
    double moonRadius = 1.0E6;
    bodySettings.at( "Moon" )->shapeModelSettings = std::make_shared< SphericalBodyShapeSettings >( moonRadius );

    // Create list of body objects
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );


    // Create ground stations
    std::pair< std::string, std::string > earthStation1 = std::pair< std::string, std::string >( "Earth", "EarthStation1" );
    std::pair< std::string, std::string > earthStation2 = std::pair< std::string, std::string >( "Earth", "EarthStation2" );
    std::pair< std::string, std::string > mslStation1 = std::pair< std::string, std::string >( "Mars", "MarsStation1" );
    std::pair< std::string, std::string > mslStation2 = std::pair< std::string, std::string >( "Mars", "MarsStation2" );
    createGroundStation( bodies.at( "Mars" ), "MarsStation1", ( Eigen::Vector3d( ) << 100.0, 0.2, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Mars" ), "MarsStation2", ( Eigen::Vector3d( ) << -2000.0, -0.4, 0.1 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "EarthStation1", ( Eigen::Vector3d( ) << 800.0, 0.12, 5.3 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "EarthStation2", ( Eigen::Vector3d( ) << 100.0, 0.15, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Define one-way range/one-way Doppler/angular position/one-way differenced range link ends
    LinkEnds oneWayLinkEnds1;
    oneWayLinkEnds1[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation1" );
    oneWayLinkEnds1[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation1" );

    LinkEnds oneWayLinkEnds2;
    oneWayLinkEnds2[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation2" );
    oneWayLinkEnds2[ receiver ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation2" );

    LinkEnds oneWayLinkEnds3;
    oneWayLinkEnds3[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation1" );
    oneWayLinkEnds3[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation2" );

    std::vector< LinkEnds > oneWayRangeLinkEnds;
    oneWayRangeLinkEnds.push_back( oneWayLinkEnds1 );
    oneWayRangeLinkEnds.push_back( oneWayLinkEnds2 );
    oneWayRangeLinkEnds.push_back( oneWayLinkEnds3 );

    // Define two-way range link ends
    LinkEnds twoWayLinkEnds1;
    twoWayLinkEnds1[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation1" );
    twoWayLinkEnds1[ reflector1 ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation1" );
    twoWayLinkEnds1[ receiver ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation1" );

    LinkEnds twoWayLinkEnds2;
    twoWayLinkEnds2[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation2" );
    twoWayLinkEnds2[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation2" );
    twoWayLinkEnds2[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation1" );

    LinkEnds twoWayLinkEnds3;
    twoWayLinkEnds3[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation2" );
    twoWayLinkEnds3[ reflector1 ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation2" );
    twoWayLinkEnds3[ receiver ] = std::make_pair< std::string, std::string >( "Mars", "MarsStation1" );
    std::vector< LinkEnds > twoWayRangeLinkEnds;
    twoWayRangeLinkEnds.push_back( twoWayLinkEnds1 );
    twoWayRangeLinkEnds.push_back( twoWayLinkEnds2 );
    twoWayRangeLinkEnds.push_back( twoWayLinkEnds3 );

    // Create list of link ends per obsevables
    std::map< ObservableType, std::vector< LinkEnds > > testLinkEndsList;
    testLinkEndsList[ one_way_range ] = oneWayRangeLinkEnds;
    testLinkEndsList[ n_way_range ] = twoWayRangeLinkEnds;

    // Define list of observation times that are to be checked: one observation every 2.5 days, over a period of 10 years.
    std::vector< double > unconstrainedObservationTimes;
    LinkEndType referenceLinkEnd = transmitter;
    double initialTime = 0.0, finalTime = 10.0 * physical_constants::JULIAN_YEAR, timeStep = physical_constants::JULIAN_DAY * 2.5;
    double currentTime = initialTime;
    while( currentTime <= finalTime )
    {
        unconstrainedObservationTimes.push_back( currentTime );
        currentTime += timeStep;
    }
    int unconstrainedNumberOfObservations = unconstrainedObservationTimes.size( );

    // Define minimum elevation angles for Earth/Mars stations
    double earthTestAngle = 4.0 * mathematical_constants::PI / 180.0;
    double marsTestAngle = 10.0 * mathematical_constants::PI / 180.0;

    // Define minimum Sun avoidance angles for Earth/Mars stations
    double earthSunAvoidanceAngle = 30.0 * mathematical_constants::PI / 180.0;
    double marsSunAvoidanceAngle = 21.0 * mathematical_constants::PI / 180.0;


    // Create observation viability settings
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                minimum_elevation_angle, std::make_pair< std::string, std::string >( "Earth", "" ), "",
                                                earthTestAngle ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                minimum_elevation_angle, std::make_pair< std::string, std::string >( "Mars", "" ), "",
                                                marsTestAngle ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_avoidance_angle, std::make_pair< std::string, std::string >( "Earth", "" ), "Sun",
                                                earthSunAvoidanceAngle ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_avoidance_angle, std::make_pair< std::string, std::string >( "Mars", "" ), "Sun",
                                                marsSunAvoidanceAngle ) );
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
                                                body_occultation, std::make_pair< std::string, std::string >( "Earth", "" ), "Moon" ) );

    for( int viabilityTest = 0; viabilityTest < 6; viabilityTest++ )
    {
        std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettingsToUse;
        if( viabilityTest < 5 )
        {
            observationViabilitySettingsToUse.push_back( observationViabilitySettings.at( viabilityTest ));
        }
        else
        {
            observationViabilitySettingsToUse = observationViabilitySettings;
        }

        // Create observation model and observation time settings for all observables
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationTimeSettings;
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationTimeSettingsConstrained;

        std::vector< std::shared_ptr< ObservationModelSettings > >  observationSettingsList;
        for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = testLinkEndsList.begin( );
             observableIterator != testLinkEndsList.end( ); observableIterator++ )
        {
            for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
            {
                LinkEnds linkEnds = observableIterator->second.at( i );
                observationSettingsList.push_back(
                            std::make_shared< ObservationModelSettings >(
                                observableIterator->first, linkEnds, std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
                observationTimeSettings.push_back(
                            std::make_shared< TabulatedObservationSimulationSettings< double > >(
                                observableIterator->first, observableIterator->second.at( i ), unconstrainedObservationTimes,
                                referenceLinkEnd ) );
                observationTimeSettingsConstrained.push_back(
                            std::make_shared< TabulatedObservationSimulationSettings< double > >(
                                observableIterator->first, observableIterator->second.at( i ), unconstrainedObservationTimes,
                                referenceLinkEnd, observationViabilitySettingsToUse ) );
            }
        }

        // Create osbervation simulatos
        std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
                createObservationSimulators( observationSettingsList, bodies );
        std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulatorsMap;
        for( unsigned int i = 0; i < observationSimulators.size( ); i++ )
        {
            observationSimulatorsMap[ observationSimulators.at( i )->getObservableType( ) ] = observationSimulators.at( i );
        }

        // Simulate observations without constraints directly from simulateObservations function
        std::shared_ptr< observation_models::ObservationCollection< > > unconstrainedSimulatedObservables =
                simulateObservations( observationTimeSettings, observationSimulators, bodies );


        // Simulate observations with viability constraints directly from simulateObservations function
        std::shared_ptr< observation_models::ObservationCollection< > > constrainedSimulatedObservables =
                simulateObservations( observationTimeSettingsConstrained, observationSimulators, bodies );

        int numberOfObservables = testLinkEndsList.size( );

        // Check consistency of simulated observations from ObservationSimulator objects/simulateObservations function
        BOOST_CHECK_EQUAL( numberOfObservables, unconstrainedSimulatedObservables->getObservationTypeStartAndSize( ).size( ) );
        BOOST_CHECK_EQUAL( numberOfObservables, constrainedSimulatedObservables->getObservationTypeStartAndSize( ).size( ) );

        // Create iterators over all simulated observations
        auto unconstrainedSortedObservations = unconstrainedSimulatedObservables->getObservationSetStartAndSize( );
        auto constrainedSortedObservations = constrainedSimulatedObservables->getObservationSetStartAndSize( );

        std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >::iterator unconstrainedIterator =
                unconstrainedSortedObservations.begin( );
        std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >::iterator constrainedIterator =
                constrainedSortedObservations.begin( );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;

        std::vector< double > unconstrainedConcatenatedTimes = unconstrainedSimulatedObservables->getConcatenatedTimeVector( );
        std::vector< double > constrainedConcatenatedTimes = constrainedSimulatedObservables->getConcatenatedTimeVector( );

        // Iterate over all observations and check viability constraints
        for( int i = 0; i < numberOfObservables; i++ )
        {
            int numberOfLinkEnds = testLinkEndsList.at( unconstrainedIterator->first ).size( );
            int currentObservableSize = getObservableSize( unconstrainedIterator->first );

            // Check consistency of simulated observations from ObservationSimulator objects/simulateObservations function
            BOOST_CHECK_EQUAL( numberOfLinkEnds, unconstrainedIterator->second.size( ) );
            BOOST_CHECK_EQUAL( numberOfLinkEnds, constrainedIterator->second.size( ) );

            // Create iterators over all simulated observations of current observable
            auto unconstrainedLinkIterator = unconstrainedIterator->second.begin( );
            auto constrainedLinkIterator = constrainedIterator->second.begin( );

            ObservableType currentObservable = unconstrainedIterator->first;

            std::shared_ptr< ObservationSimulatorBase< double, double > > currentObservationSimulator =
                    observationSimulatorsMap.at( currentObservable );

            // Iterate over all link ends of current observables.
            for( int j = 0; j < numberOfLinkEnds; j++ )
            {
                LinkEnds currentLinkEnds = unconstrainedLinkIterator->first;

                // Create manual viability calculators
                std::vector< std::shared_ptr< ObservationViabilityCalculator > > currentViabilityCalculators =
                        tudat::observation_models::createObservationViabilityCalculators(
                            bodies, currentLinkEnds, currentObservable, observationViabilitySettingsToUse );

                int unconstrainedIndex = 0;
                int constrainedIndex = 0;

                // Checks output of automatic viability calculators
                bool currentObservationWasViable = 0;

                // Checks output of manual viability calculators
                bool currentObservationIsViable = 0;
                bool isSingleViabilityConditionMet = 0;

                // Checks out[ut of manual calculation of viability
                bool computedViability = 0;

                Eigen::VectorXd currentObservation;

                // Retrieve observations for current link ends/observable type
                std::vector< std::pair< int, int > > constrainedIndices = constrainedLinkIterator->second;
                std::vector< std::pair< int, int > > unconstrainedIndices = unconstrainedLinkIterator->second;

                BOOST_CHECK_EQUAL( constrainedIndices.size( ), unconstrainedIndices.size( ) );


                for( unsigned int k = 0; k < constrainedIndices.size( ); k++ )
                {
                    // Retrieve indices from full observation vector
                    int currentUnconstrainedStartIndex = unconstrainedIndices.at( k ).first;
                    int currentUnconstrainedBlockSize = unconstrainedIndices.at( k ).second;

                    int currentConstrainedStartIndex = constrainedIndices.at( k ).first;
                    int currentConstrainedBlockSize = constrainedIndices.at( k ).second;

                    // Retrieve times from full observation vectors
                    std::vector< double > unconstrainedTimesSegment(
                                unconstrainedConcatenatedTimes.begin( ) + currentUnconstrainedStartIndex,
                                unconstrainedConcatenatedTimes.begin( ) + currentUnconstrainedStartIndex + currentUnconstrainedBlockSize );

                    std::vector< double > constrainedTimesSegment(
                                constrainedConcatenatedTimes.begin( ) + currentConstrainedStartIndex,
                                constrainedConcatenatedTimes.begin( ) + currentConstrainedStartIndex + currentConstrainedBlockSize );

                    // Iterate over all unconstrained observations for current observable/link ends
                    while( unconstrainedIndex < unconstrainedNumberOfObservations )
                    {
                        // Check if current observation was rejected by viability calculators
                        if( constrainedIndex >= static_cast< int >( constrainedTimesSegment.size( ) ) )
                        {
                            currentObservationWasViable = false;
                        }
                        else if( unconstrainedTimesSegment.at( unconstrainedIndex ) ==
                                 constrainedTimesSegment.at( constrainedIndex ) )
                        {
                            currentObservationWasViable = true;
                        }
                        else
                        {
                            currentObservationWasViable = false;
                        }

                        // Re-simulate current observation
                        if( currentObservable == angular_position )
                        {
                            currentObservation = std::dynamic_pointer_cast< ObservationSimulator< 2 > >(
                                        currentObservationSimulator )->getObservationModel( currentLinkEnds )->
                                    computeObservationsWithLinkEndData(
                                        unconstrainedObservationTimes.at( unconstrainedIndex ), referenceLinkEnd,
                                        linkEndTimes, linkEndStates );
                        }
                        else
                        {
                            currentObservation = std::dynamic_pointer_cast< ObservationSimulator< 1 > >(
                                        currentObservationSimulator )->getObservationModel( currentLinkEnds )->
                                    computeObservationsWithLinkEndData(
                                        unconstrainedObservationTimes.at( unconstrainedIndex ), referenceLinkEnd,
                                        linkEndTimes, linkEndStates );
                        }

                        // Re-compute viability according to viability calculator objects.
                        currentObservationIsViable = true;
                        for( unsigned int k = 0; k < currentViabilityCalculators.size( ); k++ )
                        {
                            isSingleViabilityConditionMet = currentViabilityCalculators.at( k )->isObservationViable(
                                        linkEndStates, linkEndTimes );
                            if( !isSingleViabilityConditionMet )
                            {
                                currentObservationIsViable = false;
                            }
                        }


                        // Manually recompute earth elevation angle condition

                        computedViability = true;
                        if( viabilityTest == 0 || viabilityTest == 5 )
                        {
                            std::vector< double > earthElevationAngles = getBodyLinkElevationAngles(
                                currentLinkEnds, currentObservable, "Earth", linkEndStates, linkEndTimes, bodies );
                            for ( unsigned int l = 0; l < earthElevationAngles.size( ); l++ )
                            {
                                if ( earthElevationAngles.at( l ) < earthTestAngle )
                                {
                                    computedViability = false;
                                }
                            }
                        }

                        // Manually recompute mars elevation angle condition
                        if( viabilityTest == 1 || viabilityTest == 5 )
                        {
                            std::vector<double> marsElevationAngles = getBodyLinkElevationAngles(
                                currentLinkEnds, currentObservable, "Mars", linkEndStates, linkEndTimes, bodies );
                            for ( unsigned int l = 0; l < marsElevationAngles.size( ); l++ )
                            {
                                if ( marsElevationAngles.at( l ) < marsTestAngle )
                                {
                                    computedViability = false;
                                }
                            }
                        }

                        if( viabilityTest == 2 || viabilityTest == 5 )
                        {
                            // Manually recompute earth-Sun avoidance angle condition
                            std::vector< double > earthSunCosineAvoidanceAngles = getBodyCosineAvoidanceAngles(
                                        currentLinkEnds, currentObservable, "Earth", "Sun", linkEndStates, linkEndTimes, bodies );
                            for( unsigned int l = 0; l < earthSunCosineAvoidanceAngles.size( ); l++ )
                            {
                                if( earthSunCosineAvoidanceAngles.at( l ) > std::cos( earthSunAvoidanceAngle ) )
                                {
                                    computedViability =  false;
                                }
                            }
                        }

                        if( viabilityTest == 3 || viabilityTest == 5 )
                        {
                            // Manually recompute mars-Sun avoidance angle condition
                            std::vector<double> marsSunCosineAvoidanceAngles = getBodyCosineAvoidanceAngles(
                                currentLinkEnds, currentObservable, "Mars", "Sun", linkEndStates, linkEndTimes,
                                bodies );
                            for ( unsigned int l = 0; l < marsSunCosineAvoidanceAngles.size( ); l++ )
                            {
                                if ( marsSunCosineAvoidanceAngles.at( l ) > std::cos( marsSunAvoidanceAngle ))
                                {
                                    computedViability = false;
                                }
                            }
                        }

                        if( viabilityTest == 4 || viabilityTest == 5 )
                        {
                            // Manually recompute occultataion consition
                            std::vector<double> moonLineOfSightDistances = getDistanceBetweenLineOfSightVectorAndPoint(
                                "Moon", linkEndStates, linkEndTimes, bodies );
                            for ( unsigned int l = 0; l < moonLineOfSightDistances.size( ); l++ )
                            {
                                //                        std::cout<<moonLineOfSightDistances.at( l )<<" "<<moonRadius<<std::endl;
                                if ( moonLineOfSightDistances.at( l ) < moonRadius )
                                {
                                    computedViability = false;
                                }
                            }
                        }

    //                    std::cout<<currentObservationIsViable<<std::endl<<std::endl;

                        // Check manual/automatic viability check
                        BOOST_CHECK_EQUAL( currentObservationIsViable, currentObservationWasViable );
                        BOOST_CHECK_EQUAL( computedViability, currentObservationWasViable );

                        if( currentObservationWasViable )
                        {
                            constrainedIndex++;
                        }
                        unconstrainedIndex++;
                    }
    //                BOOST_CHECK_EQUAL( constrainedIndex, constrainedLinkIterator->second.size( ) );

                }
                unconstrainedLinkIterator++;
                constrainedLinkIterator++;
            }

            unconstrainedIterator++;
            constrainedIterator++;
        }
    }
}

Eigen::Vector6d customStateFunction( const double time )
{
    Eigen::Vector6d currentState = Eigen::Vector6d::Zero( );
    currentState( 0 ) = ( 5.0E6 + 10.0 ) * std::sin( time / 3600.0 );
    currentState( 1 ) = ( 5.0E6 + 10.0 ) * std::cos( time / 3600.0 );
    return currentState;

}

BOOST_AUTO_TEST_CASE( testOrbiterOccultationObservationViabilityCalculators )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define environment settings
    std::vector<std::string> bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Jupiter" );

    BodyListSettings bodySettings =
        getDefaultBodySettings( bodyNames, "SSB", "ECLIPJ2000" );
    bodySettings.addSettings( "Spacecraft" );
    bodySettings.at( "Spacecraft" )->ephemerisSettings = customEphemerisSettings(  customStateFunction, "Earth", "ECLIPJ2000" );

    // Set unrealistically large radius of Moon to make iss occultation mode significant
    double earthRadius = 5.0E6;
    bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared<ApproximateJplEphemerisSettings>( "Earth", true );
    bodySettings.at( "Jupiter" )->ephemerisSettings = std::make_shared<ApproximateJplEphemerisSettings>( "Jupiter", true );
    bodySettings.at( "Sun" )->ephemerisSettings = std::make_shared<ConstantEphemerisSettings>( Eigen::Vector6d::Zero( ) );

    bodySettings.at( "Earth" )->shapeModelSettings = std::make_shared<SphericalBodyShapeSettings>( earthRadius );

    // Create list of body objects
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );


    // Define one-way range/one-way Doppler/angular position/one-way differenced range link ends
    LinkEnds oneWayLinkEnds1;
    oneWayLinkEnds1[ transmitter ] = std::make_pair<std::string, std::string>( "Jupiter", "" );
    oneWayLinkEnds1[ receiver ] = std::make_pair<std::string, std::string>( "Spacecraft", "" );

    LinkEnds oneWayLinkEnds2;
    oneWayLinkEnds2[ receiver ] = std::make_pair<std::string, std::string>( "Jupiter", "" );
    oneWayLinkEnds2[ transmitter ] = std::make_pair<std::string, std::string>( "Spacecraft", "" );

    LinkEnds twoWayLinkEnds1;
    twoWayLinkEnds1[ receiver ] = std::make_pair<std::string, std::string>( "Jupiter", "" );
    twoWayLinkEnds1[ retransmitter ] = std::make_pair<std::string, std::string>( "Spacecraft", "" );
    twoWayLinkEnds1[ transmitter ] = std::make_pair<std::string, std::string>( "Jupiter", "" );

    LinkEnds twoWayLinkEnds2;
    twoWayLinkEnds2[ receiver ] = std::make_pair<std::string, std::string>( "Spacecraft", "" );
    twoWayLinkEnds2[ retransmitter ] = std::make_pair<std::string, std::string>( "Jupiter", "" );
    twoWayLinkEnds2[ transmitter ] = std::make_pair<std::string, std::string>( "Spacecraft", "" );

    std::vector<LinkEnds> oneWayRangeLinkEnds;
    oneWayRangeLinkEnds.push_back( oneWayLinkEnds1 );
    oneWayRangeLinkEnds.push_back( oneWayLinkEnds2 );

    std::vector<LinkEnds> twoWayRangeLinkEnds;
    twoWayRangeLinkEnds.push_back( twoWayLinkEnds1 );
    twoWayRangeLinkEnds.push_back( twoWayLinkEnds2 );

    // Create list of link ends per obsevables
    std::map< ObservableType, std::vector< LinkEnds > > testLinkEndsList;
    testLinkEndsList[ one_way_range ] = oneWayRangeLinkEnds;
    testLinkEndsList[ n_way_range ] = twoWayRangeLinkEnds;

    // Define list of observation times that are to be checked: one observation every 2.5 days, over a period of 10 years.
    std::vector< double > unconstrainedObservationTimes;
    LinkEndType referenceLinkEnd = receiver;
    double initialTime = 0.0, finalTime = 10.0 * physical_constants::JULIAN_DAY, timeStep = physical_constants::JULIAN_DAY * 0.1;
    double currentTime = initialTime;
    while( currentTime <= finalTime )
    {
        unconstrainedObservationTimes.push_back( currentTime );
        currentTime += timeStep;
    }
    int unconstrainedNumberOfObservations = unconstrainedObservationTimes.size( );

    // Create observation viability settings
    std::vector< std::shared_ptr< ObservationViabilitySettings > > observationViabilitySettings;
    observationViabilitySettings.push_back( std::make_shared< ObservationViabilitySettings >(
        body_occultation, std::make_pair< std::string, std::string >( "Spacecraft", "" ), "Earth" ) );

    // Create observation model and observation time settings for all observables
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationTimeSettings;
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationTimeSettingsConstrained;

    std::vector< std::shared_ptr< ObservationModelSettings > >  observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = testLinkEndsList.begin( );
         observableIterator != testLinkEndsList.end( ); observableIterator++ )
    {
        for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
        {
            LinkEnds linkEnds = observableIterator->second.at( i );
            observationSettingsList.push_back(
                std::make_shared< ObservationModelSettings >(
                    observableIterator->first, linkEnds, std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
            observationTimeSettings.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< double > >(
                    observableIterator->first, observableIterator->second.at( i ), unconstrainedObservationTimes,
                    referenceLinkEnd ) );
            observationTimeSettingsConstrained.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< double > >(
                    observableIterator->first, observableIterator->second.at( i ), unconstrainedObservationTimes,
                    referenceLinkEnd, observationViabilitySettings ) );
        }
    }

    // Create osbervation simulatos
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
        createObservationSimulators( observationSettingsList, bodies );
    std::map< ObservableType, std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulatorsMap;
    for( unsigned int i = 0; i < observationSimulators.size( ); i++ )
    {
        observationSimulatorsMap[ observationSimulators.at( i )->getObservableType( ) ] = observationSimulators.at( i );
    }

    // Simulate observations without constraints directly from simulateObservations function
    std::shared_ptr< observation_models::ObservationCollection< > > unconstrainedSimulatedObservables =
        simulateObservations( observationTimeSettings, observationSimulators, bodies );


    // Simulate observations with viability constraints directly from simulateObservations function
    std::shared_ptr< observation_models::ObservationCollection< > > constrainedSimulatedObservables =
        simulateObservations( observationTimeSettingsConstrained, observationSimulators, bodies );

    int numberOfObservables = testLinkEndsList.size( );

    // Create iterators over all simulated observations
    auto unconstrainedSortedObservations = unconstrainedSimulatedObservables->getObservationSetStartAndSize( );
    auto constrainedSortedObservations = constrainedSimulatedObservables->getObservationSetStartAndSize( );

    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >::iterator unconstrainedIterator =
        unconstrainedSortedObservations.begin( );
    std::map< ObservableType, std::map< LinkEnds, std::vector< std::pair< int, int > > > >::iterator constrainedIterator =
        constrainedSortedObservations.begin( );

    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;

    std::vector< double > unconstrainedConcatenatedTimes = unconstrainedSimulatedObservables->getConcatenatedTimeVector( );
    std::vector< double > constrainedConcatenatedTimes = constrainedSimulatedObservables->getConcatenatedTimeVector( );

    // Iterate over all observations and check viability constraints
    for( int i = 0; i < numberOfObservables; i++ )
    {
        int numberOfLinkEnds = testLinkEndsList.at( unconstrainedIterator->first ).size( );
        int currentObservableSize = getObservableSize( unconstrainedIterator->first );

        // Create iterators over all simulated observations of current observable
        auto unconstrainedLinkIterator = unconstrainedIterator->second.begin( );
        auto constrainedLinkIterator = constrainedIterator->second.begin( );

        ObservableType currentObservable = unconstrainedIterator->first;

        std::shared_ptr<ObservationSimulatorBase<double, double> > currentObservationSimulator =
            observationSimulatorsMap.at( currentObservable );

        // Iterate over all link ends of current observables.
        for ( int j = 0; j < numberOfLinkEnds; j++ )
        {
            LinkEnds currentLinkEnds = unconstrainedLinkIterator->first;

            int unconstrainedIndex = 0;
            int constrainedIndex = 0;

            Eigen::VectorXd currentObservation;

            // Retrieve observations for current link ends/observable type
            std::vector<std::pair<int, int> > constrainedIndices = constrainedLinkIterator->second;
            std::vector<std::pair<int, int> > unconstrainedIndices = unconstrainedLinkIterator->second;

            for ( unsigned int k = 0; k < constrainedIndices.size( ); k++ )
            {
                // Retrieve indices from full observation vector
                int currentUnconstrainedStartIndex = unconstrainedIndices.at( k ).first;
                int currentUnconstrainedBlockSize = unconstrainedIndices.at( k ).second;

                int currentConstrainedStartIndex = constrainedIndices.at( k ).first;
                int currentConstrainedBlockSize = constrainedIndices.at( k ).second;

                // Retrieve times from full observation vectors
                std::vector<double> unconstrainedTimesSegment(
                    unconstrainedConcatenatedTimes.begin( ) + currentUnconstrainedStartIndex,
                    unconstrainedConcatenatedTimes.begin( ) + currentUnconstrainedStartIndex +
                    currentUnconstrainedBlockSize );

                std::vector<double> constrainedTimesSegment(
                    constrainedConcatenatedTimes.begin( ) + currentConstrainedStartIndex,
                    constrainedConcatenatedTimes.begin( ) + currentConstrainedStartIndex +
                    currentConstrainedBlockSize );

                bool currentObservationWasViable;
                // Iterate over all unconstrained observations for current observable/link ends
                while ( unconstrainedIndex < unconstrainedNumberOfObservations )
                {
                    // Check if current observation was rejected by viability calculators
                    if ( constrainedIndex >= static_cast< int >( constrainedTimesSegment.size( )))
                    {
                        currentObservationWasViable = false;
                    }
                    else if ( unconstrainedTimesSegment.at( unconstrainedIndex ) ==
                              constrainedTimesSegment.at( constrainedIndex ))
                    {
                        currentObservationWasViable = true;
                    }
                    else
                    {
                        currentObservationWasViable = false;
                    }

                    // Re-simulate current observation
                    currentObservation = std::dynamic_pointer_cast<ObservationSimulator<1> >(
                        currentObservationSimulator )->getObservationModel( currentLinkEnds )->
                        computeObservationsWithLinkEndData(
                        unconstrainedObservationTimes.at( unconstrainedIndex ), referenceLinkEnd,
                        linkEndTimes, linkEndStates );

                    int numberOfLinksToCheck = currentLinkEnds.size( ) - 1;
                    bool currentObservationIsViable = true;
                    for( int linksToCheck = 0; linksToCheck < numberOfLinksToCheck; linksToCheck++ )
                    {
                        int jupiterIndex, spacecraftIndex;
                        if( linksToCheck == 0 )
                        {
                            if ( currentLinkEnds[ transmitter ].bodyName_ == "Spacecraft" )
                            {
                                jupiterIndex = 1;
                                spacecraftIndex = 0;
                            }
                            else if ( currentLinkEnds[ transmitter ].bodyName_ == "Jupiter" )
                            {
                                jupiterIndex = 0;
                                spacecraftIndex = 1;
                            }
                            else
                            {
                                throw std::runtime_error( "Error, did not find correct link end in viability settings test" );
                            }
                        }
                        else
                        {
                            if ( currentLinkEnds[ receiver ].bodyName_ == "Spacecraft" )
                            {
                                jupiterIndex = 2;
                                spacecraftIndex = 3;
                            }
                            else if ( currentLinkEnds[ receiver ].bodyName_ == "Jupiter" )
                            {
                                jupiterIndex = 3;
                                spacecraftIndex = 2;
                            }
                            else
                            {
                                throw std::runtime_error( "Error, did not find correct link end in 2-way viability settings test" );
                            }
                        }


                        Eigen::Vector3d currentJupiterPosition =
                            bodies.getBody( "Jupiter" )->getStateInBaseFrameFromEphemeris(
                                linkEndTimes[ jupiterIndex ] ).segment( 0, 3 ) -
                            bodies.getBody( "Earth" )->getStateInBaseFrameFromEphemeris(
                                linkEndTimes[ jupiterIndex ] ).segment( 0, 3 );
                        Eigen::Vector3d currentSpacecraftPosition =
                            bodies.getBody( "Spacecraft" )->getStateInBaseFrameFromEphemeris(
                                linkEndTimes[ spacecraftIndex ] ).segment( 0, 3 ) -
                            bodies.getBody( "Earth" )->getStateInBaseFrameFromEphemeris(
                                linkEndTimes[ spacecraftIndex ] ).segment( 0, 3 );

                        double jupiterAngle = std::atan2( currentJupiterPosition( 1 ), currentJupiterPosition( 0 ));
                        Eigen::Matrix3d rotationMatrix =
                            Eigen::Matrix3d( Eigen::AngleAxisd( -jupiterAngle, Eigen::Vector3d::UnitZ( )));

                        Eigen::Vector3d rotatedJupiter = rotationMatrix * currentJupiterPosition;
                        Eigen::Vector3d rotatedSpacecraft = rotationMatrix * currentSpacecraftPosition;

                        BOOST_CHECK_CLOSE_FRACTION( rotatedJupiter( 0 ), currentJupiterPosition.norm( ),
                                                    100.0 * std::numeric_limits<double>::epsilon( ));
                        BOOST_CHECK_SMALL( std::fabs( rotatedJupiter( 1 )), 1.0E-3 );
                        BOOST_CHECK_SMALL( std::fabs( rotatedJupiter( 2 )), 1.0E-3 );


                        if ( rotatedSpacecraft( 0 ) < 0 )
                        {
                            currentObservationIsViable = false;
                        }
                    }

                    BOOST_CHECK_EQUAL( currentObservationIsViable, currentObservationWasViable );

                    if( currentObservationIsViable != currentObservationWasViable )
                    {
                        std::cout<<currentObservable<<" "<<getLinkEndsString( currentLinkEnds )<<std::endl;
                    }
                    if( currentObservationWasViable )
                    {
                        constrainedIndex++;
                    }
                    unconstrainedIndex++;
                }

            }
            unconstrainedLinkIterator++;
            constrainedLinkIterator++;
        }

        unconstrainedIterator++;
        constrainedIterator++;
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}

