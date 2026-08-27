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

#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/ground_stations/pointingAnglesCalculator.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"

using namespace tudat::simulation_setup;
using namespace tudat::basic_astrodynamics;
using namespace tudat::ground_stations;
using namespace tudat::spice_interface;
using namespace tudat::unit_conversions;

namespace tudat
{
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_ground_station_sun_visibility )

BOOST_AUTO_TEST_CASE( test_isItDarkAtGroundStation )
{
    loadStandardSpiceKernels( );

    const std::string globalFrameOrientation = "J2000";
    const std::string globalFrameOrigin = "Earth";

    // Create bodies, with an ephemeris/rotation model covering the tested epochs.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun" };
    BodyListSettings bodySettings = getDefaultBodySettings( bodiesToCreate, globalFrameOrigin, globalFrameOrientation );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create a ground station on the equator.
    const std::string stationName = "TestStation";
    createGroundStation( bodies.at( "Earth" ),
                         stationName,
                         ( Eigen::Vector3d( ) << 0.0, 0.0, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    const double testTime = 1.0E7;

    // Independently compute the Sun elevation angle at the station, using the same low-level building blocks that
    // isItDarkAtGroundStation is expected to use internally.
    std::shared_ptr< PointingAnglesCalculator > pointingAnglesCalculator =
            bodies.at( "Earth" )->getGroundStation( stationName )->getPointingAnglesCalculator( );
    std::function< Eigen::Vector6d( const double ) > groundStationStateFunction =
            getLinkEndCompleteEphemerisFunction( observation_models::LinkEndId( "Earth", stationName ), bodies );
    Eigen::Vector3d groundStationPosition = groundStationStateFunction( testTime ).segment( 0, 3 );
    Eigen::Vector3d sunPosition = bodies.at( "Sun" )->getStateInBaseFrameFromEphemeris< double, double >( testTime ).segment( 0, 3 );
    double referenceSunElevation =
            pointingAnglesCalculator->calculateElevationAngleFromInertialVector( sunPosition - groundStationPosition, testTime );

    // A threshold well above the reference elevation should be reported as dark; well below, as not dark.
    BOOST_CHECK_EQUAL( isItDarkAtGroundStation( stationName, testTime, referenceSunElevation + convertDegreesToRadians( 1.0 ), bodies ),
                       true );
    BOOST_CHECK_EQUAL( isItDarkAtGroundStation( stationName, testTime, referenceSunElevation - convertDegreesToRadians( 1.0 ), bodies ),
                       false );

    // Requesting a station that does not exist should throw.
    BOOST_CHECK_THROW( isItDarkAtGroundStation( "NonExistentStation", testTime, 0.0, bodies ), std::runtime_error );

    // The generic, state-based visibility check that isItDarkAtGroundStation is built on should agree with the reference
    // elevation: the Sun should be reported "visible" at or below its own elevation, and not visible just above it.
    BOOST_CHECK_EQUAL( isTargetVisibleFromGroundStation(
                               groundStationPosition, sunPosition, testTime, referenceSunElevation, pointingAnglesCalculator ),
                       true );
    BOOST_CHECK_EQUAL( isTargetVisibleFromGroundStation( groundStationPosition,
                                                         sunPosition,
                                                         testTime,
                                                         referenceSunElevation + convertDegreesToRadians( 1.0 ),
                                                         pointingAnglesCalculator ),
                       false );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
