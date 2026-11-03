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

#include <boost/test/included/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/astro/basic_astro/dateTime.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::basic_astrodynamics;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_default_ground_station_settings )

//! DSS-65 was built in 1987 at one location, and moved to its current location on July 3, 2005. Test that the
//! default DSN station settings reproduce this discontinuous displacement, on top of the (non-zero) linear
//! station velocity that is also part of the default DSN station motion model.
BOOST_AUTO_TEST_CASE( test_Dss65PreAndPostMovePosition )
{
    // Station position specifications (ITRF2014), as published for DSS-65 in
    // https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/stations/a_old_versions/earthstns_itrf93_050714.cmt
    const Eigen::Vector3d preMoveStationPosition( 4849336.6176, -360488.6349, 4114748.9218 );
    const Eigen::Vector3d postMoveStationPositionSpecification( 4849339.6448, -360427.6560, 4114750.7428 );

    std::shared_ptr< GroundStationSettings > dss65Settings = getDsnStationSetting( "DSS-65" );

    // Retrieve the linear velocity model that the default DSN station settings add on top of the DSS-65
    // pre-/post-move discontinuity, so that the expected computed positions below can account for it.
    std::shared_ptr< LinearGroundStationMotionSettings > linearMotionSettings;
    for( auto const& motionSettings : dss65Settings->getStationMotionSettings( ) )
    {
        linearMotionSettings = std::dynamic_pointer_cast< LinearGroundStationMotionSettings >( motionSettings );
        if( linearMotionSettings != nullptr )
        {
            break;
        }
    }
    BOOST_REQUIRE( linearMotionSettings != nullptr );
    BOOST_CHECK( linearMotionSettings->linearVelocity_.norm( ) > 0.0 );

    auto getExpectedLinearDrift = [ & ]( const double epoch ) {
        return linearMotionSettings->linearVelocity_ * ( epoch - linearMotionSettings->referenceEpoch_ );
    };

    std::shared_ptr< Body > earth = std::make_shared< Body >( );
    SystemOfBodies bodies;
    bodies.addBody( earth, "Earth" );
    createGroundStation( earth, dss65Settings );

    std::shared_ptr< ground_stations::GroundStationState > dss65State = earth->getGroundStation( "DSS-65" )->getNominalStationState( );

    // Before the move, the station should be at its (hardcoded) pre-2005 position, plus the linear drift due to
    // station velocity accumulated since the settings' reference epoch.
    double preMoveTestEpoch = DateTime( 2000, 1, 1, 0, 0, 0.0 ).epoch< double >( );
    Eigen::Vector3d preMoveComputedPosition = dss65State->getCartesianPositionInTime( preMoveTestEpoch, "Earth" );
    Eigen::Vector3d preMoveExpectedPosition = preMoveStationPosition + getExpectedLinearDrift( preMoveTestEpoch );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( preMoveComputedPosition, preMoveExpectedPosition, 1.0E-10 );

    // Exactly at, and after, the move date, the piecewise-constant displacement should vanish, leaving the
    // nominal station position used to construct the settings object, again offset by the linear drift.
    Eigen::Vector3d nominalPosition = dss65Settings->getGroundStationPosition( );

    double atMoveEpoch = DateTime( 2005, 7, 3, 0, 0, 0.0 ).epoch< double >( );
    Eigen::Vector3d atMoveComputedPosition = dss65State->getCartesianPositionInTime( atMoveEpoch, "Earth" );
    Eigen::Vector3d atMoveExpectedPosition = nominalPosition + getExpectedLinearDrift( atMoveEpoch );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( atMoveComputedPosition, atMoveExpectedPosition, 1.0E-10 );

    double postMoveTestEpoch = DateTime( 2005, 7, 4, 0, 0, 0 ).epoch< double >( );
    Eigen::Vector3d postMoveComputedPosition = dss65State->getCartesianPositionInTime( postMoveTestEpoch, "Earth" );
    Eigen::Vector3d postMoveExpectedPosition = nominalPosition + getExpectedLinearDrift( postMoveTestEpoch );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( postMoveComputedPosition, postMoveExpectedPosition, 1.0E-10 );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
