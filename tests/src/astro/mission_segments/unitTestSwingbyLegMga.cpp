/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Musegaas, P. (2012). Optimization of Space Trajectories Including Multiple Gravity Assists
 *          and Deep Space Maneuvers. MSc Thesis, Delft University of Technology, Delft,
 *          The Netherlands.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <tudat/astro/basic_astro/physicalConstants.h>
#include <tudat/basics/testMacros.h>

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/mission_segments/transferLeg.h"
#include "tudat/astro/mission_segments/transferNode.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of the swing-by leg within MGA trajectory model
BOOST_AUTO_TEST_SUITE( test_swingby_leg_mga )

//! Test delta-V computation
BOOST_AUTO_TEST_CASE( testVelocities )
{
    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
    // achieved.
    const double tolerance = 1.0e-6;

    // Expected test result based on the second leg of the ideal Cassini 1 trajectory as modelled
    // by GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 1090.64622926316;
    const Eigen::Vector3d expectedArrivalVelocity (
                37952.8844553685, -14096.9656774702, -5753.51245833761 );

    // Specify the required parameters.
    // Set the planetary positions and velocities.

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector6d planet1State =
            ( Eigen::Vector6d( ) <<
              -35554348960.8278, -102574987127.178, 648696819.780156,
              32851.224953746, -11618.7310059974, -2055.04615890989 ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris1 =
            std::make_shared< ephemerides::ConstantEphemeris >( planet1State );

    const Eigen::Vector6d planet2State =
            ( Eigen::Vector6d( ) <<
              -35568329915.7073, -102569794949.529, 650816245.825226,
              TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris2 =
            std::make_shared< ephemerides::ConstantEphemeris >( planet2State );


    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 449.385873819743 * physical_constants::JULIAN_DAY;

    // Set the gravitational parameters
    const double sunGravitationalParameter = 1.32712428e20;


    using namespace mission_segments;
    tudat::mission_segments::UnpoweredUnperturbedTransferLeg transferLeg(
                constantEphemeris1, constantEphemeris2,
                sunGravitationalParameter );
    transferLeg.updateLegParameters( ( Eigen::VectorXd( 2 )<<0.0, timeOfFlight ).finished( ) );

    // Set velocity before and after swingby body
    Eigen::Vector3d velocityBeforePlanet1 (
                34216.4827530912, -15170.1440677825, 395.792122152361 );
    Eigen::Vector3d velocityAfterPlanet1 = transferLeg.getDepartureVelocity( );

    // Set the planet gravitational parameters
    const double planet1GravitationalParameter = 3.24860e14;

    //set the minimum pericenter radius
    const double minimumRadiusPlanet1 = 6351800;

    SwingbyWithFixedIncomingFixedOutgoingVelocity transferNode(
                constantEphemeris1,
                planet1GravitationalParameter, minimumRadiusPlanet1,
                [=]( ){ return velocityBeforePlanet1; },
                [=]( ){ return velocityAfterPlanet1; } );
    transferNode.updateNodeParameters( ( Eigen::VectorXd( 1 ) << 0.0 ).finished( ) );

    BOOST_CHECK_CLOSE_FRACTION( transferLeg.getLegDeltaV( ), 0.0, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( transferNode.getNodeDeltaV( ), expectedDeltaV, tolerance );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transferLeg.getArrivalVelocity( ), expectedArrivalVelocity, tolerance );

    // Get data on 10 equispace points on trajectory
    std::map< double, Eigen::Vector6d > statesAlongTrajectory;
    transferLeg.getStatesAlongTrajectory( statesAlongTrajectory, 10 );

    // Check initial and final time on output list
    BOOST_CHECK_SMALL( statesAlongTrajectory.begin( )->first, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.rbegin( )->first, timeOfFlight, 1.0E-14 );

    // Check if Keplerian state (slow elements) is the same for each output point
    Eigen::Vector6d previousKeplerianState = Eigen::Vector6d::Constant( TUDAT_NAN );
    for( auto it : statesAlongTrajectory )
    {
        Eigen::Vector6d currentCartesianState = it.second;
        Eigen::Vector6d currentKeplerianState = tudat::orbital_element_conversions::convertCartesianToKeplerianElements(
                    currentCartesianState, sunGravitationalParameter );
        if( previousKeplerianState == previousKeplerianState )
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
                        ( previousKeplerianState.segment( 0, 5 ) ),
                        ( currentKeplerianState.segment( 0, 5 ) ),
                        1.0E-14 );

        }
        previousKeplerianState = currentKeplerianState;
    }

    // Check if output meets boundary conditions
    for( int i = 0; i < 3; i++ )
    {
        //TODO: Find out why tolerance needs to be so big for one of the legs
        BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.begin( )->second( i ) - planet1State( i ) ), 20.0E3 );
        BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.rbegin( )->second( i ) - planet2State( i ) ), 20.0E3 );
    }


    // Get data on 10 equispace points on trajectory
    std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
    transferLeg.getThrustAccelerationsAlongTrajectory( thrustAccelerationsAlongTrajectory, 10 );

    // Check initial and final time on output list
    BOOST_CHECK_SMALL( thrustAccelerationsAlongTrajectory.begin( )->first, 1.0E-14 );
    BOOST_CHECK_CLOSE_FRACTION( thrustAccelerationsAlongTrajectory.rbegin( )->first, timeOfFlight, 1.0E-14 );

    // Check if thrust acceleration is zero
    for( auto it : thrustAccelerationsAlongTrajectory )
    {
        Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
        BOOST_CHECK_SMALL( currentCartesianThrustAcceleration.norm(), 1.0E-14 );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

