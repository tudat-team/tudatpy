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

#include <cmath>
#include <limits>

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

//! Test implementation of departure leg within MGA-1DSM velocity formulation trajectory model.
BOOST_AUTO_TEST_SUITE( test_departure_leg_mga_1dsm_velocity )

//! Test delta-V computation.
BOOST_AUTO_TEST_CASE( testVelocities )
{
    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
    // achieved.
    const double tolerance = 1.0e-6;

    // Expected test result based on the first leg of the ideal Messenger trajectory as modelled by
    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedEscapeDeltaV = 1408.99421278;
    const double expectedDsmDeltaV = 910.801673341206;
    const Eigen::Vector3d expectedArrivalVelocity (
                17969.3166254715, -23543.6915939142, 6.38384671663485 );

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector6d planet1State =
            ( Eigen::Vector6d( ) <<
              -148689402143.081, 7454242895.97607, 0.0,
              -1976.48781307596, -29863.4035321021, 0.0 ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris1 =
            std::make_shared< ephemerides::ConstantEphemeris >( planet1State );

    const Eigen::Vector6d planet2State =
            ( Eigen::Vector6d( ) <<
              -128359548637.032, -78282803797.7343, 0.,
              TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris2 =
            std::make_shared< ephemerides::ConstantEphemeris >( planet2State );

    // Set the Earth gravitational parameter.
    const double earthGravitationalParameter = 3.9860119e14;

    // Set the departure parking orbit (at inifinity)
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    const double excessVelocityMagnitude = 1408.99421278;
    const double excessVelocityInPlaneAngle = 0.37992647165 * 2 * 3.14159265358979;
    const double excessVelocityOutOfPlaneAngle = std::acos(  2 * 0.498004040298 - 1 ) -
                                                 3.14159265358979 / 2;
    Eigen::VectorXd nodeFreeParameters =
            ( Eigen::VectorXd( 4 ) << 0.0, excessVelocityMagnitude,
              excessVelocityInPlaneAngle, excessVelocityOutOfPlaneAngle ).finished( );

    using namespace mission_segments;
    DepartureWithFreeOutgoingVelocityNode transferNode(
            constantEphemeris1,
            earthGravitationalParameter, semiMajorAxis, eccentricity );
    transferNode.updateNodeParameters( nodeFreeParameters );


    // Set the Sun gravitational parameter.
    const double sunGravitationalParameter = 1.32712428e20;

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 399.999999715 * physical_constants::JULIAN_DAY;

    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.234594654679;
    Eigen::VectorXd legFreeParameters =
            ( Eigen::VectorXd( 3 ) << 0.0, timeOfFlight, dsmTimeOfFlightFraction ).finished( );

    Eigen::Vector3d departureVelocity = transferNode.getOutgoingVelocity( );
    DsmVelocityBasedTransferLeg transferLeg(
            constantEphemeris1, constantEphemeris2,
                sunGravitationalParameter, [=]( ){ return departureVelocity; } );
    transferLeg.updateLegParameters( legFreeParameters );

    BOOST_CHECK_CLOSE_FRACTION( transferNode.getNodeDeltaV( ), expectedEscapeDeltaV, tolerance );
    BOOST_CHECK_CLOSE_FRACTION( transferLeg.getLegDeltaV( ), expectedDsmDeltaV, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transferLeg.getArrivalVelocity( ), expectedArrivalVelocity, tolerance );

    double dsmTime = transferLeg.getDsmTime( );
    double timeTolerance = 1.0E-9;
    Eigen::Vector3d dsmLocation = transferLeg.getDsmLocation( );

    // Test state output from departure to DSM
    {
        // Get data on 10 equispace points on trajectory
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        std::vector< double > outputTimes = utilities::linspace( 0.0, transferLeg.getDsmTime( ) - timeTolerance, 10 );
        transferLeg.getStatesAlongTrajectory( statesAlongTrajectory, outputTimes );

        // Check initial and final time on output list
        BOOST_CHECK_SMALL( statesAlongTrajectory.begin( )->first, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.rbegin( )->first, dsmTime - timeTolerance, 1.0E-14 );

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
            BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.begin( )->second( i ) - planet1State( i ) ), 1.0E-2 );
            BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.rbegin( )->second( i ) - dsmLocation( i ) ), 5.0E-2 );
        }


        // Get data on 10 equispaced points on trajectory
        std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
        transferLeg.getThrustAccelerationsAlongTrajectory( thrustAccelerationsAlongTrajectory, outputTimes );

        // Check initial and final time on output list
        BOOST_CHECK_SMALL( thrustAccelerationsAlongTrajectory.begin( )->first, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( thrustAccelerationsAlongTrajectory.rbegin( )->first, dsmTime - timeTolerance, 1.0E-14 );

        // Check if thrust acceleration is zero
        for( auto it : thrustAccelerationsAlongTrajectory )
        {
            Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
            BOOST_CHECK_SMALL( currentCartesianThrustAcceleration.norm(), 1.0E-14 );
        }

    }

    // Test state output from DSM to arrival
    {
        // Get data on 10 equispace points on trajectory
        std::map< double, Eigen::Vector6d > statesAlongTrajectory;
        std::vector< double > outputTimes = utilities::linspace( transferLeg.getDsmTime( ) + timeTolerance, timeOfFlight, 10 );
        transferLeg.getStatesAlongTrajectory( statesAlongTrajectory, outputTimes );

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION( statesAlongTrajectory.begin( )->first, dsmTime + timeTolerance, 1.0E-14 );
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
            BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.begin( )->second( i ) - dsmLocation( i ) ), 5.0E-2 );
            BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.rbegin( )->second( i ) - planet2State( i ) ), 1.0E-2 );
        }


        // Get data on 10 equispaced points on trajectory
        std::map< double, Eigen::Vector3d > thrustAccelerationsAlongTrajectory;
        transferLeg.getThrustAccelerationsAlongTrajectory( thrustAccelerationsAlongTrajectory, outputTimes );

        // Check initial and final time on output list
        BOOST_CHECK_CLOSE_FRACTION( thrustAccelerationsAlongTrajectory.begin( )->first, dsmTime + timeTolerance, 1.0E-14 );
        BOOST_CHECK_CLOSE_FRACTION( thrustAccelerationsAlongTrajectory.rbegin( )->first, timeOfFlight, 1.0E-14 );

        // Check if thrust acceleration is zero
        for( auto it : thrustAccelerationsAlongTrajectory )
        {
            Eigen::Vector3d currentCartesianThrustAcceleration = it.second;
            BOOST_CHECK_SMALL( currentCartesianThrustAcceleration.norm(), 1.0E-14 );
        }

    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
