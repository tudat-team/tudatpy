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

#include <limits>

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <tudat/astro/basic_astro/physicalConstants.h>
#include <tudat/basics/testMacros.h>

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/mission_segments/transferNode.h"
#include "tudat/astro/mission_segments/transferLeg.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of departure leg within MGA trajectory model.
BOOST_AUTO_TEST_SUITE( test_departure_leg_mga )

//! Test delta-V computation.
BOOST_AUTO_TEST_CASE( testVelocities )
{
    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
    // achieved.
    const double tolerance = 1.0e-6;

    // Expected test result based on the first leg of the ideal Cassini 1 trajectory as modelled by
    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 2754.63598732926;
    const Eigen::Vector3d expectedVelocity ( 34216.4827530912, -15170.1440677825,
                                             395.792122152361 );

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector6d planet1State =
            ( Eigen::Vector6d( ) <<
              113191651440.549, 95992973233.5064, 0.0,
              -19752.2624404406, 22607.9064746733, 0.0 ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris1 =
            std::make_shared< ephemerides::ConstantEphemeris >( planet1State );

    const Eigen::Vector6d planet2State =
            ( Eigen::Vector6d( ) <<
              -35554348960.8278, -102574987127.178, 648696819.780156,
              TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    std::shared_ptr< ephemerides::Ephemeris > constantEphemeris2 =
            std::make_shared< ephemerides::ConstantEphemeris >( planet2State );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 158.302027105278 * physical_constants::JULIAN_DAY;

    // Set the Sun gravitational parameters.
    const double sunGravitationalParameter = 1.32712428e20;

    using namespace mission_segments;
    UnpoweredUnperturbedTransferLeg transferLeg(
                constantEphemeris1, constantEphemeris2,
                sunGravitationalParameter );
    transferLeg.updateLegParameters( ( Eigen::VectorXd( 2 )<<0.0, timeOfFlight ).finished( ) );

    // Set the Earth gravitational parameters.
    const double earthGravitationalParameter = 3.9860119e14;

    // Set the departure parking orbit (at inifinity)
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    Eigen::Vector3d departureVelocity =
            transferLeg.getDepartureVelocity( );
    DepartureWithFixedOutgoingVelocityNode departureNode(
                constantEphemeris1,
                earthGravitationalParameter, semiMajorAxis, eccentricity,
                [=]( ){ return departureVelocity; } );
    departureNode.updateNodeParameters( ( Eigen::VectorXd( 1 )<< 0.0 ).finished( ) );


    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity = transferLeg.getArrivalVelocity( );
    double resultingDeltaV = departureNode.getNodeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedVelocity, resultingVelocity, tolerance );

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
        BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.begin( )->second( i ) - planet1State( i ) ), 1.0E-2 );
        BOOST_CHECK_SMALL( std::fabs( statesAlongTrajectory.rbegin( )->second( i ) - planet2State( i ) ), 1.0E-2 );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
