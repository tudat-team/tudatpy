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
                ( Eigen::VectorXd( 2 )<<0.0, timeOfFlight ).finished( ),
                sunGravitationalParameter );

    // Set the Earth gravitational parameters.
    const double earthGravitationalParameter = 3.9860119e14;

    // Set the departure parking orbit (at inifinity)
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    Eigen::Vector3d departureVelocity =
            transferLeg.getDepartureVelocity( );
    DepartureWithFixedOutgoingVelocityNode departureNode(
                constantEphemeris1,
                ( Eigen::VectorXd( 1 )<< 0.0 ).finished( ),
                earthGravitationalParameter, semiMajorAxis, eccentricity,
                [=]( ){ return departureVelocity; } );


    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity = transferLeg.getArrivalVelocity( );
    double resultingDeltaV = departureNode.getNodeDeltaV( );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedVelocity, resultingVelocity, tolerance );

}

////! Test updating the variables.
//BOOST_AUTO_TEST_CASE( testUpdatingVariables )
//{
//    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
//    // achieved.
//    const double tolerance = 1.0e-6;

//    // Expected test result based on the first leg of the ideal Cassini 1 trajectory as modelled by
//    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
//    // Astrotoolbox.
//    const double expectedDeltaV = 2754.63598732926;
//    const Eigen::Vector3d expectedVelocity ( 34216.4827530912, -15170.1440677825,
//                                             395.792122152361 );

//    // Specify the required parameters.
//    // Set the dummy positions and velocities.
//    const Eigen::Vector3d dummyPosition1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
//    const Eigen::Vector3d dummyPosition2 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
//    const Eigen::Vector3d dummyVelocity1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

//    // Set the dummy time of flight.
//    const double dummyTimeOfFlight = TUDAT_NAN;

//    // Set the gravitational parameters.
//    const double sunGravitationalParameter = 1.32712428e20;
//    const double earthGravitationalParameter = 3.9860119e14;

//    // Set the departure parking orbit (at inifinity)
//    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
//    const double eccentricity = 0.;

//    // Set test case.
//    using namespace tudat::transfer_trajectories;
//    DepartureLegMga legTest ( dummyPosition1, dummyPosition2, dummyTimeOfFlight, dummyVelocity1,
//                              sunGravitationalParameter, earthGravitationalParameter,
//                              semiMajorAxis, eccentricity );

//    // Prepare the variables for the results.
//    Eigen::Vector3d resultingVelocity;
//    double resultingDeltaV;

//    // Compute delta-V of the leg.
//    //legTest.calculateLeg( resultingVelocity, resultingDeltaV );

//    // Specify the values for the parameters that are to be updated.
//    // Set the planetary positions and velocities.
//    const Eigen::Vector3d planet1Position ( 113191651440.549, 95992973233.5064, 0. );
//    const Eigen::Vector3d planet2Position ( -35554348960.8278, -102574987127.178,
//                                            648696819.780156 );
//    const Eigen::Vector3d planet1Velocity ( -19752.2624404406, 22607.9064746733, 0. );

//    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
//    Eigen::VectorXd variableVector ( 1 );
//    variableVector[ 0 ] = 158.302027105278 * physical_constants::JULIAN_DAY;

//    // Pass both the ephemeris and the time of flight to the leg.
//    legTest.updateEphemeris( planet1Position, planet2Position, planet1Velocity );
//    legTest.updateDefiningVariables( variableVector );

//    // Compute delta-V of the leg.
//    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

//    // Test if the computed delta-V corresponds to the expected value within the specified
//    // tolerance and if the computed velocity before target planet matches the expected velocity
//    // within the specified tolerance.
//    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
//    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedVelocity, resultingVelocity, tolerance );
//}

////! Test intermediate points function.
//BOOST_AUTO_TEST_CASE( testIntermediatePoints )
//{
//    // Set tolerance.
//    const double tolerance = 1e-13;

//    // Specify the required parameters, for the first leg of Cassini1 within GTOP.
//    // Set the planetary positions and velocities.
//    const Eigen::Vector3d planet1Position ( 113191651440.549, 95992973233.5064, 0. );
//    const Eigen::Vector3d planet2Position ( -35554348960.8278, -102574987127.178,
//                                            648696819.780156 );
//    const Eigen::Vector3d planet1Velocity ( -19752.2624404406, 22607.9064746733, 0. );

//    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
//    const double timeOfFlight = 158.302027105278 * 24 * 60 * 60;

//    // Set the gravitational parameter.
//    const double sunGravitationalParameter = 1.32712428e20;
//    const double earthGravitationalParameter = 3.9860119e14;

//    // Set the departure parking orbit (at inifinity)
//    const double semiMajorAxis = 1e300;
//    const double eccentricity = 0.;

//    // Set test case.
//    using namespace tudat::transfer_trajectories;
//    DepartureLegMga legTest ( planet1Position, planet2Position, timeOfFlight, planet1Velocity,
//                              sunGravitationalParameter, earthGravitationalParameter,
//                              semiMajorAxis, eccentricity );

//    // Initiate vectors for storing the results.
//    std::vector < Eigen::Vector3d > positionVector1, positionVector2;
//    std::vector < double > timeVector1, timeVector2;

//    // Test the functionality in case the leg has not been calculated yet.
//    legTest.intermediatePoints( 80. * physical_constants::JULIAN_DAY,
//                                positionVector1, timeVector1 );

//    // Prepare the variables for calculating the leg actively.
//    Eigen::Vector3d resultingVelocity;
//    double resultingDeltaV;

//    // Calculate the leg actively.
//    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

//    // Test the functionality in case the leg has been calculated already.
//    legTest.intermediatePoints( 40. * physical_constants::JULIAN_DAY,
//                                positionVector2, timeVector2 );

//    // Test if the halfway points of the intermediate points functions match.
//    BOOST_CHECK_CLOSE_FRACTION( timeVector1[1] , timeVector2[2], tolerance );
//    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( positionVector1[1], positionVector2[2], tolerance );
//}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
