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

#include "tudat/astro/trajectory_design/departureLegMga1DsmVelocity.h"

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
    const double expectedDeltaV = 1408.99421278 + 910.801673341206;
    const Eigen::Vector3d expectedVelocity ( 17969.3166254715, -23543.6915939142,
                                             6.38384671663485 );

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -148689402143.081, 7454242895.97607, 0. );
    const Eigen::Vector3d planet2Position ( -128359548637.032, -78282803797.7343, 0. );
    const Eigen::Vector3d planet1Velocity ( -1976.48781307596, -29863.4035321021, 0. );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 399.999999715 * 24 * 60 * 60;

    // Set the gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;
    const double earthGravitationalParameter = 3.9860119e14;

    // Set the departure parking orbit (at inifinity)
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.234594654679;
    const double excessVelocityMagnitude = 1408.99421278;
    const double excessVelocityInPlaneAngle = 0.37992647165 * 2 * 3.14159265358979;
    const double excessVelocityOutOfPlaneAngle = std::acos(  2 * 0.498004040298 - 1 ) -
                                                 3.14159265358979 / 2;

    // Set test case.
    using namespace tudat::transfer_trajectories;
    DepartureLegMga1DsmVelocity legTest ( planet1Position, planet2Position, timeOfFlight,
                                          planet1Velocity, sunGravitationalParameter,
                                          earthGravitationalParameter, semiMajorAxis, eccentricity,
                                          dsmTimeOfFlightFraction, excessVelocityMagnitude,
                                          excessVelocityInPlaneAngle,
                                          excessVelocityOutOfPlaneAngle );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedVelocity, resultingVelocity, tolerance );
}

//! Test updating the variables.
BOOST_AUTO_TEST_CASE( testUpdatingVariables )
{
    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
    // achieved.
    const double tolerance = 1.0e-6;

    // Expected test result based on the first leg of the ideal Messenger trajectory as modelled by
    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 1408.99421278 + 910.801673341206;
    const Eigen::Vector3d expectedVelocity ( 17969.3166254715, -23543.6915939142,
                                             6.38384671663485 );

    // Specify the required parameters.
    // Set the dummy positions and velocities.
    const Eigen::Vector3d dummyPosition1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyPosition2 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyVelocity1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

    // Set the dummy time of flight.
    const double dummyTimeOfFlight = TUDAT_NAN;

    // Set the gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;
    const double earthGravitationalParameter = 3.9860119e14;

    // Set the departure parking orbit (at inifinity)
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    // Set the dummy DSM variables.
    const double dummyVariable1 = TUDAT_NAN;
    const double dummyVariable2 = TUDAT_NAN;
    const double dummyVariable3 = TUDAT_NAN;
    const double dummyVariable4 = TUDAT_NAN;

    // Set test case.
    using namespace tudat::transfer_trajectories;
    DepartureLegMga1DsmVelocity legTest ( dummyPosition1, dummyPosition2, dummyTimeOfFlight,
                                          dummyVelocity1, sunGravitationalParameter,
                                          earthGravitationalParameter, semiMajorAxis, eccentricity,
                                          dummyVariable1, dummyVariable2, dummyVariable3,
                                          dummyVariable4 );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg. In try/catch statement because a runtime error is expected.
    try
    {
        legTest.calculateLeg( resultingVelocity, resultingDeltaV );
    }
    catch( std::runtime_error const& )
 { }

    // Specify the values for the parameters that are to be updated.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -148689402143.081, 7454242895.97607, 0. );
    const Eigen::Vector3d planet2Position ( -128359548637.032, -78282803797.7343, 0. );
    const Eigen::Vector3d planet1Velocity ( -1976.48781307596, -29863.4035321021, 0. );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 399.999999715 * physical_constants::JULIAN_DAY;

    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.234594654679;
    const double excessVelocityMagnitude = 1408.99421278;
    const double excessVelocityInPlaneAngle = 0.37992647165 * 2 * 3.14159265358979;
    const double excessVelocityOutOfPlaneAngle = std::acos(  2 * 0.498004040298 - 1 ) -
                                                 3.14159265358979 / 2;

    // Create a variable vector containing these parameters.
    Eigen::VectorXd variableVector( 5 );
    variableVector << timeOfFlight, dsmTimeOfFlightFraction, excessVelocityMagnitude,
                      excessVelocityInPlaneAngle, excessVelocityOutOfPlaneAngle;

    // Pass both the new ephemeris and trajectory defining variables to the leg.
    legTest.updateEphemeris( planet1Position, planet2Position, planet1Velocity );
    legTest.updateDefiningVariables( variableVector );

    // Recompute the delta-V of the leg.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance and if the computed velocity before target planet matches the expected velocity
    // within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedVelocity, resultingVelocity, tolerance );
}

//! Test intermediate points function.
BOOST_AUTO_TEST_CASE( testIntermediatePoints )
{
    // Set tolerance.
    const double tolerance = 1e-13;

    // Specify the required parameters, for the first leg of Messenger within GTOP.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -148689402143.081, 7454242895.97607, 0. );
    const Eigen::Vector3d planet2Position ( -128359548637.032, -78282803797.7343, 0. );
    const Eigen::Vector3d planet1Velocity ( -1976.48781307596, -29863.4035321021, 0. );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 399.999999715 * 24 * 60 * 60;

    // Set the gravitational parameter
    const double sunGravitationalParameter = 1.32712428e20;
    const double earthGravitationalParameter = 3.9860119e14;

    // Set the departure parking orbit (at inifinity)
    const double semiMajorAxis = 1e300;
    const double eccentricity = 0.;


    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.234594654679;
    const double excessVelocityMagnitude = 1408.99421278;
    const double excessVelocityInPlaneAngle = 0.37992647165 * 2 * 3.14159265358979;
    const double excessVelocityOutOfPlaneAngle = std::acos(  2 * 0.498004040298 - 1 ) -
                                                 3.14159265358979 / 2;

    // Set test case.
    using namespace tudat::transfer_trajectories;
    DepartureLegMga1DsmVelocity legTest ( planet1Position, planet2Position, timeOfFlight,
                                          planet1Velocity, sunGravitationalParameter,
                                          earthGravitationalParameter, semiMajorAxis, eccentricity,
                                          dsmTimeOfFlightFraction, excessVelocityMagnitude,
                                          excessVelocityInPlaneAngle,
                                          excessVelocityOutOfPlaneAngle );

    // Initiate vectors for storing the results.
    std::vector < Eigen::Vector3d > positionVector1, positionVector2;
    std::vector < double > timeVector1, timeVector2;

    // Test the functionality in case the leg has not been calculated yet.
    legTest.intermediatePoints( 160. * 24 * 60 * 60, positionVector1, timeVector1 );

    // Prepare the variables for calculating the leg actively.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Calculate the leg actively.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test the functionality in case the leg has been calculated already.
    legTest.intermediatePoints( 80. * 24 * 60 * 60, positionVector2, timeVector2 );

    // Test if the halfway points in the first part of the leg match between the intermediate
    // points functions.
    BOOST_CHECK_CLOSE_FRACTION( timeVector1[1] , timeVector2[2], tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( positionVector1[1], positionVector2[2], tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
