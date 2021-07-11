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

#include "tudat/astro/trajectory_design/swingbyLegMga1DsmVelocity.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of departure leg within MGA trajectory model.
BOOST_AUTO_TEST_SUITE( test_swingby_leg_mga_1dsm_velocity )

//! Test delta-V computation for an unpowered gravity assist leg.
BOOST_AUTO_TEST_CASE( testVelocitiesUnpoweredGravityAssist )
{
    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
    // achieved.
    const double tolerance = 1.0e-6;

    // Expected test result based on the fourth leg of the ideal Messenger trajectory as modelled
    // by GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 1415.44020553569;
    const Eigen::Vector3d expectedVelocity ( -5080.63624082843, 55179.1205883319,
                                             3549.41832192081 );

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -75133023393.8197, -77873986249.455,
                                            3277461620.51787 );
    const Eigen::Vector3d planet2Position ( 53627979831.9489, -5044669560.01206,
                                            -5339232305.5444 );
    const Eigen::Vector3d planet1Velocity ( 24956.3863503886, -24481.33754925,
                                            -1774.16153112584 );

    // Set velocity before departure body.
    const Eigen::Vector3d velocityBeforePlanet1 ( 28586.0252553367, -17610.9003149933,
                                                  -1915.53135757897 );
    std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = std::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 180.510754824 * physical_constants::JULIAN_DAY;

    // Set the gravitational parameters.
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    //set model specific variables.
    const double dsmTimeOfFlightFraction = 0.317174785637;
    const double rotationAngle = 1.34317576594;
    const double pericenterRadius = 1.10000000891 * 6052000;
    const double swingbyDeltaV = 0.;

    // Set up the test leg.
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga1DsmVelocity legTest ( planet1Position, planet2Position, timeOfFlight,
                                        planet1Velocity, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, dsmTimeOfFlightFraction,
                                        rotationAngle, pericenterRadius, swingbyDeltaV );

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

//! Test delta-V computation for a powered gravity assist leg.
BOOST_AUTO_TEST_CASE( testVelocitiesPoweredGravityAssist )
{
    // Set tolerance. Due to iterative nature in the process, (much) higher accuracy cannot be
    // achieved.
    const double tolerance = 1.0e-6;

    // Expected test result based on the second leg of the ideal Cassini 1 trajectory as modelled
    // by GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 1090.64622926316;
    const Eigen::Vector3d expectedVelocity ( 37952.8844553685, -14096.9656774702,
                                             -5753.51245833761 );

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -35554348960.8278, -102574987127.178,
                                            648696819.780156 );
    const Eigen::Vector3d planet2Position ( -35568329915.7073, -102569794949.529,
                                            650816245.825226 );
    const Eigen::Vector3d planet1Velocity ( 32851.224953746, -11618.7310059974,
                                            -2055.04615890989 );

    // Set velocity before departure body
    const Eigen::Vector3d velocityBeforePlanet1 ( 34216.4827530912, -15170.1440677825,
                                                  395.792122152361 );
    std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = std::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 449.385873819743 * physical_constants::JULIAN_DAY;

    // Set the gravitational parameters
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    //set model specific variables.
    const double dsmTimeOfFlightFraction = 0.2;// This should be irrelevant for a perfect model.
    const double rotationAngle = -2.0291949514117;
    const double pericenterRadius = 6351801.04541467;
    const double swingbyDeltaV = 1090.64622870007;

    // Set up the test leg.
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga1DsmVelocity legTest ( planet1Position, planet2Position, timeOfFlight,
                                        planet1Velocity, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, dsmTimeOfFlightFraction,
                                        rotationAngle, pericenterRadius, swingbyDeltaV );

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

    // Expected test result based on the second leg of the ideal Cassini 1 trajectory as modelled
    // by GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 1090.64622926316;
    const Eigen::Vector3d expectedVelocity ( 37952.8844553685, -14096.9656774702,
                                             -5753.51245833761 );

    // Specify the required parameters.
    // Set the dummy positions and velocities.
    const Eigen::Vector3d dummyPosition1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyPosition2 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyVelocity1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

    // Set velocity before departure body
    const Eigen::Vector3d velocityBeforePlanet1 ( 34216.4827530912, -15170.1440677825,
                                                  395.792122152361 );
    std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = std::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the dummy time of flight.
    const double dummyTimeOfFlight = TUDAT_NAN;

    // Set the gravitational parameters
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    // Set the dummy DSM variables.
    const double dummyVariable1 = TUDAT_NAN;
    const double dummyVariable2 = TUDAT_NAN;
    const double dummyVariable3 = TUDAT_NAN;
    const double dummyVariable4 = TUDAT_NAN;

    // Set up the test leg.
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga1DsmVelocity legTest ( dummyPosition1, dummyPosition2, dummyTimeOfFlight,
                                        dummyVelocity1, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, dummyVariable1,
                                        dummyVariable2, dummyVariable3, dummyVariable4 );

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
    const Eigen::Vector3d planet1Position ( -35554348960.8278, -102574987127.178,
                                            648696819.780156 );
    const Eigen::Vector3d planet2Position ( -35568329915.7073, -102569794949.529,
                                            650816245.825226 );
    const Eigen::Vector3d planet1Velocity ( 32851.224953746, -11618.7310059974,
                                            -2055.04615890989 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 449.385873819743 * physical_constants::JULIAN_DAY;

    //set model specific variables.
    const double dsmTimeOfFlightFraction = 0.2;// This should be irrelevant for a perfect model.
    const double rotationAngle = -2.0291949514117;
    const double pericenterRadius = 6351801.04541467;
    const double swingbyDeltaV = 1090.64622870007;

    // Create a variable vector containing these parameters.
    Eigen::VectorXd variableVector( 5 );
    variableVector << timeOfFlight, dsmTimeOfFlightFraction, rotationAngle, pericenterRadius,
                      swingbyDeltaV;

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

    // Specify the required parameters, for the fourth leg of Messenger within GTOP.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -75133023393.8197, -77873986249.455,
                                            3277461620.51787 );
    const Eigen::Vector3d planet2Position ( 53627979831.9489, -5044669560.01206,
                                            -5339232305.5444 );
    const Eigen::Vector3d planet1Velocity ( 24956.3863503886, -24481.33754925,
                                            -1774.16153112584 );

    // Set velocity before departure body.
    const Eigen::Vector3d velocityBeforePlanet1 ( 28586.0252553367, -17610.9003149933,
                                                  -1915.53135757897 );
    std::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = std::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 180.510754824 * physical_constants::JULIAN_DAY;

    // Set the gravitational parameters.
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    //set model specific variables.
    const double dsmTimeOfFlightFraction = 0.317174785637;
    const double rotationAngle = 1.34317576594;
    const double pericenterRadius = 1.10000000891 * 6052000;
    const double swingbyDeltaV = 0;

    // Set up the test leg.
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga1DsmVelocity legTest ( planet1Position, planet2Position, timeOfFlight,
                                        planet1Velocity, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, dsmTimeOfFlightFraction,
                                        rotationAngle, pericenterRadius, swingbyDeltaV );

    // Initiate vectors for storing the results.
    std::vector < Eigen::Vector3d > positionVector1, positionVector2;
    std::vector < double > timeVector1, timeVector2;

    // Test the functionality in case the leg has not been calculated yet.
    legTest.intermediatePoints( 30. * physical_constants::JULIAN_DAY,
                                positionVector1, timeVector1 );

    // Prepare the variables for calculating the leg actively.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Calculate the leg actively.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test the functionality in case the leg has been calculated already.
    legTest.intermediatePoints( 15. * physical_constants::JULIAN_DAY,
                                positionVector2, timeVector2 );

    // Test if the halfway points in the first part of the leg match between the intermediate
    // points functions.
    BOOST_CHECK_CLOSE_FRACTION( timeVector1[1] , timeVector2[2], tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( positionVector1[1], positionVector2[2], tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
