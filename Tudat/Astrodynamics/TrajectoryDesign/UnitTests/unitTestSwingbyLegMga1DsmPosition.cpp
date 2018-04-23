/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120509    P. Musegaas       First creation of code.
 *      120611    P. Musegaas       Adaptation to new mission segments functions and update of
 *                                  of functionality.
 *      120628    P. Musegaas       Added unit test for new functionality of re-using the class.
 *
 *    References
 *      Musegaas, P. (2012). Optimization of Space Trajectories Including Multiple Gravity Assists
 *          and Deep Space Maneuvers. MSc Thesis, Delft University of Technology, Delft,
 *          The Netherlands.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include <Tudat/Basics/testMacros.h>
#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>

#include "SpaceTrajectories/swingbyLegMga1DsmPosition.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of departure leg within MGA-1DSM-PF trajectory model
BOOST_AUTO_TEST_SUITE( test_swingby_leg_mga_1dsm_position )

//! Test delta-V computation
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
    const Eigen::Vector3d planet1Velocity ( 24956.3863503886, -24481.33754925, -1774.16153112584 );

    // Set velocity before departure body.
    const Eigen::Vector3d velocityBeforePlanet1 ( 28586.0252553367, -17610.9003149933,
                                                  -1915.53135757897 );
    boost::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = boost::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 180.510754824 * physical_constants::JULIAN_DAY;

    // Set the gravitational parameters.
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    // Set the minimum pericenter radius.
    const double minimumRadiusPlanet1 = 1.1 * 6052000;

    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.317174785637;
    //?!? Eigen::Vector3d dsmLocation (  53021338819.4998, -278970556.494799, -5016845311.62288 );
    const double dimensionlessRadiusDsm = 0.491956793376231;
    const double inPlaneAngle = 2.33292208993544;
    const double outOfPlaneAngle = -0.03659314664852;

    // Set test case.
    using namespace tudat::spaceTrajectories;
    SwingbyLegMga1DsmPosition legTest ( planet1Position, planet2Position, timeOfFlight,
                                        planet1Velocity, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, minimumRadiusPlanet1,
                                        dsmTimeOfFlightFraction, //dsmLocation
                                        dimensionlessRadiusDsm, inPlaneAngle, outOfPlaneAngle );

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

    // Expected test result based on the fourth leg of the ideal Messenger trajectory as modelled
    // by GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 1415.44020553569;
    const Eigen::Vector3d expectedVelocity ( -5080.63624082843, 55179.1205883319,
                                             3549.41832192081 );

    // Specify the required parameters.
    // Set the dummy positions and velocities.
    const Eigen::Vector3d dummyPosition1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyPosition2 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyVelocity1 ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

    // Set velocity before departure body.
    const Eigen::Vector3d velocityBeforePlanet1 ( 28586.0252553367, -17610.9003149933,
                                                  -1915.53135757897 );
    boost::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = boost::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the dummy time of flight.
    const double dummyTimeOfFlight = TUDAT_NAN;

    // Set the gravitational parameters.
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    // Set the minimum pericenter radius.
    const double minimumRadiusPlanet1 = 1.1 * 6052000;

    // Set the dummy DSM variables.
    const double dummyVariable1 = TUDAT_NAN;
    const double dummyVariable2 = TUDAT_NAN;
    const double dummyVariable3 = TUDAT_NAN;
    const double dummyVariable4 = TUDAT_NAN;

    // Set test case.
    using namespace tudat::spaceTrajectories;
    SwingbyLegMga1DsmPosition legTest ( dummyPosition1, dummyPosition2, dummyTimeOfFlight,
                                        dummyVelocity1, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, minimumRadiusPlanet1,
                                        dummyVariable1,
                                        dummyVariable2, dummyVariable3, dummyVariable4 );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg.
    // legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Specify the values for the parameters that are to be updated.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -75133023393.8197, -77873986249.455,
                                            3277461620.51787 );
    const Eigen::Vector3d planet2Position ( 53627979831.9489, -5044669560.01206,
                                            -5339232305.5444 );
    const Eigen::Vector3d planet1Velocity ( 24956.3863503886, -24481.33754925, -1774.16153112584 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 180.510754824 * physical_constants::JULIAN_DAY;

    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.317174785637;
    const double dimensionlessRadiusDsm = 0.491956793376231;
    const double inPlaneAngle = 2.33292208993544;
    const double outOfPlaneAngle = -0.03659314664852;

    // Create a variable vector containing these parameters.
    Eigen::VectorXd variableVector( 5 );
    variableVector << timeOfFlight, dsmTimeOfFlightFraction, dimensionlessRadiusDsm, inPlaneAngle,
                      outOfPlaneAngle;

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
    const Eigen::Vector3d planet1Velocity ( 24956.3863503886, -24481.33754925, -1774.16153112584 );

    // Set velocity before departure body.
    const Eigen::Vector3d velocityBeforePlanet1 ( 28586.0252553367, -17610.9003149933,
                                                  -1915.53135757897 );
    boost::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet1
            = boost::make_shared< Eigen::Vector3d > ( velocityBeforePlanet1 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    const double timeOfFlight = 180.510754824 * physical_constants::JULIAN_DAY;

    // Set the gravitational parameters.
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    // Set the minimum pericenter radius.
    const double minimumRadiusPlanet1 = 1.1 * 6052000;

    // Set DSM parameters
    const double dsmTimeOfFlightFraction = 0.317174785637;
    const double dimensionlessRadiusDsm = 0.491956793376231;
    const double inPlaneAngle = 2.33292208993544;
    const double outOfPlaneAngle = -0.03659314664852;


    // Set test case.
    using namespace tudat::spaceTrajectories;
    SwingbyLegMga1DsmPosition legTest ( planet1Position, planet2Position, timeOfFlight,
                                        planet1Velocity, sunGravitationalParameter,
                                        planet1GravitationalParameter,
                                        pointerToVelocityBeforePlanet1, minimumRadiusPlanet1,
                                        dsmTimeOfFlightFraction, //dsmLocation
                                        dimensionlessRadiusDsm, inPlaneAngle, outOfPlaneAngle );

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
