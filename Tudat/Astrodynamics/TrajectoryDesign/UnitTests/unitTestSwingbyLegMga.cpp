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

#include <Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h>
#include <Tudat/Basics/testMacros.h>

#include "Tudat/Astrodynamics/TrajectoryDesign/swingbyLegMga.h"

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

    //set the minimum pericenter radius
    const double minimumRadiusPlanet1 = 6351800;

    // Set up the test leg
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga legTest ( planet1Position, planet2Position, timeOfFlight, planet1Velocity,
                            sunGravitationalParameter, planet1GravitationalParameter,
                            pointerToVelocityBeforePlanet1, minimumRadiusPlanet1 );

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

    // Set the minimum pericenter radius
    const double minimumRadiusPlanet1 = 6351800;

    // Set up the test leg
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga legTest ( dummyPosition1, dummyPosition2, dummyTimeOfFlight, dummyVelocity1,
                            sunGravitationalParameter, planet1GravitationalParameter,
                            pointerToVelocityBeforePlanet1, minimumRadiusPlanet1 );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg.
    // legTest.calculateLeg( resultingVelocity, resultingDeltaV );


    // Specify the values for the parameters that are to be updated.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planet1Position ( -35554348960.8278, -102574987127.178,
                                            648696819.780156 );
    const Eigen::Vector3d planet2Position ( -35568329915.7073, -102569794949.529,
                                            650816245.825226 );
    const Eigen::Vector3d planet1Velocity ( 32851.224953746, -11618.7310059974,
                                            -2055.04615890989 );

    // Set the time of flight, which has to be converted from JD (in GTOP) to seconds (in Tudat).
    Eigen::VectorXd variableVector ( 1 );
    variableVector[ 0 ] = 449.385873819743 * physical_constants::JULIAN_DAY;

    // Pass both the ephemeris and the time of flight to the leg.
    legTest.updateEphemeris( planet1Position, planet2Position, planet1Velocity );
    legTest.updateDefiningVariables( variableVector );

    // Compute delta-V of the leg.
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

    // Specify the required parameters, for the second leg of Cassini1 within GTOP.
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
    const double timeOfFlight = 449.385873819743 * 24 * 60 * 60;

    // Set the gravitational parameters
    const double sunGravitationalParameter = 1.32712428e20;
    const double planet1GravitationalParameter = 3.24860e14;

    //set the minimum pericenter radius
    const double minimumRadiusPlanet1 = 6351800;

    // Set up the test leg
    using namespace tudat::transfer_trajectories;
    SwingbyLegMga legTest ( planet1Position, planet2Position, timeOfFlight, planet1Velocity,
                            sunGravitationalParameter, planet1GravitationalParameter,
                            pointerToVelocityBeforePlanet1, minimumRadiusPlanet1 );

    // Initiate vectors for storing the results.
    std::vector < Eigen::Vector3d > positionVector1, positionVector2;
    std::vector < double > timeVector1, timeVector2;

    // Test the functionality in case the leg has not been calculated yet.
    legTest.intermediatePoints( 240. * physical_constants::JULIAN_DAY,
                                positionVector1, timeVector1 );

    // Prepare the variables for calculating the leg actively.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Calculate the leg actively.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test the functionality in case the leg has been calculated already.
    legTest.intermediatePoints( 120. * physical_constants::JULIAN_DAY,
                                positionVector2, timeVector2 );

    // Test if the halfway points of the intermediate points functions match.
    BOOST_CHECK_CLOSE_FRACTION( timeVector1[1] , timeVector2[2], tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( positionVector1[1], positionVector2[2], tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat

