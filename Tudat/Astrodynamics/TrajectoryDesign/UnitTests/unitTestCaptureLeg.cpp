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
 *      120611    P. Musegaas       First creation of code.
 *      120628    P. Musegaas       Added unit test for new functionality of re-using the class.
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

#include <Tudat/Basics/testMacros.h>

#include "SpaceTrajectories/captureLeg.h"

namespace tudat
{
namespace unit_tests
{

//! Test implementation of departure leg within MGA trajectory model.
BOOST_AUTO_TEST_SUITE( test_capture_leg )

//! Test delta-V computation for infinite parking orbit.
BOOST_AUTO_TEST_CASE( testVelocitiesInfiniteParkingOrbit )
{
    // Set tolerance.
    const double tolerance = 1.0e-13;

    // Expected test result based on the last part of the ideal Messenger trajectory as modelled by
    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    // Note that this concerns a capture at infinite distance and hence is not the best test case.
    const double expectedDeltaV = 4632.24252029314;

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planetPosition ( 53627979831.9492, -5044669560.01491,
                                           -5339232305.54465 );
    const Eigen::Vector3d planetVelocity ( -4857.99954791498, 50668.339570669, 4579.44661303178 );

    // Set velocity before capture body.
    const Eigen::Vector3d velocityBeforePlanet ( -5080.6362408257, 55179.1205883308,
                                                 3549.4183219232 );

    boost::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet
            = boost::make_shared< Eigen::Vector3d > ( velocityBeforePlanet );

    // Set the time of flight, which is irrelevant for the deltaV consumption.
    const double timeOfFlight = 1.0;

    // Set the gravitational parameters. The sun's is irrelevant for the deltaV consumption, but
    // required only for other functionality in the class.
    const double sunGravitationalParameter = 1.32712428e20;
    const double mercuryGravitationalParameter = 2.2321e13;

    // Set the capture orbit (at inifinity).
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    // Set test case.
    using namespace tudat::spaceTrajectories;
    CaptureLeg legTest ( planetPosition, timeOfFlight, planetVelocity, sunGravitationalParameter,
                         mercuryGravitationalParameter, pointerToVelocityBeforePlanet,
                         semiMajorAxis, eccentricity );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test delta-V computation for circular parking orbit.
BOOST_AUTO_TEST_CASE( testVelocitiesCircularParkingOrbit )
{
    // Set tolerance.
    const double tolerance = 1.0e-13;

    // Expected test result based on Table 18.2 of [Wakker, 2007]. Capture velocity to arrive at
    // a parking orbit around Mars of 1.1 times the planetary radius starting from Earth. Result
    // calculated with more precision separately using Equation 18.25 of [Wakker, 2007].
    const double expectedDeltaV = 2087.1062716740798;

    // Specify the required parameters.
    // Set the planetary positions and velocities.
    const Eigen::Vector3d planetPosition ( 227936637698.942, 0.0, 0.0 );
    const Eigen::Vector3d planetVelocity ( 0.0, 24129.4836355380, 0.0 );

    // Set velocity before capture body.
    const Eigen::Vector3d velocityBeforePlanet ( 0.0, 21480.6500358053, 0.0 );

    boost::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet
            = boost::make_shared< Eigen::Vector3d > ( velocityBeforePlanet );

    // Set the time of flight, which is irrelevant for the deltaV consumption.
    const double timeOfFlight = 1.0;

    // Set the gravitational parameters. The Sun's is irrelevant for the deltaV consumption, but
    // required only for other functionality in the class.
    const double sunGravitationalParameter = 1.32712428e20;
    const double marsGravitationalParameter = 4.2830e13;

    // Set the capture orbit.
    const double semiMajorAxis = 1.1 * 3.3895e6;;
    const double eccentricity = 0.;

    // Set test case.
    using namespace tudat::spaceTrajectories;
    CaptureLeg legTest ( planetPosition, timeOfFlight, planetVelocity, sunGravitationalParameter,
                         marsGravitationalParameter, pointerToVelocityBeforePlanet,
                         semiMajorAxis, eccentricity );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

//! Test updating the variables.
BOOST_AUTO_TEST_CASE( testUpdatingVariables )
{
    // Set tolerance.
    const double tolerance = 1.0e-13;

    // Expected test result based on the last part of the ideal Messenger trajectory as modelled by
    // GTOP software distributed and downloadable from the ESA website, or within the PaGMO
    // Astrotoolbox.
    const double expectedDeltaV = 4632.24252029314;

    // Specify the required parameters.
    // Set the dummy positions and velocities.
    const Eigen::Vector3d dummyPosition ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );
    const Eigen::Vector3d dummyVelocity ( TUDAT_NAN, TUDAT_NAN, TUDAT_NAN );

    // Set velocity before capture body.
    const Eigen::Vector3d velocityBeforePlanet ( -5080.6362408257, 55179.1205883308,
                                                 3549.4183219232 );

    boost::shared_ptr< Eigen::Vector3d > pointerToVelocityBeforePlanet
            = boost::make_shared< Eigen::Vector3d > ( velocityBeforePlanet );

    // Set the dummy time of flight.
    const double dummyTimeOfFlight = TUDAT_NAN;

    // Set the gravitational parameters. The sun's is irrelevant for the deltaV consumption, but
    // required only for other functionality in the class.
    const double sunGravitationalParameter = 1.32712428e20;
    const double mercuryGravitationalParameter = 2.2321e13;

    // Set the capture orbit (at inifinity).
    const double semiMajorAxis = std::numeric_limits< double >::infinity( );
    const double eccentricity = 0.;

    // Set test case.
    using namespace tudat::spaceTrajectories;
    CaptureLeg legTest ( dummyPosition, dummyTimeOfFlight, dummyVelocity,
                         sunGravitationalParameter, mercuryGravitationalParameter,
                         pointerToVelocityBeforePlanet, semiMajorAxis, eccentricity );

    // Prepare the variables for the results.
    Eigen::Vector3d resultingVelocity;
    double resultingDeltaV;

    // Compute delta-V of the leg.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Specify the values for the parameters that are to be updated.
    // Set the planetary positions and velocities.
    Eigen::Vector3d planetPosition ( 53627979831.9492, -5044669560.01491, -5339232305.54465 );
    Eigen::Vector3d planetVelocity ( -4857.99954791498, 50668.339570669, 4579.44661303178 );

    // Set the time of flight, which is irrelevant for the deltaV consumption.
    Eigen::VectorXd variableVector ( 1 );
    variableVector[ 0 ] = 1;

    // Pass both the ephemeris and the time of flight to the leg.
    // Note that the second position vector is unused in the capture leg, hence a dummy is passed.
    legTest.updateEphemeris( planetPosition, dummyPosition, planetVelocity );
    legTest.updateDefiningVariables( variableVector );

    // Compute delta-V of the leg.
    legTest.calculateLeg( resultingVelocity, resultingDeltaV );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedDeltaV, resultingDeltaV, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
