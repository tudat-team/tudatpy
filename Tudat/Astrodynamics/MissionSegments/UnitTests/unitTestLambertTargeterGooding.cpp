/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110113    E. Iorfida        First creation of the code.
 *      110126    E. Iorfida        Added test for velocities.
 *      110130    J. Melman         Simplified variable names, e.g., 'normOfVelocityVector' became
 *                                  'speed'. Also corrected 'tangential' to Suggested less trivial
 *                                  test case.
 *      110201    E. Iorfida        Added pointerToCelestialBody and modified variable names (from
 *                                  heliocentric, to inertial). Added elliptical case and check
 *                                  for the anti-clockwise direction of motion.
 *      110207    E. Iorfida        Added a more stric value for tolerance of semi-major axis for
 *                                  hyperbolic case.
 *      110208    E. Iorfida        Added CartesianPositionElements objects as input and
 *                                  CartesianVelocityElements objects as output.
 *      110512    K. Kumar          Updated code not to use dynamic memory allocation.
 *      110627    K. Kumar          Updated to use new predefined planets code.
 *      120326    D. Dirkx          Changed raw pointers to shared pointers.
 *      120416    T. Secretin       Boostified unit test.
 *      120623    T. Secretin       Adapted to new class implementation.
 *      120704    P. Musegaas       Removed dependency on old planet and gravity field classes.
 *                                  Changed tolerances. Minor changes.
 *      120705    T. Secretin       Added inertial velocity tests in testHyperbolicCase and
 *                                  testEllipticalCase.
 *
 *    References
 *      Mengali, G., and A.A. Quarta, Fondamenti di Meccanica del volo Spaziale.
 *      Noomen, R., Lambert targeter Excel file.
 *
 *    Notes
 *      DISCLAIMER: At the moment, the Lambert targeter only converges for about half of the cases.
 *      This is not evident from the tests below, but it was observed during simulations carried
 *      out by the author. The reason might very well be an erroneous definition of the starters.
 *
 *      The elliptical case was taken from Example 6.1, page 159-162 of ( Mengali, Quarta ). The
 *      hyperbolic case was taken from ( Noomen, R. ).
 *
 */

#define BOOST_TEST_MAIN

#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <TudatCore/Astrodynamics/BasicAstrodynamics/unitConversions.h>

#include "Tudat/Astrodynamics/MissionSegments/lambertTargeterGooding.h"

namespace tudat
{
namespace unit_tests
{

//! Test the Gooding Lambert targeting algorithm code.
BOOST_AUTO_TEST_SUITE( test_lambert_targeter_gooding )

//! Test hyperbolic case.
BOOST_AUTO_TEST_CASE( testHyperbolicCase )
{
    // Expected test result in meters.
    // Hyperbolic test case (results from excel file [1]).
    const double expectedValueOfSemiMajorAxisHyperbola = -1270129.3602e3;
    const double expectedValueOfRadialVelocityAtDepartureHyperbola = -0.74546e3;
    const double expectedValueOfRadialVelocityAtArrivalHyperbola = 0.69321e3;
    const double expectedValueOfTransverseVelocityAtDepartureHyperbola = 0.15674e3;
    const double expectedValueOfTransverseVelocityAtArrivalHyperbola = 0.10450e3;
    const Eigen::Vector3d expectedInertialVelocityAtDeparture( -745.457, 156.743, 0.0 );
    const Eigen::Vector3d expectedInertialVelocityAtArrival( 104.495, -693.209, 0.0 );

    // Tolerances.
    const double toleranceSemiMajorAxisHyperbola = 1.0e-7;
    const double toleranceVelocity = 1.0e-4;

    // Time conversions.
    const double timeOfFlightInDaysHyperbola = 100.0;
    const double timeOfFlightHyperbola = unit_conversions::convertJulianDaysToSeconds(
            timeOfFlightInDaysHyperbola );

    // Central body gravitational parameter.
    const double earthGravitationalParameter = 398600.4418e9;

    // The starting point is twice as far as L1 and L2, which is not really
    // realistic, but it is not about the case, but about the verification.
    using tudat::unit_conversions::convertAstronomicalUnitsToMeters;
    const Eigen::Vector3d positionAtDepartureHyperbola( convertAstronomicalUnitsToMeters( 0.02 ),
                                                        0.0, 0.0 ),
            positionAtArrivalHyperbola( 0.0, convertAstronomicalUnitsToMeters( -0.03 ), 0.0 );

    // Hyperbolic orbit case.
    tudat::mission_segments::LambertTargeterGooding lambertTargeterHyperbola(
                positionAtDepartureHyperbola, positionAtArrivalHyperbola, timeOfFlightHyperbola,
                earthGravitationalParameter );

    // Create local vectors for position and velocity.
    const Eigen::Vector3d positionDepartureHyperbola = positionAtDepartureHyperbola;
    const Eigen::Vector3d velocityDepartureHyperbola =
            lambertTargeterHyperbola.getInertialVelocityAtDeparture( );

    // Test if the computed semi-major axis corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getSemiMajorAxis( ),
                                expectedValueOfSemiMajorAxisHyperbola,
                                toleranceSemiMajorAxisHyperbola );

    // Test if the computed velocity components corresponds to the expected value within the
    // specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getRadialVelocityAtDeparture( ),
                                expectedValueOfRadialVelocityAtDepartureHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getRadialVelocityAtArrival( ),
                                expectedValueOfRadialVelocityAtArrivalHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getTransverseVelocityAtDeparture( ),
                                expectedValueOfTransverseVelocityAtDepartureHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getTransverseVelocityAtArrival( ),
                                expectedValueOfTransverseVelocityAtArrivalHyperbola,
                                toleranceVelocity );

    // Check that velocities match expected values within the defined tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                lambertTargeterHyperbola.getInertialVelocityAtDeparture( ).x( ),
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                lambertTargeterHyperbola.getInertialVelocityAtDeparture( ).y( ),
                                toleranceVelocity );
    BOOST_CHECK_SMALL( lambertTargeterHyperbola.getInertialVelocityAtDeparture( ).z( ),
                       toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                lambertTargeterHyperbola.getInertialVelocityAtArrival( ).x( ),
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                lambertTargeterHyperbola.getInertialVelocityAtArrival( ).y( ),
                                toleranceVelocity );
    BOOST_CHECK_SMALL( lambertTargeterHyperbola.getInertialVelocityAtArrival( ).z( ),
                       toleranceVelocity );

    // Test if the computed solution is anti-clockwise, if the z-component of the angular momentum
    // (h = r \times v) is positive.
    BOOST_CHECK_GT( positionDepartureHyperbola.cross( velocityDepartureHyperbola ).z( ),
                    std::numeric_limits< double >::epsilon( ) );
}

//! Test elliptical case.
BOOST_AUTO_TEST_CASE( testEllipticalCase )
{
    // Elliptical test case (results from example 6.1 page 159-162 [2]).
    // Set canonical units for Earth (see page 29 [2]).
    const double distanceUnit = 6.378136e6;
    const double timeUnit = 806.78;

    // Expected test result in meters.
    const double expectedValueOfSemiMajorAxisEllipse = 5.4214 * distanceUnit;
    const double expectedValueOfRadialVelocityAtDepartureEllipse = 2.73580e3;
    const double expectedValueOfRadialVelocityAtArrivalEllipse = 2.97503e3;
    const double expectedValueOfTransverseVelocityAtDepartureEllipse = 6.59430e3;
    const double expectedValueOfTransverseVelocityAtArrivalEllipse = 3.29715e3;
    const Eigen::Vector3d expectedInertialVelocityAtDeparture( 2735.8, 6594.3, 0.0 );
    const Eigen::Vector3d expectedInertialVelocityAtArrival( -1367.9, 4225.03, 0.0 );

    // Tolerance.
    const double toleranceSemiMajorAxisEllipse = 1.0e-3;
    const double toleranceVelocity = 1.0e-6;

    // Time conversions.
    const double timeOfFlightEllipse = 5.0 * timeUnit;

    // Central body gravitational parameter.
    const double earthGravitationalParameter = 398600.4418e9;

    const Eigen::Vector3d positionAtDepartureEllipse( 2.0 * distanceUnit, 0.0, 0.0 ),
            positionAtArrivalEllipse( 2.0 * distanceUnit, 2.0 * sqrt( 3.0 ) * distanceUnit, 0.0 );

    // Elliptical orbit case.
    tudat::mission_segments::LambertTargeterGooding lambertTargeterEllipse(
                positionAtDepartureEllipse, positionAtArrivalEllipse, timeOfFlightEllipse,
                earthGravitationalParameter );

    // Create local vectors for position and velocity.
    const Eigen::Vector3d positionDepartureEllipse = positionAtDepartureEllipse;
    const Eigen::Vector3d velocityDepartureEllipse =
            lambertTargeterEllipse.getInertialVelocityAtDeparture( );

    // Test if the computed semi-major axis corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getSemiMajorAxis( ),
                                expectedValueOfSemiMajorAxisEllipse,
                                toleranceSemiMajorAxisEllipse );

    // Test if the computed velocity components corresponds to the expected value within the
    // specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getRadialVelocityAtDeparture( ),
                                expectedValueOfRadialVelocityAtDepartureEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getRadialVelocityAtArrival( ),
                                expectedValueOfRadialVelocityAtArrivalEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getTransverseVelocityAtDeparture( ),
                                expectedValueOfTransverseVelocityAtDepartureEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getTransverseVelocityAtArrival( ),
                                expectedValueOfTransverseVelocityAtArrivalEllipse,
                                toleranceVelocity );

    // Check that velocities match expected values within the defined tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                lambertTargeterEllipse.getInertialVelocityAtDeparture( ).x( ),
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                lambertTargeterEllipse.getInertialVelocityAtDeparture( ).y( ),
                                toleranceVelocity );
    BOOST_CHECK_SMALL( lambertTargeterEllipse.getInertialVelocityAtDeparture( ).z( ),
                       toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                lambertTargeterEllipse.getInertialVelocityAtArrival( ).x( ),
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                lambertTargeterEllipse.getInertialVelocityAtArrival( ).y( ),
                                toleranceVelocity );
    BOOST_CHECK_SMALL( lambertTargeterEllipse.getInertialVelocityAtArrival( ).z( ),
                       toleranceVelocity );

    // Test if the computed solution is anti-clockwise, if the z-component of the angular momentum
    // (h = r \times v) is positive.
    BOOST_CHECK_GT( positionDepartureEllipse.cross( velocityDepartureEllipse ).z( ),
                    std::numeric_limits< double >::epsilon( ) );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
