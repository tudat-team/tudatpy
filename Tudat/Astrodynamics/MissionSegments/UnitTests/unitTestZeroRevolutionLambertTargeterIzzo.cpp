/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Mengali, G., and A.A. Quarta, Fondamenti di Meccanica del volo Spaziale.
 *      Noomen, R., Lambert targeter Excel file.
 *      Izzo, D., Keplerian_Toolbox.
 *
 *    Notes
 *      The elliptical case was taken from Example 6.1, page 159-162 of ( Mengali, Quarta ). The
 *      hyperbolic case was taken from ( Noomen, R. ). The retrograde and near-pi cases are verified
 *      against values found with the Lambert routine available in the Keplerian_Toolbox from
 *      ESA/ACT. It is assumed that the first two test cases are sufficient to test the radial and
 *      tangential velocity components computations. These are therefore not tested in the last two
 *      test cases.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Astrodynamics/BasicAstrodynamics/orbitalElementConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Basics/testMacros.h"

#include "Tudat/Astrodynamics/MissionSegments/zeroRevolutionLambertTargeterIzzo.h"

namespace tudat
{
namespace unit_tests
{

using namespace unit_conversions;

//! Test the Izzo Lambert targeting algorithm code.
BOOST_AUTO_TEST_SUITE( test_zero_revolution_lambert_targeter_izzo )

//! Test hyperbolic case.
BOOST_AUTO_TEST_CASE( testHyperbolicCase )
// Copied and slightly adapted from unitTestLambertTargeterIzzo (rev 466)
{
    // Expected test result in meters.
    // Hyperbolic test case (results from excel file [1]).
    const double expectedValueOfSemiMajorAxisHyperbola = -1270129.3602e3;
    const double expectedValueOfRadialSpeedAtDepartureHyperbola = -0.74546e3;
    const double expectedValueOfRadialSpeedAtArrivalHyperbola = 0.69321e3;
    const double expectedValueOfTransverseSpeedAtDepartureHyperbola = 0.15674e3;
    const double expectedValueOfTransverseSpeedAtArrivalHyperbola = 0.10450e3;
    const Eigen::Vector3d expectedInertialVelocityAtDeparture( -745.457, 156.743, 0.0 );
    const Eigen::Vector3d expectedInertialVelocityAtArrival( 104.495, -693.209, 0.0 );

    // Tolerances.
    const double toleranceSemiMajorAxisHyperbola = 1.0e-7;
    const double toleranceVelocity = 1.0e-4;

    // Time conversions.
    const double timeOfFlightInDaysHyperbola = 100.0;
    const double timeOfFlightHyperbola = convertJulianDaysToSeconds( timeOfFlightInDaysHyperbola );

    // Central body gravitational parameter.
    const double earthGravitationalParameter = 398600.4418e9;

    // The starting point is twice as far as L1 and L2, which is not really
    // realistic, but it is not about the case, but about the verification.
    const Eigen::Vector3d positionAtDepartureHyperbola( convertAstronomicalUnitsToMeters( 0.02 ),
                                                        0.0, 0.0 ),
            positionAtArrivalHyperbola( 0.0, convertAstronomicalUnitsToMeters( -0.03 ), 0.0 );

    // Compute Lambert targeting algorithms.
    mission_segments::ZeroRevolutionLambertTargeterIzzo lambertTargeterHyperbola(
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

    // Test if the computed transverse and radial velocity components corresponds to the
    // expected values within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getRadialVelocityAtDeparture( ),
                                expectedValueOfRadialSpeedAtDepartureHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getRadialVelocityAtArrival( ),
                                expectedValueOfRadialSpeedAtArrivalHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getTransverseVelocityAtDeparture( ),
                                expectedValueOfTransverseSpeedAtDepartureHyperbola,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterHyperbola.getTransverseVelocityAtArrival( ),
                                expectedValueOfTransverseSpeedAtArrivalHyperbola,
                                toleranceVelocity );

    // Check that velocities match expected values within the defined tolerance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedInertialVelocityAtDeparture,
                                       lambertTargeterHyperbola.getInertialVelocityAtDeparture( ),
                                       toleranceVelocity );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedInertialVelocityAtArrival,
                                       lambertTargeterHyperbola.getInertialVelocityAtArrival( ),
                                       toleranceVelocity );

    // Test if the computed solution is anti-clockwise, if the z-component of the angular momentum
    // (h = r \times v) is positive.
    BOOST_CHECK_GT( positionDepartureHyperbola.cross( velocityDepartureHyperbola ).z( ), 0 );
}

//! Test elliptical case.
BOOST_AUTO_TEST_CASE( testEllipticalCase )
// Copied and slightly adapted from unitTestLambertTargeterIzzo (rev 466)
{
    // Elliptical test case (results from example 6.1 page 159-162 [2]).
    // Set canonical units for Earth (see page 29 [2]).
    double distanceUnit = 6.378136e6;
    double timeUnit = 806.78;

    // Set expected test result in meters.
    double expectedValueOfSemiMajorAxisEllipse = 5.4214 * distanceUnit;
    double expectedValueOfRadialSpeedAtDepartureEllipse = 2.73580e3;
    double expectedValueOfRadialSpeedAtArrivalEllipse = 2.97503e3;
    double expectedValueOfTransverseSpeedAtDepartureEllipse = 6.59430e3;
    double expectedValueOfTransverseSpeedAtArrivalEllipse = 3.29715e3;
    const Eigen::Vector3d expectedInertialVelocityAtDeparture( 2735.8, 6594.3, 0.0 );
    const Eigen::Vector3d expectedInertialVelocityAtArrival( -1367.9, 4225.03, 0.0 );

    // Tolerance in absolute units.
    double toleranceSemiMajorAxisEllipse = 1.0e4;
    double toleranceVelocity = 1.0e-2;

    // Time conversions.
    double timeOfFlightEllipse = 5.0 * timeUnit;

    // Central body gravitational parameter.
    const double earthGravitationalParameter = 398600.4418e9;

    // Elliptical orbit case.
    const Eigen::Vector3d positionAtDepartureEllipse( 2.0 * distanceUnit, 0.0, 0.0 ),
            positionAtArrivalEllipse( 2.0 * distanceUnit, 2.0 * sqrt( 3.0 ) * distanceUnit, 0.0 );

    // Compute Lambert targeting algorithms.
    mission_segments::ZeroRevolutionLambertTargeterIzzo lambertTargeterEllipse(
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

    // Test if the computed transverse and radial velocity components corresponds to the
    // expected values within the specified tolerance.
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getRadialVelocityAtDeparture( ),
                                expectedValueOfRadialSpeedAtDepartureEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getRadialVelocityAtArrival( ),
                                expectedValueOfRadialSpeedAtArrivalEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getTransverseVelocityAtDeparture( ),
                                expectedValueOfTransverseSpeedAtDepartureEllipse,
                                toleranceVelocity );
    BOOST_CHECK_CLOSE_FRACTION( lambertTargeterEllipse.getTransverseVelocityAtArrival( ),
                                expectedValueOfTransverseSpeedAtArrivalEllipse,
                                toleranceVelocity );

    // Check that velocities match expected values within the defined tolerance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedInertialVelocityAtDeparture,
                                       lambertTargeterEllipse.getInertialVelocityAtDeparture( ),
                                       toleranceVelocity );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedInertialVelocityAtArrival,
                                       lambertTargeterEllipse.getInertialVelocityAtArrival( ),
                                       toleranceVelocity );

    // Test if the computed solution is anti-clockwise, if the z-component of the angular momentum
    // (h = r \times v) is positive.
    BOOST_CHECK_GT( positionDepartureEllipse.cross( velocityDepartureEllipse ).z( ), 0 );
}

//! Test retrograde Earth-Mars transfer.
BOOST_AUTO_TEST_CASE ( testRetrograde )
// Copied and slightly adapted from unitTestLambertTargeterIzzo (rev 466).
{
    // Set tolerance.
    const double tolerance = 1.0e-9;

    // Set positions at departure and arrival.
    /* Values taken from http://ccar.colorado.edu/~rla/lambert_j2000.html for JDi = 2456036 and
     * JDf = 2456336.
     */
    const Eigen::Vector3d positionAtDeparture( -131798187443.90068, -72114797019.4148,
                                               2343782.3918863535 ),
            positionAtArrival( 202564770723.92966, -42405023055.01754, -5861543784.413235 );

    // Set time-of-flight, coherent with initial and final positions.
    const double timeOfFlight = convertJulianDaysToSeconds( 300.0 );

    // Set central body (the Sun) gravitational parameter. Value taken from keptoolbox.
    const double solarGravitationalParameter = 1.32712428e20;

    // Set expected values for inertial velocities. Values obtained with keptoolbox.
    const Eigen::Vector3d expectedInitialVelocity( -14157.8507230353, 28751.266655828,
                                                   1395.46037631136 ),
            expectedFinalVelocity( -6609.91626743654, -22363.5220239692, -716.519714631494 );

    // Compute Lambert solution.
    mission_segments::ZeroRevolutionLambertTargeterIzzo lambertTargeterRetrograde(
                positionAtDeparture, positionAtArrival,
                timeOfFlight, solarGravitationalParameter, true );

    // Retrieve inertial velocities.
    const Eigen::Vector3d initialVelocity =
            lambertTargeterRetrograde.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d finalVelocity = lambertTargeterRetrograde.getInertialVelocityAtArrival( );

    // Check that velocities match expected values within the defined tolerance.
    TUDAT_CHECK_MATRIX_CLOSE( initialVelocity, expectedInitialVelocity, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE( finalVelocity, expectedFinalVelocity, tolerance );
}

//! Test near-pi, circular, coplanar, Earth-Mars transfer.
BOOST_AUTO_TEST_CASE ( testNearPi )
// Copied and slightly adapted from unitTestLambertTargeterIzzo (rev 466).
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Set time-of-flight, coherent with initial and final positions.
    const double timeOfFlight = convertJulianDaysToSeconds( 300.0 );

    // Set central body (the Sun) gravitational parameter. Value taken from keptoolbox.
    const double solarGravitationalParameter = 1.32712428e20;

    // Set Keplerian elements at departure and arrival.
    Eigen::Matrix< double, 6, 1 > keplerianStateAtDeparture, keplerianStateAtArrival;
    keplerianStateAtDeparture << convertAstronomicalUnitsToMeters( 1.0 ), 0.0,
            0.0, 0.0, 0.0, 0.0;
    keplerianStateAtArrival << convertAstronomicalUnitsToMeters( 1.5 ), 0.0, 0.0,
            0.0, 0.0, convertDegreesToRadians( 179.999 );

    //  Convert to Cartesian elements.
    const Eigen::Matrix< double, 6, 1 > cartesianStateAtDeparture =
            orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianStateAtDeparture, solarGravitationalParameter );
    const Eigen::Matrix< double, 6, 1 > cartesianStateAtArrival =
            orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianStateAtArrival, solarGravitationalParameter );

    // Extract positions at departure and arrival.
    const Eigen::Vector3d positionAtDeparture = cartesianStateAtDeparture.head( 3 );
    const Eigen::Vector3d positionAtArrival = cartesianStateAtArrival.head( 3 );

    // Set expected values for inertial velocities. Values obtained with keptoolbox.
    const Eigen::Vector3d expectedInitialVelocity( 3160.36638344209, 32627.4771454454, 0.0 ),
            expectedFinalVelocity( 3159.89183582648, -21751.7065841264, 0.0 );

    // Compute Lambert solution.
    mission_segments::ZeroRevolutionLambertTargeterIzzo lambertTargeterRetrograde(
                positionAtDeparture, positionAtArrival, timeOfFlight, solarGravitationalParameter );

    // Retrieve inertial velocities.
    const Eigen::Vector3d initialVelocity =
            lambertTargeterRetrograde.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d finalVelocity = lambertTargeterRetrograde.getInertialVelocityAtArrival( );

    // Check that velocities match expected values within the defined tolerance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedInitialVelocity, initialVelocity, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedFinalVelocity, finalVelocity, tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
