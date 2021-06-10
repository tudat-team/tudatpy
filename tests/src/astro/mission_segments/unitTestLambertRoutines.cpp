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
 *      Noomen, R., Lambert targeter Excel file.
 *      Mengali, G., and A.A. Quarta, Fondamenti di Meccanica del volo Spaziale.
 *      Izzo, D., Keplerian_Toolbox.
 *
 *    Notes
 *      The elliptical case was taken from Example 6.1, page 159-162 of ( Mengali, Quarta ). The
 *      hyperbolic case was taken from ( Noomen, R. ). The retrograde and near-pi cases are
 *      verified against values found with the Lambert routine available in the Keplerian_Toolbox
 *      from ESA/ACT.
 *
 *      DISCLAIMER: At the moment, the Gooding Lambert targeter only converges for about half of
 *      the cases. This is not evident from the tests below, but it was observed during simulations
 *      carried out by the author. The reason might very well be an erroneous definition of the
 *      starters.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/basics/basicTypedefs.h"

#include "tudat/astro/mission_segments/lambertRoutines.h"

namespace tudat
{
namespace unit_tests
{

using namespace unit_conversions;

//! Test the Lambert targeting routines.
BOOST_AUTO_TEST_SUITE( test_lambert_routines )

//! Test the Izzo time-of-flight computation.
BOOST_AUTO_TEST_CASE( testIzzoTimeOfFlightComputation )
{
    // Set tolerance.
    const double tolerance = 1.0e-7;

    // Set input values (taken from testEllipticalCase in test_lambert_targeter_izzo) [-].
    const double testXParameter = -0.5, testSemiPerimeter = 2.36603, testChord = 1.73205,
            testSemiMajorAxisOfTheMinimumEnergyEllipse = 1.18301;
    const bool testIsLongway = false;

    // Set expected time-of-flight [-].
    const double expectedTimeOfFlight = 9.759646;

    // Check that returned value is equal to expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedTimeOfFlight,
                                mission_segments::computeTimeOfFlightIzzo(
                                    testXParameter, testSemiPerimeter, testChord, testIsLongway,
                                    testSemiMajorAxisOfTheMinimumEnergyEllipse ),
                                tolerance );
}

//! Test the Izzo Lambert routine for an elliptical transfer.
BOOST_AUTO_TEST_CASE( testSolveLambertProblemIzzoElliptical )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Set canonical units for Earth (see page 29 [1]).
    const double distanceUnit = 6.378136e6;
    const double timeUnit = 806.78;

    // Set expected inertial vectors.
    Eigen::Vector3d expectedInertialVelocityAtDeparture( 2735.8, 6594.3, 0.0 ),
            expectedInertialVelocityAtArrival( -1367.9, 4225.03, 0.0 );

    // Time conversions.
    const double testTimeOfFlight = 5.0 * timeUnit;

    // Set central body graviational parameter.
    const double testGravitationalParameter = 398600.4418e9;

    // Set position at departure and arrival.
    const Eigen::Vector3d testCartesianPositionAtDeparture( 2.0 * distanceUnit, 0.0, 0.0 ),
            testCartesianPositionAtArrival( 2.0 * distanceUnit, 2.0 * sqrt( 3.0 ) * distanceUnit,
                                            0.0 );

    // Declare velocity vectors.
    Eigen::Vector3d testInertialVelocityAtDeparture, testInertialVelocityAtArrival;

    // Solve Lambert problem.
    mission_segments::solveLambertProblemIzzo( testCartesianPositionAtDeparture,
                                                      testCartesianPositionAtArrival,
                                                      testTimeOfFlight,
                                                      testGravitationalParameter,
                                                      testInertialVelocityAtDeparture,
                                                      testInertialVelocityAtArrival );

    // Check that returned vectors are equal to expected vectors.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                testInertialVelocityAtDeparture.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                testInertialVelocityAtDeparture.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtDeparture.z( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                testInertialVelocityAtArrival.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                testInertialVelocityAtArrival.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtArrival.z( ), tolerance );

    Eigen::Vector6d initialCartesianState;
    initialCartesianState << testCartesianPositionAtDeparture, testInertialVelocityAtDeparture;

    Eigen::Vector6d constantKeplerianState_ = orbital_element_conversions::convertCartesianToKeplerianElements(
                initialCartesianState, testGravitationalParameter );

    Eigen::Vector6d finalCartesianState = orbital_element_conversions::convertKeplerianToCartesianElements(
                orbital_element_conversions::propagateKeplerOrbit(
                    constantKeplerianState_, testTimeOfFlight, testGravitationalParameter ),
                testGravitationalParameter );

    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( finalCartesianState.segment( 0, 3 ) ), testCartesianPositionAtArrival, 1.0E-14 )
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( ( finalCartesianState.segment( 3, 3 ) ), testInertialVelocityAtArrival, 1.0E-14 )


}

//! Test the Izzo Lambert routine for a hyperbolic transfer.
BOOST_AUTO_TEST_CASE( testSolveLambertProblemIzzoHyperbolic )
{
    // Set tolerance.
    const double tolerance = 1.0e-5;

    // Set expected inertial vectors.
    Eigen::Vector3d expectedInertialVelocityAtDeparture( -745.457, 156.743, 0.0 ),
            expectedInertialVelocityAtArrival( 104.495, -693.209, 0.0 );

    // Time conversions.
    const double testTimeOfFlightInDaysHyperbola = 100.0;
    const double testTimeOfFlightHyperbola = convertJulianDaysToSeconds(
            testTimeOfFlightInDaysHyperbola );

    // Set central body graviational parameter.
    const double testGravitationalParameter = 398600.4418e9;

    // Set position at departure and arrival.
    const Eigen::Vector3d testCartesianPositionAtDeparture(
                convertAstronomicalUnitsToMeters( 0.02 ), 0.0, 0.0 ),
            testCartesianPositionAtArrival(
                0.0, convertAstronomicalUnitsToMeters( -0.03 ), 0.0 );

    // Declare velocity vectors.
    Eigen::Vector3d testInertialVelocityAtDeparture, testInertialVelocityAtArrival;

    // Solve Lambert problem.
    mission_segments::solveLambertProblemIzzo( testCartesianPositionAtDeparture,
                                                      testCartesianPositionAtArrival,
                                                      testTimeOfFlightHyperbola,
                                                      testGravitationalParameter,
                                                      testInertialVelocityAtDeparture,
                                                      testInertialVelocityAtArrival );

    // Check that returned vectors are equal to expected vectors.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                testInertialVelocityAtDeparture.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                testInertialVelocityAtDeparture.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtDeparture.z( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                testInertialVelocityAtArrival.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                testInertialVelocityAtArrival.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtArrival.z( ), tolerance );
}

//! Test the Izzo Lambert routine for retrograde orbits.
BOOST_AUTO_TEST_CASE( testSolveLambertProblemIzzoRetrograde )
{
    // Set tolerance.
    const double tolerance = 1.0e-9;

    // Set positions at departure and arrival.
    /* Values taken from http://ccar.colorado.edu/~rla/lambert_j2000.html for JDi = 2456036 and
     * JDf = 2456336.
     */
    const Eigen::Vector3d
            testPositionAtDeparture( -131798187443.90068, -72114797019.4148, 2343782.3918863535 ),
            testPositionAtArrival( 202564770723.92966, -42405023055.01754, -5861543784.413235);

    // Set time-of-flight, coherent with initial and final positions.
    const double testTimeOfFlight = convertJulianDaysToSeconds( 300.0 );

    // Set central body (the Sun) gravitational parameter. Value taken from keptoolbox.
    const double testSolarGravitationalParameter = 1.32712428e20;

    // Set expected values for inertial velocities. Values obtained with keptoolbox.
    const Eigen::Vector3d
            expectedInitialVelocity( -14157.8507230353, 28751.266655828, 1395.46037631136 ),
            expectedFinalVelocity( -6609.91626743654, -22363.5220239692, -716.519714631494 );

    // Declare initial and final velocity vectors.
    Eigen::Vector3d initialVelocity, finalVelocity;

    // Compute Lambert solution.
    mission_segments::solveLambertProblemIzzo(
                testPositionAtDeparture, testPositionAtArrival,
                testTimeOfFlight, testSolarGravitationalParameter,
                initialVelocity, finalVelocity, true );

    // Check that velocities match expected values within the defined tolerance.
    TUDAT_CHECK_MATRIX_CLOSE( expectedInitialVelocity, initialVelocity, tolerance );
    TUDAT_CHECK_MATRIX_CLOSE( expectedFinalVelocity, finalVelocity, tolerance );
}

//! Test the Izzo Lambert routine for near-pi transfers.
BOOST_AUTO_TEST_CASE( testSolveLambertProblemIzzoNearPi )
{
    // Set tolerance.
    const double tolerance = 1.0e-9;

    // Set time-of-flight, coherent with initial and final positions.
    const double testTimeOfFlight = convertJulianDaysToSeconds( 300.0 );

    // Set central body (the Sun) gravitational parameter. Value taken from keptoolbox.
    const double testSolarGravitationalParameter = 1.32712428e20;

    // Set Keplerian elements at departure and arrival.
    Eigen::Matrix< double, 6, 1 > keplerianStateAtDeparture, keplerianStateAtArrival;
    keplerianStateAtDeparture << convertAstronomicalUnitsToMeters( 1.0 ),
                                    0.0, 0.0, 0.0, 0.0, 0.0;
    keplerianStateAtArrival << convertAstronomicalUnitsToMeters( 1.5 ),
            0.0, 0.0, 0.0, 0.0, convertDegreesToRadians( 179.999 );

    //  Convert to Cartesian elements.
    const Eigen::Matrix< double, 6, 1 > cartesianStateAtDeparture =
            orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianStateAtDeparture, testSolarGravitationalParameter );
    const Eigen::Matrix< double, 6, 1 > cartesianStateAtArrival =
            orbital_element_conversions::convertKeplerianToCartesianElements(
                keplerianStateAtArrival, testSolarGravitationalParameter );

    // Extract positions at departure and arrival.
    const Eigen::Vector3d positionAtDeparture = cartesianStateAtDeparture.head( 3 );
    const Eigen::Vector3d positionAtArrival = cartesianStateAtArrival.head( 3 );

    // Set expected values for inertial velocities. Values obtained with keptoolbox.
    const Eigen::Vector3d expectedInitialVelocity( 3160.36638344209, 32627.4771454454, 0.0 ),
            expectedFinalVelocity( 3159.89183582648, -21751.7065841264, 0.0 );

    // Declare initial and final velocity vectors.
    Eigen::Vector3d initialVelocity, finalVelocity;

    // Compute Lambert solution.
    mission_segments::solveLambertProblemIzzo( positionAtDeparture, positionAtArrival,
                                               testTimeOfFlight, testSolarGravitationalParameter,
                                               initialVelocity, finalVelocity );

    // Check that velocities match expected values within the defined tolerance.
    BOOST_CHECK_CLOSE_FRACTION( expectedInitialVelocity.x( ), initialVelocity.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInitialVelocity.y( ), initialVelocity.y( ), tolerance );
    BOOST_CHECK_SMALL( initialVelocity.z( ), tolerance );

    BOOST_CHECK_CLOSE_FRACTION( expectedFinalVelocity.x( ), finalVelocity.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedFinalVelocity.y( ), finalVelocity.y( ), tolerance );
    BOOST_CHECK_SMALL( finalVelocity.z( ), tolerance );
}

//! Test the positive Gooding Lambert function.
BOOST_AUTO_TEST_CASE( testLambertFunctionPositiveGooding )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
    double testXParameter = 1.09806;
    const double testQParameter = -0.402543;
    const double testNormalizedTimeOfFlight = 0.944749;

    // Create expected value.
    const double expectedLambertFunctionPositiveValue = -0.4004214;

    // Create LambertFunctionsGooding object.
    mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                testQParameter, testNormalizedTimeOfFlight);

    // Check that the result is equal to the expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedLambertFunctionPositiveValue,
                                testLambertFunctionsGooding.lambertFunctionPositiveGooding(
                                    testXParameter ), tolerance );
}

//! Test the negative Gooding Lambert function.
BOOST_AUTO_TEST_CASE( testLambertFunctionNegativeGooding )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
    double testXParameter = 0.434564;
    const double testQParameter = -0.402543;
    const double testNormalizedTimeOfFlight = 0.944749;

    // Create expected value.
    const double expectedLambertFunctionNegativeValue = -1.1439925;

    // Create LambertFunctionsGooding object.
    mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                testQParameter, testNormalizedTimeOfFlight);

    // Check that the result is equal to the expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedLambertFunctionNegativeValue,
                                testLambertFunctionsGooding.lambertFunctionNegativeGooding(
                                    testXParameter ), tolerance );
}

//! Test the Gooding Lambert function.
BOOST_AUTO_TEST_CASE( testLambertFunctionGooding )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Test 1: Test positive case.
    {
        // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
        double testXParameter = 1.09806;
        const double testQParameter = -0.402543;
        const double testNormalizedTimeOfFlight = 0.944749;

        // Create expected value.
        const double expectedLambertFunctionValue = -0.4004214;

        // Create LambertFunctionsGooding object.
        mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                    testQParameter, testNormalizedTimeOfFlight );

        // Check that the result is equal to the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedLambertFunctionValue,
                                    testLambertFunctionsGooding.computeLambertFunctionGooding(
                                        testXParameter ), tolerance );
    }

    // Test 2: Test negative case.
    {
        // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
        double testXParameter = 0.434564;
        const double testQParameter = -0.402543;
        const double testNormalizedTimeOfFlight = 0.944749;

        // Create expected value.
        const double expectedLambertFunctionValue = -1.1439925;

        // Create LambertFunctionsGooding object.
        mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                    testQParameter, testNormalizedTimeOfFlight );

        // Check that the result is equal to the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedLambertFunctionValue,
                                    testLambertFunctionsGooding.computeLambertFunctionGooding(
                                        testXParameter ), tolerance );
    }
}

//! Test the positive Gooding Lambert first derivative function.
BOOST_AUTO_TEST_CASE( testLambertFirstDerivativeFunctionPositiveGooding )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
    double testXParameter = 1.09806;
    const double testQParameter = -0.402543;
    const double testNormalizedTimeOfFlight = 0.944749;

    // Create expected value.
    const double expectedLambertFirstDerivativeFunctionPositiveValue = 0.7261451;

    // Create LambertFunctionsGooding object.
    mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                testQParameter, testNormalizedTimeOfFlight);

    // Check that the result is equal to the expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedLambertFirstDerivativeFunctionPositiveValue,
              testLambertFunctionsGooding.lambertFirstDerivativeFunctionPositiveGooding(
                                    testXParameter ), tolerance );
}

//! Test the negative Gooding Lambert first derivative function.
BOOST_AUTO_TEST_CASE( testLambertFirstDerivativeFunctionNegativeGooding )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
    double testXParameter = 0.434564;
    const double testQParameter = -0.402543;
    const double testNormalizedTimeOfFlight = 0.944749;

    // Create expected value.
    const double expectedLambertFirstDerivativeFunctionNegativeValue = 1.72419;

    // Create LambertFunctionsGooding object.
    mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                testQParameter, testNormalizedTimeOfFlight);

    // Check that the result is equal to the expected value.
    BOOST_CHECK_CLOSE_FRACTION( expectedLambertFirstDerivativeFunctionNegativeValue,
              testLambertFunctionsGooding.lambertFirstDerivativeFunctionNegativeGooding(
                                    testXParameter ), tolerance );
}

//! Test the Gooding Lambert first derivative function.
BOOST_AUTO_TEST_CASE( testLambertFirstDerivativeFunctionGooding )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Test 1: Test the positive case.
    {
        // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
        double testXParameter = 1.09806;
        const double testQParameter = -0.402543;
        const double testNormalizedTimeOfFlight = 0.944749;

        // Create expected value.
        const double expectedLambertFirstDerivativeFunctionValue = 0.7261451;

        // Create LambertFunctionsGooding object.
        mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                    testQParameter, testNormalizedTimeOfFlight);

        // Check that the result is equal to the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedLambertFirstDerivativeFunctionValue,
                  testLambertFunctionsGooding.computeFirstDerivativeLambertFunctionGooding(
                                        testXParameter ), tolerance );
    }

    // Test 2: Test the negative case.
    {
        // Create input parameters. Values taken from validated test_lambert_targeter_gooding.
        double testXParameter = 0.434564;
        const double testQParameter = -0.402543;
        const double testNormalizedTimeOfFlight = 0.944749;

        // Create expected value.
        const double expectedLambertFirstDerivativeFunctionValue = 1.72419;

        // Create LambertFunctionsGooding object.
        mission_segments::LambertFunctionsGooding testLambertFunctionsGooding(
                    testQParameter, testNormalizedTimeOfFlight);

        // Check that the result is equal to the expected value.
        BOOST_CHECK_CLOSE_FRACTION( expectedLambertFirstDerivativeFunctionValue,
                  testLambertFunctionsGooding.computeFirstDerivativeLambertFunctionGooding(
                                        testXParameter ), tolerance );
    }
}

//! Test the Gooding Lambert routine for a hyperbolic transfer.
BOOST_AUTO_TEST_CASE( testsolveLambertProblemGoodingHyperbolic )
{
    // Set tolerance.
    const double tolerance = 1.0e-5;

    // Set expected inertial vectors.
    Eigen::Vector3d expectedInertialVelocityAtDeparture( -745.457, 156.743, 0.0 ),
            expectedInertialVelocityAtArrival( 104.495, -693.209, 0.0 );

    // Time conversions.
    const double testTimeOfFlightInDaysHyperbola = 100.0;
    const double testTimeOfFlightHyperbola = convertJulianDaysToSeconds(
            testTimeOfFlightInDaysHyperbola );

    // Set central body graviational parameter.
    const double testGravitationalParameter = 398600.4418e9;

    // Set position at departure and arrival.
    const Eigen::Vector3d testCartesianPositionAtDeparture(
                convertAstronomicalUnitsToMeters( 0.02 ), 0.0, 0.0 ),
            testCartesianPositionAtArrival(
                0.0, convertAstronomicalUnitsToMeters( -0.03 ), 0.0 );

    // Declare velocity vectors.
    Eigen::Vector3d testInertialVelocityAtDeparture, testInertialVelocityAtArrival;

    // Solve Lambert problem.
    mission_segments::solveLambertProblemGooding( testCartesianPositionAtDeparture,
                                                        testCartesianPositionAtArrival,
                                                        testTimeOfFlightHyperbola,
                                                        testGravitationalParameter,
                                                        testInertialVelocityAtDeparture,
                                                        testInertialVelocityAtArrival);

    // Check that returned vectors are equal to expected vectors.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                testInertialVelocityAtDeparture.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                testInertialVelocityAtDeparture.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtDeparture.z( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                testInertialVelocityAtArrival.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                testInertialVelocityAtArrival.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtArrival.z( ), tolerance );
}

//! Test the Gooding Lambert routine for an elliptical transfer.
BOOST_AUTO_TEST_CASE( testsolveLambertProblemGoodingElliptical )
{
    // Set tolerance.
    const double tolerance = 1.0e-6;

    // Set canonical units for Earth (see page 29 [1]).
    const double distanceUnit = 6.378136e6;
    const double timeUnit = 806.78;

    // Set expected inertial vectors.
    const Eigen::Vector3d expectedInertialVelocityAtDeparture( 2735.8, 6594.3, 0.0 ),
            expectedInertialVelocityAtArrival( -1367.9, 4225.03, 0.0 );

    // Time conversions.
    const double testTimeOfFlight = 5.0 * timeUnit;

    // Set central body graviational parameter.
    const double testGravitationalParameter = 398600.4418e9;

    // Set position at departure and arrival.
    const Eigen::Vector3d testCartesianPositionAtDeparture( 2.0 * distanceUnit, 0.0, 0.0 ),
            testCartesianPositionAtArrival( 2.0 * distanceUnit, 2.0 * sqrt( 3.0 ) * distanceUnit,
                                            0.0 );

    // Declare velocity vectors.
    Eigen::Vector3d testInertialVelocityAtDeparture, testInertialVelocityAtArrival;

    // Solve Lambert problem.
    mission_segments::solveLambertProblemGooding( testCartesianPositionAtDeparture,
                                                         testCartesianPositionAtArrival,
                                                         testTimeOfFlight,
                                                         testGravitationalParameter,
                                                         testInertialVelocityAtDeparture,
                                                         testInertialVelocityAtArrival );

    // Check that returned vectors are equal to expected vectors.
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.x( ),
                                testInertialVelocityAtDeparture.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtDeparture.y( ),
                                testInertialVelocityAtDeparture.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtDeparture.z( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.x( ),
                                testInertialVelocityAtArrival.x( ), tolerance );
    BOOST_CHECK_CLOSE_FRACTION( expectedInertialVelocityAtArrival.y( ),
                                testInertialVelocityAtArrival.y( ), tolerance );
    BOOST_CHECK_SMALL( testInertialVelocityAtArrival.z( ), tolerance );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
