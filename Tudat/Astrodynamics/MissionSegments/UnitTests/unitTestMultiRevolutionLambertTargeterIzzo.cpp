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
 *      Noomen, R., Lambert targeter Excel file.
 *      Izzo, D., Keplerian_Toolbox.
 *      JPL, HORIZONS web interface. http://ssd.jpl.nasa.gov/horizons.cgi?s_loc=1#top
 *      McConaghy, T.T., Notable Two-Synodic-Period Earth-Mars Cycler (Mar/Apr 2006). Journal of
 *          Spacecraft and Rockets, 43:2.
 *
 */

#define BOOST_TEST_MAIN

#include <boost/test/unit_test.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h"

namespace tudat
{
namespace unit_tests
{

//! Test the Izzo Lambert targeting algorithm code for multiple revolutions.
BOOST_AUTO_TEST_SUITE( test_multi_revolution_lambert_targeter_izzo )

//! Test elliptic multiple revolution case.
BOOST_AUTO_TEST_CASE( testEllipticCase )
// Tested using output from Keplerian_Toolbox PyKEP.
{
    // Required tolerance
    const double tolerance = 1.0e-6;

    // Define problem.
    const Eigen::Vector3d departurePosition( 4949101.422118526, 859402.44303969538,
                                             -151535.83799466802 );
    const Eigen::Vector3d arrivalPosition( 3648349.9884584765, 4281879.3154454567,
                                           -755010.85145052616 );
    const double timeOfFlight = 1.0307431655832210e+004;
    const double gravitationalParameter = 398600.4418e9; // From PyKEP toolbox

    // Expected test results in meters.
    int expectedMaximumNumberOfRevolutions = 4;

    // 0-rev solution.
    const Eigen::Vector3d expectedVelocityAtDeparture0Revolutions( 10096.683162831092,
                                                                   4333.9040463806396,
                                                                   -764.1842151784972 );
    const Eigen::Vector3d expectedVelocityAtArrival0Revolutions( -8110.7563754983421,
                                                                 -6018.4641067312887,
                                                                 1061.2176044421969 );
    const double expectedSemiMajorAxis0Revolutions = 10679740.431759536;

    // 1-rev solution (left branch).
    const Eigen::Vector3d expectedVelocityAtDepartureLeftBranch1Revolutions( 8918.2511158620255,
                                                                             4409.3440789101496,
                                                                             -777.48632833897398 );
    const Eigen::Vector3d expectedVelocityAtArrivalLeftBranch1Revolutions( -7506.6196898648195,
                                                                           -4929.4928888157147,
                                                                           869.20259750872378 );
    const double expectedSemiMajorAxisLeftBranch1Revolutions = 6750132.818154363;

    // 1-rev solution (right branch).
    const Eigen::Vector3d expectedVelocityAtDepartureRightBranch1Revolutions(
                -1265.9264854521089, 10660.067181950877, -1879.6574603427925 );
    const Eigen::Vector3d expectedVelocityAtArrivalRightBranch1Revolutions( -5584.601938281151,
                                                                            8204.5589196619731,
                                                                            -1446.6851023487009 );
    const double expectedSemiMajorAxisRightBranch1Revolutions = 9999999.9999999981;

    // 2-rev solution (left branch).
    const Eigen::Vector3d expectedVelocityAtDepartureLeftBranch2Revolutions( 7764.4367242290973,
                                                                             4541.6234715146611,
                                                                             -800.81075424687731 );
    const Eigen::Vector3d expectedVelocityAtArrivalLeftBranch2Revolutions( -6949.5242141064518,
                                                                           -3824.4260382745611,
                                                                           674.34949627178969 );
    const double expectedSemiMajorAxisLeftBranch2Revolutions = 5171347.6724622268;

    // 2-rev solution (right branch).
    const Eigen::Vector3d expectedVelocityAtDepartureRightBranch2Revolutions(
                -515.62716630712907, 9592.5976150608458, -1691.4337746149008 );
    const Eigen::Vector3d expectedVelocityAtArrivalRightBranch2Revolutions( -5368.5573571177865,
                                                                            6833.3233167245971,
                                                                            -1204.8992686428019 );
    const double expectedSemiMajorAxisRightBranch2Revolutions = 6278357.7897979086;

    // 3-rev solution (left branch).
    const Eigen::Vector3d expectedVelocityAtDepartureLeftBranch3Revolutions( 6500.2809323521278,
                                                                             4773.7410262238855,
                                                                             -841.73934183818676 );
    const Eigen::Vector3d expectedVelocityAtArrivalLeftBranch3Revolutions( -6390.5274190589034,
                                                                           -2555.7021700082605,
                                                                           450.63924722762863 );
    const double expectedSemiMajorAxisLeftBranch3Revolutions = 4291472.3727418669;

    // 3-rev solution (right branch).
    const Eigen::Vector3d expectedVelocityAtDepartureRightBranch3Revolutions(
                339.85968372598705, 8524.895950642951, -1503.1691638306909 );
    const Eigen::Vector3d expectedVelocityAtArrivalRightBranch3Revolutions( -5210.2694492213532,
                                                                            5369.2089601340958,
                                                                            -946.73640473328203 );
    const double expectedSemiMajorAxisRightBranch3Revolutions = 4768867.1058038883;

    // 4-rev solution (left branch).
    const Eigen::Vector3d expectedVelocityAtDepartureLeftBranch4Revolutions( 4812.3627329648789,
                                                                             5280.0203047688119,
                                                                             -931.01003841927343 );
    const Eigen::Vector3d expectedVelocityAtArrivalLeftBranch4Revolutions( -5759.8462467502432,
                                                                           -731.11592991390376,
                                                                           128.91546446958043 );
    const double expectedSemiMajorAxisLeftBranch4Revolutions = 3734713.1475614533;

    // 4-rev solution (right branch).
    const Eigen::Vector3d expectedVelocityAtDepartureRightBranch4Revolutions(
                1618.0817850471631, 7225.3760295124985, -1274.0287397672555 );
    const Eigen::Vector3d expectedVelocityAtArrivalRightBranch4Revolutions( -5148.0510872930045,
                                                                            3378.294822960549,
                                                                            -595.68452607567178 );
    const double expectedSemiMajorAxisRightBranch4Revolutions = 3900758.7239032765;


    // Constructing targeter and calculating 0-rev solution with multi-revolution class.
    mission_segments::MultiRevolutionLambertTargeterIzzo lambertTargeterEllipse(
                departurePosition, arrivalPosition, timeOfFlight, gravitationalParameter );

    // Check whether maximum number of revolutions is calculated correctly.
    const int calculatedMaximumNumberOfRevolutions
            = lambertTargeterEllipse.getMaximumNumberOfRevolutions( );
    BOOST_CHECK_EQUAL( calculatedMaximumNumberOfRevolutions, expectedMaximumNumberOfRevolutions );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDeparture0Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrival0Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxis0Revolutions = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 0 rev solutions is correct (should have been calculated upon construction)
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDeparture0Revolutions,
                                       expectedVelocityAtDeparture0Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrival0Revolutions,
                                       expectedVelocityAtArrival0Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxis0Revolutions,
                                expectedSemiMajorAxis0Revolutions,
                                tolerance );

    // Recalculating for left branch, 1-rev.
    int numberOfRevolutions = 1;
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, false );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureLeftBranch1Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalLeftBranch1Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisLeftBranch1Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 1-rev, left branch, solution is correct.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureLeftBranch1Revolutions,
                                       expectedVelocityAtDepartureLeftBranch1Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalLeftBranch1Revolutions,
                                       expectedVelocityAtArrivalLeftBranch1Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisLeftBranch1Revolutions,
                                expectedSemiMajorAxisLeftBranch1Revolutions,
                                tolerance );

    // Recalculating for right branch, 1-rev.
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, true );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureRightBranch1Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalRightBranch1Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisRightBranch1Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 1-rev, right branch, solution is correct.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureRightBranch1Revolutions,
                                       expectedVelocityAtDepartureRightBranch1Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalRightBranch1Revolutions,
                                       expectedVelocityAtArrivalRightBranch1Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisRightBranch1Revolutions,
                                expectedSemiMajorAxisRightBranch1Revolutions,
                                tolerance );

    // Recalculating for left branch, 2-rev
    numberOfRevolutions = 2;
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, false );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureLeftBranch2Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalLeftBranch2Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisLeftBranch2Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 2 rev, left branch, solution is correct.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureLeftBranch2Revolutions,
                                       expectedVelocityAtDepartureLeftBranch2Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalLeftBranch2Revolutions,
                                       expectedVelocityAtArrivalLeftBranch2Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisLeftBranch2Revolutions,
                                expectedSemiMajorAxisLeftBranch2Revolutions,
                                tolerance );

    // Recalculating for right branch, 2-rev
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, true );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureRightBranch2Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalRightBranch2Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisRightBranch2Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 2-rev, right branch, solution is correct
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureRightBranch2Revolutions,
                                       expectedVelocityAtDepartureRightBranch2Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalRightBranch2Revolutions,
                                       expectedVelocityAtArrivalRightBranch2Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisRightBranch2Revolutions,
                                expectedSemiMajorAxisRightBranch2Revolutions,
                                tolerance );


    // Recalculating for left branch, 3-rev.
    numberOfRevolutions = 3;
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, false );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureLeftBranch3Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalLeftBranch3Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisLeftBranch3Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 3-rev, left branch, solution is correct.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureLeftBranch3Revolutions,
                                       expectedVelocityAtDepartureLeftBranch3Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalLeftBranch3Revolutions,
                                       expectedVelocityAtArrivalLeftBranch3Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisLeftBranch3Revolutions,
                                expectedSemiMajorAxisLeftBranch3Revolutions,
                                tolerance );

    // Recalculating for right branch, 3-rev.
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, true );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureRightBranch3Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalRightBranch3Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisRightBranch3Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 3-rev, right branch, solution is correct.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureRightBranch3Revolutions,
                                       expectedVelocityAtDepartureRightBranch3Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalRightBranch3Revolutions,
                                       expectedVelocityAtArrivalRightBranch3Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisRightBranch3Revolutions,
                                expectedSemiMajorAxisRightBranch3Revolutions,
                                tolerance );

    // Recalculating for left branch, 4-rev
    numberOfRevolutions = 4;
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, false );

    // Extracting velocities from targeter object
    const Eigen::Vector3d calculatedVelocityAtDepartureLeftBranch4Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalLeftBranch4Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisLeftBranch4Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 4-rev, left branch, solution is correct
    // (should have been calculated upon construction)
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureLeftBranch4Revolutions,
                                       expectedVelocityAtDepartureLeftBranch4Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalLeftBranch4Revolutions,
                                       expectedVelocityAtArrivalLeftBranch4Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisLeftBranch4Revolutions,
                                expectedSemiMajorAxisLeftBranch4Revolutions,
                                tolerance );

    // Recalculating for right branch, 4-rev.
    lambertTargeterEllipse.computeForRevolutionsAndBranch( numberOfRevolutions, true );

    // Extracting velocities from targeter object.
    const Eigen::Vector3d calculatedVelocityAtDepartureRightBranch4Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtDeparture( );
    const Eigen::Vector3d calculatedVelocityAtArrivalRightBranch4Revolutions
            = lambertTargeterEllipse.getInertialVelocityAtArrival( );
    const double calculatedSemiMajorAxisRightBranch4Revolutions
            = lambertTargeterEllipse.getSemiMajorAxis( );

    // Check whether 4-rev, right branch, solution is correct.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtDepartureRightBranch4Revolutions,
                                       expectedVelocityAtDepartureRightBranch4Revolutions,
                                       tolerance );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedVelocityAtArrivalRightBranch4Revolutions,
                                       expectedVelocityAtArrivalRightBranch4Revolutions,
                                       tolerance );
    BOOST_CHECK_CLOSE_FRACTION( calculatedSemiMajorAxisRightBranch4Revolutions,
                                expectedSemiMajorAxisRightBranch4Revolutions,
                                tolerance );
}

//! Test hyperbolic case.
BOOST_AUTO_TEST_CASE( testHyperbolicCase )
// Copied and slightly adapted from unitTestLambertTargeterIzzo (rev 466). While also being tested
// in the zero revolution unit test, its application to the multi revolution case proves that the
// multi revolution case deals properly with (obvious) zero revolution cases.
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
    const double timeOfFlightHyperbola = unit_conversions::
            convertJulianDaysToSeconds( timeOfFlightInDaysHyperbola );

    // Central body gravitational parameter.
    const double earthGravitationalParameter = 398600.4418e9;

    // The starting point is twice as far as L1 and L2, which is not really
    // realistic, but it is not about the case, but about the verification.
    using unit_conversions::convertAstronomicalUnitsToMeters;
    const Eigen::Vector3d positionAtDepartureHyperbola( convertAstronomicalUnitsToMeters( 0.02 ),
                                                        0.0, 0.0 ),
            positionAtArrivalHyperbola( 0.0, convertAstronomicalUnitsToMeters( -0.03 ), 0.0 );

    // Compute Lambert targeting algorithms.
    mission_segments::MultiRevolutionLambertTargeterIzzo lambertTargeterHyperbola(
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

//! Test hyperbolic case.
BOOST_AUTO_TEST_CASE( testMaximumNumberOfRevolutions )
// Tested using reference output from Keplerian_Toolbox PyKEP. Positions are derived from cycler
// itineraries in McConaghy (2006), transfers Earth-1 to Mars-2 and Mars-2 to Earth-3, table 2.
{
    /* The problem was that for certain combinations of parameters (but not all), the computed
       maximum number of revolutions was one too low. As a result, not all possible solutions could
       be found resulting in a wrong result if one was looking for the best solution amongst all
       possible solutions. This test shows that this is no longer the case, as it uses two of those
       'border' cases.
      */

    // Expected results
    int expectedMaximumNumberOfRevolutions1 = 1;
    int expectedMaximumNumberOfRevolutions2 = 1;

    // Departure position (Cartesian, m)
    Eigen::Vector3d departurePosition1( -231427903073.9245, 81487782204.595016,
                                        43627778250.279877 );
    Eigen::Vector3d departurePosition2( 59996317860.871559, -128191631368.53741,
                                        -55574089046.065445 );

    // Arrival position (Cartesian, m)
    Eigen::Vector3d arrivalPosition1 ( 59996317860.871559, -128191631368.53741,
                                       -55574089046.065445 );
    Eigen::Vector3d arrivalPosition2 ( -45319765987.585777, 128402601206.42529,
                                       55665242520.139145 );

    // Time of flight (s)
    double timeOfFlight1 = 68169600.998501450;
    double timeOfFlight2 = 46828800.000441134;

    // Sun gravitational parameter
    const double sunGravitationalParameter = 1.32712440018e20; // m/s, HORIZONS

    // Constructing Lambert problem
    mission_segments::MultiRevolutionLambertTargeterIzzo lambertTargeter1(
                departurePosition1, arrivalPosition1, timeOfFlight1, sunGravitationalParameter );
    mission_segments::MultiRevolutionLambertTargeterIzzo lambertTargeter2(
                departurePosition2, arrivalPosition2, timeOfFlight2, sunGravitationalParameter );

    // Extracting number of maximum revolutions
    int computedMaximumNumberOfRevolutions1 = lambertTargeter1.getMaximumNumberOfRevolutions();
    int computedMaximumNumberOfRevolutions2 = lambertTargeter2.getMaximumNumberOfRevolutions();

    // Comparing
    BOOST_CHECK_EQUAL( computedMaximumNumberOfRevolutions1, expectedMaximumNumberOfRevolutions1 );
    BOOST_CHECK_EQUAL( computedMaximumNumberOfRevolutions2, expectedMaximumNumberOfRevolutions2 );
}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
