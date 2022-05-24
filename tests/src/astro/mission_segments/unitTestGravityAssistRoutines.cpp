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
 *      Izzo, D. and Vinko, T. ACT - Informatics - GTOP Database, ESA Advanced Concept Team, last
 *          accessed on 2012-01-12. http://www.esa.int/gsp/ACT/inf/op/globopt.htm.
 *      Musegaas, P. Gravity Assist calculation Verification.xlsx, last accessed: 3 December 2012,
 *          http://tudat.tudelft.nl/projects/tudat/wiki/Unit_tests, 2012.
 *
 *    Notes
 *      Three main functions are tested in these unit tests.
 *        Regarding the deltaV calculation gravity assist method:
 *          There is a complicated if-statement in this method. Hence many unit test are performed
 *          to test the functionality. Also various limit cases failed previously, hence many tests
 *          for this are also included:
 *              Case 1: required bending angle > maximum bending angle:
 *                  Two tests were written. In the first one no velocity effect is needed. This
 *                  test has a low accuracy, which should be replaced one day (it still relies on
 *                  hand calculator calculations done in 2011). In the second one a combination of
 *                  bending-effect deltaV and velocity-effect deltaV is calculated. This test has
 *                  been calculated using Tudat, and was verified using Excel.
 *                  Could definitely be improved.
 *              Case 2: no assist is required:
 *                  One test was written.
 *              Case 3: velocity effect deltaV only, using eccentricity iteration scheme:
 *                  Four tests were written. The first one calculates a case from Cassini-1 of GTOP
 *                  with high precision. The other three test limit cases: low incoming, high
 *                  outgoing velocity; high incoming, low outgoing velocity; low incoming, low
 *                  outgoing velocity. These tests were calculated using Tudat, but verified in
 *                  Excel to be exactly correct.
 *              Case 4: velocity effect deltaV only, using pericenter radius iteration scheme:
 *                  The same four tests as for case 3 were used.
 *        Regarding the unpowered gravity assist propagator:
 *          One test was written, based on GTOP. This should be a satisfactory test.
 *        Regarding the powered gravity assist propagator:
 *          Two tests were written. The first one is similar to the unpowered gravity assist
 *          propagator. The second one is reverse engineered from the Cassini-1 test, similar to
 *          the one in the deltaV calculation test.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <cmath>
#include <Eigen/Dense>

#include <boost/test/tools/floating_point_comparison.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/basics/testMacros.h"
#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/astro/mission_segments/gravityAssist.h"

namespace tudat
{
namespace unit_tests
{

//! Test of gravity assist code.
BOOST_AUTO_TEST_SUITE( test_gravity_assist )

//! Test bending angle Delta-V computation.
BOOST_AUTO_TEST_CASE( testBendingAngleDeltaV )
{
    // Tolerance, determined primarily by the accuracy of the hand calculations for this test case.
    const double velocityTolerance = 0.0002;

    // In the first test case, the incoming and outgoing inertial velocities are defined such that
    // the hyperbolic excess velocities are equal. In that way, a delta-V is only needed to rotate
    // the velocity vector, which has been calculated by hand.
    // Expected delta-V for a powered swing-by around Mars.
    const double expectedDeltaV = 3.652e3;

    // Define swingby body gravitational parameter.
    const double marsGravitationalParameter = 4.2828018915e13;

    // Define Sun gravitational parameter.
    const double gravitationalParameterSun = 1.32712440018e20;

    // Define planet-Sun distance.
    const double distanceMarsToSun = unit_conversions::
            convertAstronomicalUnitsToMeters( 1.5 );

    // Define smallest periapsis distance factor.
    const double marsSmallestPeriapsisDistance = 3656248.0;

    // Define planet heliocentric velocity vector. The orbit is considered to be circular.
    const Eigen::Vector3d marsVelocity( 0.0,
                                        std::sqrt( gravitationalParameterSun / distanceMarsToSun ),
                                        0.0 );

    // Define satellite incoming vector.
    using mathematical_constants::PI;
    const Eigen::Vector3d incomingVelocity( -25.0e3 * std::sin( PI / 6.0 ),
                                            25.0e3 * std::cos( PI / 6.0 ),
                                            0.0 );

    // Define satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( incomingVelocity( 0 ),
                                            2.0 * marsVelocity( 1 ) - incomingVelocity( 1 ),
                                            0.0 );
    
    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV(
                marsGravitationalParameter,
                marsVelocity, incomingVelocity,
                outgoingVelocity,
                marsSmallestPeriapsisDistance );
    
    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, velocityTolerance );
}

//! Test case for both bending angle and velocity effect delta V.
BOOST_AUTO_TEST_CASE( testBendingAngleAndVelocityEffectDeltaVPericenter )
{
    // Tolerance.
    const double tolerance = 1e-12;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 183.8481861944;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 36000.0 , 0.0 , 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 34500.0, 0.0, 0.0 );

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV(
                venusGravitationalParameter,
                venusVelocity,incomingVelocity,
                outgoingVelocity,
                venusSmallestPeriapsisDistance );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test case in which no assist is required.
BOOST_AUTO_TEST_CASE( testNoDeltaVRequired )
{
    // Tolerance.
    const double tolerance = std::numeric_limits< double >::epsilon( );

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 0.0;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 36000.0 , 0.0 , 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.0, 1000.0, 0.0 );

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test velocity effect Delta-V computation using the eccentricity iteration scheme.
BOOST_AUTO_TEST_CASE( testVelocityEffectDeltaVEccentricity )
{
    // Tolerance. Benchmark obtained directly from the GTOP code, based on the first swing-by for
    // the ideal Cassini-1 trajectory. Values were obtained with a 15-digit accuracy from GTOP,
    // resulting in an accuracy of 6e-14 in the final results.
    const double tolerance = 1.0e-13;

    // Expected deltaV cost, as obtained from GTOP code.
    const double expectedDeltaV = 1090.64622870007;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 32851.224953746, -11618.7310059974, -2055.04615890989 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 34216.4827530912, -15170.1440677825,
                                            395.792122152361 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 37954.2431376052, -14093.0467234774,
                                            -5753.53728279429 );

    // Set flag to use iteration scheme on eccentricity.
    const bool useEccentricity = true;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test limit case for deltaV computation with low incoming velocity with eccentricity iteration.
BOOST_AUTO_TEST_CASE( testLimitCaseDeltaVLowIncomingVelocityEccentricity )
{
    // Tolerance.
    const double tolerance = 1.0e-11;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 966.37867363;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 35000.01, 0.0, 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.0 , 1000.0 , 0.0 );

    // Set flag to use iteration scheme on eccentricity.
    const bool useEccentricity = true;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test limit case for deltaV computation with low outgoing velocity with eccentricity iteration.
BOOST_AUTO_TEST_CASE( testLimitCaseDeltaVLowOutgoingVelocityEccentricity )
{
    // Tolerance.
    const double tolerance = 1.0e-11;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 966.37867363;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 35000.0 , 1000.0 , 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.01, 0.0, 0.0 );

    // Set flag to use iteration scheme on eccentricity.
    const bool useEccentricity = true;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity, incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test limit case for deltaV computation with low velocities with eccentricity iteration.
BOOST_AUTO_TEST_CASE( testLimitCaseDeltaVLowVelocitiesEccentricity )
{
    // Tolerance.
    const double tolerance = 1.0e-9;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 0.004260780473;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 35000.0 , 0.02 , 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.01, 0.0, 0.0 );

    // Set flag to use iteration scheme on eccentricity.
    const bool useEccentricity = true;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test velocity effect Delta-V computation using the pericenter iteration scheme.
BOOST_AUTO_TEST_CASE( testVelocityEffectDeltaVPericenter )
{
    // Tolerance. Benchmark obtained directly from the GTOP code, based on the first swing-by for
    // the ideal Cassini-1 trajectory. Values were obtained with a 15-digit accuracy from GTOP,
    // resulting in an accuracy of 6e-14 in the final results.
    const double tolerance = 1.0e-13;

    // Expected deltaV cost, as obtained from GTOP code.
    const double expectedDeltaV = 1090.64622870007;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 32851.224953746, -11618.7310059974, -2055.04615890989 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 34216.4827530912, -15170.1440677825,
                                            395.792122152361 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 37954.2431376052, -14093.0467234774,
                                            -5753.53728279429 );

    // Set flag to use iteration scheme on pericenter.
    const bool useEccentricity = false;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test limit case for deltaV computation with low incoming velocity with pericenter iteration.
BOOST_AUTO_TEST_CASE( testLimitCaseDeltaVLowIncomingVelocityPericenter )
{
    // Tolerance.
    const double tolerance = 1.0e-11;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 966.37867363;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 35000.01, 0.0, 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.0 , 1000.0 , 0.0 );

    // Set flag to use iteration scheme on eccentricity.
    const bool useEccentricity = false;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test limit case for deltaV computation with low outgoing velocity with pericenter iteration.
BOOST_AUTO_TEST_CASE( testLimitCaseDeltaVLowOutgoingVelocityPericenter )
{
    // Tolerance.
    const double tolerance = 1.0e-11;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 966.37867363;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 35000.0 , 1000.0 , 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.01, 0.0, 0.0 );

    // Set flag to use iteration scheme on pericenter.
    const bool useEccentricity = false;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test limit case for deltaV computation with low velocities with pericenter iteration.
BOOST_AUTO_TEST_CASE( testLimitCaseDeltaVLowVelocitiesPericenter )
{
    // Tolerance.
    const double tolerance = 1.0e-9;

    // Expected deltaV cost, as obtained from this code, verified in Excel (Musegaas, 2012).
    const double expectedDeltaV = 0.004260780473;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define smallest periapsis distance factor.
    const double venusSmallestPeriapsisDistance = 6351800.0;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 35000.0, 0.0 , 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 35000.0 , 0.02 , 0.0 );

    // Define heliocentric satellite outgoing vector.
    const Eigen::Vector3d outgoingVelocity( 35000.01, 0.0, 0.0 );

    // Set flag to use iteration scheme on pericenter.
    const bool useEccentricity = false;

    // Perform the gravity assist.
    const double deltaV = mission_segments::calculateGravityAssistDeltaV( venusGravitationalParameter,
                                                           venusVelocity,incomingVelocity,
                                                           outgoingVelocity,
                                                           venusSmallestPeriapsisDistance,
                                                           useEccentricity );

    // Test if the computed delta-V corresponds to the expected value within the specified
    // tolerance.
    BOOST_CHECK_CLOSE_FRACTION( deltaV, expectedDeltaV, tolerance );
}

//! Test unpowered gravity assist propagation.
BOOST_AUTO_TEST_CASE( testUnpoweredGravityAssistPropagation )
{
    // Tolerance. Benchmark obtained directly from the GTOP code, based on the first swing-by for
    // the ideal Messenger trajectory. Values were obtained with a 15-digit accuracy from GTOP,
    // resulting in an accuracy of 2e-14 in the final results.
    const double tolerance = 1.0e-13;

    // Expected deltaV cost, as obtained from GTOP code.
    const Eigen::Vector3d expectedOutgoingVelocity( 12868.5248737923, -22821.444560174,
                                                    -775.698475033994 );

    // Define swingby body gravitational parameter.
    const double earthGravitationalParameter = 3.9860119e14;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d earthVelocity( 15025.522196446, -25544.3782752036, 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 17969.3166254716, -23543.691593914, 6.38384671663496 );

    // Define rotation angle.
    const double rotationAngle = 1.35077257078;

    // Define pericenter radius.
    const double pericenterRadius = 1.80629232251 * 6378000.0;

    // Perform the gravity assist.
    const Eigen::Vector3d outgoingVelocity = mission_segments::calculateUnpoweredGravityAssistOutgoingVelocity(
                                                    earthGravitationalParameter, earthVelocity,
                                                    incomingVelocity, rotationAngle,
                                                    pericenterRadius );

    // Test if the computed outgoing velocity corresponds to the expected velocity within the
    // specified tolerance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedOutgoingVelocity, outgoingVelocity, tolerance );
}

//! Test powered gravity assist propagation function for unpowered gravity assist.
BOOST_AUTO_TEST_CASE( testPoweredGravityAssistPropagationForUnpoweredGravityAssist )
{
    // Tolerance. Benchmark obtained directly from the GTOP code, based on the first swing-by for
    // the ideal Messenger trajectory. Values were obtained with a 15-digit accuracy from GTOP,
    // resulting in an accuracy of 2e-14 in the final results.
    const double tolerance = 1.0e-13;

    // Expected deltaV cost, as obtained from GTOP code.
    const Eigen::Vector3d expectedOutgoingVelocity( 12868.5248737923, -22821.444560174,
                                                    -775.698475033994 );

    // Define swingby body gravitational parameter.
    const double earthGravitationalParameter = 3.9860119e14;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d earthVelocity( 15025.522196446, -25544.3782752036, 0.0 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 17969.3166254716, -23543.691593914, 6.38384671663496 );

    // Define rotation angle.
    const double rotationAngle = 1.35077257078;

    // Define pericenter radius.
    const double pericenterRadius = 1.80629232251 * 6378000.0;

    // Define deltaV.
    const double deltaV = 0.0;

    // Perform the gravity assist.
    const Eigen::Vector3d outgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                                                    earthGravitationalParameter, earthVelocity,
                                                    incomingVelocity, rotationAngle,
                                                    pericenterRadius, deltaV );

    // Test if the computed outgoing velocity corresponds to the expected velocity within the
    // specified tolerance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedOutgoingVelocity, outgoingVelocity, tolerance );
}

//! Test powered gravity assist propagation function for reverse engineered test case.
BOOST_AUTO_TEST_CASE( testPoweredGravityAssistPropagationReverseEngineered )
{
    // Tolerance. Benchmark obtained by reverse engineering the GTOP code, based on the first
    // swing-by for the ideal Cassini-1 trajectory. Values were obtained with a 15-digit accuracy
    // from GTOP, resulting in an accuracy of 1e-14 in the final results.
    const double tolerance = 1.0e-14;

    // Expected deltaV cost, as obtained from GTOP code.
    const Eigen::Vector3d expectedOutgoingVelocity( 37954.2431376052, -14093.0467234774,
                                                    -5753.53728279429 );

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;

    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 32851.224953746, -11618.7310059974, -2055.04615890989 );

    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 34216.4827530912, -15170.1440677825,
                                            395.792122152361 );

    // Define rotation angle.
    const double rotationAngle = -2.0291949514117;

    // Define pericenter radius.
    const double pericenterRadius = 6351801.04541467;

    // Define deltaV.
    const double deltaV = 1090.64622870007;

    // Perform the gravity assist.
    const Eigen::Vector3d outgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                                                    venusGravitationalParameter, venusVelocity,
                                                    incomingVelocity, rotationAngle,
                                                    pericenterRadius, deltaV );

    // Test if the computed outgoing velocity corresponds to the expected velocity within the
    // specified tolerance.
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedOutgoingVelocity, outgoingVelocity, tolerance );
}

//! Test powered gravity assist backward propagation function by checking consistenty with respect to forward
//! propagation. Used the data from the reverse engineered GTOP code.
BOOST_AUTO_TEST_CASE( testPoweredGravityAssistBackwardPropagationConsistencyCheck )
{
    // Benchmark obtained by reverse engineering the GTOP code, based on the first
    // swing-by for the ideal Cassini-1 trajectory. Values were obtained with a 15-digit accuracy
    // from GTOP, resulting in an accuracy of 1e-14 in the final results.
    const double tolerance = 1.0e-14;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;
    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 32851.224953746, -11618.7310059974, -2055.04615890989 );
    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 34216.4827530912, -15170.1440677825,
                                            395.792122152361 );
    // Define ougoing velocity rotation angle.
    const double outgoingVelocityRotationAngle = -2.0291949514117;
    // Define pericenter radius.
    const double pericenterRadius = 6351801.04541467;
    // Define deltaV.
    const double deltaV = 1090.64622870007;

    // Define dummy value for incoming velocity rotation angle (not defined in the same way as ougoing velocity rotation angle!)
    const double incomingVelocityRotationAngle = 0.0;

    // Perform the gravity assist.
    const Eigen::Vector3d outgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
            venusGravitationalParameter, venusVelocity,
            incomingVelocity, outgoingVelocityRotationAngle,
            pericenterRadius, deltaV );

    const Eigen::Vector3d calculatedIncomingVelocity = mission_segments::calculatePoweredGravityAssistIncomingVelocity(
            venusGravitationalParameter, venusVelocity,
            outgoingVelocity, incomingVelocityRotationAngle,
            pericenterRadius, deltaV );

    // Check if deltaV applied according to the calculatedIncomingVelocity is the selected one
    const bool useEccentricity = true;
    const double calculatedDeltaV = mission_segments::calculateGravityAssistDeltaV(
            venusGravitationalParameter, venusVelocity, calculatedIncomingVelocity, outgoingVelocity, pericenterRadius,
            useEccentricity );
    BOOST_CHECK_CLOSE_FRACTION( deltaV, calculatedDeltaV, tolerance );

    // Check if bending angle between the outgoing velocity and the incoming velocity is the expected. To do that, first need
    // to calculate the bending angle
    const Eigen::Vector3d relativeIncomingVelocity = incomingVelocity - venusVelocity;
    const double absoluteRelativeIncomingVelocity = relativeIncomingVelocity.norm( );
    const double incomingEccentricity = 1.0 + pericenterRadius /
            venusGravitationalParameter * absoluteRelativeIncomingVelocity * absoluteRelativeIncomingVelocity;
    const double incomingBendingAngle = std::asin ( 1.0 / incomingEccentricity );
    const double incomingPericenterVelocity = std::sqrt( absoluteRelativeIncomingVelocity * absoluteRelativeIncomingVelocity *
            ( incomingEccentricity + 1.0 ) / ( incomingEccentricity - 1.0 ) );
    const double outgoingPericenterVelocity = incomingPericenterVelocity + deltaV;
    const double absoluteRelativeOutgoingVelocity = std::sqrt( outgoingPericenterVelocity * outgoingPericenterVelocity -
            2.0 * venusGravitationalParameter / pericenterRadius );
    const double outgoingBendingAngle = std::asin ( 1.0 / ( 1.0 + absoluteRelativeOutgoingVelocity *
            absoluteRelativeOutgoingVelocity * pericenterRadius / venusGravitationalParameter ) );
    const double bendingAngle = incomingBendingAngle + outgoingBendingAngle;

    const Eigen::Vector3d planetocentricOutgoingVelocity = outgoingVelocity - venusVelocity;
    const Eigen::Vector3d planetocentricCalculatedIncomingVelocity = calculatedIncomingVelocity - venusVelocity;

    const double calculatedBendingAngle = std::acos(
            planetocentricOutgoingVelocity.dot(planetocentricCalculatedIncomingVelocity) / planetocentricOutgoingVelocity.norm() /
                    planetocentricCalculatedIncomingVelocity.norm() );
    BOOST_CHECK_CLOSE_FRACTION( calculatedBendingAngle, bendingAngle, tolerance );
}

//! Test powered gravity assist backward propagation function: use forward propagation function to compute the outgoing
//! velocity from the incoming velocity. Based on those two, extract the incoming velocity rotation angle. Using that,
//! backward propagate the outgoing velocity, and check if the result is consistent with the initial incoming velocity.
//! Data used for the initial forward propagation comes from the reverse engineered GTOP code.
BOOST_AUTO_TEST_CASE( testPoweredGravityAssistBackwardPropagation )
{
    // Benchmark obtained by reverse engineering the GTOP code, based on the first
    // swing-by for the ideal Cassini-1 trajectory. Values were obtained with a 15-digit accuracy
    // from GTOP, resulting in an accuracy of 1e-14 in the final results.
    const double tolerance = 1.0e-14;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;
    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 32851.224953746, -11618.7310059974, -2055.04615890989 );
    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 34216.4827530912, -15170.1440677825,
                                            395.792122152361 );
    // Define ougoing velocity rotation angle.
    const double outgoingVelocityRotationAngle = -2.0291949514117;
    // Define pericenter radius.
    const double pericenterRadius = 6351801.04541467;
    // Define deltaV.
    const double deltaV = 1090.64622870007;

    // Perform the gravity assist.
    const Eigen::Vector3d outgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
            venusGravitationalParameter, venusVelocity,
            incomingVelocity, outgoingVelocityRotationAngle,
            pericenterRadius, deltaV );

    const Eigen::Vector3d relativeIncomingVelocity = incomingVelocity - venusVelocity;
    const Eigen::Vector3d relativeOutgoingVelocity = outgoingVelocity - venusVelocity;

    // Compute the unit vectors used in the forward computation of the gravity assist
    const Eigen::Vector3d unitVectorIncoming1 = relativeIncomingVelocity.normalized( );
    const Eigen::Vector3d unitVectorIncoming2 = unitVectorIncoming1.cross( venusVelocity ).normalized( );
    const Eigen::Vector3d unitVectorIncoming3 = unitVectorIncoming1.cross( unitVectorIncoming2 );

    // Compute the unit vectors used in the backward computation of the gravity assist
    const Eigen::Vector3d unitVectorOutgoing1 = relativeOutgoingVelocity.normalized();
    const Eigen::Vector3d unitVectorOutgoing2 = unitVectorOutgoing1.cross( venusVelocity ).normalized( );
    const Eigen::Vector3d unitVectorOutgoing3 = unitVectorOutgoing1.cross( unitVectorOutgoing2 );

    // Compute rotation matrix between base frame and incoming-velocity frame
    Eigen::MatrixXd A(3,3);
    A.row(0) = unitVectorIncoming1;
    A.row(1) = unitVectorIncoming2;
    A.row(2) = unitVectorIncoming3;

    // Compute rotation matrix between base frame and outgoing-velocity frame
    Eigen::MatrixXd B(3,3);
    B.row(0) = unitVectorOutgoing1;
    B.row(1) = unitVectorOutgoing2;
    B.row(2) = unitVectorOutgoing3;

    // Rotation matrix between incoming-velocity frame and outgoing-velocity frame
    Eigen::MatrixXd C = A * B.inverse();

    // Extract bending angle and incoming velocity rotation angle
    const double bendingAngle = std::acos(C(0,0));
    const double incomingVelocityRotationAngleCosine = C(0,1) / sin(bendingAngle);
    const double incomingVelocityRotationAngleSine = C(0,2) / sin(bendingAngle);
    double incomingVelocityRotationAngle;

    if (incomingVelocityRotationAngleCosine > 0)
    {
        incomingVelocityRotationAngle = std::asin(incomingVelocityRotationAngleSine);
    }
    else if ( incomingVelocityRotationAngleSine > 0)
    {
        incomingVelocityRotationAngle = std::acos(incomingVelocityRotationAngleCosine);
    }
    else
    {
        incomingVelocityRotationAngle = mathematical_constants::PI - std::acos(incomingVelocityRotationAngleCosine);
    }

    // Backward propagate gravity assist to get incoming velocity
    const Eigen::Vector3d calculatedIncomingVelocity = mission_segments::calculatePoweredGravityAssistIncomingVelocity(
            venusGravitationalParameter, venusVelocity,
            outgoingVelocity, incomingVelocityRotationAngle,
            pericenterRadius, deltaV );

    // Check if backward propagated incoming velocity coincides with original incoming velocity
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION( calculatedIncomingVelocity, incomingVelocity, tolerance );
}

//! Test forward and backward powered gravity assist when the periapsis is infinity
BOOST_AUTO_TEST_CASE( testGravityInfinitePeriapsis )
{

    // Benchmark obtained by reverse engineering the GTOP code, based on the first
    // swing-by for the ideal Cassini-1 trajectory. Values were obtained with a 15-digit accuracy
    // from GTOP, resulting in an accuracy of 1e-14 in the final results.
    const double tolerance = 1.0e-14;

    // Define swingby body gravitational parameter.
    const double venusGravitationalParameter = 3.24860e14;
    // Define heliocentric planet velocity vector.
    const Eigen::Vector3d venusVelocity( 32851.224953746, -11618.7310059974, -2055.04615890989 );
    // Define heliocentric satellite incoming vector.
    const Eigen::Vector3d incomingVelocity( 34216.4827530912, -15170.1440677825,
                                            395.792122152361 );
    // Define ougoing velocity rotation angle.
    const double outgoingVelocityRotationAngle = -2.0291949514117;
    // Define pericenter radius.
    const double pericenterRadius = std::numeric_limits< double >::infinity();

    for ( unsigned int deltaVCase : {0,1})
    {
        double deltaV;
        Eigen::Vector3d expectedOutgoingVelocity;
        
        if ( deltaVCase == 0 )
        {
            // Define deltaV.
            deltaV = 0.0;
            expectedOutgoingVelocity = incomingVelocity;
        }
        else
        {
            deltaV = 1090.64622870007;
            const Eigen::Vector3d relativeIncomingVelocity = incomingVelocity - venusVelocity;
            expectedOutgoingVelocity = incomingVelocity + deltaV * relativeIncomingVelocity / relativeIncomingVelocity.norm();
        }

        // Perform the gravity assist.
        const Eigen::Vector3d computedOutgoingVelocity = mission_segments::calculatePoweredGravityAssistOutgoingVelocity(
                venusGravitationalParameter, venusVelocity,
                incomingVelocity, outgoingVelocityRotationAngle,
                pericenterRadius, deltaV );

        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( expectedOutgoingVelocity, computedOutgoingVelocity, tolerance );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests
} // namespace tudat
