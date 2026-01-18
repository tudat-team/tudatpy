/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Ephemerides WASM tests - Ported from native tests.
 *    Based on: unitTestSimpleRotationalEphemeris.cpp, unitTestKeplerEphemeris.cpp
 */

#include "wasmTestFramework.h"

#include <memory>
#include <map>

// Mathematics
#include "tudat/math/basic/mathematicalConstants.h"

// Basic astrodynamics
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/orbitalElementConversions.h"

// Ephemerides
#include "tudat/astro/ephemerides/simpleRotationalEphemeris.h"
#include "tudat/astro/ephemerides/keplerEphemeris.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/astro/ephemerides/constantEphemeris.h"

// Reference frames
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

// Interpolators
#include "tudat/math/interpolators/cubicSplineInterpolator.h"

using namespace tudat;
using namespace tudat::ephemerides;
using namespace tudat::unit_conversions;
using namespace tudat::orbital_element_conversions;

// =============================================================================
// Kepler Propagator Benchmark Data (from keplerPropagatorTestData.h)
// =============================================================================

typedef std::map<double, Eigen::Vector6d> PropagationHistory;

double getMelmanEarthGravitationalParameter()
{
    return 3.986004415e14;  // m^3/s^2
}

PropagationHistory getODTBXBenchmarkData()
{
    PropagationHistory benchmarkPropagationHistory;

    Eigen::Vector6d stateInKeplerianElements;
    stateInKeplerianElements << 42165.3431351313e3, 0.26248354351331, 0.30281462522101,
                                4.71463172847351, 4.85569272927819, 2.37248926702153;

    double timeStep = 8640.0;

    for (unsigned int i = 0; i < 11; i++) {
        benchmarkPropagationHistory[static_cast<double>(i) * timeStep] = stateInKeplerianElements;
    }

    benchmarkPropagationHistory[1.0 * timeStep](5) = 2.79722436211144;
    benchmarkPropagationHistory[2.0 * timeStep](5) = 3.18337407409023;
    benchmarkPropagationHistory[3.0 * timeStep](5) = 3.57400974200765;
    benchmarkPropagationHistory[4.0 * timeStep](5) = 4.01425565759545;
    benchmarkPropagationHistory[5.0 * timeStep](5) = 4.57232665706546;
    benchmarkPropagationHistory[6.0 * timeStep](5) = 5.35956850972672;
    benchmarkPropagationHistory[7.0 * timeStep](5) = 0.137251905665217;
    benchmarkPropagationHistory[8.0 * timeStep](5) = 1.14521863765007;
    benchmarkPropagationHistory[9.0 * timeStep](5) = 1.86433634881636;
    benchmarkPropagationHistory[10.0 * timeStep](5) = 2.38486787064101;

    return benchmarkPropagationHistory;
}

double getGTOPGravitationalParameter()
{
    return 1.327e20;  // Sun GM, m^3/s^2
}

PropagationHistory getGTOPBenchmarkData()
{
    PropagationHistory benchmarkPropagationHistory;

    Eigen::Vector6d stateInCartesianElements, stateInKeplerianElements;

    stateInCartesianElements << 1.5e11, 0.0, 0.0, 0.0, 6.0e4, 0.0;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
        stateInCartesianElements, getGTOPGravitationalParameter());
    benchmarkPropagationHistory[0.0] = stateInKeplerianElements;

    stateInCartesianElements << 50369576778.98602, 453006898372.5074, 2.156946592732799e-005,
                                -14654.13750690802, 46884.94068619227, 4.665334803219454e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
        stateInCartesianElements, getGTOPGravitationalParameter());
    benchmarkPropagationHistory[8640000.0] = stateInKeplerianElements;

    stateInCartesianElements << -76810236076.38216, 842661848023.4473, 6.100268297443444e-005,
                                -14683.57015580225, 43917.12010513522, 4.48721854707566e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
        stateInCartesianElements, getGTOPGravitationalParameter());
    benchmarkPropagationHistory[8640000.0 * 2.0] = stateInKeplerianElements;

    stateInCartesianElements << -203052258817.9893, 1216808019495.603, 9.937145023346651e-005,
                                -14543.34378775917, 42828.66589049961, 4.403399939593385e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
        stateInCartesianElements, getGTOPGravitationalParameter());
    benchmarkPropagationHistory[8640000.0 * 3.0] = stateInKeplerianElements;

    stateInCartesianElements << -328225472457.8796, 1584186440047.591, 0.0001371949389119038,
                                -14437.813524927732, 42264.20425643964, 4.355914471377053e-012;
    stateInKeplerianElements = convertCartesianToKeplerianElements(
        stateInCartesianElements, getGTOPGravitationalParameter());
    benchmarkPropagationHistory[8640000.0 * 4.0] = stateInKeplerianElements;

    return benchmarkPropagationHistory;
}

// =============================================================================
// Simple Rotational Ephemeris Tests
// =============================================================================

void testSimpleRotationalEphemeris()
{
    std::cout << "\n=== Simple Rotational Ephemeris (Venus Benchmark) ===" << std::endl;

    // Reference: unitTestSimpleRotationalEphemeris.cpp
    // Data from pck00010.tpc SPICE kernel for Venus rotation

    const double venusPoleRightAscension = convertDegreesToRadians(272.76);
    const double venusPoleDeclination = convertDegreesToRadians(67.16);
    const double venusPrimeMeridianAtJ2000 = convertDegreesToRadians(160.20);
    const double venusRotationRate = convertDegreesToRadians(-1.4813688) / physical_constants::JULIAN_DAY;

    const std::string baseFrame = "J2000";
    const std::string targetFrame = "IAU_VENUS";

    // Calculate initial rotation quaternion to frame
    const Eigen::Quaterniond initialRotationToTargetFrame =
        reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
            venusPoleDeclination, venusPoleRightAscension, venusPrimeMeridianAtJ2000);

    // Create rotational ephemeris from angles
    SimpleRotationalEphemeris venusRotationalEphemerisFromAngles(
        venusPoleRightAscension, venusPoleDeclination, venusPrimeMeridianAtJ2000,
        venusRotationRate, 0.0, baseFrame, targetFrame);

    // Create rotational ephemeris from initial state
    SimpleRotationalEphemeris venusRotationalEphemerisFromInitialState(
        initialRotationToTargetFrame, venusRotationRate, 0.0, baseFrame, targetFrame);

    // Test 1: Initial rotation matrix against SPICE benchmark
    {
        // Benchmark data from SPICE (pck00010.tpc)
        Eigen::Matrix3d spiceInitialRotationToTargetFrameMatrix;
        spiceInitialRotationToTargetFrameMatrix <<
            -0.9548214974296336, 0.2665104385944917, 0.1314841974018291,
            -0.296591573568662, -0.882413772579987, -0.3652114078848295,
            0.01869081416890202, -0.3877088083617989, 0.9215923900425705;

        // Check calculated initial state against SPICE benchmark
        Eigen::Matrix3d calculatedInitialRotation = initialRotationToTargetFrame.toRotationMatrix();

        double maxError = (calculatedInitialRotation - spiceInitialRotationToTargetFrameMatrix).cwiseAbs().maxCoeff();
        checkTrue("Venus initial rotation vs SPICE (max err < 1e-15)", maxError < 1e-15);

        // Check ephemeris calculations at t=0
        Eigen::Matrix3d ephemerisRotation = venusRotationalEphemerisFromAngles.getRotationToTargetFrame(0.0).toRotationMatrix();
        double ephemerisError = (ephemerisRotation - calculatedInitialRotation).cwiseAbs().maxCoeff();
        checkTrue("Venus ephemeris at t=0 vs initial (max err < 1e-15)", ephemerisError < 1e-15);
    }

    // Test 2: Rotation at specified time (1e6 seconds since J2000)
    const double secondsSinceJ2000 = 1.0E6;
    {
        // Benchmark rotation matrix from SPICE
        Eigen::Matrix3d spiceRotationMatrix;
        spiceRotationMatrix <<
            -0.8249537745726603, 0.5148010526833556, 0.2333048348715243,
            -0.5648910720519699, -0.7646317780963481, -0.3102197940834743,
            0.01869081416890206, -0.3877088083617987, 0.9215923900425707;

        // Benchmark rotation matrix derivative from SPICE
        Eigen::Matrix3d spiceRotationMatrixDerivative;
        spiceRotationMatrixDerivative <<
            1.690407961416589e-07, 2.288121921543265e-07, 9.283170431475241e-08,
            -2.468632444964533e-07, 1.540516111965609e-07, 6.981529179974795e-08,
            0.0, 0.0, 0.0;

        Eigen::Quaterniond ephemerisRotation = venusRotationalEphemerisFromAngles.getRotationToTargetFrame(secondsSinceJ2000);
        Eigen::Matrix3d ephemerisRotationDerivative =
            venusRotationalEphemerisFromAngles.getDerivativeOfRotationToTargetFrame(secondsSinceJ2000);

        // Check rotation matrix
        double maxRotError = (ephemerisRotation.toRotationMatrix() - spiceRotationMatrix).cwiseAbs().maxCoeff();
        checkTrue("Venus rotation at 1e6s vs SPICE (max err < 2e-15)", maxRotError < 2e-15);

        // Check rotation derivative
        double maxDerivError = (ephemerisRotationDerivative - spiceRotationMatrixDerivative).cwiseAbs().maxCoeff();
        checkTrue("Venus rotation deriv at 1e6s vs SPICE (max err < 2e-22)", maxDerivError < 2e-22);
    }

    // Test 3: Rotation matrix derivative via finite differences
    {
        double timeStep = 10.0;

        Eigen::Quaterniond ephemerisRotation1 = venusRotationalEphemerisFromAngles.getRotationToTargetFrame(secondsSinceJ2000);
        Eigen::Quaterniond ephemerisRotation2 = venusRotationalEphemerisFromAngles.getRotationToTargetFrame(secondsSinceJ2000 + timeStep);

        Eigen::Matrix3d numericalDerivative =
            (ephemerisRotation2.toRotationMatrix() - ephemerisRotation1.toRotationMatrix()) / timeStep;

        Eigen::Matrix3d analyticalDerivative =
            venusRotationalEphemerisFromAngles.getDerivativeOfRotationToTargetFrame(secondsSinceJ2000 + timeStep / 2.0);

        double maxError = (numericalDerivative - analyticalDerivative).cwiseAbs().maxCoeff();
        checkTrue("Venus rotation deriv (numerical vs analytical, err < 2e-16)", maxError < 2e-16);
    }

    // Test 4: Rotational velocity vector
    {
        // Benchmark from SPICE
        Eigen::Vector3d spiceRotationalVelocityVectorInBaseFrame;
        spiceRotationalVelocityVectorInBaseFrame << -5.593131603532092e-09, 1.160198999048488e-07, -2.75781861386115e-07;

        Eigen::Vector3d computedVelocity =
            venusRotationalEphemerisFromAngles.getRotationalVelocityVectorInBaseFrame(secondsSinceJ2000);

        double maxError = (computedVelocity - spiceRotationalVelocityVectorInBaseFrame).cwiseAbs().maxCoeff();
        checkTrue("Venus rotational velocity vs SPICE (max err < 1e-21)", maxError < 1e-21);
    }

    // Test 5: Orthonormality of rotation matrices
    {
        Eigen::Matrix3d productOfOppositeRotations =
            venusRotationalEphemerisFromAngles.getRotationToTargetFrame(secondsSinceJ2000).toRotationMatrix() *
            venusRotationalEphemerisFromAngles.getRotationToBaseFrame(secondsSinceJ2000).toRotationMatrix();

        double maxError = (productOfOppositeRotations - Eigen::Matrix3d::Identity()).cwiseAbs().maxCoeff();
        checkTrue("Venus rotation orthonormality (max err < 1e-15)", maxError < 1e-15);
    }
}

// =============================================================================
// Kepler Ephemeris Tests
// =============================================================================

void testKeplerEphemerisElliptical()
{
    std::cout << "\n=== Kepler Ephemeris (Elliptical - ODTBX Benchmark) ===" << std::endl;

    // Reference: unitTestKeplerEphemeris.cpp - testKeplerEphemerisElliptical
    // Benchmark data from ODTBX (Orbit Determination Toolbox)

    const double earthGravitationalParameter = 398600.4415e9;  // m^3/s^2

    PropagationHistory expectedPropagationHistory = getODTBXBenchmarkData();

    // Create Kepler ephemeris from initial state
    KeplerEphemeris keplerEphemeris(expectedPropagationHistory[0.0], 0.0, earthGravitationalParameter);

    int testCount = 0;
    int passCount = 0;

    for (auto& stateIterator : expectedPropagationHistory)
    {
        Eigen::Vector6d computedState = convertCartesianToKeplerianElements(
            keplerEphemeris.getCartesianState(stateIterator.first),
            earthGravitationalParameter);

        // Check true anomaly (element 5)
        double relError = std::abs(computedState(5) - stateIterator.second(5)) / std::abs(stateIterator.second(5));
        testCount++;
        if (relError < 2.5e-14) {
            passCount++;
        }
    }

    checkTrue("Kepler elliptical true anomaly (all " + std::to_string(testCount) + " points pass)",
              passCount == testCount);

    std::cout << "[INFO] Tested " << testCount << " time points, " << passCount << " passed" << std::endl;
}

void testKeplerEphemerisHyperbolic()
{
    std::cout << "\n=== Kepler Ephemeris (Hyperbolic - GTOP Benchmark) ===" << std::endl;

    // Reference: unitTestKeplerEphemeris.cpp - testKeplerEphemerisHyperbolic
    // Benchmark data from GTOP (Global Trajectory Optimisation Portal)

    PropagationHistory expectedPropagationHistory = getGTOPBenchmarkData();

    // Create Kepler ephemeris from initial state
    KeplerEphemeris keplerEphemeris(expectedPropagationHistory[0.0], 0.0, getGTOPGravitationalParameter());

    int testCount = 0;
    int passCount = 0;

    for (auto& stateIterator : expectedPropagationHistory)
    {
        Eigen::Vector6d computedState = convertCartesianToKeplerianElements(
            keplerEphemeris.getCartesianState(stateIterator.first),
            getGTOPGravitationalParameter());

        // Check true anomaly (element 5)
        // Note: At t=0 the expected true anomaly is 0, so use absolute error
        // For non-zero expected values, use relative error with 1e-13 tolerance
        testCount++;
        double expected = stateIterator.second(5);
        double computed = computedState(5);
        bool passed;

        if (std::abs(expected) < 1e-10) {
            // Absolute error for near-zero expected values
            passed = std::abs(computed - expected) < 1e-10;
        } else {
            // Relative error for non-zero expected values
            double relError = std::abs(computed - expected) / std::abs(expected);
            passed = relError < 1e-13;
        }

        if (passed) {
            passCount++;
        }
    }

    checkTrue("Kepler hyperbolic true anomaly (all " + std::to_string(testCount) + " points pass)",
              passCount == testCount);

    std::cout << "[INFO] Tested " << testCount << " time points, " << passCount << " passed" << std::endl;
}

// =============================================================================
// Tabulated Ephemeris Tests
// =============================================================================

void testTabulatedEphemeris()
{
    std::cout << "\n=== Tabulated Ephemeris ===" << std::endl;

    // Reference: unitTestTabulatedEphemeris.cpp
    // Tests the tabulated ephemeris using interpolation

    // Create a simple analytical ephemeris (circular orbit)
    const double earthGM = 3.986004418e14;
    const double semiMajorAxis = 7000.0e3;  // 7000 km
    const double orbitalPeriod = 2.0 * mathematical_constants::PI * std::sqrt(
        semiMajorAxis * semiMajorAxis * semiMajorAxis / earthGM);

    // Generate state history for a circular orbit
    std::map<double, Eigen::Vector6d> stateHistoryMap;

    double startTime = 0.0;
    double finalTime = orbitalPeriod;
    double timeStep = 60.0;  // 1 minute steps

    for (double t = startTime; t <= finalTime; t += timeStep)
    {
        double meanAnomaly = 2.0 * mathematical_constants::PI * t / orbitalPeriod;
        Eigen::Vector6d keplerState;
        keplerState << semiMajorAxis, 0.0, 0.0, 0.0, 0.0, meanAnomaly;
        stateHistoryMap[t] = convertKeplerianToCartesianElements(keplerState, earthGM);
    }

    // Create interpolator
    auto stateInterpolator = std::make_shared<interpolators::CubicSplineInterpolator<double, Eigen::Vector6d>>(
        stateHistoryMap);

    // Create tabulated ephemeris
    TabulatedCartesianEphemeris<> tabulatedEphemeris(stateInterpolator, "Earth", "J2000");

    // Test at an interpolated time point
    double testTime = orbitalPeriod / 4.0;  // Quarter orbit

    Eigen::Vector6d interpolatedState = stateInterpolator->interpolate(testTime);
    Eigen::Vector6d ephemerisState = tabulatedEphemeris.getCartesianState(testTime);

    // Should match exactly (same interpolator)
    double stateError = (interpolatedState - ephemerisState).norm();
    checkTrue("Tabulated ephemeris matches interpolator (err = 0)", stateError < 1e-15);

    // Test position magnitude (should be approximately semi-major axis for circular orbit)
    double positionMagnitude = ephemerisState.head<3>().norm();
    double positionRelError = std::abs(positionMagnitude - semiMajorAxis) / semiMajorAxis;
    checkTrue("Tabulated ephemeris position magnitude (rel err < 1e-3)", positionRelError < 1e-3);

    // Test velocity magnitude (should be sqrt(GM/a) for circular orbit)
    double expectedVelocity = std::sqrt(earthGM / semiMajorAxis);
    double velocityMagnitude = ephemerisState.tail<3>().norm();
    double velocityRelError = std::abs(velocityMagnitude - expectedVelocity) / expectedVelocity;
    checkTrue("Tabulated ephemeris velocity magnitude (rel err < 1e-3)", velocityRelError < 1e-3);

    std::cout << "[INFO] Position: " << positionMagnitude / 1e3 << " km (expected: "
              << semiMajorAxis / 1e3 << " km)" << std::endl;
    std::cout << "[INFO] Velocity: " << velocityMagnitude / 1e3 << " km/s (expected: "
              << expectedVelocity / 1e3 << " km/s)" << std::endl;
}

void testConstantEphemeris()
{
    std::cout << "\n=== Constant Ephemeris ===" << std::endl;

    // Test a simple constant state ephemeris
    Eigen::Vector6d constantState;
    constantState << 7000.0e3, 0.0, 0.0, 0.0, 7.5e3, 0.0;

    ConstantEphemeris constantEphemeris(constantState, "Earth", "J2000");

    // Test at various times
    Eigen::Vector6d state0 = constantEphemeris.getCartesianState(0.0);
    Eigen::Vector6d state1 = constantEphemeris.getCartesianState(1000.0);
    Eigen::Vector6d state2 = constantEphemeris.getCartesianState(-1000.0);

    checkVector6dClose("Constant ephemeris at t=0", state0, constantState, 1e-15);
    checkVector6dClose("Constant ephemeris at t=1000", state1, constantState, 1e-15);
    checkVector6dClose("Constant ephemeris at t=-1000", state2, constantState, 1e-15);
}
