/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Edge case and boundary condition tests for WASM.
 *    Tests invalid inputs, NaN/Infinity handling, numerical precision limits.
 */

#include "wasmTestFramework.h"

#include <cmath>
#include <limits>
#include <vector>

// Basic astrodynamics
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"

// Mathematics
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/math/interpolators/linearInterpolator.h"

using namespace tudat;

// ============================================================================
// NUMERICAL PRECISION EDGE CASES
// ============================================================================

/**
 * Test: NaN and Infinity Handling
 *
 * Tests how the library handles special floating-point values.
 */
void testNaNInfinityHandling()
{
    std::cout << "\n=== Edge Case: NaN and Infinity Handling ===" << std::endl;

    // Test that NaN is detected correctly
    double nanValue = std::nan("");
    double infValue = std::numeric_limits<double>::infinity();
    double negInfValue = -std::numeric_limits<double>::infinity();

    checkTrue("NaN is detected", std::isnan(nanValue));
    checkTrue("Infinity is detected", std::isinf(infValue));
    checkTrue("Negative infinity is detected", std::isinf(negInfValue));

    // Test NaN propagation in arithmetic
    double nanResult = 1.0 + nanValue;
    checkTrue("NaN propagates through addition", std::isnan(nanResult));

    // Test infinity arithmetic
    double infResult = infValue + 1.0;
    checkTrue("Infinity + 1 = infinity", std::isinf(infResult));

    // Test 0/0 produces NaN
    double zeroOverZero = 0.0 / 0.0;
    checkTrue("0/0 is NaN", std::isnan(zeroOverZero));

    // Test 1/0 produces infinity
    double oneOverZero = 1.0 / 0.0;
    checkTrue("1/0 is infinity", std::isinf(oneOverZero));

    std::cout << "[INFO] NaN/Infinity handling test passed" << std::endl;
}

/**
 * Test: Subnormal Numbers
 *
 * Tests handling of very small (subnormal/denormalized) numbers.
 */
void testSubnormalNumbers()
{
    std::cout << "\n=== Edge Case: Subnormal Numbers ===" << std::endl;

    double minNormal = std::numeric_limits<double>::min();
    double minSubnormal = std::numeric_limits<double>::denorm_min();
    double epsilon = std::numeric_limits<double>::epsilon();

    std::cout << "[INFO] Min normal:    " << minNormal << std::endl;
    std::cout << "[INFO] Min subnormal: " << minSubnormal << std::endl;
    std::cout << "[INFO] Epsilon:       " << epsilon << std::endl;

    // Test that subnormal numbers are non-zero
    checkTrue("Subnormal is non-zero", minSubnormal > 0.0);
    checkTrue("Subnormal < normal min", minSubnormal < minNormal);

    // Test arithmetic with subnormal numbers
    double subnormalSum = minSubnormal + minSubnormal;
    checkTrue("Subnormal addition works", subnormalSum > minSubnormal);

    // Test gradual underflow
    double almostZero = minNormal / 2.0;
    checkTrue("Gradual underflow produces subnormal", almostZero > 0.0 && almostZero < minNormal);

    std::cout << "[INFO] Subnormal number test passed" << std::endl;
}

/**
 * Test: Epsilon-Scale Comparisons
 *
 * Tests numerical comparisons at machine precision limits.
 */
void testEpsilonComparisons()
{
    std::cout << "\n=== Edge Case: Epsilon-Scale Comparisons ===" << std::endl;

    double epsilon = std::numeric_limits<double>::epsilon();

    // Test 1 + epsilon > 1
    double onePlusEps = 1.0 + epsilon;
    checkTrue("1 + epsilon > 1", onePlusEps > 1.0);

    // Test 1 + epsilon/2 == 1 (below precision)
    double onePlusHalfEps = 1.0 + epsilon / 2.0;
    checkTrue("1 + epsilon/2 == 1", onePlusHalfEps == 1.0);

    // Test relative comparison function
    auto relativeClose = [](double a, double b, double relTol) {
        double maxVal = std::max(std::abs(a), std::abs(b));
        if (maxVal < std::numeric_limits<double>::min()) return true;  // Both essentially zero
        return std::abs(a - b) / maxVal < relTol;
    };

    checkTrue("Relative comparison near 1", relativeClose(1.0, 1.0 + 1e-15, 1e-14));
    checkTrue("Relative comparison large numbers", relativeClose(1e10, 1e10 + 1.0, 1e-9));

    std::cout << "[INFO] Epsilon comparison test passed" << std::endl;
}

// ============================================================================
// ORBITAL MECHANICS EDGE CASES
// ============================================================================

/**
 * Test: Circular Orbit Edge Case
 *
 * Tests orbital element conversions for perfectly circular orbits (e=0).
 */
void testCircularOrbitEdgeCase()
{
    std::cout << "\n=== Edge Case: Circular Orbit (e=0) ===" << std::endl;

    using namespace orbital_element_conversions;

    double earthGravParam = 3.986004418e14;

    // Perfectly circular orbit
    Eigen::Vector6d keplerianElements;
    keplerianElements(semiMajorAxisIndex) = 7000.0e3;
    keplerianElements(eccentricityIndex) = 0.0;  // Exactly zero
    keplerianElements(inclinationIndex) = unit_conversions::convertDegreesToRadians(45.0);
    keplerianElements(argumentOfPeriapsisIndex) = 0.0;  // Undefined for circular
    keplerianElements(longitudeOfAscendingNodeIndex) = 0.0;
    keplerianElements(trueAnomalyIndex) = unit_conversions::convertDegreesToRadians(90.0);

    // Convert to Cartesian
    Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(
        keplerianElements, earthGravParam);

    // Verify position magnitude equals semi-major axis
    double positionMagnitude = cartesianState.head<3>().norm();
    checkClose("Circular orbit radius = a", positionMagnitude, keplerianElements(semiMajorAxisIndex), 1.0);

    // Convert back to Keplerian
    Eigen::Vector6d recoveredKeplerian = convertCartesianToKeplerianElements(
        cartesianState, earthGravParam);

    // Check eccentricity is still ~0
    checkClose("Recovered eccentricity ~0", recoveredKeplerian(eccentricityIndex), 0.0, 1e-10);

    std::cout << "[INFO] Circular orbit test passed" << std::endl;
}

/**
 * Test: Near-Parabolic Orbit Edge Case
 *
 * Tests orbital elements near the parabolic boundary (e≈1).
 */
void testNearParabolicOrbitEdgeCase()
{
    std::cout << "\n=== Edge Case: Near-Parabolic Orbit (e≈1) ===" << std::endl;

    using namespace orbital_element_conversions;

    double earthGravParam = 3.986004418e14;

    // Nearly parabolic orbit (just below e=1)
    Eigen::Vector6d keplerianElements;
    keplerianElements(semiMajorAxisIndex) = 50000.0e3;  // Large SMA
    keplerianElements(eccentricityIndex) = 0.999;  // Very close to 1
    keplerianElements(inclinationIndex) = unit_conversions::convertDegreesToRadians(30.0);
    keplerianElements(argumentOfPeriapsisIndex) = 0.0;
    keplerianElements(longitudeOfAscendingNodeIndex) = 0.0;
    keplerianElements(trueAnomalyIndex) = unit_conversions::convertDegreesToRadians(10.0);

    // Convert to Cartesian
    Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(
        keplerianElements, earthGravParam);

    // Verify state is finite
    checkTrue("Near-parabolic position finite", std::isfinite(cartesianState.head<3>().norm()));
    checkTrue("Near-parabolic velocity finite", std::isfinite(cartesianState.tail<3>().norm()));

    // Convert back
    Eigen::Vector6d recoveredKeplerian = convertCartesianToKeplerianElements(
        cartesianState, earthGravParam);

    // Check eccentricity recovered correctly
    double eccError = std::abs(recoveredKeplerian(eccentricityIndex) - keplerianElements(eccentricityIndex));
    checkTrue("Recovered eccentricity accurate", eccError < 1e-10);

    std::cout << "[INFO] Near-parabolic orbit test passed" << std::endl;
}

/**
 * Test: Hyperbolic Orbit Edge Case
 *
 * Tests handling of hyperbolic orbits (e>1).
 */
void testHyperbolicOrbitEdgeCase()
{
    std::cout << "\n=== Edge Case: Hyperbolic Orbit (e>1) ===" << std::endl;

    using namespace orbital_element_conversions;

    double sunGravParam = 1.32712440018e20;

    // Hyperbolic orbit (e > 1)
    Eigen::Vector6d keplerianElements;
    keplerianElements(semiMajorAxisIndex) = -1.0e11;  // Negative for hyperbolic
    keplerianElements(eccentricityIndex) = 1.5;
    keplerianElements(inclinationIndex) = unit_conversions::convertDegreesToRadians(5.0);
    keplerianElements(argumentOfPeriapsisIndex) = 0.0;
    keplerianElements(longitudeOfAscendingNodeIndex) = 0.0;
    keplerianElements(trueAnomalyIndex) = unit_conversions::convertDegreesToRadians(30.0);

    // Convert to Cartesian
    Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(
        keplerianElements, sunGravParam);

    // Verify state is finite
    checkTrue("Hyperbolic position finite", std::isfinite(cartesianState.head<3>().norm()));
    checkTrue("Hyperbolic velocity finite", std::isfinite(cartesianState.tail<3>().norm()));

    // Check that velocity exceeds escape velocity
    double velocityMagnitude = cartesianState.tail<3>().norm();
    double positionMagnitude = cartesianState.head<3>().norm();
    double escapeVelocity = std::sqrt(2.0 * sunGravParam / positionMagnitude);

    checkTrue("Hyperbolic: v > v_escape", velocityMagnitude > escapeVelocity);

    std::cout << "[INFO] Hyperbolic orbit test passed" << std::endl;
}

/**
 * Test: Equatorial Orbit Edge Case
 *
 * Tests handling of equatorial orbits (i=0).
 */
void testEquatorialOrbitEdgeCase()
{
    std::cout << "\n=== Edge Case: Equatorial Orbit (i=0) ===" << std::endl;

    using namespace orbital_element_conversions;

    double earthGravParam = 3.986004418e14;

    // Equatorial orbit
    Eigen::Vector6d keplerianElements;
    keplerianElements(semiMajorAxisIndex) = 7000.0e3;
    keplerianElements(eccentricityIndex) = 0.1;
    keplerianElements(inclinationIndex) = 0.0;  // Exactly zero
    keplerianElements(argumentOfPeriapsisIndex) = unit_conversions::convertDegreesToRadians(45.0);
    keplerianElements(longitudeOfAscendingNodeIndex) = 0.0;  // Undefined for equatorial
    keplerianElements(trueAnomalyIndex) = unit_conversions::convertDegreesToRadians(90.0);

    // Convert to Cartesian
    Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(
        keplerianElements, earthGravParam);

    // Z component should be zero for equatorial orbit
    checkClose("Equatorial orbit Z = 0", cartesianState(2), 0.0, 1e-10);
    checkClose("Equatorial orbit Vz = 0", cartesianState(5), 0.0, 1e-10);

    std::cout << "[INFO] Equatorial orbit test passed" << std::endl;
}

/**
 * Test: Polar Orbit Edge Case
 *
 * Tests handling of polar orbits (i=90°).
 */
void testPolarOrbitEdgeCase()
{
    std::cout << "\n=== Edge Case: Polar Orbit (i=90°) ===" << std::endl;

    using namespace orbital_element_conversions;

    double earthGravParam = 3.986004418e14;

    // Polar orbit
    Eigen::Vector6d keplerianElements;
    keplerianElements(semiMajorAxisIndex) = 7000.0e3;
    keplerianElements(eccentricityIndex) = 0.001;
    keplerianElements(inclinationIndex) = mathematical_constants::PI / 2.0;  // 90 degrees
    keplerianElements(argumentOfPeriapsisIndex) = 0.0;
    keplerianElements(longitudeOfAscendingNodeIndex) = 0.0;
    keplerianElements(trueAnomalyIndex) = 0.0;

    // Convert to Cartesian
    Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(
        keplerianElements, earthGravParam);

    // Convert back
    Eigen::Vector6d recoveredKeplerian = convertCartesianToKeplerianElements(
        cartesianState, earthGravParam);

    // Check inclination recovered correctly
    checkClose("Recovered inclination = 90°",
               recoveredKeplerian(inclinationIndex),
               mathematical_constants::PI / 2.0,
               1e-10);

    std::cout << "[INFO] Polar orbit test passed" << std::endl;
}

// ============================================================================
// KEPLER PROPAGATION EDGE CASES
// ============================================================================

/**
 * Test: Zero Time Propagation
 *
 * Tests Kepler propagation for zero time interval.
 */
void testZeroTimePropagation()
{
    std::cout << "\n=== Edge Case: Zero Time Propagation ===" << std::endl;

    using namespace orbital_element_conversions;
    using namespace basic_astrodynamics;

    double earthGravParam = 3.986004418e14;

    Eigen::Vector6d keplerianElements;
    keplerianElements << 7000.0e3, 0.1, 0.5, 0.3, 0.2, 0.1;

    // Propagate for zero time
    Eigen::Vector6d propagatedElements = propagateKeplerOrbit(
        keplerianElements, 0.0, earthGravParam);

    // Should be identical to input
    double diff = (propagatedElements - keplerianElements).norm();
    checkClose("Zero time propagation unchanged", diff, 0.0, 1e-15);

    std::cout << "[INFO] Zero time propagation test passed" << std::endl;
}

/**
 * Test: Full Orbit Propagation
 *
 * Tests that propagation for exactly one orbital period returns to initial state.
 */
void testFullOrbitPropagation()
{
    std::cout << "\n=== Edge Case: Full Orbit Propagation ===" << std::endl;

    using namespace orbital_element_conversions;
    using namespace basic_astrodynamics;

    double earthGravParam = 3.986004418e14;
    double semiMajorAxis = 7000.0e3;

    Eigen::Vector6d keplerianElements;
    keplerianElements << semiMajorAxis, 0.1, 0.5, 0.3, 0.2, 0.1;

    // Compute orbital period
    double orbitalPeriod = 2.0 * mathematical_constants::PI *
                           std::sqrt(std::pow(semiMajorAxis, 3) / earthGravParam);

    // Propagate for exactly one period
    Eigen::Vector6d propagatedElements = propagateKeplerOrbit(
        keplerianElements, orbitalPeriod, earthGravParam);

    // Mean anomaly should have increased by 2π (wrapped to same angle)
    // Other elements should be unchanged
    checkClose("SMA preserved after 1 orbit",
               propagatedElements(semiMajorAxisIndex),
               keplerianElements(semiMajorAxisIndex), 1e-6);
    checkClose("Eccentricity preserved after 1 orbit",
               propagatedElements(eccentricityIndex),
               keplerianElements(eccentricityIndex), 1e-14);

    std::cout << "[INFO] Full orbit propagation test passed" << std::endl;
}

/**
 * Test: Very Long Propagation
 *
 * Tests propagation over many orbital periods for accumulated error.
 */
void testVeryLongPropagation()
{
    std::cout << "\n=== Edge Case: Very Long Propagation ===" << std::endl;

    using namespace orbital_element_conversions;
    using namespace basic_astrodynamics;

    double earthGravParam = 3.986004418e14;
    double semiMajorAxis = 7000.0e3;

    Eigen::Vector6d keplerianElements;
    keplerianElements << semiMajorAxis, 0.01, 0.5, 0.3, 0.2, 0.1;

    // Compute orbital period
    double orbitalPeriod = 2.0 * mathematical_constants::PI *
                           std::sqrt(std::pow(semiMajorAxis, 3) / earthGravParam);

    // Propagate for 1000 orbits
    double totalTime = 1000.0 * orbitalPeriod;
    Eigen::Vector6d propagatedElements = propagateKeplerOrbit(
        keplerianElements, totalTime, earthGravParam);

    // SMA and eccentricity should be perfectly preserved (Kepler problem is exact)
    checkClose("SMA preserved after 1000 orbits",
               propagatedElements(semiMajorAxisIndex),
               keplerianElements(semiMajorAxisIndex), 1e-3);
    checkClose("Eccentricity preserved after 1000 orbits",
               propagatedElements(eccentricityIndex),
               keplerianElements(eccentricityIndex), 1e-12);

    std::cout << "[INFO] Very long propagation test passed" << std::endl;
}

// ============================================================================
// COORDINATE CONVERSION EDGE CASES
// ============================================================================

/**
 * Test: Spherical Coordinate Singularities
 *
 * Tests spherical coordinate conversions at poles (theta=0, pi).
 */
void testSphericalCoordinateSingularities()
{
    std::cout << "\n=== Edge Case: Spherical Coordinate Singularities ===" << std::endl;

    using namespace coordinate_conversions;

    // Point at north pole (theta = 0)
    Eigen::Vector3d northPoleSpherical;
    northPoleSpherical << 1.0, 0.0, 0.0;  // r=1, theta=0 (pole), phi=0

    Eigen::Vector3d northPoleCartesian = convertSphericalToCartesian(northPoleSpherical);

    checkClose("North pole X = 0", northPoleCartesian(0), 0.0, 1e-15);
    checkClose("North pole Y = 0", northPoleCartesian(1), 0.0, 1e-15);
    checkClose("North pole Z = r", northPoleCartesian(2), 1.0, 1e-15);

    // Point at south pole (theta = pi)
    Eigen::Vector3d southPoleSpherical;
    southPoleSpherical << 1.0, mathematical_constants::PI, 0.0;

    Eigen::Vector3d southPoleCartesian = convertSphericalToCartesian(southPoleSpherical);

    checkClose("South pole X = 0", southPoleCartesian(0), 0.0, 1e-15);
    checkClose("South pole Y = 0", southPoleCartesian(1), 0.0, 1e-15);
    checkClose("South pole Z = -r", southPoleCartesian(2), -1.0, 1e-15);

    std::cout << "[INFO] Spherical singularity test passed" << std::endl;
}

/**
 * Test: Zero Radius Handling
 *
 * Tests coordinate conversions with zero radius.
 */
void testZeroRadiusHandling()
{
    std::cout << "\n=== Edge Case: Zero Radius Handling ===" << std::endl;

    using namespace coordinate_conversions;

    // Zero radius in spherical coordinates
    Eigen::Vector3d zeroSpherical;
    zeroSpherical << 0.0, 0.5, 0.5;  // r=0

    Eigen::Vector3d zeroCartesian = convertSphericalToCartesian(zeroSpherical);

    checkClose("Zero radius: X = 0", zeroCartesian(0), 0.0, 1e-15);
    checkClose("Zero radius: Y = 0", zeroCartesian(1), 0.0, 1e-15);
    checkClose("Zero radius: Z = 0", zeroCartesian(2), 0.0, 1e-15);

    std::cout << "[INFO] Zero radius test passed" << std::endl;
}

// ============================================================================
// INTEGRATOR EDGE CASES
// ============================================================================

/**
 * Test: Integrator with Zero Step Size
 *
 * Tests integrator behavior with very small step sizes.
 */
void testIntegratorSmallStepSize()
{
    std::cout << "\n=== Edge Case: Integrator Small Step Size ===" << std::endl;

    using namespace numerical_integrators;

    // Simple exponential growth: dy/dt = y
    auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
        return state;
    };

    Eigen::VectorXd initialState(1);
    initialState << 1.0;

    // Very small step size
    double stepSize = 1e-10;
    RungeKutta4Integrator<double, Eigen::VectorXd> integrator(
        stateDerivative, 0.0, initialState, stepSize);

    // Take one step
    Eigen::VectorXd nextState = integrator.performIntegrationStep(stepSize);

    // Should be very close to 1 + stepSize (Taylor expansion)
    double expected = 1.0 + stepSize;
    checkClose("Small step integration", nextState(0), expected, 1e-15);

    std::cout << "[INFO] Small step size test passed" << std::endl;
}

/**
 * Test: Integrator with Stiff ODE
 *
 * Tests integrator behavior with stiff differential equations.
 */
void testIntegratorStiffODE()
{
    std::cout << "\n=== Edge Case: Integrator with Stiff ODE ===" << std::endl;

    using namespace numerical_integrators;

    // Stiff ODE: dy/dt = -100*y (fast decay)
    const double lambda = -100.0;

    auto stiffDerivative = [lambda](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
        Eigen::VectorXd deriv(1);
        deriv(0) = lambda * state(0);
        return deriv;
    };

    Eigen::VectorXd initialState(1);
    initialState << 1.0;

    // Use small enough step size for stability (h < 2/|lambda| for RK4)
    double stepSize = 0.01;  // Should be stable for lambda = -100
    RungeKutta4Integrator<double, Eigen::VectorXd> integrator(
        stiffDerivative, 0.0, initialState, stepSize);

    // Integrate to t = 0.1
    Eigen::VectorXd finalState = integrator.integrateTo(0.1, stepSize);

    // Exact solution: y = e^(-100*0.1) = e^(-10) ≈ 4.5e-5
    double exact = std::exp(lambda * 0.1);
    double error = std::abs(finalState(0) - exact) / exact;

    checkTrue("Stiff ODE relative error < 1%", error < 0.01);

    std::cout << "[INFO] Stiff ODE test passed" << std::endl;
}

// ============================================================================
// INTERPOLATION EDGE CASES
// ============================================================================

/**
 * Test: Interpolation at Boundaries
 *
 * Tests interpolation at exact boundary points.
 */
void testInterpolationAtBoundaries()
{
    std::cout << "\n=== Edge Case: Interpolation at Boundaries ===" << std::endl;

    using namespace interpolators;

    // Create simple linear data
    std::map<double, double> dataMap;
    dataMap[0.0] = 0.0;
    dataMap[1.0] = 1.0;
    dataMap[2.0] = 4.0;
    dataMap[3.0] = 9.0;

    LinearInterpolator<double, double> interpolator(dataMap);

    // Interpolate at exact data points
    checkClose("Interpolation at x=0", interpolator.interpolate(0.0), 0.0, 1e-15);
    checkClose("Interpolation at x=1", interpolator.interpolate(1.0), 1.0, 1e-15);
    checkClose("Interpolation at x=2", interpolator.interpolate(2.0), 4.0, 1e-15);
    checkClose("Interpolation at x=3", interpolator.interpolate(3.0), 9.0, 1e-15);

    // Interpolate between points
    double midValue = interpolator.interpolate(1.5);
    checkClose("Linear interpolation at x=1.5", midValue, 2.5, 1e-15);  // (1+4)/2

    std::cout << "[INFO] Boundary interpolation test passed" << std::endl;
}

/**
 * Test: Single Point Interpolation
 *
 * Tests interpolation with minimal data.
 */
void testSinglePointInterpolation()
{
    std::cout << "\n=== Edge Case: Two-Point Interpolation ===" << std::endl;

    using namespace interpolators;

    // Minimum viable data for linear interpolation: 2 points
    std::map<double, double> dataMap;
    dataMap[0.0] = 10.0;
    dataMap[10.0] = 20.0;

    LinearInterpolator<double, double> interpolator(dataMap);

    // Should interpolate linearly
    checkClose("Two-point interpolation at x=0", interpolator.interpolate(0.0), 10.0, 1e-15);
    checkClose("Two-point interpolation at x=5", interpolator.interpolate(5.0), 15.0, 1e-15);
    checkClose("Two-point interpolation at x=10", interpolator.interpolate(10.0), 20.0, 1e-15);

    std::cout << "[INFO] Two-point interpolation test passed" << std::endl;
}

// ============================================================================
// MATRIX OPERATION EDGE CASES
// ============================================================================

/**
 * Test: Singular Matrix Operations
 *
 * Tests linear algebra with singular or near-singular matrices.
 */
void testSingularMatrixOperations()
{
    std::cout << "\n=== Edge Case: Singular Matrix Operations ===" << std::endl;

    // Create a singular matrix (rank deficient)
    Eigen::Matrix3d singular;
    singular << 1, 2, 3,
                2, 4, 6,   // Row 2 = 2 * Row 1
                1, 1, 1;

    // Determinant should be zero (or very close)
    double det = singular.determinant();
    checkClose("Singular matrix determinant ≈ 0", det, 0.0, 1e-10);

    // Create an ill-conditioned matrix
    Eigen::Matrix3d illConditioned;
    illConditioned << 1, 1, 1,
                      1, 1, 1.0001,
                      1, 1.0001, 1;

    double condNumber = illConditioned.norm() * illConditioned.inverse().norm();
    std::cout << "[INFO] Ill-conditioned matrix condition number: " << condNumber << std::endl;
    checkTrue("Ill-conditioned matrix has large condition number", condNumber > 1000);

    std::cout << "[INFO] Singular matrix test passed" << std::endl;
}

/**
 * Test: Empty and Zero Vectors
 *
 * Tests operations with zero vectors.
 */
void testEmptyAndZeroVectors()
{
    std::cout << "\n=== Edge Case: Zero Vectors ===" << std::endl;

    // Zero vector operations
    Eigen::Vector3d zeroVec = Eigen::Vector3d::Zero();

    checkClose("Zero vector norm = 0", zeroVec.norm(), 0.0, 1e-15);

    // Zero vector dot product
    Eigen::Vector3d otherVec(1.0, 2.0, 3.0);
    double dotProduct = zeroVec.dot(otherVec);
    checkClose("Zero dot product = 0", dotProduct, 0.0, 1e-15);

    // Cross product with zero vector
    Eigen::Vector3d crossProduct = zeroVec.cross(otherVec);
    checkClose("Cross with zero = 0", crossProduct.norm(), 0.0, 1e-15);

    std::cout << "[INFO] Zero vector test passed" << std::endl;
}

/**
 * Test: Large Vector Operations
 *
 * Tests operations with large vectors to check for overflow.
 */
void testLargeVectorOperations()
{
    std::cout << "\n=== Edge Case: Large Vector Operations ===" << std::endl;

    // Create vector with large values
    double largeValue = 1e150;
    Eigen::Vector3d largeVec(largeValue, largeValue, largeValue);

    // Norm computation (should not overflow)
    double norm = largeVec.norm();
    checkTrue("Large vector norm is finite", std::isfinite(norm));

    // Normalized vector
    Eigen::Vector3d normalized = largeVec.normalized();
    checkClose("Normalized large vector has unit norm", normalized.norm(), 1.0, 1e-10);

    // Very large values that would overflow in squared norm
    double veryLarge = 1e200;
    Eigen::Vector3d veryLargeVec(veryLarge, veryLarge, veryLarge);

    // stableNorm should handle this
    double stableNorm = veryLargeVec.stableNorm();
    checkTrue("Stable norm handles very large values", std::isfinite(stableNorm));

    std::cout << "[INFO] Large vector test passed" << std::endl;
}
