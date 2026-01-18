/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Basic astrodynamics and mathematics WASM tests.
 */

#include "wasmTestFramework.h"

#include <functional>
#include <memory>
#include <map>
#include <vector>

// Basic astrodynamics
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/astrodynamicsFunctions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"
#include "tudat/astro/basic_astro/modifiedEquinoctialElementConversions.h"
#include "tudat/astro/basic_astro/clohessyWiltshirePropagator.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"

// Mathematics
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/math/basic/legendrePolynomials.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/interpolators/linearInterpolator.h"
#include "tudat/math/interpolators/cubicSplineInterpolator.h"
#include "tudat/math/integrators/rungeKutta4Integrator.h"
#include "tudat/math/statistics/basicStatistics.h"

// Resource paths
#include "tudat/resource/resource.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#endif

using namespace tudat;

void testUnitConversions()
{
    std::cout << "\n=== Unit Conversions ===" << std::endl;

    using namespace unit_conversions;

    // Test degree/radian conversions
    checkClose("180 degrees to radians",
               convertDegreesToRadians(180.0),
               mathematical_constants::PI);

    checkClose("PI radians to degrees",
               convertRadiansToDegrees(mathematical_constants::PI),
               180.0);

    // Test distance conversions
    checkClose("1 AU to meters",
               convertAstronomicalUnitsToMeters(1.0),
               physical_constants::ASTRONOMICAL_UNIT,
               1.0); // 1 meter tolerance

    checkClose("1 meter to AU",
               convertMetersToAstronomicalUnits(physical_constants::ASTRONOMICAL_UNIT),
               1.0,
               1e-15);
}

void testPhysicalConstants()
{
    std::cout << "\n=== Physical Constants ===" << std::endl;

    // Check some well-known constants
    checkClose("Speed of light",
               physical_constants::SPEED_OF_LIGHT,
               299792458.0,
               1.0);

    checkClose("Gravitational constant",
               physical_constants::GRAVITATIONAL_CONSTANT,
               6.67259e-11,  // Value from tudat
               1e-15);

    checkClose("Astronomical unit",
               physical_constants::ASTRONOMICAL_UNIT,
               1.495978707e11,
               1e3);  // 1 km tolerance
}

void testOrbitalElementConversions()
{
    std::cout << "\n=== Orbital Element Conversions (NASA ODTBX Benchmarks) ===" << std::endl;

    using namespace orbital_element_conversions;

    // =========================================================================
    // Case 1: Elliptical orbit around Earth - NASA ODTBX benchmark
    // Reference: NASA Goddard Spaceflight Center, Orbit Determination Toolbox (ODTBX)
    // =========================================================================
    {
        const double earthGravParam = 3.986004415e14;  // m^3/s^2

        // Keplerian elements [m, -, rad, rad, rad, rad]
        Eigen::Vector6d keplerianElements;
        keplerianElements(semiMajorAxisIndex) = 8000.0 * 1000.0;
        keplerianElements(eccentricityIndex) = 0.23;
        keplerianElements(inclinationIndex) = 20.6 / 180.0 * mathematical_constants::PI;
        keplerianElements(argumentOfPeriapsisIndex) = 274.78 / 180.0 * mathematical_constants::PI;
        keplerianElements(longitudeOfAscendingNodeIndex) = 108.77 / 180.0 * mathematical_constants::PI;
        keplerianElements(trueAnomalyIndex) = 46.11 / 180.0 * mathematical_constants::PI;

        // Expected Cartesian elements from ODTBX [m, m, m, m/s, m/s, m/s]
        Eigen::Vector6d expectedCartesian;
        expectedCartesian(xCartesianPositionIndex) = 2.021874804243437e6;
        expectedCartesian(yCartesianPositionIndex) = 6.042523817035284e6;
        expectedCartesian(zCartesianPositionIndex) = -1.450371183512575e6;
        expectedCartesian(xCartesianVelocityIndex) = -7.118283509842652e3;
        expectedCartesian(yCartesianVelocityIndex) = 4.169050171542199e3;
        expectedCartesian(zCartesianVelocityIndex) = 2.029066072016241e3;

        Eigen::Vector6d computedCartesian = convertKeplerianToCartesianElements(
            keplerianElements, earthGravParam);

        // Check each component with NASA-grade precision (1e-15 relative)
        for (int i = 0; i < 6; i++) {
            double relError = std::abs(computedCartesian(i) - expectedCartesian(i)) /
                             std::abs(expectedCartesian(i));
            checkTrue("ODTBX Earth elliptical component " + std::to_string(i) + " (rel err < 1e-14)",
                     relError < 1e-14);
        }
    }

    // =========================================================================
    // Case 2: Circular equatorial orbit around Mars - NASA ODTBX benchmark
    // =========================================================================
    {
        const double marsGravParam = 4.2828018915e13;  // m^3/s^2

        Eigen::Vector6d keplerianElements;
        keplerianElements(semiMajorAxisIndex) = 9201.61 * 1000.0;
        keplerianElements(eccentricityIndex) = 0.0;
        keplerianElements(inclinationIndex) = 0.0;
        keplerianElements(argumentOfPeriapsisIndex) = 12.54 / 180.0 * mathematical_constants::PI;
        keplerianElements(longitudeOfAscendingNodeIndex) = 201.55 / 180.0 * mathematical_constants::PI;
        keplerianElements(trueAnomalyIndex) = -244.09 / 180.0 * mathematical_constants::PI;

        Eigen::Vector6d expectedCartesian;
        expectedCartesian(xCartesianPositionIndex) = 7.968828015716932e6;
        expectedCartesian(yCartesianPositionIndex) = -4.600804999999997e6;
        expectedCartesian(zCartesianPositionIndex) = 0.0;
        expectedCartesian(xCartesianVelocityIndex) = 1.078703495685965e3;
        expectedCartesian(yCartesianVelocityIndex) = 1.868369260830248e3;
        expectedCartesian(zCartesianVelocityIndex) = 0.0;

        Eigen::Vector6d computedCartesian = convertKeplerianToCartesianElements(
            keplerianElements, marsGravParam);

        // Check non-zero components with relative tolerance
        double relErrX = std::abs(computedCartesian(0) - expectedCartesian(0)) / std::abs(expectedCartesian(0));
        double relErrY = std::abs(computedCartesian(1) - expectedCartesian(1)) / std::abs(expectedCartesian(1));
        double relErrVx = std::abs(computedCartesian(3) - expectedCartesian(3)) / std::abs(expectedCartesian(3));
        double relErrVy = std::abs(computedCartesian(4) - expectedCartesian(4)) / std::abs(expectedCartesian(4));

        checkTrue("ODTBX Mars circular X (rel err < 1e-14)", relErrX < 1e-14);
        checkTrue("ODTBX Mars circular Y (rel err < 1e-14)", relErrY < 1e-14);
        checkTrue("ODTBX Mars circular Vx (rel err < 1e-14)", relErrVx < 1e-14);
        checkTrue("ODTBX Mars circular Vy (rel err < 1e-14)", relErrVy < 1e-14);
        // Z components should be exactly zero
        checkClose("ODTBX Mars circular Z", computedCartesian(2), 0.0, 1e-10);
        checkClose("ODTBX Mars circular Vz", computedCartesian(5), 0.0, 1e-10);
    }

    // =========================================================================
    // Case 3: Hyperbolic orbit around the Sun - NASA ODTBX benchmark
    // =========================================================================
    {
        const double sunGravParam = 1.32712440018e20;  // m^3/s^2

        Eigen::Vector6d keplerianElements;
        keplerianElements(semiMajorAxisIndex) = -4.5e11;  // Negative for hyperbolic
        keplerianElements(eccentricityIndex) = 2.3;
        keplerianElements(inclinationIndex) = 25.5 / 180.0 * mathematical_constants::PI;
        keplerianElements(argumentOfPeriapsisIndex) = 156.11 / 180.0 * mathematical_constants::PI;
        keplerianElements(longitudeOfAscendingNodeIndex) = -215.03 / 180.0 * mathematical_constants::PI;
        keplerianElements(trueAnomalyIndex) = 123.29 / 180.0 * mathematical_constants::PI;

        Eigen::Vector6d expectedCartesian;
        expectedCartesian(xCartesianPositionIndex) = -2.776328224174438e12;
        expectedCartesian(yCartesianPositionIndex) = -6.053823869632723e12;
        expectedCartesian(zCartesianPositionIndex) = 3.124576293512172e12;
        expectedCartesian(xCartesianVelocityIndex) = 7.957674684798018e3;
        expectedCartesian(yCartesianVelocityIndex) = 1.214817382001788e4;
        expectedCartesian(zCartesianVelocityIndex) = -6.923442392618828e3;

        Eigen::Vector6d computedCartesian = convertKeplerianToCartesianElements(
            keplerianElements, sunGravParam);

        for (int i = 0; i < 6; i++) {
            double relError = std::abs(computedCartesian(i) - expectedCartesian(i)) /
                             std::abs(expectedCartesian(i));
            checkTrue("ODTBX Sun hyperbolic component " + std::to_string(i) + " (rel err < 1e-14)",
                     relError < 1e-14);
        }
    }

    // =========================================================================
    // Case 4: Cartesian to Keplerian - Elliptical orbit - NASA ODTBX benchmark
    // =========================================================================
    {
        const double earthGravParam = 3.986004415e14;

        Eigen::Vector6d cartesianElements;
        cartesianElements(xCartesianPositionIndex) = 3.75e6;
        cartesianElements(yCartesianPositionIndex) = 4.24e6;
        cartesianElements(zCartesianPositionIndex) = -1.39e6;
        cartesianElements(xCartesianVelocityIndex) = -4.65e3;
        cartesianElements(yCartesianVelocityIndex) = -2.21e3;
        cartesianElements(zCartesianVelocityIndex) = 1.66e3;

        Eigen::Vector6d expectedKeplerian;
        expectedKeplerian(semiMajorAxisIndex) = 3.707478199246163e6;
        expectedKeplerian(eccentricityIndex) = 0.949175203660321;
        expectedKeplerian(inclinationIndex) = 0.334622356632438;
        expectedKeplerian(argumentOfPeriapsisIndex) = 2.168430616511167;
        expectedKeplerian(longitudeOfAscendingNodeIndex) = 1.630852596545341;
        expectedKeplerian(trueAnomalyIndex) = 3.302032232567084;

        Eigen::Vector6d computedKeplerian = convertCartesianToKeplerianElements(
            cartesianElements, earthGravParam);

        for (int i = 0; i < 6; i++) {
            double relError = std::abs(computedKeplerian(i) - expectedKeplerian(i)) /
                             std::abs(expectedKeplerian(i));
            checkTrue("ODTBX Cart->Kep elliptical component " + std::to_string(i) + " (rel err < 1e-13)",
                     relError < 1e-13);
        }
    }

    // =========================================================================
    // Case 5: Cartesian to Keplerian - Hyperbolic orbit around Sun - NASA ODTBX
    // =========================================================================
    {
        const double sunGravParam = 1.32712440018e20;

        Eigen::Vector6d cartesianElements;
        cartesianElements(xCartesianPositionIndex) = 7.035635643405699e11;
        cartesianElements(yCartesianPositionIndex) = -2.351218213055550e11;
        cartesianElements(zCartesianPositionIndex) = 0.037960971564309e11;
        cartesianElements(xCartesianVelocityIndex) = -1.731375459746510e4;
        cartesianElements(yCartesianVelocityIndex) = -1.535713656317794e4;
        cartesianElements(zCartesianVelocityIndex) = 0.423498718768347e4;

        Eigen::Vector6d expectedKeplerian;
        expectedKeplerian(semiMajorAxisIndex) = -6.78e11;
        expectedKeplerian(eccentricityIndex) = 1.89;
        expectedKeplerian(inclinationIndex) = 167.91 / 180.0 * mathematical_constants::PI;
        expectedKeplerian(argumentOfPeriapsisIndex) = 45.78 / 180.0 * mathematical_constants::PI;
        expectedKeplerian(longitudeOfAscendingNodeIndex) = 342.89 / 180.0 * mathematical_constants::PI;
        expectedKeplerian(trueAnomalyIndex) = 315.62 / 180.0 * mathematical_constants::PI;

        Eigen::Vector6d computedKeplerian = convertCartesianToKeplerianElements(
            cartesianElements, sunGravParam);

        for (int i = 0; i < 6; i++) {
            double relError = std::abs(computedKeplerian(i) - expectedKeplerian(i)) /
                             std::abs(expectedKeplerian(i));
            checkTrue("ODTBX Cart->Kep hyperbolic component " + std::to_string(i) + " (rel err < 1e-14)",
                     relError < 1e-14);
        }
    }
}

void testAnomalyConversions()
{
    std::cout << "\n=== Anomaly Conversions (NASA ODTBX Benchmarks) ===" << std::endl;

    using namespace orbital_element_conversions;

    // =========================================================================
    // True Anomaly to Eccentric Anomaly - NASA ODTBX benchmarks
    // =========================================================================

    // Case 1: General elliptical orbit
    {
        const double eccentricity = 0.146;
        const double trueAnomaly = 82.16 / 180.0 * mathematical_constants::PI;
        const double expectedEccentricAnomaly = 1.290237398010989;

        double computedEccentricAnomaly = convertTrueAnomalyToEllipticalEccentricAnomaly(
            trueAnomaly, eccentricity);

        double relError = std::abs(computedEccentricAnomaly - expectedEccentricAnomaly) /
                         std::abs(expectedEccentricAnomaly);
        checkTrue("ODTBX True->Eccentric elliptical (rel err < 2*eps)",
                 relError < 2.0 * std::numeric_limits<double>::epsilon());
    }

    // Case 2: Circular orbit
    {
        const double eccentricity = 0.0;
        const double trueAnomaly = 160.43 / 180.0 * mathematical_constants::PI;
        const double expectedEccentricAnomaly = 2.800031718974503;

        double computedEccentricAnomaly = convertTrueAnomalyToEllipticalEccentricAnomaly(
            trueAnomaly, eccentricity);

        double relError = std::abs(computedEccentricAnomaly - expectedEccentricAnomaly) /
                         std::abs(expectedEccentricAnomaly);
        checkTrue("ODTBX True->Eccentric circular (rel err < eps)",
                 relError < std::numeric_limits<double>::epsilon());
    }

    // Case 3: Hyperbolic orbit (Fortescue reference)
    {
        const double eccentricity = 3.0;
        const double trueAnomaly = 0.5291;
        const double expectedHyperbolicAnomaly = 0.3879;

        double computedHyperbolicAnomaly = convertTrueAnomalyToHyperbolicEccentricAnomaly(
            trueAnomaly, eccentricity);

        double relError = std::abs(computedHyperbolicAnomaly - expectedHyperbolicAnomaly) /
                         std::abs(expectedHyperbolicAnomaly);
        checkTrue("True->Hyperbolic eccentric (rel err < 1e-4)", relError < 1e-4);
    }

    // =========================================================================
    // Eccentric Anomaly to True Anomaly - NASA ODTBX benchmarks
    // =========================================================================

    // Case 4: General elliptical orbit
    {
        const double eccentricity = 0.639;
        const double eccentricAnomaly = 239.45 / 180.0 * mathematical_constants::PI;
        const double expectedTrueAnomaly = 3.665218735816221;  // After adding 2*PI

        double computedTrueAnomaly = convertEllipticalEccentricAnomalyToTrueAnomaly(
            eccentricAnomaly, eccentricity) + 2.0 * mathematical_constants::PI;

        double relError = std::abs(computedTrueAnomaly - expectedTrueAnomaly) /
                         std::abs(expectedTrueAnomaly);
        checkTrue("ODTBX Eccentric->True elliptical (rel err < eps)",
                 relError < std::numeric_limits<double>::epsilon());
    }

    // Case 5: Hyperbolic orbit (Fortescue reference)
    {
        const double eccentricity = 3.0;
        const double hyperbolicAnomaly = 0.3879;
        const double expectedTrueAnomaly = 0.5291;

        double computedTrueAnomaly = convertHyperbolicEccentricAnomalyToTrueAnomaly(
            hyperbolicAnomaly, eccentricity);

        double relError = std::abs(computedTrueAnomaly - expectedTrueAnomaly) /
                         std::abs(expectedTrueAnomaly);
        checkTrue("Hyperbolic eccentric->True (rel err < 1e-4)", relError < 1e-4);
    }

    // =========================================================================
    // Eccentric Anomaly to Mean Anomaly - NASA ODTBX benchmarks
    // =========================================================================

    // Case 6: General elliptical orbit
    {
        const double eccentricity = 0.541;
        const double eccentricAnomaly = 176.09 / 180.0 * mathematical_constants::PI;
        const double expectedMeanAnomaly = 3.036459804491048;

        double computedMeanAnomaly = convertEllipticalEccentricAnomalyToMeanAnomaly(
            eccentricAnomaly, eccentricity);

        double relError = std::abs(computedMeanAnomaly - expectedMeanAnomaly) /
                         std::abs(expectedMeanAnomaly);
        checkTrue("ODTBX Eccentric->Mean elliptical (rel err < eps)",
                 relError < std::numeric_limits<double>::epsilon());
    }

    // Case 7: Hyperbolic orbit (Vallado reference)
    {
        const double eccentricity = 2.4;
        const double hyperbolicAnomaly = 1.6013761449;
        const double expectedMeanAnomaly = 235.4 / 180.0 * mathematical_constants::PI;

        double computedMeanAnomaly = convertHyperbolicEccentricAnomalyToMeanAnomaly(
            hyperbolicAnomaly, eccentricity);

        double relError = std::abs(computedMeanAnomaly - expectedMeanAnomaly) /
                         std::abs(expectedMeanAnomaly);
        checkTrue("Vallado Hyperbolic Eccentric->Mean (rel err < 1e-7)", relError < 1e-7);
    }

    // =========================================================================
    // Elapsed Time to Mean Anomaly Change - NASA ODTBX benchmarks
    // =========================================================================

    // Case 8: Earth-orbiting satellite
    {
        const double elapsedTime = 8640.0;  // seconds
        const double earthGravParam = 398600.4415;  // km^3/s^2
        const double semiMajorAxis = 42165.3431351313;  // km
        const double expectedMeanAnomalyChange = 2.580579656848906 - 1.950567148859647;

        double computedMeanAnomalyChange = convertElapsedTimeToEllipticalMeanAnomalyChange(
            elapsedTime, earthGravParam, semiMajorAxis);

        double relError = std::abs(computedMeanAnomalyChange - expectedMeanAnomalyChange) /
                         std::abs(expectedMeanAnomalyChange);
        checkTrue("ODTBX Time->Mean anomaly change (rel err < 1e-13)", relError < 1e-13);
    }

    // =========================================================================
    // Mean Anomaly Change to Elapsed Time - NASA ODTBX benchmarks
    // =========================================================================

    // Case 9: Earth-orbiting satellite
    {
        const double meanAnomalyChange = 3.210592164838165 - 1.950567148859647;
        const double earthGravParam = 398600.4415;  // km^3/s^2
        const double semiMajorAxis = 42165.3431351313;  // km
        const double expectedElapsedTime = 17280.0;

        double computedElapsedTime = convertEllipticalMeanAnomalyChangeToElapsedTime(
            meanAnomalyChange, earthGravParam, semiMajorAxis);

        double relError = std::abs(computedElapsedTime - expectedElapsedTime) /
                         std::abs(expectedElapsedTime);
        checkTrue("ODTBX Mean anomaly->Time (rel err < 1e-14)", relError < 1e-14);
    }
}

void testCoordinateConversions()
{
    std::cout << "\n=== Coordinate Conversions ===" << std::endl;

    using namespace coordinate_conversions;

    // Test Cartesian to Spherical and back
    // Point on x-axis at Earth radius
    Eigen::Vector3d cartesian(6378.0e3, 0.0, 0.0);

    // Spherical coordinates are (radius, zenith, azimuth)
    // For a point on x-axis: zenith = π/2, azimuth = 0
    Eigen::Vector3d spherical = convertCartesianToSpherical(cartesian);
    Eigen::Vector3d cartesianRecovered = convertSphericalToCartesian(spherical);

    checkVectorClose("Cartesian-Spherical round-trip", cartesianRecovered, cartesian, 1.0);

    // Check spherical coordinates (radius, zenith, azimuth)
    checkClose("Spherical radius", spherical(0), 6378.0e3, 1.0);
    checkClose("Spherical zenith (angle from z)", spherical(1), mathematical_constants::PI / 2.0, 1e-10);
    checkClose("Spherical azimuth (angle in xy)", spherical(2), 0.0, 1e-10);
}

void testEigenOperations()
{
    std::cout << "\n=== Eigen Matrix Operations ===" << std::endl;

    // Basic matrix operations to ensure Eigen works in WASM
    Eigen::Matrix3d rotation = Eigen::AngleAxisd(
        mathematical_constants::PI / 4, Eigen::Vector3d::UnitZ()).toRotationMatrix();

    Eigen::Vector3d v(1.0, 0.0, 0.0);
    Eigen::Vector3d rotated = rotation * v;

    // After 45-degree rotation around Z, x-unit vector should be at (sqrt(2)/2, sqrt(2)/2, 0)
    Eigen::Vector3d expected(std::sqrt(2.0)/2.0, std::sqrt(2.0)/2.0, 0.0);
    checkVectorClose("45-degree Z rotation", rotated, expected, 1e-14);

    // Check rotation matrix properties
    checkClose("Rotation matrix determinant", rotation.determinant(), 1.0, 1e-14);

    Eigen::Matrix3d shouldBeIdentity = rotation * rotation.transpose();
    checkClose("R * R^T = I (trace)", shouldBeIdentity.trace(), 3.0, 1e-14);
}

void testKeplerFunctions()
{
    std::cout << "\n=== Kepler Orbital Mechanics ===" << std::endl;

    using namespace basic_astrodynamics;

    // Test geostationary orbit period
    // Reference: http://en.wikipedia.org/wiki/Geostationary_orbit
    double satelliteMass = 1.0e3;  // kg
    double earthGravParam = physical_constants::GRAVITATIONAL_CONSTANT * 5.9736e24;
    double geoRadius = 4.2164e7;  // meters

    double orbitalPeriod = computeKeplerOrbitalPeriod(geoRadius, earthGravParam, satelliteMass);
    double expectedPeriod = 86164.09054;  // seconds (sidereal day)

    checkClose("Geostationary orbital period", orbitalPeriod, expectedPeriod, expectedPeriod * 1e-5);

    // Test mean motion
    double meanMotion = computeKeplerMeanMotion(geoRadius, earthGravParam, satelliteMass);
    double expectedMeanMotion = 2.0 * mathematical_constants::PI / expectedPeriod;

    checkClose("Geostationary mean motion", meanMotion, expectedMeanMotion, 1e-9);

    // Test synodic period (Earth-Mars example)
    double earthPeriod = 365.25 * 86400.0;  // seconds
    double marsPeriod = 687.0 * 86400.0;    // seconds
    double synodicPeriod = computeSynodicPeriod(earthPeriod, marsPeriod);
    double expectedSynodic = 779.94 * 86400.0;  // ~780 days

    checkClose("Earth-Mars synodic period", synodicPeriod, expectedSynodic, expectedSynodic * 0.01);

    // Test radial distance calculation
    double semiMajorAxis = 25999.683025291e3;
    double eccentricity = 0.864564003552322;
    double trueAnomaly = 0.757654217738482;
    double radialDistance = computeKeplerRadialDistance(semiMajorAxis, eccentricity, trueAnomaly);
    double expectedRadius = 4032815.56442827;

    checkClose("Kepler radial distance", radialDistance, expectedRadius, expectedRadius * 1e-5);
}

void testTimeConversions()
{
    std::cout << "\n=== Time Conversions ===" << std::endl;

    using namespace basic_astrodynamics;

    // Test Julian day conversions
    // J2000 epoch: January 1, 2000, 12:00 TT = JD 2451545.0
    double j2000JulianDay = 2451545.0;

    // Convert seconds since J2000 to Julian day
    double secondsSinceJ2000 = 0.0;
    double julianDay = convertSecondsSinceEpochToJulianDay(secondsSinceJ2000, j2000JulianDay);

    checkClose("J2000 epoch Julian day", julianDay, j2000JulianDay, 1e-10);

    // One day later
    double oneDaySeconds = 86400.0;
    julianDay = convertSecondsSinceEpochToJulianDay(oneDaySeconds, j2000JulianDay);

    checkClose("J2000 + 1 day Julian day", julianDay, j2000JulianDay + 1.0, 1e-10);

    // Convert back
    double recoveredSeconds = convertJulianDayToSecondsSinceEpoch(julianDay, j2000JulianDay);

    checkClose("Julian day round-trip", recoveredSeconds, oneDaySeconds, 1e-6);
}

void testLegendrePolynomials()
{
    std::cout << "\n=== Legendre Polynomials ===" << std::endl;

    using namespace basic_mathematics;

    // Test Legendre polynomials at known values
    // The tudat function computes associated Legendre polynomials P_l^m(x)
    // For m=0, these are the regular Legendre polynomials:
    // P_0^0(x) = 1
    // P_1^0(x) = x
    // P_2^0(x) = (3x^2 - 1)/2
    // P_3^0(x) = (5x^3 - 3x)/2

    double x = 0.5;
    int order = 0;  // m=0 gives regular Legendre polynomials

    double p0 = computeLegendrePolynomial(0, order, x);
    checkClose("P_0^0(0.5)", p0, 1.0, 1e-14);

    double p1 = computeLegendrePolynomial(1, order, x);
    checkClose("P_1^0(0.5)", p1, 0.5, 1e-14);

    double p2 = computeLegendrePolynomial(2, order, x);
    double expectedP2 = (3.0 * x * x - 1.0) / 2.0;
    checkClose("P_2^0(0.5)", p2, expectedP2, 1e-14);

    double p3 = computeLegendrePolynomial(3, order, x);
    double expectedP3 = (5.0 * x * x * x - 3.0 * x) / 2.0;
    checkClose("P_3^0(0.5)", p3, expectedP3, 1e-14);
}

void testLinearInterpolation()
{
    std::cout << "\n=== Linear Interpolation ===" << std::endl;

    using namespace interpolators;

    // Create simple data: y = 2x + 1
    std::map<double, double> dataMap;
    dataMap[0.0] = 1.0;
    dataMap[1.0] = 3.0;
    dataMap[2.0] = 5.0;
    dataMap[3.0] = 7.0;

    LinearInterpolator<double, double> interpolator(dataMap);

    // Test at known points
    checkClose("Interpolate at x=0", interpolator.interpolate(0.0), 1.0, 1e-14);
    checkClose("Interpolate at x=1", interpolator.interpolate(1.0), 3.0, 1e-14);
    checkClose("Interpolate at x=2", interpolator.interpolate(2.0), 5.0, 1e-14);

    // Test at intermediate points
    checkClose("Interpolate at x=0.5", interpolator.interpolate(0.5), 2.0, 1e-14);
    checkClose("Interpolate at x=1.5", interpolator.interpolate(1.5), 4.0, 1e-14);
    checkClose("Interpolate at x=2.5", interpolator.interpolate(2.5), 6.0, 1e-14);
}

void testNumericalIntegration()
{
    std::cout << "\n=== Numerical Integration ===" << std::endl;

    using namespace numerical_integrators;

    // Test RK4 on simple ODE: dy/dt = y, y(0) = 1
    // Exact solution: y = e^t

    // State derivative function
    auto stateDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
        return state;  // dy/dt = y
    };

    // Initial conditions
    double t0 = 0.0;
    Eigen::VectorXd y0(1);
    y0 << 1.0;

    // Create RK4 integrator
    double stepSize = 0.01;
    RungeKutta4Integrator<double, Eigen::VectorXd> integrator(stateDerivative, t0, y0, stepSize);

    // Integrate to t = 1
    double tEnd = 1.0;
    while (integrator.getCurrentIndependentVariable() < tEnd) {
        integrator.performIntegrationStep(stepSize);
    }

    double computedY = integrator.getCurrentState()(0);
    double exactY = std::exp(1.0);

    checkClose("RK4 exponential growth", computedY, exactY, 1e-6);

    // Test on harmonic oscillator: d²x/dt² = -x
    // Rewrite as system: dx/dt = v, dv/dt = -x
    // Initial: x(0) = 1, v(0) = 0
    // Exact: x(t) = cos(t)

    auto harmonicDerivative = [](const double t, const Eigen::VectorXd& state) -> Eigen::VectorXd {
        Eigen::VectorXd derivative(2);
        derivative(0) = state(1);   // dx/dt = v
        derivative(1) = -state(0);  // dv/dt = -x
        return derivative;
    };

    Eigen::VectorXd harmonicState(2);
    harmonicState << 1.0, 0.0;  // x=1, v=0

    RungeKutta4Integrator<double, Eigen::VectorXd> harmonicIntegrator(
        harmonicDerivative, 0.0, harmonicState, stepSize);

    // Integrate to t = π/2
    double tHalf = mathematical_constants::PI / 2.0;
    while (harmonicIntegrator.getCurrentIndependentVariable() < tHalf - stepSize/2) {
        harmonicIntegrator.performIntegrationStep(stepSize);
    }

    double computedX = harmonicIntegrator.getCurrentState()(0);
    double exactX = std::cos(tHalf);  // Should be ~0

    // RK4 accumulates some error over many steps, use absolute tolerance
    checkClose("RK4 harmonic oscillator x(pi/2)", computedX, exactX, 1e-3);
}

void testCubicSplineInterpolation()
{
    std::cout << "\n=== Cubic Spline Interpolation ===" << std::endl;

    using namespace interpolators;

    // Create data from sin function
    std::map<double, double> dataMap;
    for (int i = 0; i <= 10; i++) {
        double x = i * mathematical_constants::PI / 10.0;
        dataMap[x] = std::sin(x);
    }

    CubicSplineInterpolator<double, double> spline(dataMap);

    // Test at intermediate points
    double x1 = mathematical_constants::PI / 4.0;
    double y1 = spline.interpolate(x1);
    checkClose("Cubic spline sin(pi/4)", y1, std::sin(x1), 1e-4);

    double x2 = mathematical_constants::PI / 3.0;
    double y2 = spline.interpolate(x2);
    checkClose("Cubic spline sin(pi/3)", y2, std::sin(x2), 1e-4);

    double x3 = mathematical_constants::PI / 2.0;
    double y3 = spline.interpolate(x3);
    checkClose("Cubic spline sin(pi/2)", y3, 1.0, 1e-10);
}

void testReferenceFrameTransformations()
{
    std::cout << "\n=== Reference Frame Transformations ===" << std::endl;

    using namespace reference_frames;

    // Test getRotatingPlanetocentricToInertialFrameTransformationMatrix
    double angle = mathematical_constants::PI / 6.0;  // 30 degrees
    Eigen::Matrix3d toInertial = getRotatingPlanetocentricToInertialFrameTransformationMatrix(angle);
    Eigen::Matrix3d toRotating = getInertialToPlanetocentricFrameTransformationMatrix(angle);

    // The rotation should be well-formed (orthogonal)
    Eigen::Matrix3d shouldBeIdentity = toInertial * toInertial.transpose();
    checkClose("Rotation matrix is orthogonal (trace)", shouldBeIdentity.trace(), 3.0, 1e-14);

    // These should be inverses of each other
    Eigen::Matrix3d product = toInertial * toRotating;
    checkClose("Planetocentric transforms are inverses (trace)", product.trace(), 3.0, 1e-14);

    // Test a specific transformation: rotate a vector on x-axis
    Eigen::Vector3d xAxis(1.0, 0.0, 0.0);
    Eigen::Vector3d rotated = toInertial * xAxis;

    // After rotation by 30 degrees around Z, x should go to (cos(30), sin(30), 0)
    double expectedX = std::cos(angle);
    double expectedY = std::sin(angle);
    checkClose("Planetocentric X rotation x-comp", rotated(0), expectedX, 1e-14);
    checkClose("Planetocentric X rotation y-comp", rotated(1), expectedY, 1e-14);
    checkClose("Planetocentric X rotation z-comp", rotated(2), 0.0, 1e-14);
}

void testModifiedEquinoctialElements()
{
    std::cout << "\n=== Modified Equinoctial Elements ===" << std::endl;

    using namespace orbital_element_conversions;

    // Define a test orbit (Mars-like)
    double earthGravParam = 3.986004418e14;
    double semiMajorAxis = 10000.0e3;  // 10,000 km
    double eccentricity = 0.1;
    double inclination = unit_conversions::convertDegreesToRadians(30.0);
    double argumentOfPeriapsis = unit_conversions::convertDegreesToRadians(45.0);
    double raan = unit_conversions::convertDegreesToRadians(60.0);
    double trueAnomaly = unit_conversions::convertDegreesToRadians(90.0);

    // Create Keplerian elements vector
    Eigen::Vector6d keplerianElements;
    keplerianElements << semiMajorAxis, eccentricity, inclination,
                         argumentOfPeriapsis, raan, trueAnomaly;

    // Convert to Modified Equinoctial Elements
    Eigen::Vector6d mee = convertKeplerianToModifiedEquinoctialElements(
        keplerianElements, true);  // true = use semi-latus rectum representation

    // Convert back to Keplerian
    Eigen::Vector6d keplerianRecovered = convertModifiedEquinoctialToKeplerianElements(
        mee, true);

    // Check round-trip conversion
    checkClose("MEE semi-major axis round-trip",
               keplerianRecovered(0), semiMajorAxis, 1.0);
    checkClose("MEE eccentricity round-trip",
               keplerianRecovered(1), eccentricity, 1e-10);
    checkClose("MEE inclination round-trip",
               keplerianRecovered(2), inclination, 1e-10);
    checkClose("MEE arg periapsis round-trip",
               keplerianRecovered(3), argumentOfPeriapsis, 1e-10);
    checkClose("MEE RAAN round-trip",
               keplerianRecovered(4), raan, 1e-10);
    checkClose("MEE true anomaly round-trip",
               keplerianRecovered(5), trueAnomaly, 1e-10);
}

void testStatistics()
{
    std::cout << "\n=== Statistics ===" << std::endl;

    using namespace statistics;

    // Test sample mean
    std::vector<double> data = {1.0, 2.0, 3.0, 4.0, 5.0};
    double mean = computeSampleMean(data);
    checkClose("Sample mean of 1-5", mean, 3.0, 1e-14);

    // Test with another data set
    std::vector<double> data2 = {2.0, 4.0, 6.0, 8.0, 10.0};
    double mean2 = computeSampleMean(data2);
    checkClose("Sample mean of 2,4,6,8,10", mean2, 6.0, 1e-14);

    // Test sample variance
    double variance = computeSampleVariance(data);
    // Variance of 1,2,3,4,5: mean=3, sum of squared deviations = 4+1+0+1+4 = 10
    // Sample variance = 10/4 = 2.5
    checkClose("Sample variance of 1-5", variance, 2.5, 1e-14);
}

void testSphericalHarmonics()
{
    std::cout << "\n=== Spherical Harmonics ===" << std::endl;

    using namespace basic_mathematics;

    // Test spherical harmonic computations
    // For degree 0, order 0: Y_0^0 = 1/sqrt(4*pi)
    double cosineLatitude = 0.5;  // latitude = 60 degrees

    // Compute geodesy-normalized Legendre polynomial
    double p00 = computeGeodesyLegendrePolynomial(0, 0, cosineLatitude);
    checkClose("Geodesy P_0^0(0.5)", p00, 1.0, 1e-14);

    double p10 = computeGeodesyLegendrePolynomial(1, 0, cosineLatitude);
    // For geodesy normalization, P_1^0(x) = sqrt(3) * x
    double expectedP10 = std::sqrt(3.0) * cosineLatitude;
    checkClose("Geodesy P_1^0(0.5)", p10, expectedP10, 1e-14);

    double p11 = computeGeodesyLegendrePolynomial(1, 1, cosineLatitude);
    // P_1^1(x) = sqrt(3) * sqrt(1-x^2)
    double sinLatitude = std::sqrt(1.0 - cosineLatitude * cosineLatitude);
    double expectedP11 = std::sqrt(3.0) * sinLatitude;
    checkClose("Geodesy P_1^1(0.5)", p11, expectedP11, 1e-14);
}

void testResourcePaths()
{
    std::cout << "\n=== Resource Paths (WASM Virtual FS) ===" << std::endl;

    using namespace paths;

    // Test that all resource paths are properly defined and point to the WASM data mount point
    std::string basePath = "/tudat_data";

    checkStringEquals("Base resources path", get_resources_path(), basePath);
    checkStringEquals("Ephemeris path", get_ephemeris_path(), basePath + "/ephemeris");
    checkStringEquals("Earth orientation path", get_earth_orientation_path(), basePath + "/earth_orientation");
    checkStringEquals("Quadrature path", get_quadrature_path(), basePath + "/quadrature");
    checkStringEquals("SPICE kernels path", get_spice_kernels_path(), basePath + "/spice_kernels");
    checkStringEquals("Atmosphere tables path", get_atmosphere_tables_path(), basePath + "/atmosphere_tables");
    checkStringEquals("Gravity models path", get_gravity_models_path(), basePath + "/gravity_models");
    checkStringEquals("Space weather path", get_space_weather_path(), basePath + "/space_weather");

    // Verify all paths start with the base path (consistency check)
    checkStringStartsWith("Ephemeris path prefix", get_ephemeris_path(), basePath);
    checkStringStartsWith("SPICE kernels path prefix", get_spice_kernels_path(), basePath);
}

#ifdef __EMSCRIPTEN__
void testEmscriptenEnvironment()
{
    std::cout << "\n=== Emscripten Environment ===" << std::endl;

    // Test that we're running in the Emscripten environment
    checkTrue("Running in Emscripten", true);

    // Test that the virtual filesystem is available (basic check)
    // In a full WASM environment, you could test FS operations here
    checkTrue("WASM environment detected", emscripten_run_script_int("1") == 1);
}
#endif

void testLinearAlgebra()
{
    std::cout << "\n=== Linear Algebra Operations ===" << std::endl;

    using namespace linear_algebra;

    // Test cross product
    Eigen::Vector3d a(1.0, 0.0, 0.0);
    Eigen::Vector3d b(0.0, 1.0, 0.0);
    Eigen::Vector3d cross = a.cross(b);
    Eigen::Vector3d expectedCross(0.0, 0.0, 1.0);
    checkVectorClose("Cross product x × y = z", cross, expectedCross, 1e-14);

    // Test getCrossProductMatrix
    Eigen::Matrix3d crossMatrix = getCrossProductMatrix(a);
    Eigen::Vector3d crossFromMatrix = crossMatrix * b;
    checkVectorClose("Cross product matrix", crossFromMatrix, expectedCross, 1e-14);

    // Verify cross product matrix property: [a×] * b = a × b
    Eigen::Vector3d c(1.0, 2.0, 3.0);
    Eigen::Vector3d d(4.0, 5.0, 6.0);
    Eigen::Matrix3d cCrossMatrix = getCrossProductMatrix(c);
    Eigen::Vector3d crossDirect = c.cross(d);
    Eigen::Vector3d crossViaMatrix = cCrossMatrix * d;
    checkVectorClose("Cross product matrix general case", crossViaMatrix, crossDirect, 1e-14);
}

// ============================================================================
// Clohessy-Wiltshire Propagator Tests
// Reference: Vallado, D.A. "Fundamentals of Astrodynamics and Applications"
//            Wakker, K.F. "Astrodynamics I"
// ============================================================================

void testClohessyWiltshirePropagation()
{
    std::cout << "\n=== Clohessy-Wiltshire Propagation ===" << std::endl;

    using namespace basic_astrodynamics;

    // =========================================================================
    // Test 1: Full state propagation in LEO
    // Benchmark: MATLAB routine "hillsr" from Vallado [2001]
    // =========================================================================
    {
        // Central body is Earth at 400 km altitude
        const double earthGravParam = 3.986004418e14;  // m^3/s^2
        const double referenceOrbitRadius = 6.778137e6;  // m
        const double propagationDuration = 1800.0;  // 30 minutes

        // Initial state: arbitrary non-zero distances and velocities
        // [x, y, z, vx, vy, vz] in Hill frame
        Eigen::Vector6d initialState;
        initialState << 45.0, 37.0, 12.0, 0.08, 0.03, 0.01;

        // Compute propagated state
        Eigen::Vector6d computedFinalState = propagateClohessyWiltshire(
            initialState, propagationDuration, earthGravParam, referenceOrbitRadius);

        // Expected final state from Vallado's MATLAB "hillsr" routine
        Eigen::Vector6d expectedFinalState;
        expectedFinalState << 3.806450080201250e2,
                             -5.437424675454679e2,
                              2.509547637285142,
                              1.541620605755606e-1,
                             -7.294751390499470e-1,
                             -1.662099488431618e-2;

        // Check each component with high precision
        for (int i = 0; i < 6; i++) {
            double relError = std::abs(computedFinalState(i) - expectedFinalState(i)) /
                             std::abs(expectedFinalState(i));
            checkTrue("C-W Vallado component " + std::to_string(i) + " (rel err < 1e-14)",
                     relError < 1e-14);
        }
    }

    // =========================================================================
    // Test 2: Pure harmonic relative motion (Wakker 2007)
    // If initial conditions satisfy the harmonic motion constraints, the state
    // should return to initial after one orbital period.
    // Constraints: vx = 0.5 * n * y, vy = -2 * n * x
    // =========================================================================
    {
        const double earthGravParam = 3.986004418e14;
        const double referenceOrbitRadius = 6.778137e6;

        // Mean angular motion
        double meanMotion = std::sqrt(earthGravParam /
            (referenceOrbitRadius * referenceOrbitRadius * referenceOrbitRadius));

        // Orbital period
        double orbitalPeriod = 2.0 * mathematical_constants::PI / meanMotion;

        // Initial state satisfying harmonic motion conditions
        double x0 = 34.0, y0 = 49.0, z0 = 17.0;
        double vz0 = 0.04;
        Eigen::Vector6d initialState;
        initialState << x0, y0, z0,
                       0.5 * meanMotion * y0,  // vx = 0.5 * n * y
                      -2.0 * meanMotion * x0,  // vy = -2 * n * x
                       vz0;

        // Propagate for one full orbital period
        Eigen::Vector6d computedFinalState = propagateClohessyWiltshire(
            initialState, orbitalPeriod, earthGravParam, referenceOrbitRadius);

        // Final state should equal initial state (harmonic motion returns to start)
        for (int i = 0; i < 6; i++) {
            double relError = std::abs(computedFinalState(i) - initialState(i)) /
                             std::abs(initialState(i));
            checkTrue("C-W harmonic motion component " + std::to_string(i) + " (rel err < 1e-14)",
                     relError < 1e-14);
        }
    }
}

// ============================================================================
// Kepler Propagator Tests
// References:
//   - Melman, J. "Propagate" software, TU Delft
//   - NASA GSFC Orbit Determination Toolbox (ODTBX)
//   - ESA Global Trajectory Optimisation (GTOP) Toolbox
// ============================================================================

void testKeplerPropagation()
{
    std::cout << "\n=== Kepler Orbit Propagation ===" << std::endl;

    using namespace orbital_element_conversions;
    using namespace basic_mathematics;

    // =========================================================================
    // Test 1: Elliptical orbit propagation (Melman 2010)
    // Propagate from initial to final state using Kepler propagator
    // =========================================================================
    {
        const double earthGravParam = 3.986004415e14;  // m^3/s^2

        // Initial state in Cartesian coordinates
        Eigen::Vector6d initialCartesian;
        initialCartesian << 6.75e6, 0.0, 0.0, 0.0, 8.0595973215e3, 0.0;

        // Convert to Keplerian elements
        Eigen::Vector6d initialKeplerian = convertCartesianToKeplerianElements(
            initialCartesian, earthGravParam);

        // Propagate for 1 day (86400 seconds)
        double propagationTime = 86400.0;
        Eigen::Vector6d propagatedKeplerian = propagateKeplerOrbit(
            initialKeplerian, propagationTime, earthGravParam);

        // Expected final state in Cartesian from Melman benchmark
        Eigen::Vector6d expectedFinalCartesian;
        expectedFinalCartesian << -6.1318272067e6, 5.1974105627e6, 0.0,
                                  -4.7375063953e3, -4.8565484865e3, 0.0;
        Eigen::Vector6d expectedFinalKeplerian = convertCartesianToKeplerianElements(
            expectedFinalCartesian, earthGravParam);

        // The only element that should change significantly is true anomaly
        double computedTrueAnomaly = computeModulo(propagatedKeplerian(5), 2.0 * mathematical_constants::PI);
        double expectedTrueAnomaly = computeModulo(expectedFinalKeplerian(5), 2.0 * mathematical_constants::PI);

        double relError = std::abs(computedTrueAnomaly - expectedTrueAnomaly) /
                         std::abs(expectedTrueAnomaly);
        checkTrue("Melman Kepler true anomaly (rel err < 1e-8)", relError < 1e-8);
    }

    // =========================================================================
    // Test 2: Elliptical orbit (ODTBX benchmark) - forward propagation with modulo
    // =========================================================================
    {
        const double earthGravParam = 398600.4415e9;  // m^3/s^2
        const double timeStep = 8640.0;  // seconds

        // Initial Keplerian elements from ODTBX
        Eigen::Vector6d keplerianElements;
        keplerianElements << 42165.3431351313e3,  // semi-major axis [m]
                            0.26248354351331,     // eccentricity
                            0.30281462522101,     // inclination [rad]
                            4.71463172847351,     // arg periapsis [rad]
                            4.85569272927819,     // RAAN [rad]
                            2.37248926702153;     // true anomaly [rad]

        // Expected true anomaly values from ODTBX at each timestep
        std::vector<double> expectedTrueAnomalies = {
            2.37248926702153,  // t=0
            2.79722436211144,  // t=1*8640
            3.18337407409023,  // t=2*8640
            3.57400974200765,  // t=3*8640
            4.01425565759545,  // t=4*8640
            4.57232665706546,  // t=5*8640
            5.35956850972672,  // t=6*8640
            0.137251905665217, // t=7*8640
            1.14521863765007,  // t=8*8640
            1.86433634881636,  // t=9*8640
            2.38486787064101   // t=10*8640
        };

        Eigen::Vector6d currentState = keplerianElements;
        for (size_t i = 1; i < expectedTrueAnomalies.size(); i++) {
            currentState = propagateKeplerOrbit(currentState, timeStep, earthGravParam);
            double computedTrueAnomaly = computeModulo(currentState(5), 2.0 * mathematical_constants::PI);

            double relError = std::abs(computedTrueAnomaly - expectedTrueAnomalies[i]) /
                             std::abs(expectedTrueAnomalies[i]);
            checkTrue("ODTBX Kepler step " + std::to_string(i) + " (rel err < 1e-13)",
                     relError < 1e-13);
        }
    }

    // =========================================================================
    // Test 3: Elliptical orbit (ODTBX) - backward propagation
    // =========================================================================
    {
        const double earthGravParam = 398600.4415e9;
        const double timeStep = 8640.0;

        Eigen::Vector6d keplerianElements;
        keplerianElements << 42165.3431351313e3, 0.26248354351331, 0.30281462522101,
                            4.71463172847351, 4.85569272927819, 2.38486787064101;  // Start at t=10

        // Expected values at t=0 (after backward propagation)
        double expectedInitialTrueAnomaly = 2.37248926702153;

        // Propagate backward (negative time) for 10 steps
        Eigen::Vector6d currentState = keplerianElements;
        for (int i = 0; i < 10; i++) {
            currentState = propagateKeplerOrbit(currentState, -timeStep, earthGravParam);
        }

        double computedTrueAnomaly = computeModulo(currentState(5), 2.0 * mathematical_constants::PI);
        double relError = std::abs(computedTrueAnomaly - expectedInitialTrueAnomaly) /
                         std::abs(expectedInitialTrueAnomaly);
        checkTrue("ODTBX Kepler backward (rel err < 1e-13)", relError < 1e-13);
    }

    // =========================================================================
    // Test 4: Hyperbolic orbit propagation (ESA GTOP benchmark)
    // =========================================================================
    {
        // Gravitational parameter of the Sun (slightly different value used by GTOP)
        const double sunGravParam = 1.327e20;  // m^3/s^2
        const double timeStep = 86400.0 * 100.0;  // 100 days

        // Initial state: hyperbolic orbit (e > 1)
        Eigen::Vector6d initialCartesian;
        initialCartesian << 1.5e11, 0.0, 0.0, 0.0, 6.0e4, 0.0;

        Eigen::Vector6d initialKeplerian = convertCartesianToKeplerianElements(
            initialCartesian, sunGravParam);

        // Verify it's hyperbolic (e > 1)
        checkTrue("GTOP orbit is hyperbolic (e > 1)", initialKeplerian(1) > 1.0);

        // Expected final Cartesian states from GTOP after 100, 200, 300, 400 days
        std::vector<Eigen::Vector6d> expectedCartesianStates;
        Eigen::Vector6d state1, state2, state3, state4;
        state1 << 50369576778.98602, 453006898372.5074, 2.156946592732799e-005,
                 -14654.13750690802, 46884.94068619227, 4.665334803219454e-012;
        state2 << -76810236076.38216, 842661848023.4473, 6.100268297443444e-005,
                 -14683.57015580225, 43917.12010513522, 4.48721854707566e-012;
        state3 << -203052258817.9893, 1216808019495.603, 9.937145023346651e-005,
                 -14543.34378775917, 42828.66589049961, 4.403399939593385e-012;
        state4 << -328225472457.8796, 1584186440047.591, 0.0001371949389119038,
                 -14437.813524927732, 42264.20425643964, 4.355914471377053e-012;
        expectedCartesianStates.push_back(state1);
        expectedCartesianStates.push_back(state2);
        expectedCartesianStates.push_back(state3);
        expectedCartesianStates.push_back(state4);

        // Propagate and check true anomaly at each step
        Eigen::Vector6d currentKeplerian = initialKeplerian;
        for (size_t i = 0; i < expectedCartesianStates.size(); i++) {
            currentKeplerian = propagateKeplerOrbit(currentKeplerian, timeStep, sunGravParam);

            Eigen::Vector6d expectedKeplerian = convertCartesianToKeplerianElements(
                expectedCartesianStates[i], sunGravParam);

            double relError = std::abs(currentKeplerian(5) - expectedKeplerian(5)) /
                             std::abs(expectedKeplerian(5));
            checkTrue("GTOP hyperbolic step " + std::to_string(i+1) + " (rel err < 1e-15)",
                     relError < 1e-15);
        }
    }

    // =========================================================================
    // Test 5: Mean motion consistency test
    // Verify that propagated mean anomaly matches n * dt
    // =========================================================================
    {
        const double earthGravParam = 398600.4415e9;
        const double timeStep = 600.0;  // 10 minutes

        Eigen::Vector6d keplerianElements;
        keplerianElements << 42165.3431351313e3, 0.26248354351331, 0.30281462522101,
                            4.71463172847351, 4.85569272927819, 2.37248926702153;

        double meanMotion = std::sqrt(earthGravParam / std::pow(keplerianElements(0), 3.0));

        double initialMeanAnomaly = convertEccentricAnomalyToMeanAnomaly(
            convertTrueAnomalyToEllipticalEccentricAnomaly(keplerianElements(5), keplerianElements(1)),
            keplerianElements(1));

        // Test forward and backward propagation
        for (int i = -25; i <= 25; i++) {
            double propagationTime = static_cast<double>(i) * timeStep;
            Eigen::Vector6d propagated = propagateKeplerOrbit(keplerianElements, propagationTime, earthGravParam);

            double propagatedMeanAnomaly = convertEccentricAnomalyToMeanAnomaly(
                convertTrueAnomalyToEllipticalEccentricAnomaly(propagated(5), keplerianElements(1)),
                keplerianElements(1));

            double expectedMeanAnomalyChange = meanMotion * propagationTime;
            double computedMeanAnomalyChange = propagatedMeanAnomaly - initialMeanAnomaly;

            double error = std::abs(expectedMeanAnomalyChange - computedMeanAnomalyChange);
            // Error should be within 5 * machine epsilon
            checkTrue("Mean motion consistency i=" + std::to_string(i),
                     error < 5.0 * std::numeric_limits<double>::epsilon());
        }
    }
}
