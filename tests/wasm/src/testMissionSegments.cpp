/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Mission segments WASM tests - Lambert targeting.
 */

#include "wasmTestFramework.h"

// Basic astrodynamics
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/physicalConstants.h"

// Mission segments
#include "tudat/astro/mission_segments/lambertTargeterIzzo.h"
#include "tudat/astro/mission_segments/gravityAssist.h"
#include "tudat/astro/mission_segments/escapeAndCapture.h"
#include "tudat/astro/basic_astro/unitConversions.h"

using namespace tudat;

void testLambertTargetingIzzo()
{
    std::cout << "\n=== Lambert Targeting (Izzo Algorithm) ===" << std::endl;

    using namespace mission_segments;
    using namespace orbital_element_conversions;

    // Test 1: Elliptical case (Earth orbit)
    // Reference: Mengali & Quarta "Fondamenti di Meccanica del volo Spaziale", Example 6.1
    {
        double distanceUnit = 6.378136e6;  // Earth radius [m]
        double timeUnit = 806.78;  // Canonical time unit [s]
        double gravitationalParameter = 398600.4418e9;  // Earth GM [m^3/s^2]

        Eigen::Vector3d departurePosition(2.0 * distanceUnit, 0.0, 0.0);
        Eigen::Vector3d arrivalPosition(2.0 * distanceUnit, 2.0 * std::sqrt(3.0) * distanceUnit, 0.0);
        double timeOfFlight = 5.0 * timeUnit;

        LambertTargeterIzzo lambertTargeter(departurePosition, arrivalPosition,
                                            timeOfFlight, gravitationalParameter);

        Eigen::Vector3d departureVelocity = lambertTargeter.getInertialVelocityAtDeparture();
        Eigen::Vector3d arrivalVelocity = lambertTargeter.getInertialVelocityAtArrival();

        // Expected velocities from reference
        Eigen::Vector3d expectedDepartureVelocity(2735.8, 6594.3, 0.0);
        Eigen::Vector3d expectedArrivalVelocity(-1367.9, 4225.03, 0.0);

        // Check departure velocity
        for (int i = 0; i < 3; i++) {
            if (std::abs(expectedDepartureVelocity(i)) > 1e-10) {
                checkClose("Lambert elliptical departure V" + std::to_string(i),
                          departureVelocity(i), expectedDepartureVelocity(i), 1.0);
            }
        }

        // Check arrival velocity
        for (int i = 0; i < 3; i++) {
            if (std::abs(expectedArrivalVelocity(i)) > 1e-10) {
                checkClose("Lambert elliptical arrival V" + std::to_string(i),
                          arrivalVelocity(i), expectedArrivalVelocity(i), 1.0);
            }
        }

        // Verify orbit is prograde (positive angular momentum in z)
        Eigen::Vector3d angularMomentum = departurePosition.cross(departureVelocity);
        checkTrue("Lambert elliptical orbit is prograde", angularMomentum(2) > 0.0);
    }

    // Test 2: Hyperbolic case
    // Reference: unitTestLambertTargeterIzzo.cpp
    {
        double AU = physical_constants::ASTRONOMICAL_UNIT;
        double gravitationalParameter = 398600.4418e9;  // Earth GM

        Eigen::Vector3d departurePosition(0.02 * AU, 0.0, 0.0);
        Eigen::Vector3d arrivalPosition(0.0, -0.03 * AU, 0.0);
        double timeOfFlight = 100.0 * physical_constants::JULIAN_DAY;

        LambertTargeterIzzo lambertTargeter(departurePosition, arrivalPosition,
                                            timeOfFlight, gravitationalParameter);

        Eigen::Vector3d departureVelocity = lambertTargeter.getInertialVelocityAtDeparture();
        Eigen::Vector3d arrivalVelocity = lambertTargeter.getInertialVelocityAtArrival();

        // Expected radial and transverse velocities from reference
        double expectedRadialDep = -745.457;
        double expectedTransverseDep = 156.743;

        // For hyperbolic orbit, velocities should have these magnitudes
        checkClose("Lambert hyperbolic departure Vx", departureVelocity(0),
                  expectedRadialDep, std::abs(expectedRadialDep) * 1e-3);
        checkClose("Lambert hyperbolic departure Vy", departureVelocity(1),
                  expectedTransverseDep, std::abs(expectedTransverseDep) * 1e-3);
    }
}

void testGravityAssistRoutines()
{
    std::cout << "\n=== Gravity Assist Delta-V Computation ===" << std::endl;

    using namespace mission_segments;
    using namespace unit_conversions;
    using mathematical_constants::PI;

    // Test 1: Bending angle Delta-V (powered swing-by around Mars)
    // Reference: Hand calculations, Musegaas (2012)
    {
        const double expectedDeltaV = 3.652e3;
        const double tolerance = 0.0002;

        // Mars parameters
        const double marsGravitationalParameter = 4.2828018915e13;
        const double gravitationalParameterSun = 1.32712440018e20;
        const double distanceMarsToSun = convertAstronomicalUnitsToMeters(1.5);
        const double marsSmallestPeriapsisDistance = 3656248.0;

        // Mars heliocentric velocity (circular orbit)
        const Eigen::Vector3d marsVelocity(0.0, std::sqrt(gravitationalParameterSun / distanceMarsToSun), 0.0);

        // Satellite velocities
        const Eigen::Vector3d incomingVelocity(-25.0e3 * std::sin(PI / 6.0), 25.0e3 * std::cos(PI / 6.0), 0.0);
        const Eigen::Vector3d outgoingVelocity(incomingVelocity(0), 2.0 * marsVelocity(1) - incomingVelocity(1), 0.0);

        const double deltaV = calculateGravityAssistDeltaV(
            marsGravitationalParameter, marsVelocity, incomingVelocity, outgoingVelocity, marsSmallestPeriapsisDistance);

        double relError = std::abs(deltaV - expectedDeltaV) / expectedDeltaV;
        checkTrue("Mars gravity assist bending deltaV (rel err)", relError < tolerance);
    }

    // Test 2: Bending angle and velocity effect combined (Venus swing-by)
    // Reference: Musegaas (2012), verified in Excel
    {
        const double expectedDeltaV = 183.8481861944;
        const double tolerance = 1e-12;

        const double venusGravitationalParameter = 3.24860e14;
        const double venusSmallestPeriapsisDistance = 6351800.0;

        const Eigen::Vector3d venusVelocity(35000.0, 0.0, 0.0);
        const Eigen::Vector3d incomingVelocity(36000.0, 0.0, 0.0);
        const Eigen::Vector3d outgoingVelocity(34500.0, 0.0, 0.0);

        const double deltaV = calculateGravityAssistDeltaV(
            venusGravitationalParameter, venusVelocity, incomingVelocity, outgoingVelocity, venusSmallestPeriapsisDistance);

        double relError = std::abs(deltaV - expectedDeltaV) / expectedDeltaV;
        checkTrue("Venus gravity assist combined deltaV (rel err < 1e-12)", relError < tolerance);
    }

    // Test 3: No assist required
    {
        const double expectedDeltaV = 0.0;

        const double venusGravitationalParameter = 3.24860e14;
        const double venusSmallestPeriapsisDistance = 6351800.0;

        const Eigen::Vector3d venusVelocity(35000.0, 0.0, 0.0);
        const Eigen::Vector3d incomingVelocity(36000.0, 0.0, 0.0);
        const Eigen::Vector3d outgoingVelocity(35000.0, 1000.0, 0.0);

        const double deltaV = calculateGravityAssistDeltaV(
            venusGravitationalParameter, venusVelocity, incomingVelocity, outgoingVelocity, venusSmallestPeriapsisDistance);

        checkClose("No gravity assist required", deltaV, expectedDeltaV, 1e-20);
    }

    // Test 4: Velocity effect Delta-V (Cassini-1 trajectory from GTOP)
    // Reference: ESA GTOP code, 15-digit precision
    {
        const double expectedDeltaV = 1090.64622870007;
        const double tolerance = 1.0e-13;

        const double venusGravitationalParameter = 3.24860e14;
        const double venusSmallestPeriapsisDistance = 6351800.0;

        const Eigen::Vector3d venusVelocity(32851.224953746, -11618.7310059974, -2055.04615890989);
        const Eigen::Vector3d incomingVelocity(34216.4827530912, -15170.1440677825, 395.792122152361);
        const Eigen::Vector3d outgoingVelocity(37954.2431376052, -14093.0467234774, -5753.53728279429);

        // Test with eccentricity iteration scheme
        const double deltaVEcc = calculateGravityAssistDeltaV(
            venusGravitationalParameter, venusVelocity, incomingVelocity, outgoingVelocity,
            venusSmallestPeriapsisDistance, true);

        double relErrorEcc = std::abs(deltaVEcc - expectedDeltaV) / expectedDeltaV;
        checkTrue("Cassini-1 Venus gravity assist (eccentricity iteration)", relErrorEcc < tolerance);

        // Test with pericenter iteration scheme
        const double deltaVPeri = calculateGravityAssistDeltaV(
            venusGravitationalParameter, venusVelocity, incomingVelocity, outgoingVelocity,
            venusSmallestPeriapsisDistance, false);

        double relErrorPeri = std::abs(deltaVPeri - expectedDeltaV) / expectedDeltaV;
        checkTrue("Cassini-1 Venus gravity assist (pericenter iteration)", relErrorPeri < tolerance);
    }
}

void testUnpoweredGravityAssistPropagation()
{
    std::cout << "\n=== Unpowered Gravity Assist Propagation ===" << std::endl;

    using namespace mission_segments;

    // Reference: GTOP code, Messenger trajectory first swing-by
    // 15-digit precision benchmark
    const double tolerance = 1.0e-13;

    const Eigen::Vector3d expectedOutgoingVelocity(12868.5248737923, -22821.444560174, -775.698475033994);

    const double earthGravitationalParameter = 3.9860119e14;
    const Eigen::Vector3d earthVelocity(15025.522196446, -25544.3782752036, 0.0);
    const Eigen::Vector3d incomingVelocity(17969.3166254716, -23543.691593914, 6.38384671663496);

    const double rotationAngle = 1.35077257078;
    const double pericenterRadius = 1.80629232251 * 6378000.0;

    const Eigen::Vector3d outgoingVelocity = calculateUnpoweredGravityAssistOutgoingVelocity(
        earthGravitationalParameter, earthVelocity, incomingVelocity, rotationAngle, pericenterRadius);

    for (int i = 0; i < 3; i++) {
        double relError = std::abs(outgoingVelocity(i) - expectedOutgoingVelocity(i)) / std::abs(expectedOutgoingVelocity(i));
        checkTrue("Unpowered gravity assist outgoing V" + std::to_string(i), relError < tolerance);
    }
}

void testPoweredGravityAssistPropagation()
{
    std::cout << "\n=== Powered Gravity Assist Propagation ===" << std::endl;

    using namespace mission_segments;

    // Test 1: Powered GA function for unpowered case (deltaV = 0)
    // Reference: GTOP code, Messenger trajectory
    {
        const double tolerance = 1.0e-13;
        const Eigen::Vector3d expectedOutgoingVelocity(12868.5248737923, -22821.444560174, -775.698475033994);

        const double earthGravitationalParameter = 3.9860119e14;
        const Eigen::Vector3d earthVelocity(15025.522196446, -25544.3782752036, 0.0);
        const Eigen::Vector3d incomingVelocity(17969.3166254716, -23543.691593914, 6.38384671663496);

        const double rotationAngle = 1.35077257078;
        const double pericenterRadius = 1.80629232251 * 6378000.0;
        const double deltaV = 0.0;

        const Eigen::Vector3d outgoingVelocity = calculatePoweredGravityAssistOutgoingVelocity(
            earthGravitationalParameter, earthVelocity, incomingVelocity, rotationAngle, pericenterRadius, deltaV);

        for (int i = 0; i < 3; i++) {
            double relError = std::abs(outgoingVelocity(i) - expectedOutgoingVelocity(i)) / std::abs(expectedOutgoingVelocity(i));
            checkTrue("Powered GA (dV=0) outgoing V" + std::to_string(i), relError < tolerance);
        }
    }

    // Test 2: Powered GA with non-zero deltaV (Cassini-1 reverse engineered)
    // Reference: GTOP code, 15-digit precision
    {
        const double tolerance = 1.0e-14;
        const Eigen::Vector3d expectedOutgoingVelocity(37954.2431376052, -14093.0467234774, -5753.53728279429);

        const double venusGravitationalParameter = 3.24860e14;
        const Eigen::Vector3d venusVelocity(32851.224953746, -11618.7310059974, -2055.04615890989);
        const Eigen::Vector3d incomingVelocity(34216.4827530912, -15170.1440677825, 395.792122152361);

        const double rotationAngle = -2.0291949514117;
        const double pericenterRadius = 6351801.04541467;
        const double deltaV = 1090.64622870007;

        const Eigen::Vector3d outgoingVelocity = calculatePoweredGravityAssistOutgoingVelocity(
            venusGravitationalParameter, venusVelocity, incomingVelocity, rotationAngle, pericenterRadius, deltaV);

        for (int i = 0; i < 3; i++) {
            double relError = std::abs(outgoingVelocity(i) - expectedOutgoingVelocity(i)) / std::abs(expectedOutgoingVelocity(i));
            checkTrue("Cassini-1 powered GA outgoing V" + std::to_string(i), relError < tolerance);
        }
    }
}

void testEscapeAndCapture()
{
    std::cout << "\n=== Escape and Capture Delta-V ===" << std::endl;

    using namespace mission_segments;

    // Reference: Wakker (2007), Lecture Notes astro II, Table 18.2

    // Test 1: Delta-V for escape from Earth parking orbit to Mars transfer
    {
        const double tolerance = 1.0e-15;
        const double expectedDeltaVEscape = 3614.64460281887;

        const double gravitationalParameterEarth = 3.986e14;
        const double semiMajorAxis = 6.378e6 + 1.85e5;  // 185 km altitude
        const double eccentricity = 0.0;
        const double excessVelocity = 2944.61246668719;

        const double deltaVEscape = computeEscapeOrCaptureDeltaV(
            gravitationalParameterEarth, semiMajorAxis, eccentricity, excessVelocity);

        double relError = std::abs(deltaVEscape - expectedDeltaVEscape) / expectedDeltaVEscape;
        checkTrue("Earth escape deltaV (rel err < 1e-15)", relError < tolerance);
    }

    // Test 2: Delta-V for capture at Mars from Earth transfer
    {
        const double tolerance = 1.0e-15;
        const double expectedDeltaVCapture = 2087.1062716740798;

        const double gravitationalParameterMars = 4.2830e13;
        const double semiMajorAxis = 1.1 * 3.3895e6;  // 1.1 Mars radii
        const double eccentricity = 0.0;
        const double excessVelocity = 2648.83359973278;

        const double deltaVCapture = computeEscapeOrCaptureDeltaV(
            gravitationalParameterMars, semiMajorAxis, eccentricity, excessVelocity);

        double relError = std::abs(deltaVCapture - expectedDeltaVCapture) / expectedDeltaVCapture;
        checkTrue("Mars capture deltaV (rel err < 1e-15)", relError < tolerance);
    }

    // Test 3: Delta-V for escape from eccentric parking orbit
    {
        const double tolerance = 1.0e-15;
        const double expectedDeltaVEscape = 3200.2178506657729;

        const double gravitationalParameterEarth = 3.986e14;
        const double semiMajorAxis = 6.378e6 + 1.0e6;  // 200 km perigee, 1800 km apogee
        const double eccentricity = 0.108430468961778;
        const double excessVelocity = 2944.61246668719;

        const double deltaVEscape = computeEscapeOrCaptureDeltaV(
            gravitationalParameterEarth, semiMajorAxis, eccentricity, excessVelocity);

        double relError = std::abs(deltaVEscape - expectedDeltaVEscape) / expectedDeltaVEscape;
        checkTrue("Elliptical orbit escape deltaV (rel err < 1e-15)", relError < tolerance);
    }
}
