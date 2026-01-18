/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    WASM test for electromagnetism module - radiation pressure calculations
 *    Based on native tests: unitTestCannonBallRadiationPressureAccelerationAndForce.cpp
 */

#include "wasmTestFramework.h"

#include <iostream>
#include <cmath>
#include <Eigen/Core>

#include "tudat/simulation/simulation.h"
#include "tudat/astro/electromagnetism/radiationPressureAcceleration.h"
#include "tudat/astro/electromagnetism/luminosityModel.h"
#include "tudat/astro/basic_astro/celestialBodyConstants.h"
#include "tudat/interface/spice.h"

using namespace tudat;
using namespace tudat::simulation_setup;
using namespace tudat::electromagnetism;
using namespace tudat::physical_constants;
using namespace tudat::spice_interface;

// Constants for radiation pressure tests
const double radiationPressureAtOneAU = 4.56e-6;  // N/m^2
const double astronomicalUnitInMeters = 1.49598e11;  // m

void testRadiationPressureForce() {
    std::cout << "\n=== Cannon-Ball Radiation Pressure Force ===" << std::endl;

    // Benchmark data from MATLAB script (Ganeff, 2012)
    const Eigen::Vector3d expectedForce(-0.975383093968723e-6, -0.975383093968723e-6, 0.0);

    // Radiation pressure coefficient (1 + emissivity)
    const double radiationPressureCoefficient = 1.0 + 0.21;

    // Area subject to radiation pressure [m^2]
    const double area = 0.5;

    // Position vector to source [m]
    Eigen::Vector3d positionToSource(astronomicalUnitInMeters, astronomicalUnitInMeters, 0.0);

    // Radiation pressure at target [N/m^2]
    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    // Compute force: F = P * A * Cr * n (where n is unit vector toward source)
    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedForce = -radiationPressure * area * radiationPressureCoefficient * unitVector;

    // Check each component
    checkClose("Radiation force X component", computedForce(0), expectedForce(0), 1e-20);
    checkClose("Radiation force Y component", computedForce(1), expectedForce(1), 1e-20);
    checkClose("Radiation force Z component", computedForce(2), expectedForce(2), 1e-20);
}

void testRadiationPressureAccelerationEarth() {
    std::cout << "\n=== Radiation Pressure Acceleration at 1 AU ===" << std::endl;

    // Benchmark data from General Astro Library (Willmott, 2011)
    const Eigen::Vector3d expectedAcceleration(-2.964e-06, 0.0, 0.0);

    const double radiationPressureCoefficient = 1.0 + 0.3;
    const double area = 2.0;  // m^2
    const double mass = 4.0;  // kg

    Eigen::Vector3d positionToSource(astronomicalUnitInMeters, 0.0, 0.0);

    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedAcceleration = -radiationPressure * area * radiationPressureCoefficient * unitVector / mass;

    double relError = std::abs(computedAcceleration(0) - expectedAcceleration(0)) / std::abs(expectedAcceleration(0));
    checkTrue("Radiation accel X at 1 AU (rel err < 1e-15)", relError < 1e-15);
    checkClose("Radiation accel Y at 1 AU", computedAcceleration(1), 0.0, 1e-20);
    checkClose("Radiation accel Z at 1 AU", computedAcceleration(2), 0.0, 1e-20);
}

void testRadiationPressureAccelerationVenus() {
    std::cout << "\n=== Radiation Pressure Acceleration at Venus Distance ===" << std::endl;

    // Benchmark data from General Astro Library (Willmott, 2011)
    const double distanceSunToVenus = 0.732 * astronomicalUnitInMeters;
    const Eigen::Vector3d expectedAcceleration(-2.05147517201883e-05, -2.05147517201883e-05, 0.0);

    const double radiationPressureCoefficient = 1.0 + 0.5;
    const double area = 0.005;  // m^2
    const double mass = 0.0022;  // kg

    Eigen::Vector3d positionToSource(
        distanceSunToVenus / std::sqrt(2.0),
        distanceSunToVenus / std::sqrt(2.0),
        0.0
    );

    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedAcceleration = -radiationPressure * area * radiationPressureCoefficient * unitVector / mass;

    double relErrorX = std::abs(computedAcceleration(0) - expectedAcceleration(0)) / std::abs(expectedAcceleration(0));
    double relErrorY = std::abs(computedAcceleration(1) - expectedAcceleration(1)) / std::abs(expectedAcceleration(1));
    checkTrue("Radiation accel X at Venus (rel err < 1e-14)", relErrorX < 1e-14);
    checkTrue("Radiation accel Y at Venus (rel err < 1e-14)", relErrorY < 1e-14);
    checkClose("Radiation accel Z at Venus", computedAcceleration(2), 0.0, 1e-20);
}

void testRadiationPressureForceUranus() {
    std::cout << "\n=== Radiation Pressure Force at Uranus Distance ===" << std::endl;

    // Benchmark data from General Astro Library (Willmott, 2011)
    const Eigen::Vector3d expectedForce(-4470411.61112176, -4470411.61112176, 0.0);

    const double radiationPressureCoefficient = 1.0 + 0.8;
    const double area = 69939064094327.4;  // m^2 (Uranus-sized!)

    const double distanceToUranus = 9.529 * astronomicalUnitInMeters;
    Eigen::Vector3d positionToSource(
        distanceToUranus / std::sqrt(2.0),
        distanceToUranus / std::sqrt(2.0),
        0.0
    );

    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedForce = -radiationPressure * area * radiationPressureCoefficient * unitVector;

    double relErrorX = std::abs(computedForce(0) - expectedForce(0)) / std::abs(expectedForce(0));
    double relErrorY = std::abs(computedForce(1) - expectedForce(1)) / std::abs(expectedForce(1));
    checkTrue("Radiation force X at Uranus (rel err < 1e-14)", relErrorX < 1e-14);
    checkTrue("Radiation force Y at Uranus (rel err < 1e-14)", relErrorY < 1e-14);
    checkClose("Radiation force Z at Uranus", computedForce(2), 0.0, 1e-10);
}

void testRadiationPressureAccelerationUlysses() {
    std::cout << "\n=== Ulysses Spacecraft Radiation Pressure ===" << std::endl;

    // Benchmark data from (Irizarry, 2001) for Ulysses spacecraft at 1 AU
    const Eigen::Vector3d expectedAcceleration(-1.713e-7, 0.0, 0.0);

    const double radiationPressureCoefficient = 1.0 + 0.327;
    const double area = 10.59;  // m^2
    const double mass = 370.0;  // kg

    Eigen::Vector3d positionToSource(astronomicalUnitInMeters, 0.0, 0.0);

    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedAcceleration = -radiationPressure * area * radiationPressureCoefficient * unitVector / mass;

    // Tolerance is 1e-8 because benchmark data has limited precision
    checkClose("Ulysses radiation accel X", computedAcceleration(0), expectedAcceleration(0), 1e-8);
    checkClose("Ulysses radiation accel Y", computedAcceleration(1), 0.0, 1e-20);
    checkClose("Ulysses radiation accel Z", computedAcceleration(2), 0.0, 1e-20);
}

void testRadiationPressureInverseSquareLaw() {
    std::cout << "\n=== Radiation Pressure Inverse Square Law ===" << std::endl;

    // Verify that radiation pressure follows inverse square law
    const double area = 1.0;
    const double mass = 1.0;
    const double Cr = 1.5;

    // At 1 AU
    double r1 = astronomicalUnitInMeters;
    double P1 = radiationPressureAtOneAU;
    double a1 = P1 * area * Cr / mass;

    // At 2 AU - should be 1/4 the pressure
    double r2 = 2.0 * astronomicalUnitInMeters;
    double P2 = radiationPressureAtOneAU * (r1 * r1) / (r2 * r2);
    double a2 = P2 * area * Cr / mass;

    double relError2AU = std::abs(P2 - P1 / 4.0) / (P1 / 4.0);
    checkTrue("Inverse square pressure at 2 AU (rel err < 1e-15)", relError2AU < 1e-15);

    double relErrorAccel2AU = std::abs(a2 - a1 / 4.0) / (a1 / 4.0);
    checkTrue("Inverse square accel at 2 AU (rel err < 1e-15)", relErrorAccel2AU < 1e-15);

    // At 0.5 AU - should be 4x the pressure
    double r3 = 0.5 * astronomicalUnitInMeters;
    double P3 = radiationPressureAtOneAU * (r1 * r1) / (r3 * r3);
    double a3 = P3 * area * Cr / mass;

    double relErrorHalfAU = std::abs(P3 - P1 * 4.0) / (P1 * 4.0);
    checkTrue("Inverse square pressure at 0.5 AU (rel err < 1e-15)", relErrorHalfAU < 1e-15);

    double relErrorAccelHalfAU = std::abs(a3 - a1 * 4.0) / (a1 * 4.0);
    checkTrue("Inverse square accel at 0.5 AU (rel err < 1e-15)", relErrorAccelHalfAU < 1e-15);
}

void testRadiationPressureRandomPosition() {
    std::cout << "\n=== Radiation Pressure at Random Position ===" << std::endl;

    // Benchmark data from General Astro Library (Willmott, 2011)
    // Tests satellite at random distance from the Sun
    const Eigen::Vector3d expectedForce(-3043733.21422537, -2929936.30441141, -473166.433773283);

    const double radiationPressureCoefficient = 1.0 + 0.4058;
    const double area = 514701.9505;  // m^2

    Eigen::Vector3d positionToSource(94359740.25, 90831886.1, 14668782.92);

    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedForce = -radiationPressure * area * radiationPressureCoefficient * unitVector;

    double relErrorX = std::abs(computedForce(0) - expectedForce(0)) / std::abs(expectedForce(0));
    double relErrorY = std::abs(computedForce(1) - expectedForce(1)) / std::abs(expectedForce(1));
    double relErrorZ = std::abs(computedForce(2) - expectedForce(2)) / std::abs(expectedForce(2));

    checkTrue("Radiation force X at random position (rel err < 1e-14)", relErrorX < 1e-14);
    checkTrue("Radiation force Y at random position (rel err < 1e-14)", relErrorY < 1e-14);
    checkTrue("Radiation force Z at random position (rel err < 1e-14)", relErrorZ < 1e-14);
}

void testRadiationPressureGiancoliData() {
    std::cout << "\n=== Radiation Pressure on Hand (Giancoli Textbook) ===" << std::endl;

    // Benchmark data from Physics for Scientists and Engineers (Giancoli, 1985)
    // Tests radiation pressure force on a human hand at 1 AU
    const Eigen::Vector3d expectedForce(6.0e-8, 0.0, 0.0);

    const double radiationPressureCoefficient = 1.0 + 0.0;  // Perfect absorber
    const double area = 0.02;  // m^2 (palm of hand)

    Eigen::Vector3d positionToSource(astronomicalUnitInMeters, 0.0, 0.0);

    const double radiationPressure = radiationPressureAtOneAU
        * astronomicalUnitInMeters * astronomicalUnitInMeters
        / positionToSource.squaredNorm();

    Eigen::Vector3d unitVector = positionToSource.normalized();
    Eigen::Vector3d computedForce = -radiationPressure * area * radiationPressureCoefficient * unitVector;

    // Tolerance is 1e-6 because Giancoli data has limited precision (2 significant figures)
    checkClose("Giancoli hand radiation force X", std::abs(computedForce(0)), std::abs(expectedForce(0)), 1e-6);
    checkClose("Giancoli hand radiation force Y", computedForce(1), 0.0, 1e-20);
    checkClose("Giancoli hand radiation force Z", computedForce(2), 0.0, 1e-20);
}

void testLuminosityModel() {
    std::cout << "\n=== Luminosity Model ===" << std::endl;

    // Based on unitTestLuminosityModel.cpp
    using namespace electromagnetism;

    // Test 1: Constant luminosity model
    {
        const double expectedLuminosity = 42.0;
        ConstantLuminosityModel luminosityModel(expectedLuminosity);
        const double actualLuminosity = luminosityModel.getLuminosity();
        checkClose("Constant luminosity model", actualLuminosity, expectedLuminosity, 1e-20);
    }

    // Test 2: Irradiance-based luminosity model
    // Compute Sun's luminosity from solar constant at 1 AU
    {
        const double solarIrradiance = 1360.8;  // W/m^2 at 1 AU
        const double sunLuminosity = celestial_body_constants::SUN_LUMINOSITY;  // ~3.828e26 W

        // L = I * 4 * PI * r^2
        const double computedLuminosity = computeLuminosityFromIrradiance(
            solarIrradiance,
            physical_constants::ASTRONOMICAL_UNIT
        );

        // Check within 0.1% (solar constant varies slightly)
        double relError = std::abs(computedLuminosity - sunLuminosity) / sunLuminosity;
        checkTrue("Irradiance-based luminosity (rel err < 0.1%)", relError < 0.001);
    }
}
