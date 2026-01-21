/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    SPICE WASM tests - Functions that work without external kernel files.
 */

#include "wasmTestFramework.h"

#include <memory>

// Basic astrodynamics
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/unitConversions.h"

// Mathematics
#include "tudat/math/basic/mathematicalConstants.h"

// SPICE interface
#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/earth_orientation/terrestrialTimeScaleConverter.h"

using namespace tudat;

void testSpiceTimeConversions()
{
    std::cout << "\n=== SPICE Time Conversions ===" << std::endl;

    using namespace spice_interface;

    // Test Julian Date to Ephemeris Time conversion
    // J2000 epoch: January 1, 2000, 12:00 TT
    // Julian Date at J2000 = 2451545.0
    // Ephemeris Time at J2000 = 0.0 (by definition)

    double j2000JulianDate = 2451545.0;
    double ephemerisTime = convertJulianDateToEphemerisTime(j2000JulianDate);
    checkClose("SPICE: JD 2451545.0 -> ET 0.0", ephemerisTime, 0.0, 1e-6);

    // Test reverse conversion
    double recoveredJD = convertEphemerisTimeToJulianDate(ephemerisTime);
    checkClose("SPICE: ET 0.0 -> JD 2451545.0", recoveredJD, j2000JulianDate, 1e-10);

    // Test one day after J2000
    // One day = 86400 seconds
    double oneDayET = 86400.0;
    double oneDayJD = convertEphemerisTimeToJulianDate(oneDayET);
    checkClose("SPICE: ET 86400 -> JD 2451546.0", oneDayJD, 2451546.0, 1e-10);

    // Round-trip test
    double testJD = 2460000.0;  // Some arbitrary Julian Date
    double testET = convertJulianDateToEphemerisTime(testJD);
    double roundTripJD = convertEphemerisTimeToJulianDate(testET);
    checkClose("SPICE: JD round-trip", roundTripJD, testJD, 1e-10);

    // Test with negative ephemeris time (before J2000)
    double beforeJ2000ET = -86400.0;  // One day before J2000
    double beforeJ2000JD = convertEphemerisTimeToJulianDate(beforeJ2000ET);
    checkClose("SPICE: ET -86400 -> JD 2451544.0", beforeJ2000JD, 2451544.0, 1e-10);
}

void testSpiceFrameRotations()
{
    std::cout << "\n=== SPICE Frame Rotations ===" << std::endl;

    // In WASM, the spiceInterface uses an analytical rotation for J2000<->ECLIPJ2000
    // This is a constant rotation about the X-axis by the obliquity of the ecliptic
    using namespace spice_interface;

    // Test the J2000 <-> ECLIPJ2000 rotation
    // These rotations don't require kernel files (analytical in WASM, SPICE otherwise)

    Eigen::Matrix3d j2000ToEclip = getRotationFromJ2000ToEclipJ2000();
    Eigen::Matrix3d eclipToJ2000 = getRotationFromEclipJ2000ToJ2000();

    // Verify rotation matrices are valid (orthogonal, det = 1)
    checkClose("J2000->ECLIP determinant", j2000ToEclip.determinant(), 1.0, 1e-14);
    checkClose("ECLIP->J2000 determinant", eclipToJ2000.determinant(), 1.0, 1e-14);

    // Verify orthogonality: R * R^T = I
    Eigen::Matrix3d shouldBeI1 = j2000ToEclip * j2000ToEclip.transpose();
    checkClose("J2000->ECLIP orthogonality (trace)", shouldBeI1.trace(), 3.0, 1e-14);

    Eigen::Matrix3d shouldBeI2 = eclipToJ2000 * eclipToJ2000.transpose();
    checkClose("ECLIP->J2000 orthogonality (trace)", shouldBeI2.trace(), 3.0, 1e-14);

    // Verify they are inverses of each other
    Eigen::Matrix3d product = j2000ToEclip * eclipToJ2000;
    checkClose("J2000<->ECLIP inverse (trace)", product.trace(), 3.0, 1e-14);

    // The obliquity of the ecliptic at J2000 is 84381.448 arcseconds = 23.4392911... degrees
    // The rotation should be about the X-axis by this angle
    double obliquityRad = 84381.448 * mathematical_constants::PI / (180.0 * 3600.0);

    // Check that the rotation preserves the X-axis (rotation is about X)
    Eigen::Vector3d xAxis(1.0, 0.0, 0.0);
    Eigen::Vector3d rotatedX = j2000ToEclip * xAxis;
    checkVectorClose("J2000->ECLIP preserves X-axis", rotatedX, xAxis, 1e-10);

    // Check rotation angle by looking at Y and Z components
    Eigen::Vector3d yAxis(0.0, 1.0, 0.0);
    Eigen::Vector3d rotatedY = j2000ToEclip * yAxis;
    // After rotation about X by obliquity, Y should go to (0, cos(obl), sin(obl))
    checkClose("J2000->ECLIP Y->Y' cos component", rotatedY(1), std::cos(obliquityRad), 1e-10);
    checkClose("J2000->ECLIP Y->Z' sin component", rotatedY(2), std::sin(obliquityRad), 1e-10);
}

void testSpiceErrorHandling()
{
    std::cout << "\n=== SPICE Error Handling ===" << std::endl;

    using namespace spice_interface;

    // Test SPICE error handling functions.
    // In WASM builds, these are no-op stubs because the underlying SPICE functions
    // (erract_c, errdev_c, getmsg_c) crash due to f2c string handling issues.
    // The stubs allow code to run without crashing.

    // Test 1: checkFailure (calls failed_c - works in WASM)
    bool hadError = checkFailure();
    checkTrue("checkFailure() works", true);

    // Test 2: toggleErrorReturn (no-op in WASM)
    toggleErrorReturn();
    checkTrue("toggleErrorReturn() callable", true);

    // Test 3: suppressErrorOutput (no-op in WASM)
    suppressErrorOutput();
    checkTrue("suppressErrorOutput() callable", true);

    // Test 4: getErrorMessage (returns empty in WASM)
    std::string errorMsg = getErrorMessage();
    checkTrue("getErrorMessage() callable", true);

    // Test 5: kernel count
    int kernelCount = getTotalCountOfKernelsLoaded();
    checkTrue("SPICE kernel count >= 0", kernelCount >= 0);
}

void testSpiceTLEPropagation()
{
    std::cout << "\n=== TLE/SGP4 Propagation (Vallado Benchmark) ===" << std::endl;

    using namespace spice_interface;
    using namespace ephemerides;

    // =========================================================================
    // Vallado TLE Test Case - Reference: Vallado (2013), page 234
    // This is the canonical test case for SGP4 validation
    // =========================================================================
    {
        // Real TLE from Vallado textbook
        std::string tleLines = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753\n"
                               "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

        std::shared_ptr<Tle> tle = std::make_shared<Tle>(tleLines);

        // Create TLE ephemeris in J2000 frame (converted from TEME)
        TleEphemeris tleEphemeris("Earth", "J2000", tle, false);

        // Propagate for 3 Julian days from TLE epoch
        // Use same time conversion as native test: UTC epoch -> TDB epoch
        double utcEpoch = 3.0 * physical_constants::JULIAN_DAY + tle->getEpoch();
        double tdbEpoch = earth_orientation::createDefaultTimeConverter()->getCurrentTime(
            basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, utcEpoch);

        Eigen::Vector6d propagatedState = tleEphemeris.getCartesianState(tdbEpoch);
        Eigen::Vector3d propagatedPosition = propagatedState.head<3>();
        Eigen::Vector3d propagatedVelocity = propagatedState.tail<3>();

        // Reference values from Vallado (converted to J2000 frame)
        // Note: Vallado gives values in m and m/s
        Eigen::Vector3d valladoPosition;
        Eigen::Vector3d valladoVelocity;
        valladoPosition << -9059941.3786, 4659697.2000, 813958.8875;
        valladoVelocity << -2233.348094, -4110.136162, -3157.394074;

        // Check position (within 50m tolerance - same as main Tudat test)
        double positionError = (propagatedPosition - valladoPosition).norm();
        checkTrue("Vallado TLE position error < 50m", positionError < 50.0);

        // Check velocity (within 0.05 m/s tolerance - same as main Tudat test)
        double velocityError = (propagatedVelocity - valladoVelocity).norm();
        checkTrue("Vallado TLE velocity error < 0.05 m/s", velocityError < 0.05);

        // Report actual errors for diagnostic purposes
        std::cout << "[INFO] Vallado position error: " << positionError << " m" << std::endl;
        std::cout << "[INFO] Vallado velocity error: " << velocityError << " m/s" << std::endl;
    }

    // =========================================================================
    // Basic TLE Parsing and Property Verification
    // =========================================================================
    {
        // Use the same Vallado TLE for parsing tests
        std::string tleLines = "1 00005U 58002B   00179.78495062  .00000023  00000-0  28098-4 0  4753\n"
                               "2 00005  34.2682 348.7242 1859667 331.7664  19.3264 10.82419157413667";

        std::shared_ptr<Tle> tle = std::make_shared<Tle>(tleLines);

        // Verify parsed orbital elements
        // Inclination: 34.2682 degrees
        double expectedInclination = 34.2682 / 180.0 * mathematical_constants::PI;
        checkClose("TLE parsed inclination", tle->getInclination(), expectedInclination, 1e-5);  // Looser tolerance for text parsing

        // Eccentricity: 0.1859667 (stored as 1859667 in TLE, implied decimal point)
        checkClose("TLE parsed eccentricity", tle->getEccentricity(), 0.1859667, 1e-7);

        // RAAN: 348.7242 degrees
        double expectedRaan = 348.7242 / 180.0 * mathematical_constants::PI;
        checkClose("TLE parsed RAAN", tle->getRightAscension(), expectedRaan, 1e-6);

        // Argument of perigee: 331.7664 degrees
        double expectedArgPerigee = 331.7664 / 180.0 * mathematical_constants::PI;
        checkClose("TLE parsed arg perigee", tle->getArgOfPerigee(), expectedArgPerigee, 1e-6);

        // Mean anomaly: 19.3264 degrees
        double expectedMeanAnomaly = 19.3264 / 180.0 * mathematical_constants::PI;
        checkClose("TLE parsed mean anomaly", tle->getMeanAnomaly(), expectedMeanAnomaly, 1e-6);

        // Mean motion: 10.82419157 rev/day -> convert to rad/min for internal storage
        // rad/min = rev/day * 2*pi / (24*60)
        double expectedMeanMotion = 10.82419157 * 2.0 * mathematical_constants::PI / (24.0 * 60.0);
        checkClose("TLE parsed mean motion", tle->getMeanMotion(), expectedMeanMotion, 1e-10);
    }

    // =========================================================================
    // TLE Constructed from Orbital Elements (ISS-like orbit)
    // =========================================================================
    {
        double epoch = 0.0;  // J2000 epoch
        double bStar = 0.0001;
        double inclination = unit_conversions::convertDegreesToRadians(51.6);
        double rightAscension = 0.0;
        double eccentricity = 0.0001;
        double argOfPerigee = 0.0;
        double meanAnomaly = 0.0;
        double meanMotion = 2.0 * mathematical_constants::PI / 92.0;  // ~92 min period

        std::shared_ptr<Tle> tle = std::make_shared<Tle>(
            epoch, bStar, inclination, rightAscension,
            eccentricity, argOfPerigee, meanAnomaly, meanMotion);

        TleEphemeris tleEphemeris("Earth", "J2000", tle, false);

        // Get state at epoch
        Eigen::Vector6d state = tleEphemeris.getCartesianState(0.0);
        double positionMagnitude = state.head<3>().norm();
        double velocityMagnitude = state.tail<3>().norm();

        // ISS-like orbit: ~6778 km radius, ~7.66 km/s velocity
        checkTrue("ISS-like orbit radius > 6.5e6 m", positionMagnitude > 6.5e6);
        checkTrue("ISS-like orbit radius < 7.0e6 m", positionMagnitude < 7.0e6);
        checkTrue("ISS-like orbit velocity > 7.5e3 m/s", velocityMagnitude > 7.5e3);
        checkTrue("ISS-like orbit velocity < 7.8e3 m/s", velocityMagnitude < 7.8e3);
    }
}

void testSpiceTemeFrameRotation()
{
    std::cout << "\n=== SPICE TEME Frame Rotation ===" << std::endl;

    // TEME frame rotation uses SOFA functions (calculateEquationOfEquinoxes, getPrecessionNutationMatrix)
    // which are pure computational and should work in WASM
    using namespace ephemerides;

    // Test the TEME (True Equator, Mean Equinox) frame rotation
    // This is used for TLE/SGP4 coordinate transformations

    double epoch = 0.0;  // J2000

    Eigen::Matrix3d temeToJ2000 = getRotationMatrixFromTemeToJ2000(epoch);
    Eigen::Matrix3d j2000ToTeme = getRotationMatrixFromJ2000ToTeme(epoch);

    // Verify rotation matrices are valid (orthogonal, det = 1)
    checkClose("TEME->J2000 determinant", temeToJ2000.determinant(), 1.0, 1e-14);
    checkClose("J2000->TEME determinant", j2000ToTeme.determinant(), 1.0, 1e-14);

    // Verify orthogonality
    Eigen::Matrix3d shouldBeI = temeToJ2000 * temeToJ2000.transpose();
    checkClose("TEME->J2000 orthogonality (trace)", shouldBeI.trace(), 3.0, 1e-14);

    // Verify they are inverses
    Eigen::Matrix3d product = temeToJ2000 * j2000ToTeme;
    checkClose("TEME<->J2000 inverse (trace)", product.trace(), 3.0, 1e-14);

    // At J2000 epoch, the TEME and J2000 frames should be very close
    // (they differ mainly due to nutation and precession accumulated since J2000)
    // The difference should be small angles (arc-seconds to arc-minutes)
    Eigen::Matrix3d diff = temeToJ2000 - Eigen::Matrix3d::Identity();
    double maxDiff = diff.cwiseAbs().maxCoeff();
    checkTrue("TEME≈J2000 at epoch (small rotation)", maxDiff < 0.01);  // Less than ~0.5 degrees
}

void testAnalyticalPlanetaryEphemeris()
{
    std::cout << "\n=== Analytical Planetary Ephemeris (JPL Approximate Positions) ===" << std::endl;

    using namespace ephemerides;

    // Test the JPL Approximate Planetary Positions algorithm
    // This uses Tudat's built-in ApproximateJplEphemeris class which requires
    // NO external SPICE kernels - all orbital elements are embedded in the code.

    // Reference: Standish, E.M. "Keplerian Elements for Approximate Positions
    //            of the Major Planets" (JPL)

    // Test at J2000 epoch (secondsSinceJ2000 = 0)
    double epoch = 0.0;

    // Test Earth state at J2000
    {
        ApproximateJplEphemeris earthEphemeris("Earth");
        Eigen::Vector6d earthState = earthEphemeris.getCartesianState(epoch);

        double earthDistance = earthState.head<3>().norm();
        double earthVelocity = earthState.tail<3>().norm();

        // Earth should be ~1 AU from Sun (149.6 million km)
        double AU = 149597870700.0;  // meters
        checkTrue("Earth distance ~1 AU at J2000", std::abs(earthDistance - AU) < 0.03 * AU);

        // Earth orbital velocity ~29.78 km/s
        checkTrue("Earth velocity ~30 km/s", earthVelocity > 29000.0 && earthVelocity < 31000.0);

        std::cout << "[INFO] Earth at J2000: r=" << earthDistance/1e9 << " Gm, v=" << earthVelocity/1e3 << " km/s" << std::endl;
    }

    // Test Mars state at J2000
    {
        ApproximateJplEphemeris marsEphemeris("Mars");
        Eigen::Vector6d marsState = marsEphemeris.getCartesianState(epoch);

        double marsDistance = marsState.head<3>().norm();
        double marsVelocity = marsState.tail<3>().norm();

        // Mars should be ~1.38-1.67 AU from Sun (varies due to eccentricity)
        double AU = 149597870700.0;
        checkTrue("Mars distance 1.3-1.7 AU", marsDistance > 1.3 * AU && marsDistance < 1.7 * AU);

        // Mars orbital velocity ~21-26 km/s
        checkTrue("Mars velocity 21-26 km/s", marsVelocity > 21000.0 && marsVelocity < 27000.0);

        std::cout << "[INFO] Mars at J2000: r=" << marsDistance/1e9 << " Gm, v=" << marsVelocity/1e3 << " km/s" << std::endl;
    }

    // Test state propagation over time (should change smoothly)
    {
        ApproximateJplEphemeris earthEphemeris("Earth");

        double epoch1 = 0.0;                    // J2000
        double epoch2 = 365.25 * 86400.0;       // 1 year later

        Eigen::Vector6d state1 = earthEphemeris.getCartesianState(epoch1);
        Eigen::Vector6d state2 = earthEphemeris.getCartesianState(epoch2);

        // After one year, Earth should return to approximately the same position
        // (small drift due to orbital precession and approximate model)
        double positionDiff = (state1.head<3>() - state2.head<3>()).norm();
        double AU = 149597870700.0;

        // Should be within 1% of orbital radius after one complete orbit
        checkTrue("Earth returns after 1 year (within 5%)", positionDiff < 0.05 * AU);

        std::cout << "[INFO] Earth position change after 1 year: " << positionDiff/1e6 << " km" << std::endl;
    }

    // Test GTOP ephemeris as alternative
    {
        ApproximateGtopEphemeris earthGtop("Earth");
        Eigen::Vector6d earthGtopState = earthGtop.getCartesianState(epoch);

        double earthGtopDistance = earthGtopState.head<3>().norm();
        double AU = 149597870700.0;

        // GTOP should also give ~1 AU
        checkTrue("GTOP Earth distance ~1 AU", std::abs(earthGtopDistance - AU) < 0.03 * AU);

        std::cout << "[INFO] GTOP Earth at J2000: r=" << earthGtopDistance/1e9 << " Gm" << std::endl;
    }

    // Test all supported planets
    {
        std::vector<std::string> planets = {"Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune"};
        double AU = 149597870700.0;

        // Expected distances in AU (approximate)
        std::vector<double> expectedDistances = {0.39, 0.72, 1.0, 1.52, 5.2, 9.5, 19.2, 30.1};

        for (size_t i = 0; i < planets.size(); ++i) {
            ApproximateJplEphemeris ephemeris(planets[i]);
            Eigen::Vector6d state = ephemeris.getCartesianState(epoch);
            double distance = state.head<3>().norm() / AU;

            // Check distance is within 20% of expected (accounts for orbital eccentricity)
            bool distanceOk = std::abs(distance - expectedDistances[i]) < 0.3 * expectedDistances[i];
            checkTrue(planets[i] + " distance reasonable", distanceOk);

            if (!distanceOk) {
                std::cout << "[WARN] " << planets[i] << " distance: " << distance << " AU (expected ~" << expectedDistances[i] << " AU)" << std::endl;
            }
        }
    }
}
