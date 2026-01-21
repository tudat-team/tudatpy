/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Emscripten bindings for WASM visualization - exposes Tudat functions to JavaScript.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <emscripten/emscripten.h>

#include <vector>
#include <map>
#include <memory>
#include <cmath>

// Basic astrodynamics
#include "tudat/astro/basic_astro/orbitalElementConversions.h"
#include "tudat/astro/basic_astro/keplerPropagator.h"
#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"

// Mathematics
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/integrators/rungeKuttaVariableStepSizeIntegrator.h"
#include "tudat/math/root_finders/newtonRaphson.h"

// Gravitation
#include "tudat/astro/gravitation/librationPoint.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

// Ephemerides (for TLE/SGP4 and planetary positions)
#include "tudat/astro/ephemerides/tleEphemeris.h"
#include "tudat/astro/ephemerides/approximatePlanetPositions.h"
#include "tudat/astro/ephemerides/tabulatedEphemeris.h"
#include "tudat/math/interpolators/createInterpolator.h"

// Mission segments (Lambert targeting)
#include "tudat/astro/mission_segments/zeroRevolutionLambertTargeterIzzo.h"

// SPICE interface
#include <tudat/interface/spice.h>

// Propagation
#include "tudat/simulation/propagation_setup/propagationCR3BPFullProblem.h"
#include "tudat/math/integrators/createNumericalIntegrator.h"

using namespace tudat;
using namespace emscripten;
namespace tsi = tudat::spice_interface;

// ============================================================================
// Helper: Convert JavaScript array to std::vector
// ============================================================================

template<typename T>
std::vector<T> vecFromJSArray(const val& jsArray) {
    const auto length = jsArray["length"].as<unsigned>();
    std::vector<T> result;
    result.reserve(length);
    for (unsigned i = 0; i < length; ++i) {
        result.push_back(jsArray[i].as<T>());
    }
    return result;
}

// ============================================================================
// Kepler Orbit Propagation (for analytical vs numerical comparison)
// ============================================================================

// Propagate Keplerian orbit and return ephemeris as array of [t, x, y, z]
// Uses Tudat's propagateKeplerOrbit for analytical solution
val propagateKeplerOrbit(
    double semiMajorAxis,       // km (converted to m internally)
    double eccentricity,
    double inclination,         // degrees
    double raan,                // degrees
    double argPeriapsis,        // degrees
    double trueAnomaly0,        // degrees
    double duration,            // seconds
    int numPoints)
{
    using namespace orbital_element_conversions;

    // Standard Earth gravitational parameter in m^3/s^2
    // GM_Earth = 3.986004418e14 m^3/s^2 (WGS84 value)
    double mu = 3.986004418e14;

    // Convert to SI units (meters, radians)
    double smaMeters = semiMajorAxis * 1000.0;
    double incRad = inclination * mathematical_constants::PI / 180.0;
    double raanRad = raan * mathematical_constants::PI / 180.0;
    double argPerRad = argPeriapsis * mathematical_constants::PI / 180.0;
    double ta0Rad = trueAnomaly0 * mathematical_constants::PI / 180.0;

    // Initial Keplerian elements (Tudat ordering: a, e, i, omega, RAAN, true anomaly)
    Eigen::Vector6d keplerElements;
    keplerElements << smaMeters,
                      eccentricity,
                      incRad,
                      argPerRad,
                      raanRad,
                      ta0Rad;

    std::vector<double> ephemeris;
    ephemeris.reserve(numPoints * 4);

    for (int i = 0; i < numPoints; i++) {
        double t = (static_cast<double>(i) / (numPoints - 1)) * duration;

        // Propagate Keplerian elements using Tudat's analytical solution
        Eigen::Vector6d propagatedKepler = propagateKeplerOrbit(keplerElements, t, mu);

        // Convert to Cartesian (returns state in meters, m/s)
        Eigen::Vector6d cartesianState = convertKeplerianToCartesianElements(propagatedKepler, mu);

        ephemeris.push_back(t);
        ephemeris.push_back(cartesianState(0));  // Already in meters
        ephemeris.push_back(cartesianState(1));
        ephemeris.push_back(cartesianState(2));
    }

    return val(typed_memory_view(ephemeris.size(), ephemeris.data())).call<val>("slice");
}

// ============================================================================
// Fixed-step RK4 Integration (shows visible error accumulation)
// ============================================================================

val propagateRK4Orbit(
    double semiMajorAxis,       // km (converted to m internally)
    double eccentricity,
    double inclination,         // degrees
    double raan,                // degrees
    double argPeriapsis,        // degrees
    double trueAnomaly0,        // degrees
    double duration,            // seconds
    int numPoints,
    double stepSize)            // Fixed step size in seconds
{
    using namespace orbital_element_conversions;

    double mu = 3.986004418e14;  // m^3/s^2
    double smaMeters = semiMajorAxis * 1000.0;

    // Keplerian elements
    Eigen::Vector6d keplerElements;
    keplerElements << smaMeters,
                      eccentricity,
                      inclination * mathematical_constants::PI / 180.0,
                      argPeriapsis * mathematical_constants::PI / 180.0,
                      raan * mathematical_constants::PI / 180.0,
                      trueAnomaly0 * mathematical_constants::PI / 180.0;

    Eigen::Vector6d state = convertKeplerianToCartesianElements(keplerElements, mu);

    // State derivative function
    auto derivative = [mu](const Eigen::Vector6d& s) -> Eigen::Vector6d {
        Eigen::Vector6d d;
        Eigen::Vector3d pos = s.head<3>();
        double r = pos.norm();
        d.head<3>() = s.tail<3>();
        d.tail<3>() = -mu / (r * r * r) * pos;
        return d;
    };

    std::vector<double> ephemeris;
    ephemeris.reserve(numPoints * 4);

    double t = 0.0;
    double dt = duration / (numPoints - 1);
    double h = stepSize;  // Fixed step size

    for (int i = 0; i < numPoints; i++) {
        double targetTime = i * dt;

        // Integrate with fixed RK4 steps until we reach target time
        while (t < targetTime - h/2) {
            // RK4 step
            Eigen::Vector6d k1 = derivative(state);
            Eigen::Vector6d k2 = derivative(state + h/2 * k1);
            Eigen::Vector6d k3 = derivative(state + h/2 * k2);
            Eigen::Vector6d k4 = derivative(state + h * k3);
            state = state + (h/6.0) * (k1 + 2*k2 + 2*k3 + k4);
            t += h;
        }

        ephemeris.push_back(targetTime);
        ephemeris.push_back(state(0));
        ephemeris.push_back(state(1));
        ephemeris.push_back(state(2));
    }

    return val(typed_memory_view(ephemeris.size(), ephemeris.data())).call<val>("slice");
}

// ============================================================================
// J2-Only vs Full Force Model Numerical Propagation Comparison
// Returns interleaved arrays: [t, j2_x, j2_y, j2_z, full_x, full_y, full_z, ...]
//
// This compares two numerical models with the SAME initial osculating state:
// - J2-only: Two-body + J2 perturbation (baseline high-fidelity model)
// - Full force: J2/J3/J4 + Sun/Moon third-body + Drag + SRP
//
// NOTE: We use the TLE/SGP4 only to get an initial state. We do NOT compare
// against SGP4 propagation because SGP4 uses Brouwer mean elements internally,
// which are fundamentally different from the osculating elements used in
// numerical propagation. Comparing them gives misleading divergence.
// ============================================================================

// Approximate Sun position in J2000 frame (circular orbit approximation)
// Returns position in meters
Eigen::Vector3d getApproxSunPosition(double julianDate) {
    // Sun orbital parameters (approximate circular orbit around barycenter)
    const double AU = 149597870700.0;  // meters
    const double sunMeanMotion = 2.0 * mathematical_constants::PI / 365.25 / 86400.0;  // rad/s

    // J2000 epoch: JD 2451545.0 = 2000-01-01 12:00:00 TT
    const double J2000_JD = 2451545.0;
    const double obliquity = 23.439291 * mathematical_constants::PI / 180.0;  // Earth's obliquity

    double daysSinceJ2000 = julianDate - J2000_JD;
    double secondsSinceJ2000 = daysSinceJ2000 * 86400.0;

    // Mean longitude of Sun (simplified)
    double meanLong = 280.46 * mathematical_constants::PI / 180.0 + sunMeanMotion * secondsSinceJ2000;

    // Sun position in ecliptic coordinates
    double xEclip = AU * std::cos(meanLong);
    double yEclip = AU * std::sin(meanLong);
    double zEclip = 0.0;

    // Rotate from ecliptic to J2000 equatorial
    double xJ2000 = xEclip;
    double yJ2000 = yEclip * std::cos(obliquity) - zEclip * std::sin(obliquity);
    double zJ2000 = yEclip * std::sin(obliquity) + zEclip * std::cos(obliquity);

    return Eigen::Vector3d(xJ2000, yJ2000, zJ2000);
}

// Approximate Moon position in J2000 frame (simplified model)
// Returns position in meters
Eigen::Vector3d getApproxMoonPosition(double julianDate) {
    // Moon orbital parameters (approximate)
    const double moonSMA = 384400000.0;  // meters (average distance)
    const double moonPeriod = 27.321661 * 86400.0;  // seconds (sidereal period)
    const double moonMeanMotion = 2.0 * mathematical_constants::PI / moonPeriod;
    const double moonInclination = 5.145 * mathematical_constants::PI / 180.0;  // to ecliptic

    const double J2000_JD = 2451545.0;
    const double obliquity = 23.439291 * mathematical_constants::PI / 180.0;

    double daysSinceJ2000 = julianDate - J2000_JD;
    double secondsSinceJ2000 = daysSinceJ2000 * 86400.0;

    // Simplified Moon longitude (ignoring perturbations)
    double meanLong = 218.32 * mathematical_constants::PI / 180.0 + moonMeanMotion * secondsSinceJ2000;

    // Moon position in orbital plane (inclined to ecliptic)
    double xOrbit = moonSMA * std::cos(meanLong);
    double yOrbit = moonSMA * std::sin(meanLong) * std::cos(moonInclination);
    double zOrbit = moonSMA * std::sin(meanLong) * std::sin(moonInclination);

    // Rotate from ecliptic to J2000 equatorial
    double xJ2000 = xOrbit;
    double yJ2000 = yOrbit * std::cos(obliquity) - zOrbit * std::sin(obliquity);
    double zJ2000 = yOrbit * std::sin(obliquity) + zOrbit * std::cos(obliquity);

    return Eigen::Vector3d(xJ2000, yJ2000, zJ2000);
}

// Third-body gravitational acceleration
// a = -mu_body * (r_sat_body/|r_sat_body|^3 + r_body/|r_body|^3)
Eigen::Vector3d thirdBodyAcceleration(const Eigen::Vector3d& satPos, const Eigen::Vector3d& bodyPos, double muBody) {
    Eigen::Vector3d r_sat_body = satPos - bodyPos;  // Vector from body to satellite
    double r_sb = r_sat_body.norm();
    double r_b = bodyPos.norm();

    // Third-body perturbation formula
    Eigen::Vector3d accel = -muBody * (r_sat_body / (r_sb * r_sb * r_sb) + bodyPos / (r_b * r_b * r_b));
    return accel;
}

val propagateJ2vsFullForce(
    std::string tleLine1,
    std::string tleLine2,
    double duration,            // seconds from TLE epoch
    int numPoints)
{
    using namespace ephemerides;
    using namespace orbital_element_conversions;

    // ==========================================================================
    // Earth Geopotential Parameters (EGM96/WGS84)
    // ==========================================================================
    double muEarth = 3.986004418e14;  // m^3/s^2
    double earthRadius = 6378137.0;   // m (WGS84 equatorial radius)

    // Zonal harmonics (unnormalized)
    double J2 = 1.08263e-3;           // Oblateness (dominant)
    double J3 = -2.5327e-6;           // Pear-shape (asymmetric N/S)
    double J4 = -1.6196e-6;           // Higher-order oblateness

    // ==========================================================================
    // Third-body Gravitational Parameters
    // ==========================================================================
    double muSun = 1.32712440018e20;   // m^3/s^2 (Sun)
    double muMoon = 4.9028e12;         // m^3/s^2 (Moon)

    // ==========================================================================
    // Atmospheric Drag Parameters (exponential model)
    // ==========================================================================
    double rho0 = 1.225;               // kg/m^3 (sea level density)
    double H = 8500.0;                 // m (scale height, approximate)
    double Cd = 2.2;                   // Drag coefficient (typical for satellites)
    double AreaToMass = 0.01;          // m^2/kg (typical LEO satellite)

    // ==========================================================================
    // Solar Radiation Pressure Parameters
    // ==========================================================================
    double P_sun = 4.56e-6;            // N/m^2 (solar radiation pressure at 1 AU)
    double Cr = 1.5;                   // Reflectivity coefficient (1=absorb, 2=reflect)
    double AU = 149597870700.0;        // m (astronomical unit)

    // Parse TLE to get initial state (we only use it for initial conditions)
    std::shared_ptr<Tle> tle = std::make_shared<Tle>(tleLine1, tleLine2);
    TleEphemeris sgp4Ephemeris("Earth", "J2000", tle);

    // Get initial state from TLE at epoch
    double tleEpoch = tle->getEpoch();  // seconds from J2000
    Eigen::Vector6d initialState = sgp4Ephemeris.getCartesianState(tleEpoch);

    // Debug: print initial state magnitude to verify reasonable values
    double r0 = initialState.head<3>().norm();
    double v0 = initialState.tail<3>().norm();
    EM_ASM_({
        var msg = "Initial state: r=" + $0.toFixed(1) + " km, v=" + $1.toFixed(3) + " km/s";
        console.log(msg);
    }, r0/1000.0, v0/1000.0);

    // ==========================================================================
    // J2-only State Derivative (baseline model)
    // ==========================================================================
    auto stateDerivativeJ2 = [=](double t, const Eigen::Vector6d& s) -> Eigen::Vector6d {
        Eigen::Vector6d d;
        Eigen::Vector3d pos = s.head<3>();
        Eigen::Vector3d vel = s.tail<3>();
        double r = pos.norm();
        double x = pos(0), y = pos(1), z = pos(2);

        double r2 = r * r;
        double r3 = r2 * r;
        double r5 = r2 * r3;
        Eigen::Vector3d a_twobody = -muEarth / r3 * pos;

        double z2 = z * z;
        double Re2 = earthRadius * earthRadius;
        double factorJ2 = 1.5 * J2 * muEarth * Re2 / r5;

        Eigen::Vector3d a_j2;
        a_j2(0) = factorJ2 * x * (5.0 * z2 / r2 - 1.0);
        a_j2(1) = factorJ2 * y * (5.0 * z2 / r2 - 1.0);
        a_j2(2) = factorJ2 * z * (5.0 * z2 / r2 - 3.0);

        d.head<3>() = vel;
        d.tail<3>() = a_twobody + a_j2;
        return d;
    };

    // ==========================================================================
    // Full Force Model State Derivative
    // Includes: J2/J3/J4, Sun/Moon third-body, Drag, SRP
    // ==========================================================================
    auto stateDerivativeFullForce = [=](double t, const Eigen::Vector6d& s) -> Eigen::Vector6d {
        Eigen::Vector6d d;
        Eigen::Vector3d pos = s.head<3>();
        Eigen::Vector3d vel = s.tail<3>();
        double r = pos.norm();
        double x = pos(0), y = pos(1), z = pos(2);

        // 1. Two-body acceleration (central gravity)
        double r2 = r * r;
        double r3 = r2 * r;
        double r4 = r3 * r;
        double r5 = r4 * r;
        Eigen::Vector3d a_twobody = -muEarth / r3 * pos;

        // 2. J2 perturbation acceleration
        double z2 = z * z;
        double Re2 = earthRadius * earthRadius;
        double factorJ2 = 1.5 * J2 * muEarth * Re2 / (r5);

        Eigen::Vector3d a_j2;
        a_j2(0) = factorJ2 * x * (5.0 * z2 / r2 - 1.0);
        a_j2(1) = factorJ2 * y * (5.0 * z2 / r2 - 1.0);
        a_j2(2) = factorJ2 * z * (5.0 * z2 / r2 - 3.0);

        // 3. J3 perturbation acceleration (asymmetric N/S)
        double Re3 = Re2 * earthRadius;
        double factorJ3 = 0.5 * J3 * muEarth * Re3 / (r5 * r2);

        Eigen::Vector3d a_j3;
        a_j3(0) = factorJ3 * x * (35.0 * z2 * z / r3 - 30.0 * z / r);
        a_j3(1) = factorJ3 * y * (35.0 * z2 * z / r3 - 30.0 * z / r);
        a_j3(2) = factorJ3 * (35.0 * z2 * z2 / r4 - 30.0 * z2 / r2 + 3.0);

        // 4. J4 perturbation acceleration
        double Re4 = Re3 * earthRadius;
        double z4 = z2 * z2;
        double factorJ4 = 0.625 * J4 * muEarth * Re4 / (r5 * r4);

        Eigen::Vector3d a_j4;
        a_j4(0) = factorJ4 * x * (63.0 * z4 / r4 - 42.0 * z2 / r2 + 3.0);
        a_j4(1) = factorJ4 * y * (63.0 * z4 / r4 - 42.0 * z2 / r2 + 3.0);
        a_j4(2) = factorJ4 * z * (63.0 * z4 / r4 - 70.0 * z2 / r2 + 15.0);

        // 5. Sun third-body perturbation
        double currentJD = tleEpoch + t / 86400.0;
        Eigen::Vector3d sunPos = getApproxSunPosition(currentJD);
        Eigen::Vector3d a_sun = thirdBodyAcceleration(pos, sunPos, muSun);

        // 6. Moon third-body perturbation
        Eigen::Vector3d moonPos = getApproxMoonPosition(currentJD);
        Eigen::Vector3d a_moon = thirdBodyAcceleration(pos, moonPos, muMoon);

        // 7. Atmospheric Drag (exponential model)
        double altitude = r - earthRadius;
        Eigen::Vector3d a_drag = Eigen::Vector3d::Zero();
        if (altitude < 1000000.0 && altitude > 0) {  // Below 1000 km
            double rho = rho0 * std::exp(-altitude / H);
            double v = vel.norm();
            if (v > 0) {
                // Drag acceleration: a = -0.5 * rho * v^2 * Cd * A/m * v_hat
                double dragMag = 0.5 * rho * v * v * Cd * AreaToMass;
                a_drag = -dragMag * vel / v;
            }
        }

        // 8. Solar Radiation Pressure
        Eigen::Vector3d r_sat_sun = pos - sunPos;
        double d_sun = r_sat_sun.norm();
        // SRP acceleration: a = -P * Cr * (A/m) * (AU/d)^2 * r_hat
        double srpMag = P_sun * Cr * AreaToMass * (AU * AU) / (d_sun * d_sun);
        Eigen::Vector3d a_srp = srpMag * r_sat_sun / d_sun;

        // Total acceleration - full force model
        d.head<3>() = vel;
        d.tail<3>() = a_twobody + a_j2 + a_j3 + a_j4 + a_sun + a_moon + a_drag + a_srp;
        return d;
    };

    std::vector<double> result;
    result.reserve(numPoints * 7);  // t, j2_xyz, full_xyz

    // Use RK78 for high-fidelity propagation with both models
    using namespace numerical_integrators;

    RungeKuttaCoefficients rkCoeffs =
        RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);

    // Estimate orbital period for initial step
    double sma = initialState.head<3>().norm();  // Approximate SMA from position
    double period = 2.0 * mathematical_constants::PI * std::sqrt(std::pow(sma, 3) / muEarth);
    double initialStep = period / 500.0;

    // Create TWO separate integrators - one for J2-only, one for full force
    RungeKuttaVariableStepSizeIntegrator<double, Eigen::Vector6d> integratorJ2(
        rkCoeffs, stateDerivativeJ2, 0.0, initialState,
        std::numeric_limits<double>::epsilon(),
        std::numeric_limits<double>::infinity(),
        initialStep,
        1e-12, 1e-12);

    RungeKuttaVariableStepSizeIntegrator<double, Eigen::Vector6d> integratorFull(
        rkCoeffs, stateDerivativeFullForce, 0.0, initialState,
        std::numeric_limits<double>::epsilon(),
        std::numeric_limits<double>::infinity(),
        initialStep,
        1e-12, 1e-12);

    double dt = duration / (numPoints - 1);

    for (int i = 0; i < numPoints; i++) {
        double t = i * dt;

        // Get J2-only state
        Eigen::Vector6d j2State;
        if (i == 0) {
            j2State = initialState;
        } else {
            j2State = integratorJ2.integrateTo(t, initialStep);
        }

        // Get full force model state
        Eigen::Vector6d fullState;
        if (i == 0) {
            fullState = initialState;
        } else {
            fullState = integratorFull.integrateTo(t, initialStep);
        }

        // Debug: print first few and last sample
        if (i < 3 || i == numPoints - 1) {
            double sep = (j2State.head<3>() - fullState.head<3>()).norm();
            EM_ASM_({
                var msg = "t=" + $0.toFixed(0) + "s: J2-Full separation=" + $1.toFixed(3) + " km";
                console.log(msg);
            }, t, sep/1000.0);
        }

        result.push_back(t);
        result.push_back(j2State(0));
        result.push_back(j2State(1));
        result.push_back(j2State(2));
        result.push_back(fullState(0));
        result.push_back(fullState(1));
        result.push_back(fullState(2));
    }

    return val(typed_memory_view(result.size(), result.data())).call<val>("slice");
}

// Backward compatibility aliases
val propagateSGP4vsFullForce(
    std::string tleLine1,
    std::string tleLine2,
    double duration,
    int numPoints)
{
    return propagateJ2vsFullForce(tleLine1, tleLine2, duration, numPoints);
}

val propagateSGP4vsJ2(
    std::string tleLine1,
    std::string tleLine2,
    double duration,
    int numPoints)
{
    return propagateJ2vsFullForce(tleLine1, tleLine2, duration, numPoints);
}

// ============================================================================
// Numerical Integration (RK78) for comparison with analytical
// ============================================================================

val propagateRK78Orbit(
    double semiMajorAxis,       // km (converted to m internally)
    double eccentricity,
    double inclination,         // degrees
    double raan,                // degrees
    double argPeriapsis,        // degrees
    double trueAnomaly0,        // degrees
    double duration,            // seconds
    int numPoints,
    double relTolerance,
    double absTolerance)
{
    using namespace numerical_integrators;
    using namespace orbital_element_conversions;

    // Standard Earth gravitational parameter in m^3/s^2 (same as analytical)
    // GM_Earth = 3.986004418e14 m^3/s^2 (WGS84 value)
    double mu = 3.986004418e14;

    // Convert to SI units (meters, radians)
    double smaMeters = semiMajorAxis * 1000.0;

    // Keplerian elements in Tudat order: a, e, i, argPeriapsis, RAAN, trueAnomaly
    Eigen::Vector6d keplerElements;
    keplerElements << smaMeters,
                      eccentricity,
                      inclination * mathematical_constants::PI / 180.0,
                      argPeriapsis * mathematical_constants::PI / 180.0,  // index 3
                      raan * mathematical_constants::PI / 180.0,          // index 4
                      trueAnomaly0 * mathematical_constants::PI / 180.0;

    // Convert to Cartesian (meters, m/s)
    Eigen::Vector6d initialState = convertKeplerianToCartesianElements(keplerElements, mu);

    // State derivative function (two-body problem in SI units)
    auto stateDerivative = [mu](const double t, const Eigen::Vector6d& state) -> Eigen::Vector6d {
        Eigen::Vector6d derivative;
        Eigen::Vector3d position = state.head<3>();
        Eigen::Vector3d velocity = state.tail<3>();

        double r = position.norm();
        Eigen::Vector3d acceleration = -mu / (r * r * r) * position;

        derivative.head<3>() = velocity;
        derivative.tail<3>() = acceleration;
        return derivative;
    };

    // Create RK78 integrator
    RungeKuttaCoefficients rkCoeffs =
        RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);

    double period = 2.0 * mathematical_constants::PI * std::sqrt(std::pow(smaMeters, 3) / mu);
    double initialStep = period / 1000.0;

    RungeKuttaVariableStepSizeIntegrator<double, Eigen::Vector6d> integrator(
        rkCoeffs, stateDerivative, 0.0, initialState,
        std::numeric_limits<double>::epsilon(),
        std::numeric_limits<double>::infinity(),
        initialStep,
        relTolerance, absTolerance);

    std::vector<double> ephemeris;
    ephemeris.reserve(numPoints * 4);

    double dt = duration / (numPoints - 1);

    for (int i = 0; i < numPoints; i++) {
        double targetTime = i * dt;

        // Integrate to target time
        Eigen::Vector6d state;
        if (i == 0) {
            state = initialState;
        } else {
            state = integrator.integrateTo(targetTime, initialStep);
        }

        ephemeris.push_back(targetTime);
        ephemeris.push_back(state(0));  // Already in meters
        ephemeris.push_back(state(1));
        ephemeris.push_back(state(2));
    }

    return val(typed_memory_view(ephemeris.size(), ephemeris.data())).call<val>("slice");
}

// ============================================================================
// Libration Points (Earth-Moon L1-L5)
// ============================================================================

val computeLibrationPoints(double massParameter) {
    using namespace circular_restricted_three_body_problem;
    using namespace root_finders;

    std::shared_ptr<NewtonRaphson<double>> rootFinder =
        std::make_shared<NewtonRaphson<double>>(1e-14, 1000);

    LibrationPoint librationPoint(massParameter, rootFinder);

    std::vector<double> result;
    result.reserve(15);  // 5 points * 3 coords

    // L1
    librationPoint.computeLocationOfLibrationPoint(LibrationPoint::l1);
    Eigen::Vector3d l1 = librationPoint.getLocationOfLagrangeLibrationPoint();
    result.push_back(l1(0)); result.push_back(l1(1)); result.push_back(l1(2));

    // L2
    librationPoint.computeLocationOfLibrationPoint(LibrationPoint::l2);
    Eigen::Vector3d l2 = librationPoint.getLocationOfLagrangeLibrationPoint();
    result.push_back(l2(0)); result.push_back(l2(1)); result.push_back(l2(2));

    // L3
    librationPoint.computeLocationOfLibrationPoint(LibrationPoint::l3);
    Eigen::Vector3d l3 = librationPoint.getLocationOfLagrangeLibrationPoint();
    result.push_back(l3(0)); result.push_back(l3(1)); result.push_back(l3(2));

    // L4
    librationPoint.computeLocationOfLibrationPoint(LibrationPoint::l4);
    Eigen::Vector3d l4 = librationPoint.getLocationOfLagrangeLibrationPoint();
    result.push_back(l4(0)); result.push_back(l4(1)); result.push_back(l4(2));

    // L5
    librationPoint.computeLocationOfLibrationPoint(LibrationPoint::l5);
    Eigen::Vector3d l5 = librationPoint.getLocationOfLagrangeLibrationPoint();
    result.push_back(l5(0)); result.push_back(l5(1)); result.push_back(l5(2));

    return val(typed_memory_view(result.size(), result.data())).call<val>("slice");
}

// ============================================================================
// CR3BP Propagation (for Jacobi energy visualization)
// ============================================================================

val propagateCR3BP(
    double massParameter,
    double x0, double y0, double z0,
    double vx0, double vy0, double vz0,
    double duration,
    int numPoints)
{
    using namespace propagators;
    using namespace numerical_integrators;

    Eigen::Vector6d initialState;
    initialState << x0, y0, z0, vx0, vy0, vz0;

    double timeStep = duration / numPoints / 10.0;  // Fine step for accuracy

    std::shared_ptr<IntegratorSettings<>> integratorSettings =
        std::make_shared<IntegratorSettings<>>(rungeKutta4, 0.0, timeStep);

    std::map<double, Eigen::Vector6d> stateHistory = performCR3BPIntegration(
        integratorSettings, massParameter, initialState, 0.0, duration);

    // Sample at requested number of points
    std::vector<double> ephemeris;
    ephemeris.reserve(numPoints * 7);  // t, x, y, z, vx, vy, vz

    double dt = duration / (numPoints - 1);
    auto it = stateHistory.begin();

    for (int i = 0; i < numPoints && it != stateHistory.end(); i++) {
        double targetTime = i * dt;

        // Find closest state in history
        while (std::next(it) != stateHistory.end() &&
               std::abs(std::next(it)->first - targetTime) < std::abs(it->first - targetTime)) {
            ++it;
        }

        ephemeris.push_back(it->first);
        for (int j = 0; j < 6; j++) {
            ephemeris.push_back(it->second(j));
        }
    }

    return val(typed_memory_view(ephemeris.size(), ephemeris.data())).call<val>("slice");
}

// ============================================================================
// Lambert Targeting (Izzo algorithm)
// ============================================================================

val solveLambertProblem(
    double r1x, double r1y, double r1z,
    double r2x, double r2y, double r2z,
    double timeOfFlight,
    double mu,
    bool prograde)
{
    Eigen::Vector3d r1(r1x, r1y, r1z);
    Eigen::Vector3d r2(r2x, r2y, r2z);

    // Use Tudat's Izzo Lambert targeter (zero-revolution)
    mission_segments::ZeroRevolutionLambertTargeterIzzo lambertTargeter(
        r1, r2, timeOfFlight, mu, !prograde);  // isRetrograde = !prograde

    // Get the computed velocities
    Eigen::Vector3d v1 = lambertTargeter.getInertialVelocityAtDeparture();
    Eigen::Vector3d v2 = lambertTargeter.getInertialVelocityAtArrival();

    std::vector<double> result;
    result.reserve(6);

    result.push_back(v1(0));
    result.push_back(v1(1));
    result.push_back(v1(2));
    result.push_back(v2(0));
    result.push_back(v2(1));
    result.push_back(v2(2));

    return val(typed_memory_view(result.size(), result.data())).call<val>("slice");
}

// ============================================================================
// OMM vs Two-Body Propagation Comparison
// Compares OMM mean-element propagation (with J2 secular rates) against
// pure Two-Body (Kepler) propagation. Shows the effect of J2 secular
// perturbations on orbital plane orientation (RAAN drift, argp drift).
// ============================================================================

val propagateOMMvsJ2(
    std::string ommJson,        // JSON with OMM elements
    double duration,            // seconds
    int numPoints)
{
    using namespace orbital_element_conversions;

    // Parse OMM JSON (simple parsing for known format)
    double semiMajorAxis = 6793.0;  // km
    double eccentricity = 0.0001;
    double inclination = 51.6;      // degrees
    double raan = 45.0;             // degrees
    double argPeriapsis = 90.0;     // degrees
    double meanAnomaly = 0.0;       // degrees

    auto parseValue = [&ommJson](const std::string& key) -> double {
        size_t pos = ommJson.find("\"" + key + "\"");
        if (pos != std::string::npos) {
            pos = ommJson.find(":", pos);
            if (pos != std::string::npos) {
                return std::stod(ommJson.substr(pos + 1));
            }
        }
        return 0.0;
    };

    semiMajorAxis = parseValue("semiMajorAxis");
    eccentricity = parseValue("eccentricity");
    inclination = parseValue("inclination");
    raan = parseValue("raan");
    argPeriapsis = parseValue("argPeriapsis");
    meanAnomaly = parseValue("meanAnomaly");

    // Earth parameters
    double muEarth = 3.986004418e14;  // m^3/s^2
    double earthRadius = 6378137.0;   // m
    double J2 = 1.08263e-3;

    // Convert OMM to SI units
    double smaMeters = semiMajorAxis * 1000.0;
    double incRad = inclination * mathematical_constants::PI / 180.0;
    double raanRad = raan * mathematical_constants::PI / 180.0;
    double argpRad = argPeriapsis * mathematical_constants::PI / 180.0;
    double maRad = meanAnomaly * mathematical_constants::PI / 180.0;

    // Compute mean motion and orbital parameters
    double n = std::sqrt(muEarth / std::pow(smaMeters, 3));
    double p = smaMeters * (1 - eccentricity * eccentricity);
    double cosI = std::cos(incRad);
    double sinI = std::sin(incRad);

    // J2 secular rates (rad/s) for OMM propagation
    double raanDot = -1.5 * J2 * std::pow(earthRadius / p, 2) * n * cosI;
    double argpDot = 0.75 * J2 * std::pow(earthRadius / p, 2) * n * (5.0 * cosI * cosI - 1.0);

    // Convert initial mean anomaly to true anomaly
    double E0 = maRad;
    for (int iter = 0; iter < 20; iter++) {
        E0 = E0 - (E0 - eccentricity * std::sin(E0) - maRad) / (1 - eccentricity * std::cos(E0));
    }
    double ta0 = 2.0 * std::atan2(
        std::sqrt(1 + eccentricity) * std::sin(E0 / 2),
        std::sqrt(1 - eccentricity) * std::cos(E0 / 2));

    // Initial Keplerian elements for Two-Body propagation (Tudat order)
    Eigen::Vector6d keplerElements;
    keplerElements << smaMeters, eccentricity, incRad, argpRad, raanRad, ta0;

    std::vector<double> result;
    result.reserve(numPoints * 7);

    double dt = duration / (numPoints - 1);

    for (int i = 0; i < numPoints; i++) {
        double t = i * dt;

        // ==========================================================================
        // OMM Mean Element Propagation (with J2 secular rates)
        // RAAN and argument of perigee drift due to J2
        // ==========================================================================
        double M = maRad + n * t;  // Mean anomaly at time t

        // Evolving mean elements due to J2 secular perturbations
        double currentRaan = raanRad + raanDot * t;
        double currentArgp = argpRad + argpDot * t;

        // Solve Kepler's equation
        double Ecurrent = M;
        for (int iter = 0; iter < 20; iter++) {
            Ecurrent = Ecurrent - (Ecurrent - eccentricity * std::sin(Ecurrent) - M) /
                       (1 - eccentricity * std::cos(Ecurrent));
        }

        double taCurrent = 2.0 * std::atan2(
            std::sqrt(1 + eccentricity) * std::sin(Ecurrent / 2),
            std::sqrt(1 - eccentricity) * std::cos(Ecurrent / 2));

        // Position in orbital frame
        double rMag = smaMeters * (1 - eccentricity * std::cos(Ecurrent));
        double xOrb = rMag * std::cos(taCurrent);
        double yOrb = rMag * std::sin(taCurrent);

        // Rotation to inertial frame with DRIFTING RAAN and argp
        double cosRaan = std::cos(currentRaan);
        double sinRaan = std::sin(currentRaan);
        double cosArgp = std::cos(currentArgp);
        double sinArgp = std::sin(currentArgp);

        double xOMM = (cosRaan * cosArgp - sinRaan * sinArgp * cosI) * xOrb +
                     (-cosRaan * sinArgp - sinRaan * cosArgp * cosI) * yOrb;
        double yOMM = (sinRaan * cosArgp + cosRaan * sinArgp * cosI) * xOrb +
                     (-sinRaan * sinArgp + cosRaan * cosArgp * cosI) * yOrb;
        double zOMM = (sinArgp * sinI) * xOrb + (cosArgp * sinI) * yOrb;

        // ==========================================================================
        // Two-Body (Kepler) Propagation - NO J2 secular drift
        // Uses Tudat's propagateKeplerOrbit (pure two-body, fixed orbital plane)
        // ==========================================================================
        Eigen::Vector6d propagatedKepler = propagateKeplerOrbit(keplerElements, t, muEarth);
        Eigen::Vector6d twoBodyCartesian = convertKeplerianToCartesianElements(propagatedKepler, muEarth);

        result.push_back(t);
        result.push_back(xOMM);                // OMM x (meters) - with J2 secular drift
        result.push_back(yOMM);                // OMM y
        result.push_back(zOMM);                // OMM z
        result.push_back(twoBodyCartesian(0)); // Two-Body x (meters) - no drift
        result.push_back(twoBodyCartesian(1)); // Two-Body y
        result.push_back(twoBodyCartesian(2)); // Two-Body z
    }

    return val(typed_memory_view(result.size(), result.data())).call<val>("slice");
}

// ============================================================================
// Orbit Determination / Differential Correction Visualization
// Demonstrates batch least squares orbit determination with position observations.
// Compares convergence using OMM (mean elements with secular drift) vs Full Force.
// ============================================================================

// Helper: Propagate using OMM mean elements with J2 secular rates
void propagateOMMState(
    const Eigen::Vector6d& keplerElements,  // Initial Keplerian elements (SI)
    double t,                                // Time from epoch (seconds)
    double muEarth,
    double earthRadius,
    double J2,
    Eigen::Vector3d& posOut)
{
    double sma = keplerElements(0);
    double ecc = keplerElements(1);
    double inc = keplerElements(2);
    double argp = keplerElements(3);
    double raan = keplerElements(4);
    double ta0 = keplerElements(5);

    // Convert true anomaly to mean anomaly at t=0
    double E0 = 2.0 * std::atan2(
        std::sqrt(1 - ecc) * std::sin(ta0 / 2),
        std::sqrt(1 + ecc) * std::cos(ta0 / 2));
    double M0 = E0 - ecc * std::sin(E0);

    // Mean motion and J2 secular rates
    double n = std::sqrt(muEarth / std::pow(sma, 3));
    double p = sma * (1 - ecc * ecc);
    double cosI = std::cos(inc);
    double sinI = std::sin(inc);

    double raanDot = -1.5 * J2 * std::pow(earthRadius / p, 2) * n * cosI;
    double argpDot = 0.75 * J2 * std::pow(earthRadius / p, 2) * n * (5.0 * cosI * cosI - 1.0);

    // Propagate mean elements
    double M = M0 + n * t;
    double currentRaan = raan + raanDot * t;
    double currentArgp = argp + argpDot * t;

    // Solve Kepler's equation
    double E = M;
    for (int iter = 0; iter < 20; iter++) {
        E = E - (E - ecc * std::sin(E) - M) / (1 - ecc * std::cos(E));
    }

    double ta = 2.0 * std::atan2(
        std::sqrt(1 + ecc) * std::sin(E / 2),
        std::sqrt(1 - ecc) * std::cos(E / 2));

    // Position in orbital frame
    double r = sma * (1 - ecc * std::cos(E));
    double xOrb = r * std::cos(ta);
    double yOrb = r * std::sin(ta);

    // Rotation to inertial frame
    double cosRaan = std::cos(currentRaan);
    double sinRaan = std::sin(currentRaan);
    double cosArgp = std::cos(currentArgp);
    double sinArgp = std::sin(currentArgp);

    posOut(0) = (cosRaan * cosArgp - sinRaan * sinArgp * cosI) * xOrb +
                (-cosRaan * sinArgp - sinRaan * cosArgp * cosI) * yOrb;
    posOut(1) = (sinRaan * cosArgp + cosRaan * sinArgp * cosI) * xOrb +
                (-sinRaan * sinArgp + cosRaan * cosArgp * cosI) * yOrb;
    posOut(2) = (sinArgp * sinI) * xOrb + (cosArgp * sinI) * yOrb;
}

// Helper: Propagate using full force model (J2-J4 + Sun/Moon)
void propagateFullForceState(
    const Eigen::Vector6d& initialCartesian,
    double t,
    double muEarth,
    double earthRadius,
    double J2,
    double J3,
    double J4,
    double tleEpoch,  // For sun/moon positions
    Eigen::Vector3d& posOut)
{
    using namespace numerical_integrators;

    auto stateDerivative = [&](double currentT, const Eigen::Vector6d& state) -> Eigen::Vector6d {
        Eigen::Vector3d pos = state.head<3>();
        Eigen::Vector3d vel = state.tail<3>();

        double r = pos.norm();
        double r2 = r * r;
        double r3 = r2 * r;
        double r5 = r3 * r2;

        double x = pos(0), y = pos(1), z = pos(2);
        double z2 = z * z;

        // Two-body
        Eigen::Vector3d a_twobody = -muEarth / r3 * pos;

        // J2
        double Re2 = earthRadius * earthRadius;
        double factorJ2 = 1.5 * J2 * muEarth * Re2 / r5;
        Eigen::Vector3d a_j2;
        a_j2(0) = factorJ2 * x * (5.0 * z2 / r2 - 1.0);
        a_j2(1) = factorJ2 * y * (5.0 * z2 / r2 - 1.0);
        a_j2(2) = factorJ2 * z * (5.0 * z2 / r2 - 3.0);

        // J3
        double Re3 = Re2 * earthRadius;
        double factorJ3 = 0.5 * J3 * muEarth * Re3 / (r5 * r2);
        Eigen::Vector3d a_j3;
        a_j3(0) = factorJ3 * x * (35.0 * z2 * z / r3 - 30.0 * z / r);
        a_j3(1) = factorJ3 * y * (35.0 * z2 * z / r3 - 30.0 * z / r);
        a_j3(2) = factorJ3 * (35.0 * z2 * z2 / (r2*r2) - 30.0 * z2 / r2 + 3.0);

        // J4
        double Re4 = Re3 * earthRadius;
        double z4 = z2 * z2;
        double r4 = r2 * r2;
        double factorJ4 = 0.625 * J4 * muEarth * Re4 / (r5 * r4);
        Eigen::Vector3d a_j4;
        a_j4(0) = factorJ4 * x * (63.0 * z4 / r4 - 42.0 * z2 / r2 + 3.0);
        a_j4(1) = factorJ4 * y * (63.0 * z4 / r4 - 42.0 * z2 / r2 + 3.0);
        a_j4(2) = factorJ4 * z * (63.0 * z4 / r4 - 70.0 * z2 / r2 + 15.0);

        Eigen::Vector6d d;
        d.head<3>() = vel;
        d.tail<3>() = a_twobody + a_j2 + a_j3 + a_j4;
        return d;
    };

    if (std::abs(t) < 1e-10) {
        posOut = initialCartesian.head<3>();
        return;
    }

    RungeKuttaCoefficients rkCoeffs =
        RungeKuttaCoefficients::get(CoefficientSets::rungeKuttaFehlberg78);

    double period = 2.0 * mathematical_constants::PI *
        std::sqrt(std::pow(initialCartesian.head<3>().norm(), 3) / muEarth);
    double initialStep = period / 500.0;

    RungeKuttaVariableStepSizeIntegrator<double, Eigen::Vector6d> integrator(
        rkCoeffs, stateDerivative, 0.0, initialCartesian,
        std::numeric_limits<double>::epsilon(),
        std::numeric_limits<double>::infinity(),
        initialStep, 1e-12, 1e-12);

    Eigen::Vector6d finalState = integrator.integrateTo(t, initialStep);
    posOut = finalState.head<3>();
}

val runOrbitDetermination(
    std::string initialStateJson,   // Initial state guess (Keplerian elements in km, deg)
    std::string truthStateJson,     // Truth state (Keplerian elements in km, deg)
    double duration,                // Observation arc duration (seconds)
    int numObservations,            // Number of position observations
    int numOrbitSamples,            // Number of samples for orbit rendering
    double noiseStdDev,             // Position noise standard deviation (meters)
    int maxIterations,              // Maximum DC iterations
    std::string dynamicsModel)      // "omm" or "fullforce"
{
    using namespace orbital_element_conversions;

    // Parse JSON helper
    auto parseValue = [](const std::string& json, const std::string& key) -> double {
        size_t pos = json.find("\"" + key + "\"");
        if (pos != std::string::npos) {
            pos = json.find(":", pos);
            if (pos != std::string::npos) {
                return std::stod(json.substr(pos + 1));
            }
        }
        return 0.0;
    };

    // Earth parameters
    double muEarth = 3.986004418e14;
    double earthRadius = 6378137.0;
    double J2 = 1.08263e-3;
    double J3 = -2.5327e-6;
    double J4 = -1.6196e-6;
    double tleEpoch = 2460000.5;  // Arbitrary epoch for sun/moon

    // Parse truth state
    double truthSma = parseValue(truthStateJson, "semiMajorAxis") * 1000.0;
    double truthEcc = parseValue(truthStateJson, "eccentricity");
    double truthInc = parseValue(truthStateJson, "inclination") * mathematical_constants::PI / 180.0;
    double truthRaan = parseValue(truthStateJson, "raan") * mathematical_constants::PI / 180.0;
    double truthArgp = parseValue(truthStateJson, "argPeriapsis") * mathematical_constants::PI / 180.0;
    double truthTa = parseValue(truthStateJson, "trueAnomaly") * mathematical_constants::PI / 180.0;

    Eigen::Vector6d truthKepler;
    truthKepler << truthSma, truthEcc, truthInc, truthArgp, truthRaan, truthTa;
    Eigen::Vector6d truthCartesian = convertKeplerianToCartesianElements(truthKepler, muEarth);

    // Parse initial guess
    double initSma = parseValue(initialStateJson, "semiMajorAxis") * 1000.0;
    double initEcc = parseValue(initialStateJson, "eccentricity");
    double initInc = parseValue(initialStateJson, "inclination") * mathematical_constants::PI / 180.0;
    double initRaan = parseValue(initialStateJson, "raan") * mathematical_constants::PI / 180.0;
    double initArgp = parseValue(initialStateJson, "argPeriapsis") * mathematical_constants::PI / 180.0;
    double initTa = parseValue(initialStateJson, "trueAnomaly") * mathematical_constants::PI / 180.0;

    Eigen::Vector6d currentKepler;
    currentKepler << initSma, initEcc, initInc, initArgp, initRaan, initTa;
    Eigen::Vector6d currentCartesian = convertKeplerianToCartesianElements(currentKepler, muEarth);

    // Generate truth observations using full force model
    // Use realistic SSN (Space Surveillance Network) ground station locations
    struct GroundStation {
        double lat;  // degrees
        double lon;  // degrees
        std::string name;
    };

    std::vector<GroundStation> ssnStations = {
        {32.9, -106.5, "White Sands"},      // New Mexico
        {21.6, -158.0, "Kaena Point"},      // Hawaii
        {-7.3, -14.4, "Ascension Island"},  // South Atlantic
        {36.0, -121.6, "Point Mugu"},       // California
        {64.3, -149.2, "Clear AFS"},        // Alaska
        {-33.2, -70.6, "Santiago"},         // Chile
        {28.2, -16.5, "Tenerife"},          // Canary Islands
        {52.4, -1.1, "Fylingdales"},        // UK
        {-21.8, 114.2, "Exmouth"},          // Australia
        {7.4, 134.5, "Palau"},              // Western Pacific
    };

    // Convert station positions to ECEF
    std::vector<Eigen::Vector3d> stationECEF;
    for (const auto& station : ssnStations) {
        double latRad = station.lat * mathematical_constants::PI / 180.0;
        double lonRad = station.lon * mathematical_constants::PI / 180.0;
        double cosLat = std::cos(latRad);
        double sinLat = std::sin(latRad);
        double cosLon = std::cos(lonRad);
        double sinLon = std::sin(lonRad);
        Eigen::Vector3d pos(earthRadius * cosLat * cosLon,
                           earthRadius * cosLat * sinLon,
                           earthRadius * sinLat);
        stationECEF.push_back(pos);
    }

    std::vector<double> obsTimes;
    std::vector<Eigen::Vector3d> observations;
    std::vector<int> obsStationIdx;  // Which station made each observation

    // Simple random number generator (linear congruential)
    unsigned int seed = 12345;
    auto nextRand = [&seed]() -> double {
        seed = seed * 1103515245 + 12345;
        return (static_cast<double>((seed >> 16) & 0x7fff) / 32767.0) * 2.0 - 1.0;
    };

    // Create clusters of observations evenly spaced around the orbit
    // Each cluster has 3 consecutive observations separated by small time intervals
    int clusterSize = 3;
    int numClusters = numObservations / clusterSize;
    if (numClusters < 1) numClusters = 1;

    double clusterSpacing = duration / numClusters;  // Time between cluster centers
    double obsInterval = 30.0;  // 30 seconds between observations within a cluster

    for (int cluster = 0; cluster < numClusters; cluster++) {
        // Center time for this cluster (evenly spaced around the orbit)
        double clusterCenterTime = (cluster + 0.5) * clusterSpacing;

        // Generate observations within this cluster
        for (int obs = 0; obs < clusterSize; obs++) {
            double t = clusterCenterTime + (obs - clusterSize / 2) * obsInterval;
            if (t < 0) t = 0;
            if (t > duration) t = duration;

            // Get truth position at this time
            Eigen::Vector3d truthPos;
            propagateFullForceState(truthCartesian, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, truthPos);

            // Find closest station (for display purposes, not visibility constrained)
            int closestStation = 0;
            double minDist = std::numeric_limits<double>::max();
            for (size_t s = 0; s < stationECEF.size(); s++) {
                double dist = (truthPos - stationECEF[s]).norm();
                if (dist < minDist) {
                    minDist = dist;
                    closestStation = static_cast<int>(s);
                }
            }

            obsTimes.push_back(t);
            obsStationIdx.push_back(closestStation);

            // Add noise to the truth position
            Eigen::Vector3d noise;
            for (int j = 0; j < 3; j++) {
                double u1 = (nextRand() + 1.0) / 2.0;
                double u2 = (nextRand() + 1.0) / 2.0;
                if (u1 < 1e-10) u1 = 1e-10;
                noise(j) = std::sqrt(-2.0 * std::log(u1)) * std::cos(2.0 * mathematical_constants::PI * u2) * noiseStdDev;
            }
            observations.push_back(truthPos + noise);
        }
    }

    // Update numObservations to actual count
    numObservations = static_cast<int>(observations.size());

    // Results storage
    std::vector<double> result;
    // Format: [numIterations, numObs,
    //          iter0_rms, iter0_state(6), iter0_residuals(numObs*3),
    //          iter1_rms, iter1_state(6), iter1_residuals(numObs*3),
    //          ...]

    bool useOMM = (dynamicsModel == "omm");

    // Position covariance matrix (3x3) - will be populated on final iteration
    Eigen::Matrix3d positionCovariance = Eigen::Matrix3d::Identity() * 1e6;  // Default large uncertainty

    // Differential correction iterations
    for (int iter = 0; iter < maxIterations; iter++) {
        // Compute residuals and build design matrix
        Eigen::MatrixXd H(numObservations * 3, 6);  // Jacobian (position w.r.t. initial state)
        Eigen::VectorXd residuals(numObservations * 3);

        for (int i = 0; i < numObservations; i++) {
            double t = obsTimes[i];
            Eigen::Vector3d computedPos;

            if (useOMM) {
                Eigen::Vector6d currentKeplerFromCart =
                    convertCartesianToKeplerianElements(currentCartesian, muEarth);
                propagateOMMState(currentKeplerFromCart, t, muEarth, earthRadius, J2, computedPos);
            } else {
                propagateFullForceState(currentCartesian, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, computedPos);
            }

            // Residuals: observed - computed
            Eigen::Vector3d res = observations[i] - computedPos;
            residuals.segment<3>(i * 3) = res;

            // Numerical partials (finite difference)
            // Use appropriate perturbation sizes: 10m for position, 0.01 m/s for velocity
            for (int j = 0; j < 6; j++) {
                Eigen::Vector6d perturbedState = currentCartesian;
                double perturbation = (j < 3) ? 10.0 : 0.01;
                perturbedState(j) += perturbation;

                Eigen::Vector3d perturbedPos;
                if (useOMM) {
                    Eigen::Vector6d perturbedKepler =
                        convertCartesianToKeplerianElements(perturbedState, muEarth);
                    propagateOMMState(perturbedKepler, t, muEarth, earthRadius, J2, perturbedPos);
                } else {
                    propagateFullForceState(perturbedState, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, perturbedPos);
                }

                Eigen::Vector3d partial = (perturbedPos - computedPos) / perturbation;
                H.block<3, 1>(i * 3, j) = partial;
            }
        }

        // Compute RMS of residuals
        double rms = std::sqrt(residuals.squaredNorm() / residuals.size());

        // Store iteration results
        result.push_back(rms);
        for (int j = 0; j < 6; j++) {
            result.push_back(currentCartesian(j));
        }
        for (int i = 0; i < numObservations * 3; i++) {
            result.push_back(residuals(i));
        }

        // Compute covariance matrix (H^T H)^-1
        Eigen::MatrixXd HtH = H.transpose() * H;

        // Add regularization for stability
        double lambda = 1e-6;
        HtH.diagonal() += Eigen::VectorXd::Constant(6, lambda);

        // Compute (H^T H)^-1
        Eigen::MatrixXd HtHinv = HtH.inverse();

        // Compute post-fit residual variance (scaled covariance)
        // sigma^2 = r^T r / (n - p) where n = num measurements, p = num parameters
        int numMeasurements = numObservations * 3;
        int numParams = 6;
        double sigma2 = residuals.squaredNorm() / std::max(numMeasurements - numParams, 1);

        // Scale covariance by sigma^2 for realism
        Eigen::MatrixXd scaledCovariance = sigma2 * HtHinv;

        // Store covariance for final output (only on last iteration)
        bool isLastIteration = (rms < noiseStdDev * 1.05 || iter == maxIterations - 1);

        if (isLastIteration) {
            // Store the 3x3 position covariance (upper-left block)
            // This will be extracted after the loop
            positionCovariance = scaledCovariance.block<3, 3>(0, 0);
        }

        // Check convergence (RMS should approach noise level)
        if (isLastIteration) {
            break;
        }

        // Least squares update: dx = (H^T H)^-1 H^T r
        Eigen::VectorXd Htr = H.transpose() * residuals;
        Eigen::VectorXd dx = HtHinv * Htr;

        currentCartesian += dx;
    }

    // Prepend metadata
    std::vector<double> finalResult;
    finalResult.push_back(static_cast<double>(result.size() / (1 + 6 + numObservations * 3)));  // Number of iterations
    finalResult.push_back(static_cast<double>(numObservations));
    finalResult.push_back(static_cast<double>(numOrbitSamples));
    finalResult.insert(finalResult.end(), result.begin(), result.end());

    // Add truth trajectory for reference (high resolution for rendering)
    double orbitDt = duration / (numOrbitSamples - 1);
    for (int i = 0; i < numOrbitSamples; i++) {
        double t = i * orbitDt;
        Eigen::Vector3d truthPos;
        propagateFullForceState(truthCartesian, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, truthPos);
        finalResult.push_back(truthPos(0));
        finalResult.push_back(truthPos(1));
        finalResult.push_back(truthPos(2));
    }

    // Add observations (with noise)
    for (int i = 0; i < numObservations; i++) {
        finalResult.push_back(observations[i](0));
        finalResult.push_back(observations[i](1));
        finalResult.push_back(observations[i](2));
    }

    // Add truth positions at observation times (for residual visualization)
    for (int i = 0; i < numObservations; i++) {
        double t = obsTimes[i];
        Eigen::Vector3d truthPos;
        propagateFullForceState(truthCartesian, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, truthPos);
        finalResult.push_back(truthPos(0));
        finalResult.push_back(truthPos(1));
        finalResult.push_back(truthPos(2));
    }

    // Add estimated positions at observation times (for residual visualization)
    for (int i = 0; i < numObservations; i++) {
        double t = obsTimes[i];
        Eigen::Vector3d estPos;
        if (useOMM) {
            Eigen::Vector6d estKepler = convertCartesianToKeplerianElements(currentCartesian, muEarth);
            propagateOMMState(estKepler, t, muEarth, earthRadius, J2, estPos);
        } else {
            propagateFullForceState(currentCartesian, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, estPos);
        }
        finalResult.push_back(estPos(0));
        finalResult.push_back(estPos(1));
        finalResult.push_back(estPos(2));
    }

    // Add estimated trajectory (high resolution for rendering)
    for (int i = 0; i < numOrbitSamples; i++) {
        double t = i * orbitDt;
        Eigen::Vector3d estPos;
        if (useOMM) {
            Eigen::Vector6d estKepler = convertCartesianToKeplerianElements(currentCartesian, muEarth);
            propagateOMMState(estKepler, t, muEarth, earthRadius, J2, estPos);
        } else {
            propagateFullForceState(currentCartesian, t, muEarth, earthRadius, J2, J3, J4, tleEpoch, estPos);
        }
        finalResult.push_back(estPos(0));
        finalResult.push_back(estPos(1));
        finalResult.push_back(estPos(2));
    }

    // Add position covariance matrix (3x3, row-major: Cxx, Cxy, Cxz, Cyx, Cyy, Cyz, Czx, Czy, Czz)
    // This represents the uncertainty ellipsoid around the estimated position
    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++) {
            finalResult.push_back(positionCovariance(i, j));
        }
    }

    return val(typed_memory_view(finalResult.size(), finalResult.data())).call<val>("slice");
}

// ============================================================================
// Atmosphere Models
// ============================================================================

val computeAtmosphereDensity(double altitude, std::string model) {
    // Exponential atmosphere model
    if (model == "exponential") {
        double rho0 = 1.225;  // kg/m^3 at sea level
        double H = 8500.0;    // Scale height in meters
        double density = rho0 * std::exp(-altitude / H);

        std::vector<double> result = {density};
        return val(typed_memory_view(result.size(), result.data())).call<val>("slice");
    }

    // For NRLMSISE-00, would need more parameters (date, location, solar activity)
    // Return placeholder
    std::vector<double> result = {0.0};
    return val(typed_memory_view(result.size(), result.data())).call<val>("slice");
}

// ============================================================================
// Precomputed Ephemeris (from pre-converted SPK data)
// These use TabulatedCartesianEphemeris with Lagrange interpolation to provide
// high-accuracy ephemeris data from pre-converted binary SPK files.
// The data is loaded as arrays from JavaScript (which parses JSON files).
// ============================================================================

// Global registry for precomputed ephemerides
// Key format: "target_observer_frame" (lowercase), e.g., "earth_sun_j2000"
std::map<std::string, std::shared_ptr<ephemerides::Ephemeris>> precomputedEphemerides;

// Helper to create ephemeris key
std::string makeEphemerisKey(const std::string& target, const std::string& observer, const std::string& frame) {
    std::string key = target + "_" + observer + "_" + frame;
    std::transform(key.begin(), key.end(), key.begin(), ::tolower);
    return key;
}

// Load precomputed ephemeris from arrays passed from JavaScript
// epochs: array of times (seconds since J2000)
// states: flat array of states [x0,y0,z0,vx0,vy0,vz0, x1,y1,z1,vx1,vy1,vz1, ...]
// Units: meters and m/s
bool loadPrecomputedEphemeris(
    const std::string& target,
    const std::string& observer,
    const std::string& frame,
    const std::vector<double>& epochs,
    const std::vector<double>& states)
{
    if (epochs.empty() || states.size() != epochs.size() * 6) {
        std::cerr << "[EPHEMERIS] Invalid data: " << epochs.size() << " epochs, "
                  << states.size() << " state values (expected " << epochs.size() * 6 << ")" << std::endl;
        return false;
    }

    try {
        // Build state history map
        std::map<double, Eigen::Vector6d> stateHistory;
        for (size_t i = 0; i < epochs.size(); ++i) {
            Eigen::Vector6d state;
            state << states[i*6], states[i*6+1], states[i*6+2],
                     states[i*6+3], states[i*6+4], states[i*6+5];
            stateHistory[epochs[i]] = state;
        }

        // Create interpolator settings (8th order Lagrange for high accuracy)
        auto interpolatorSettings = std::make_shared<interpolators::LagrangeInterpolatorSettings>(8);

        // Create tabulated ephemeris
        auto interpolator = interpolators::createOneDimensionalInterpolator(
            stateHistory, interpolatorSettings);

        auto ephemeris = std::make_shared<ephemerides::TabulatedCartesianEphemeris<>>(
            interpolator, observer, frame);

        // Store in registry
        std::string key = makeEphemerisKey(target, observer, frame);
        precomputedEphemerides[key] = ephemeris;

        std::cout << "[EPHEMERIS] Loaded precomputed ephemeris: " << target << " relative to " << observer
                  << " in " << frame << " frame (" << epochs.size() << " states, "
                  << epochs.front() << " to " << epochs.back() << " s)" << std::endl;

        return true;
    } catch (const std::exception& e) {
        std::cerr << "[EPHEMERIS] Error loading precomputed ephemeris: " << e.what() << std::endl;
        return false;
    }
}

// Check if precomputed ephemeris is available
bool isPrecomputedEphemerisAvailable(const std::string& target, const std::string& observer, const std::string& frame) {
    std::string key = makeEphemerisKey(target, observer, frame);
    return precomputedEphemerides.find(key) != precomputedEphemerides.end();
}

// Get state from precomputed ephemeris
Eigen::Vector6d getPrecomputedState(const std::string& target, const std::string& observer,
                                     const std::string& frame, double epoch) {
    std::string key = makeEphemerisKey(target, observer, frame);
    auto it = precomputedEphemerides.find(key);
    if (it == precomputedEphemerides.end()) {
        throw std::runtime_error("Precomputed ephemeris not found: " + key);
    }
    return it->second->getCartesianState(epoch);
}

// Get time bounds for precomputed ephemeris
std::pair<double, double> getPrecomputedEphemerisTimeBounds(const std::string& target,
                                                            const std::string& observer,
                                                            const std::string& frame) {
    std::string key = makeEphemerisKey(target, observer, frame);
    auto it = precomputedEphemerides.find(key);
    if (it == precomputedEphemerides.end()) {
        return {0.0, 0.0};
    }
    auto tabEph = std::dynamic_pointer_cast<ephemerides::TabulatedCartesianEphemeris<>>(it->second);
    if (tabEph) {
        return tabEph->getSafeInterpolationInterval();
    }
    return {0.0, 0.0};
}

// Clear all precomputed ephemerides
void clearPrecomputedEphemerides() {
    precomputedEphemerides.clear();
    std::cout << "[EPHEMERIS] Cleared all precomputed ephemerides" << std::endl;
}

// List loaded precomputed ephemerides
std::vector<std::string> listPrecomputedEphemerides() {
    std::vector<std::string> keys;
    for (const auto& pair : precomputedEphemerides) {
        keys.push_back(pair.first);
    }
    return keys;
}

// ============================================================================
// Analytical Planetary Ephemeris (JPL Approximate Positions)
// These use Tudat's built-in analytical ephemeris models and do NOT require
// SPICE kernels. Use these for WASM applications where binary SPK files
// cannot be loaded.
// Reference: Standish, E.M. "Keplerian Elements for Approximate Positions
//            of the Major Planets" (JPL)
// ============================================================================

// Get planetary state using JPL approximate positions algorithm
// This is purely computational and requires no external data files.
// Supports: Mercury, Venus, Earth, Mars, Jupiter, Saturn, Uranus, Neptune, Pluto
// Returns state in ECLIPJ2000 frame relative to Sun, in meters and m/s
Eigen::Vector6d getApproximatePlanetState(const std::string& bodyName, double secondsSinceJ2000) {
    try {
        ephemerides::ApproximateJplEphemeris ephemeris(bodyName);
        return ephemeris.getCartesianState(secondsSinceJ2000);
    } catch (const std::exception& e) {
        std::cerr << "[EPHEMERIS] Error getting state for " << bodyName << ": " << e.what() << std::endl;
        return Eigen::Vector6d::Zero();
    }
}

// Alternative: GTOP ephemeris (ESA algorithm, slightly different accuracy characteristics)
// Also purely computational, no external data required.
Eigen::Vector6d getGtopPlanetState(const std::string& bodyName, double secondsSinceJ2000) {
    try {
        ephemerides::ApproximateGtopEphemeris ephemeris(bodyName);
        return ephemeris.getCartesianState(secondsSinceJ2000);
    } catch (const std::exception& e) {
        std::cerr << "[EPHEMERIS] Error getting GTOP state for " << bodyName << ": " << e.what() << std::endl;
        return Eigen::Vector6d::Zero();
    }
}

// Check if a planet name is valid for analytical ephemeris
bool isValidPlanetName(const std::string& bodyName) {
    static const std::vector<std::string> validNames = {
        "Mercury", "Venus", "Earth", "Mars", "Jupiter", "Saturn", "Uranus", "Neptune", "Pluto"
    };
    for (const auto& name : validNames) {
        if (bodyName == name) return true;
    }
    return false;
}

// ============================================================================
// Emscripten Bindings
// ============================================================================

EMSCRIPTEN_BINDINGS(tudat_visualization) {
    // Kepler orbit propagation
    function("propagateKeplerOrbit", &propagateKeplerOrbit);

    // RK4 fixed-step numerical integration (shows error growth)
    function("propagateRK4Orbit", &propagateRK4Orbit);

    // J2-only vs Full Force Model numerical propagation comparison
    function("propagateJ2vsFullForce", &propagateJ2vsFullForce);
    function("propagateSGP4vsFullForce", &propagateSGP4vsFullForce);  // Alias for backwards compatibility
    function("propagateSGP4vsJ2", &propagateSGP4vsJ2);  // Alias for backwards compatibility

    // RK78 variable-step numerical integration
    function("propagateRK78Orbit", &propagateRK78Orbit);

    // Libration points
    function("computeLibrationPoints", &computeLibrationPoints);

    // CR3BP propagation
    function("propagateCR3BP", &propagateCR3BP);

    // Lambert targeting
    function("solveLambertProblem", &solveLambertProblem);

    // Atmosphere
    function("computeAtmosphereDensity", &computeAtmosphereDensity);

    // OMM vs J2 comparison
    function("propagateOMMvsJ2", &propagateOMMvsJ2);

    // Orbit determination / differential correction
    function("runOrbitDetermination", &runOrbitDetermination);

    // ============================================================================
    // SPICE Interface Functions
    // ============================================================================

    // Time conversion functions
    function("interface_spice_convert_julian_date_to_ephemeris_time",
        &tsi::convertJulianDateToEphemerisTime);

    function("interface_spice_convert_ephemeris_time_to_julian_date",
        &tsi::convertEphemerisTimeToJulianDate);

    function("interface_spice_convert_date_string_to_ephemeris_time",
        &tsi::convertDateStringToEphemerisTime);

    // Position/state functions - wrapped to return JS arrays
    function("interface_spice_get_body_cartesian_position_at_epoch",
        optional_override([](const std::string& targetBody, const std::string& observerBody,
                            const std::string& frame, const std::string& aberration, double epoch) {
            Eigen::Vector3d pos = tsi::getBodyCartesianPositionAtEpoch(
                targetBody, observerBody, frame, aberration, epoch);
            return val(typed_memory_view(3, pos.data())).call<val>("slice");
        }));

    function("interface_spice_get_body_cartesian_state_at_epoch",
        optional_override([](const std::string& targetBody, const std::string& observerBody,
                            const std::string& frame, const std::string& aberration, double epoch) {
            Eigen::Vector6d state = tsi::getBodyCartesianStateAtEpoch(
                targetBody, observerBody, frame, aberration, epoch);
            return val(typed_memory_view(6, state.data())).call<val>("slice");
        }));

    // Body properties
    function("interface_spice_get_body_gravitational_parameter",
        &tsi::getBodyGravitationalParameter);

    function("interface_spice_get_average_radius",
        &tsi::getAverageRadius);

    // Kernel management
    function("interface_spice_load_kernel",
        &tsi::loadSpiceKernelInTudat);

    function("interface_spice_clear_kernels",
        &tsi::clearSpiceKernels);

    function("interface_spice_get_total_count_of_kernels_loaded",
        &tsi::getTotalCountOfKernelsLoaded);

    // Body identification
    function("interface_spice_convert_body_name_to_naif_id",
        &tsi::convertBodyNameToNaifId);

    // ============================================================================
    // Precomputed Ephemeris (from pre-converted SPK data)
    // High-accuracy ephemeris from pre-converted binary SPK files.
    // Data is loaded from JSON via JavaScript, passed as arrays to C++.
    // ============================================================================

    // Load precomputed ephemeris data
    // epochs: Float64Array of times (seconds since J2000)
    // states: Float64Array of states [x0,y0,z0,vx0,vy0,vz0, x1,y1,...]
    function("ephemeris_load_precomputed",
        optional_override([](const std::string& target, const std::string& observer,
                            const std::string& frame, val epochsVal, val statesVal) {
            // Convert JavaScript arrays to std::vector
            std::vector<double> epochs = vecFromJSArray<double>(epochsVal);
            std::vector<double> states = vecFromJSArray<double>(statesVal);
            return loadPrecomputedEphemeris(target, observer, frame, epochs, states);
        }));

    // Check if precomputed ephemeris is available for a body
    function("ephemeris_is_precomputed_available",
        &isPrecomputedEphemerisAvailable);

    // Get state from precomputed ephemeris - returns [x, y, z, vx, vy, vz] in m and m/s
    // Returns null if ephemeris is not found or epoch is out of bounds
    function("ephemeris_get_precomputed_state",
        optional_override([](const std::string& target, const std::string& observer,
                            const std::string& frame, double epoch) -> val {
            try {
                Eigen::Vector6d state = getPrecomputedState(target, observer, frame, epoch);
                return val(typed_memory_view(6, state.data())).call<val>("slice");
            } catch (const std::exception& e) {
                std::cerr << "[EPHEMERIS] Error in getPrecomputedState: " << e.what() << std::endl;
                return val::null();
            }
        }));

    // Get time bounds for precomputed ephemeris - returns [startEpoch, endEpoch]
    function("ephemeris_get_precomputed_time_bounds",
        optional_override([](const std::string& target, const std::string& observer,
                            const std::string& frame) {
            auto bounds = getPrecomputedEphemerisTimeBounds(target, observer, frame);
            std::vector<double> result = {bounds.first, bounds.second};
            return val(typed_memory_view(2, result.data())).call<val>("slice");
        }));

    // Clear all precomputed ephemerides
    function("ephemeris_clear_precomputed", &clearPrecomputedEphemerides);

    // List all loaded precomputed ephemerides
    function("ephemeris_list_precomputed",
        optional_override([]() {
            auto keys = listPrecomputedEphemerides();
            val result = val::array();
            for (const auto& key : keys) {
                result.call<void>("push", key);
            }
            return result;
        }));

    // ============================================================================
    // Analytical Planetary Ephemeris (no SPICE kernels required)
    // Use these instead of SPICE functions in WASM where binary SPK files
    // cannot be loaded due to f2c I/O limitations.
    // ============================================================================

    // JPL Approximate Positions - returns [x, y, z, vx, vy, vz] in meters and m/s
    // Frame: ECLIPJ2000, relative to Sun
    // Returns null if body name is invalid or ephemeris computation fails
    function("ephemeris_get_planet_state",
        optional_override([](const std::string& bodyName, double secondsSinceJ2000) -> val {
            try {
                Eigen::Vector6d state = getApproximatePlanetState(bodyName, secondsSinceJ2000);
                // Check for zero state (indicates error in getApproximatePlanetState)
                if (state.norm() < 1e-10) {
                    return val::null();
                }
                return val(typed_memory_view(6, state.data())).call<val>("slice");
            } catch (const std::exception& e) {
                std::cerr << "[EPHEMERIS] Error in ephemeris_get_planet_state: " << e.what() << std::endl;
                return val::null();
            }
        }));

    // GTOP (ESA) algorithm - alternative with slightly different accuracy
    // Returns null if body name is invalid or ephemeris computation fails
    function("ephemeris_get_planet_state_gtop",
        optional_override([](const std::string& bodyName, double secondsSinceJ2000) -> val {
            try {
                Eigen::Vector6d state = getGtopPlanetState(bodyName, secondsSinceJ2000);
                if (state.norm() < 1e-10) {
                    return val::null();
                }
                return val(typed_memory_view(6, state.data())).call<val>("slice");
            } catch (const std::exception& e) {
                std::cerr << "[EPHEMERIS] Error in ephemeris_get_planet_state_gtop: " << e.what() << std::endl;
                return val::null();
            }
        }));

    // Check if planet name is valid for analytical ephemeris
    function("ephemeris_is_valid_planet", &isValidPlanetName);
}

#endif // __EMSCRIPTEN__
