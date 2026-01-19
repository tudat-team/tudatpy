/**
 * Example 23: Linear Sensitivity Analysis
 *
 * Ported from: examples/tudatpy/propagation/linear_sensitivity_analysis.py
 *
 * This example demonstrates how to analyze the sensitivity of an orbit
 * to changes in initial conditions and parameters.
 *
 * Run with: node 23_linear_sensitivity.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Linear Sensitivity Analysis Example ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const earthGM = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;  // m
    const J2 = 1.08263e-3;  // Earth's J2 coefficient

    // Nominal orbit parameters (LEO sun-synchronous-like)
    const altitude = 600.0e3;  // 600 km
    const semiMajorAxis = earthRadius + altitude;
    const eccentricity = 0.001;
    const inclination = 97.8 * Math.PI / 180;  // Sun-sync inclination

    // Compute orbital characteristics
    const n = Math.sqrt(earthGM / Math.pow(semiMajorAxis, 3));  // mean motion
    const period = 2 * Math.PI / n;

    console.log('Nominal Orbit:');
    console.log(`  Altitude: ${(altitude / 1000).toFixed(1)} km`);
    console.log(`  Semi-major axis: ${(semiMajorAxis / 1000).toFixed(2)} km`);
    console.log(`  Eccentricity: ${eccentricity}`);
    console.log(`  Inclination: ${(inclination * 180 / Math.PI).toFixed(2)} deg`);
    console.log(`  Orbital period: ${(period / 60).toFixed(2)} minutes`);

    // Perturbation sizes for finite differencing
    const da = 100.0;  // m (semi-major axis perturbation)
    const de = 1e-5;   // eccentricity perturbation
    const di = 0.01 * Math.PI / 180;  // inclination perturbation (0.01 deg)

    console.log('\nPerturbation sizes:');
    console.log(`  da: ${da.toFixed(0)} m`);
    console.log(`  de: ${de.toExponential(1)}`);
    console.log(`  di: ${(di * 180 / Math.PI * 3600).toFixed(1)} arcsec`);

    // Compute J2 secular drift rates
    // dOmega/dt = -3/2 * n * J2 * (Re/a)^2 * cos(i)
    // domega/dt = 3/4 * n * J2 * (Re/a)^2 * (5*cos^2(i) - 1)
    // dM/dt = n + 3/4 * n * J2 * (Re/a)^2 * sqrt(1-e^2) * (3*cos^2(i) - 1)

    const aRatio = earthRadius / semiMajorAxis;
    const aRatio2 = aRatio * aRatio;
    const cosI = Math.cos(inclination);
    const cosI2 = cosI * cosI;
    const sqrtE = Math.sqrt(1 - eccentricity * eccentricity);

    const OmegaDot = -1.5 * n * J2 * aRatio2 * cosI;
    const omegaDot = 0.75 * n * J2 * aRatio2 * (5 * cosI2 - 1);
    const MDot = n + 0.75 * n * J2 * aRatio2 * sqrtE * (3 * cosI2 - 1);

    console.log('\nJ2 Secular Drift Rates:');
    console.log(`  RAAN drift: ${(OmegaDot * 180 / Math.PI * 86400).toFixed(4)} deg/day`);
    console.log(`  Arg perigee drift: ${(omegaDot * 180 / Math.PI * 86400).toFixed(4)} deg/day`);
    console.log(`  Mean motion correction: ${((MDot - n) / n * 1e6).toFixed(2)} ppm`);

    // Sensitivity of RAAN drift to inclination
    // d(OmegaDot)/di = 3/2 * n * J2 * (Re/a)^2 * sin(i)
    const sinI = Math.sin(inclination);
    const dOmegaDot_di = 1.5 * n * J2 * aRatio2 * sinI;

    console.log('\nSensitivity Analysis:');
    console.log('  Sensitivity of RAAN drift to inclination:');
    console.log(`    d(OmegaDot)/di = ${(dOmegaDot_di * 180 / Math.PI * 86400 / (Math.PI/180)).toFixed(4)} (deg/day)/(deg)`);

    // Change in RAAN drift for 0.01 deg inclination change
    const deltaOmegaDot = dOmegaDot_di * di;
    console.log(`    For di = ${(di * 180 / Math.PI * 3600).toFixed(1)} arcsec:`);
    console.log(`      Delta RAAN drift = ${(deltaOmegaDot * 180 / Math.PI * 86400 * 1000).toFixed(3)} mdeg/day`);

    // Sensitivity of period to semi-major axis
    // T = 2*pi*sqrt(a^3/GM)
    // dT/da = 3*pi*sqrt(a/GM) = 3*T/(2*a)
    const dT_da = 1.5 * period / semiMajorAxis;
    console.log('\n  Sensitivity of orbital period to semi-major axis:');
    console.log(`    dT/da = ${(dT_da * 1000).toFixed(4)} ms/m`);
    console.log(`    For da = ${da.toFixed(0)} m: dT = ${(dT_da * da).toFixed(3)} s`);

    // Ground track repeat sensitivity
    // For exact repeat: N * T = k * T_earth (N orbits in k days)
    const T_earth = 86164.0905;  // sidereal day in seconds
    const orbitsPerDay = T_earth / period;
    console.log('\n  Ground track repeat orbit:');
    console.log(`    Orbits per sidereal day: ${orbitsPerDay.toFixed(4)}`);

    // Find nearest repeat orbit (e.g., 15 orbits per day)
    const targetOrbits = Math.round(orbitsPerDay);
    const targetPeriod = T_earth / targetOrbits;
    const targetSMA = Math.pow(earthGM * Math.pow(targetPeriod / (2 * Math.PI), 2), 1/3);

    console.log(`    For ${targetOrbits}-orbit repeat:`);
    console.log(`      Required period: ${(targetPeriod / 60).toFixed(4)} min`);
    console.log(`      Required SMA: ${(targetSMA / 1000).toFixed(3)} km`);
    console.log(`      SMA adjustment: ${((targetSMA - semiMajorAxis) / 1000).toFixed(3)} km`);

    // Monte Carlo-style uncertainty propagation (simplified)
    console.log('\n--- Monte Carlo Uncertainty Analysis ---');

    const numSamples = 1000;
    const positionErrors = [];
    const rng = () => {
        // Box-Muller transform for normal distribution
        const u1 = Math.random();
        const u2 = Math.random();
        return Math.sqrt(-2 * Math.log(u1)) * Math.cos(2 * Math.PI * u2);
    };

    // Initial position uncertainty: 10 m
    const posSigma = 10.0;

    for (let i = 0; i < numSamples; i++) {
        // Random initial position error
        const dx = rng() * posSigma;
        const dy = rng() * posSigma;
        const dz = rng() * posSigma;

        // After one orbit, along-track error grows due to period difference
        // da ~ dr (radial) leads to dT/T ~ 3/2 * da/a
        // Along-track error after one orbit: ~3*pi*da
        const dr = Math.sqrt(dx*dx + dy*dy + dz*dz);
        const alongTrackError = 3 * Math.PI * dr;

        positionErrors.push(alongTrackError);
    }

    // Compute statistics
    const mean = positionErrors.reduce((a, b) => a + b) / numSamples;
    const variance = positionErrors.reduce((a, b) => a + Math.pow(b - mean, 2), 0) / numSamples;
    const stdDev = Math.sqrt(variance);

    console.log(`  Initial position sigma: ${posSigma.toFixed(1)} m`);
    console.log(`  After 1 orbit along-track error:`);
    console.log(`    Mean: ${mean.toFixed(1)} m`);
    console.log(`    Std dev: ${stdDev.toFixed(1)} m`);
    console.log(`    Expected (3*pi*sigma): ${(3 * Math.PI * posSigma).toFixed(1)} m`);

    console.log('\n=== Sensitivity analysis complete ===');
}

main().catch(console.error);
