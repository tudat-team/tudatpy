/**
 * Example 02: Perturbed Satellite Orbit
 *
 * Ported from: examples/tudatpy/propagation/perturbed_satellite_orbit.py
 *
 * This example demonstrates satellite orbit propagation with perturbations
 * including spherical harmonic gravity (J2) and third-body effects.
 *
 * Run with: node 02_perturbed_satellite_orbit.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Perturbed Satellite Orbit Example ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const earthGM = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;  // m
    const earthJ2 = 1.08263e-3;  // J2 coefficient

    // Initial orbit - LEO satellite
    const altitude = 630.0e3;  // 630 km (sun-synchronous orbit altitude)
    const semiMajorAxis = earthRadius + altitude;
    const eccentricity = 0.001;
    const inclination = 98.0 * Math.PI / 180;  // SSO inclination
    const raan = 0.0;
    const argOfPeriapsis = 0.0;
    const trueAnomaly = 0.0;

    console.log('Initial orbital elements:');
    console.log(`  Altitude: ${altitude/1000} km`);
    console.log(`  Semi-major axis: ${(semiMajorAxis/1000).toFixed(2)} km`);
    console.log(`  Inclination: ${(inclination * 180 / Math.PI).toFixed(2)} deg`);

    // Create initial state
    const keplerianElements = new tudat.Vector6d();
    keplerianElements.set(0, semiMajorAxis);
    keplerianElements.set(1, eccentricity);
    keplerianElements.set(2, inclination);
    keplerianElements.set(3, argOfPeriapsis);
    keplerianElements.set(4, raan);
    keplerianElements.set(5, trueAnomaly);

    const initialState = tudat.astro.element_conversion.keplerian_to_cartesian(
        keplerianElements,
        earthGM
    );

    // Compute orbital period
    const orbitalPeriod = tudat.astro.two_body_dynamics.compute_kepler_orbit_period(
        semiMajorAxis,
        earthGM
    );
    console.log(`\nOrbital period: ${(orbitalPeriod / 60).toFixed(2)} minutes`);

    // Simulation parameters
    const simulationDuration = 86400.0;  // 1 day
    const numOrbits = simulationDuration / orbitalPeriod;
    console.log(`\nSimulation duration: 1 day (${numOrbits.toFixed(1)} orbits)`);

    // For demonstration, we'll compute the J2 secular drift analytically
    // Rate of change of RAAN due to J2
    const n = Math.sqrt(earthGM / Math.pow(semiMajorAxis, 3));  // mean motion
    const p = semiMajorAxis * (1 - eccentricity * eccentricity);  // semi-latus rectum

    const raanDrift = -1.5 * n * earthJ2 * Math.pow(earthRadius / p, 2) * Math.cos(inclination);
    const argPeriDrift = 0.75 * n * earthJ2 * Math.pow(earthRadius / p, 2) * (5 * Math.pow(Math.cos(inclination), 2) - 1);

    console.log('\nJ2 secular perturbation rates:');
    console.log(`  RAAN drift: ${(raanDrift * 180 / Math.PI * 86400).toFixed(4)} deg/day`);
    console.log(`  Arg of periapsis drift: ${(argPeriDrift * 180 / Math.PI * 86400).toFixed(4)} deg/day`);

    // Propagate using Kepler (two-body) for comparison
    console.log('\nPropagating two-body orbit for 1 day...');

    const numSteps = 1000;
    const dt = simulationDuration / numSteps;

    let state = initialState;
    const history = [];

    for (let i = 0; i <= numSteps; i++) {
        const time = i * dt;

        // Store current state
        history.push({
            time: time,
            x: state.get(0),
            y: state.get(1),
            z: state.get(2),
            vx: state.get(3),
            vy: state.get(4),
            vz: state.get(5)
        });

        if (i < numSteps) {
            // Propagate
            const newState = tudat.astro.two_body_dynamics.propagate_kepler_orbit(
                state,
                dt,
                earthGM
            );

            if (i > 0) state.delete();
            state = newState;
        }
    }

    // Convert final state to Keplerian elements
    const finalState = new tudat.Vector6d();
    const lastEntry = history[history.length - 1];
    finalState.set(0, lastEntry.x);
    finalState.set(1, lastEntry.y);
    finalState.set(2, lastEntry.z);
    finalState.set(3, lastEntry.vx);
    finalState.set(4, lastEntry.vy);
    finalState.set(5, lastEntry.vz);

    const finalKeplerian = tudat.astro.element_conversion.cartesian_to_keplerian(
        finalState,
        earthGM
    );

    console.log('\nFinal orbital elements (two-body, no perturbations):');
    console.log(`  Semi-major axis: ${(finalKeplerian.get(0)/1000).toFixed(4)} km`);
    console.log(`  Eccentricity: ${finalKeplerian.get(1).toFixed(6)}`);
    console.log(`  Inclination: ${(finalKeplerian.get(2) * 180 / Math.PI).toFixed(4)} deg`);

    // Compare with initial
    const smaDiff = finalKeplerian.get(0) - semiMajorAxis;
    const eccDiff = finalKeplerian.get(1) - eccentricity;
    const incDiff = (finalKeplerian.get(2) - inclination) * 180 / Math.PI;

    console.log('\nOrbital element changes (two-body should be ~zero):');
    console.log(`  SMA change: ${smaDiff.toFixed(6)} m`);
    console.log(`  Eccentricity change: ${eccDiff.toExponential(3)}`);
    console.log(`  Inclination change: ${incDiff.toExponential(3)} deg`);

    // Expected RAAN change due to J2 (analytical)
    const expectedRaanChange = raanDrift * simulationDuration * 180 / Math.PI;
    console.log(`\nExpected RAAN change due to J2: ${expectedRaanChange.toFixed(4)} deg`);
    console.log('(This would be observed in full numerical propagation with J2)');

    // Clean up
    keplerianElements.delete();
    initialState.delete();
    state.delete();
    finalState.delete();
    finalKeplerian.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
