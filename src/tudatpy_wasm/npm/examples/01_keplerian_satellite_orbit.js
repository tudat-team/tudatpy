/**
 * Example 01: Keplerian Satellite Orbit
 *
 * Ported from: examples/tudatpy/propagation/keplerian_satellite_orbit.py
 *
 * This example demonstrates basic two-body orbit propagation using Tudat WASM.
 * A satellite is propagated around Earth using only point-mass gravity.
 *
 * Run with: node 01_keplerian_satellite_orbit.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Keplerian Satellite Orbit Example ===\n');

    // Initialize Tudat WASM module
    const tudat = await createTudatModule();

    // Physical constants
    const earthGravitationalParameter = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;  // m

    // Define initial Keplerian elements
    const semiMajorAxis = earthRadius + 400.0e3;  // 400 km altitude (ISS-like)
    const eccentricity = 0.001;
    const inclination = 51.6 * Math.PI / 180;  // ISS inclination
    const argumentOfPeriapsis = 0.0;
    const raan = 0.0;  // Right Ascension of Ascending Node
    const trueAnomaly = 0.0;

    // Create Keplerian elements vector
    const keplerianElements = new tudat.Vector6d();
    keplerianElements.set(0, semiMajorAxis);
    keplerianElements.set(1, eccentricity);
    keplerianElements.set(2, inclination);
    keplerianElements.set(3, argumentOfPeriapsis);
    keplerianElements.set(4, raan);
    keplerianElements.set(5, trueAnomaly);

    console.log('Initial Keplerian elements:');
    console.log(`  Semi-major axis: ${(semiMajorAxis / 1000).toFixed(2)} km`);
    console.log(`  Eccentricity: ${eccentricity}`);
    console.log(`  Inclination: ${(inclination * 180 / Math.PI).toFixed(2)} deg`);
    console.log(`  Arg of periapsis: ${argumentOfPeriapsis} rad`);
    console.log(`  RAAN: ${raan} rad`);
    console.log(`  True anomaly: ${trueAnomaly} rad`);

    // Convert to Cartesian state
    const initialState = tudat.astro.element_conversion.keplerian_to_cartesian(
        keplerianElements,
        earthGravitationalParameter
    );

    console.log('\nInitial Cartesian state:');
    console.log(`  Position: [${(initialState.get(0)/1000).toFixed(3)}, ${(initialState.get(1)/1000).toFixed(3)}, ${(initialState.get(2)/1000).toFixed(3)}] km`);
    console.log(`  Velocity: [${initialState.get(3).toFixed(3)}, ${initialState.get(4).toFixed(3)}, ${initialState.get(5).toFixed(3)}] m/s`);

    // Compute orbital period
    const orbitalPeriod = tudat.astro.two_body_dynamics.compute_kepler_orbit_period(
        semiMajorAxis,
        earthGravitationalParameter
    );
    console.log(`\nOrbital period: ${(orbitalPeriod / 60).toFixed(2)} minutes`);

    // Propagate for one complete orbit
    console.log('\nPropagating orbit...');

    const numSteps = 100;
    const dt = orbitalPeriod / numSteps;

    const positions = [];
    let state = initialState;

    for (let i = 0; i <= numSteps; i++) {
        const time = i * dt;

        // Propagate Kepler orbit
        const newState = tudat.astro.two_body_dynamics.propagate_kepler_orbit(
            state,
            dt,
            earthGravitationalParameter
        );

        // Store position
        positions.push({
            time: time,
            x: newState.get(0),
            y: newState.get(1),
            z: newState.get(2)
        });

        // Clean up old state and update
        if (i > 0) state.delete();
        state = newState;
    }

    // Verify orbit closure (final position should match initial)
    const finalPosition = positions[positions.length - 1];
    const initialPosition = positions[0];

    const positionError = Math.sqrt(
        Math.pow(finalPosition.x - initialPosition.x, 2) +
        Math.pow(finalPosition.y - initialPosition.y, 2) +
        Math.pow(finalPosition.z - initialPosition.z, 2)
    );

    console.log(`\nOrbit closure error: ${positionError.toFixed(6)} m`);
    console.log(`(Should be < 1 m for well-propagated Keplerian orbit)`);

    // Print some sample positions
    console.log('\nSample positions during orbit:');
    const sampleIndices = [0, 25, 50, 75, 100];
    for (const i of sampleIndices) {
        const pos = positions[i];
        const altitude = Math.sqrt(pos.x*pos.x + pos.y*pos.y + pos.z*pos.z) - earthRadius;
        console.log(`  t=${(pos.time/60).toFixed(1)} min: altitude = ${(altitude/1000).toFixed(2)} km`);
    }

    // Clean up
    keplerianElements.delete();
    initialState.delete();
    state.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
