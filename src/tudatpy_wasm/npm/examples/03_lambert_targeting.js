/**
 * Example 03: Lambert Targeting for Interplanetary Transfers
 *
 * Ported from: examples/tudatpy/propagation/lambert_targeting_example.py
 *
 * This example demonstrates computing interplanetary transfer trajectories
 * using Lambert's problem solver (Izzo algorithm).
 *
 * Run with: node 03_lambert_targeting.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Lambert Targeting Example ===\n');

    const tudat = await createTudatModule();

    // Sun gravitational parameter
    const sunGM = 1.32712440018e20;  // m^3/s^2

    // Approximate planetary semi-major axes
    const earthSMA = 1.496e11;  // m (1 AU)
    const marsSMA = 2.279e11;   // m (~1.52 AU)

    // Earth departure position (at 0 degrees in heliocentric frame)
    const earthPos = new tudat.Vector3d();
    earthPos.set(0, earthSMA);
    earthPos.set(1, 0);
    earthPos.set(2, 0);

    // Mars arrival position (assuming ~44 degree offset for Hohmann-like transfer)
    const marsAngle = 44 * Math.PI / 180;  // Mars ahead of Earth
    const marsPos = new tudat.Vector3d();
    marsPos.set(0, marsSMA * Math.cos(Math.PI + marsAngle));
    marsPos.set(1, marsSMA * Math.sin(Math.PI + marsAngle));
    marsPos.set(2, 0);

    console.log('Earth-Mars Transfer Setup:');
    console.log(`  Earth position: [${(earthPos.get(0)/1e9).toFixed(3)}, ${(earthPos.get(1)/1e9).toFixed(3)}, 0] x10^9 m`);
    console.log(`  Mars position: [${(marsPos.get(0)/1e9).toFixed(3)}, ${(marsPos.get(1)/1e9).toFixed(3)}, 0] x10^9 m`);

    // Time of flight for Earth-Mars transfer (approximately 8-9 months)
    const timeOfFlight = 259 * 86400;  // 259 days in seconds
    console.log(`  Time of flight: ${(timeOfFlight / 86400).toFixed(0)} days`);

    // Solve Lambert's problem
    console.log('\nSolving Lambert\'s problem...');

    const lambertResult = tudat.trajectory_design.transfer_trajectory.solve_lambert_problem(
        earthPos,
        marsPos,
        timeOfFlight,
        sunGM,
        false  // prograde transfer
    );

    const departureVelocity = lambertResult.departure_velocity;
    const arrivalVelocity = lambertResult.arrival_velocity;

    console.log('\nTransfer trajectory:');
    console.log(`  Departure velocity: [${(departureVelocity.get(0)/1000).toFixed(3)}, ${(departureVelocity.get(1)/1000).toFixed(3)}, ${(departureVelocity.get(2)/1000).toFixed(3)}] km/s`);
    console.log(`  Arrival velocity: [${(arrivalVelocity.get(0)/1000).toFixed(3)}, ${(arrivalVelocity.get(1)/1000).toFixed(3)}, ${(arrivalVelocity.get(2)/1000).toFixed(3)}] km/s`);

    // Compute Earth's circular velocity
    const earthVelocity = new tudat.Vector3d();
    const earthCircularVel = Math.sqrt(sunGM / earthSMA);
    earthVelocity.set(0, 0);
    earthVelocity.set(1, earthCircularVel);
    earthVelocity.set(2, 0);

    // Compute Mars's circular velocity
    const marsVelocity = new tudat.Vector3d();
    const marsCircularVel = Math.sqrt(sunGM / marsSMA);
    marsVelocity.set(0, -marsCircularVel * Math.sin(Math.PI + marsAngle));
    marsVelocity.set(1, marsCircularVel * Math.cos(Math.PI + marsAngle));
    marsVelocity.set(2, 0);

    // Compute delta-V at departure (Earth escape)
    const departureDeltaV = Math.sqrt(
        Math.pow(departureVelocity.get(0) - earthVelocity.get(0), 2) +
        Math.pow(departureVelocity.get(1) - earthVelocity.get(1), 2) +
        Math.pow(departureVelocity.get(2) - earthVelocity.get(2), 2)
    );

    // Compute delta-V at arrival (Mars capture)
    const arrivalDeltaV = Math.sqrt(
        Math.pow(arrivalVelocity.get(0) - marsVelocity.get(0), 2) +
        Math.pow(arrivalVelocity.get(1) - marsVelocity.get(1), 2) +
        Math.pow(arrivalVelocity.get(2) - marsVelocity.get(2), 2)
    );

    const totalDeltaV = departureDeltaV + arrivalDeltaV;

    console.log('\nDelta-V budget:');
    console.log(`  Earth departure: ${(departureDeltaV/1000).toFixed(3)} km/s`);
    console.log(`  Mars arrival: ${(arrivalDeltaV/1000).toFixed(3)} km/s`);
    console.log(`  Total: ${(totalDeltaV/1000).toFixed(3)} km/s`);

    // Compute transfer orbit characteristics
    const r1 = earthSMA;
    const r2 = marsSMA;

    // Hohmann transfer for comparison
    const hohmannSMA = (r1 + r2) / 2;
    const hohmannPeriod = 2 * Math.PI * Math.sqrt(Math.pow(hohmannSMA, 3) / sunGM);
    const hohmannTOF = hohmannPeriod / 2;

    console.log('\nComparison with Hohmann transfer:');
    console.log(`  Hohmann TOF: ${(hohmannTOF / 86400).toFixed(0)} days`);
    console.log(`  Our TOF: ${(timeOfFlight / 86400).toFixed(0)} days`);

    // Hohmann delta-V calculation
    const v1_circ = Math.sqrt(sunGM / r1);
    const v1_trans = Math.sqrt(sunGM * (2/r1 - 1/hohmannSMA));
    const v2_trans = Math.sqrt(sunGM * (2/r2 - 1/hohmannSMA));
    const v2_circ = Math.sqrt(sunGM / r2);

    const hohmannDV1 = Math.abs(v1_trans - v1_circ);
    const hohmannDV2 = Math.abs(v2_circ - v2_trans);
    const hohmannTotal = hohmannDV1 + hohmannDV2;

    console.log(`  Hohmann total delta-V: ${(hohmannTotal/1000).toFixed(3)} km/s`);
    console.log(`  Our total delta-V: ${(totalDeltaV/1000).toFixed(3)} km/s`);

    // Clean up
    earthPos.delete();
    marsPos.delete();
    departureVelocity.delete();
    arrivalVelocity.delete();
    earthVelocity.delete();
    marsVelocity.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
