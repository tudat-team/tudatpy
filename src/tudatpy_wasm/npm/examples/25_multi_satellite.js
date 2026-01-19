/**
 * Example 25: Multi-Satellite Propagation
 *
 * Ported from: examples/tudatpy/propagation/separation_satellites_diff_drag.py
 *
 * This example demonstrates propagating multiple satellites simultaneously
 * and analyzing their relative motion and differential drag effects.
 *
 * Run with: node 25_multi_satellite.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Multi-Satellite Propagation Example ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const earthGM = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;  // m
    const J2 = 1.08263e-3;

    // Chief satellite orbit (ISS-like)
    const altitude = 400.0e3;  // 400 km
    const a = earthRadius + altitude;
    const e = 0.0005;
    const i = 51.6 * Math.PI / 180;
    const n = Math.sqrt(earthGM / (a * a * a));
    const period = 2 * Math.PI / n;

    console.log('Chief Satellite Orbit:');
    console.log(`  Altitude: ${(altitude / 1000).toFixed(1)} km`);
    console.log(`  Eccentricity: ${e}`);
    console.log(`  Inclination: ${(i * 180 / Math.PI).toFixed(2)} deg`);
    console.log(`  Period: ${(period / 60).toFixed(2)} min`);

    // Deputy satellites with different ballistic coefficients
    // Ballistic coefficient = mass / (Cd * A)
    const satellites = [
        { name: 'Chief', bc: 50.0, dx0: 0, dy0: 0, dz0: 0, dvx0: 0, dvy0: 0, dvz0: 0 },
        { name: 'Deputy-1 (ahead)', bc: 50.0, dx0: 0, dy0: 1000, dz0: 0, dvx0: 0, dvy0: 0, dvz0: 0 },  // 1 km ahead
        { name: 'Deputy-2 (above)', bc: 50.0, dx0: 100, dy0: 0, dz0: 0, dvx0: 0, dvy0: 0, dvz0: 0 },   // 100 m higher
        { name: 'Deputy-3 (drag)', bc: 25.0, dx0: 0, dy0: 500, dz0: 0, dvx0: 0, dvy0: 0, dvz0: 0 },   // Higher drag
    ];

    console.log('\nSatellite Configuration:');
    for (const sat of satellites) {
        console.log(`  ${sat.name}: BC=${sat.bc} kg/m^2, offset=(${sat.dx0}, ${sat.dy0}, ${sat.dz0}) m`);
    }

    // Atmospheric density at 400 km (exponential model)
    const rho0 = 1.225;  // kg/m^3
    const H = 50000;  // Scale height at this altitude
    const rho = rho0 * Math.exp(-altitude / H) * 1e6;  // Approximate density at 400 km

    // Simulation parameters
    const dt = 10.0;  // 10 second time step
    const totalTime = 24 * 3600;  // 24 hours
    const numSteps = Math.floor(totalTime / dt);

    // Initialize relative states (Hill/CW frame: x=radial, y=along-track, z=cross-track)
    const states = satellites.map(sat => ({
        name: sat.name,
        bc: sat.bc,
        x: sat.dx0,
        y: sat.dy0,
        z: sat.dz0,
        vx: sat.dvx0,
        vy: sat.dvy0,
        vz: sat.dvz0
    }));

    console.log('\nPropagating relative motion with differential drag...\n');
    console.log('Time (hr) | Deputy-1 y (m) | Deputy-2 x (m) | Deputy-3 y (m)');
    console.log('----------|----------------|----------------|---------------');

    // Propagate using linearized CW equations with drag
    for (let step = 0; step <= numSteps; step++) {
        const t = step * dt;

        // Print every hour
        if (step % 360 === 0) {
            console.log(`${(t / 3600).toFixed(1).padStart(9)} | ${states[1].y.toFixed(1).padStart(14)} | ${states[2].x.toFixed(1).padStart(14)} | ${states[3].y.toFixed(1).padStart(13)}`);
        }

        // Update each satellite
        for (let i = 1; i < states.length; i++) {  // Skip chief (reference)
            const s = states[i];

            // Differential drag acceleration (relative to chief)
            // a_drag = -0.5 * rho * v^2 / BC
            const v = n * a;  // Orbital velocity
            const dragChief = 0.5 * rho * v * v / satellites[0].bc;
            const dragDeputy = 0.5 * rho * v * v / s.bc;
            const diffDrag = dragDeputy - dragChief;  // Along-track deceleration

            // CW equations with drag
            // x'' - 2n*y' - 3n^2*x = 0
            // y'' + 2n*x' = -diffDrag
            // z'' + n^2*z = 0

            // Simple Euler integration
            const ax = 3 * n * n * s.x + 2 * n * s.vy;
            const ay = -2 * n * s.vx - diffDrag;
            const az = -n * n * s.z;

            s.x += s.vx * dt;
            s.y += s.vy * dt;
            s.z += s.vz * dt;
            s.vx += ax * dt;
            s.vy += ay * dt;
            s.vz += az * dt;
        }
    }

    // Final separation analysis
    console.log('\n=== Final Relative States (after 24 hours) ===');
    console.log('\nSatellite     | x (radial) | y (along-track) | z (cross-track) | Range (m)');
    console.log('--------------|------------|-----------------|-----------------|----------');

    for (let i = 1; i < states.length; i++) {
        const s = states[i];
        const range = Math.sqrt(s.x*s.x + s.y*s.y + s.z*s.z);
        console.log(`${s.name.padEnd(13)} | ${s.x.toFixed(1).padStart(10)} | ${s.y.toFixed(1).padStart(15)} | ${s.z.toFixed(1).padStart(15)} | ${range.toFixed(1).padStart(8)}`);
    }

    // Drift rate analysis
    console.log('\n=== Drift Rate Analysis ===');

    // Deputy-1: Pure along-track offset - should stay relatively stable
    const drift1 = (states[1].y - satellites[1].dy0) / 24;  // m/hour
    console.log(`\nDeputy-1 (same BC, ahead):`);
    console.log(`  Along-track drift: ${drift1.toFixed(2)} m/hour`);
    console.log(`  Note: Minimal drift due to same ballistic coefficient`);

    // Deputy-2: Radial offset leads to along-track drift
    const drift2y = (states[2].y - satellites[2].dy0) / 24;
    console.log(`\nDeputy-2 (same BC, higher):`);
    console.log(`  Along-track drift: ${drift2y.toFixed(2)} m/hour`);
    console.log(`  Note: Higher altitude = longer period = drifts behind`);

    // Deputy-3: Different drag coefficient
    const drift3 = (states[3].y - satellites[3].dy0) / 24;
    console.log(`\nDeputy-3 (higher drag):`);
    console.log(`  Along-track drift: ${drift3.toFixed(2)} m/hour`);
    console.log(`  Note: Higher drag = loses energy = drifts ahead`);

    // Differential drag formula
    const v = n * a;
    const dragDiff = 0.5 * rho * v * v * (1/satellites[3].bc - 1/satellites[0].bc);
    const expectedDrift = 0.5 * dragDiff * 3600 * 3600;  // 1/2 * a * t^2 per hour
    console.log(`\n  Differential drag acceleration: ${(dragDiff * 1e6).toFixed(4)} micro-m/s^2`);
    console.log(`  Theoretical drift rate: ${(3 * dragDiff * 3600 / n).toFixed(2)} m/hour`);

    console.log('\n=== Multi-satellite propagation complete ===');
}

main().catch(console.error);
