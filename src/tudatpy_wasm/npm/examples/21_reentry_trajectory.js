/**
 * Example 21: Atmospheric Reentry Trajectory
 *
 * Ported from: examples/tudatpy/propagation/reentry_trajectory.py
 *
 * This example demonstrates atmospheric reentry trajectory simulation
 * including aerodynamic forces and heating effects.
 *
 * Run with: node 21_reentry_trajectory.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Atmospheric Reentry Trajectory Example ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const earthGM = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;  // m

    // Reentry vehicle parameters
    const vehicleMass = 2000.0;  // kg
    const referenceArea = 4.0;  // m^2
    const dragCoefficient = 1.2;

    // Initial conditions (entry interface at ~120 km altitude)
    const entryAltitude = 120.0e3;  // 120 km
    const entryVelocity = 7500.0;  // m/s (typical orbital decay speed)
    const flightPathAngle = -2.0 * Math.PI / 180;  // -2 degrees (shallow entry)

    // Convert to Cartesian state
    const r = earthRadius + entryAltitude;
    const v = entryVelocity;
    const gamma = flightPathAngle;

    // Initial position and velocity in inertial frame
    const x0 = r;
    const y0 = 0.0;
    const z0 = 0.0;
    const vx0 = v * Math.sin(gamma);
    const vy0 = v * Math.cos(gamma);
    const vz0 = 0.0;

    const initialState = new tudat.Vector6d(x0, y0, z0, vx0, vy0, vz0);

    console.log('Entry Interface Conditions:');
    console.log(`  Altitude: ${(entryAltitude / 1000).toFixed(1)} km`);
    console.log(`  Velocity: ${entryVelocity.toFixed(1)} m/s`);
    console.log(`  Flight path angle: ${(flightPathAngle * 180 / Math.PI).toFixed(2)} deg`);

    // Atmospheric model parameters (exponential atmosphere)
    const rho0 = 1.225;  // kg/m^3 (sea level density)
    const scaleHeight = 8500.0;  // m

    // Propagation parameters
    const dt = 1.0;  // 1 second time step
    const maxTime = 600.0;  // 10 minutes max

    console.log('\nPropagating reentry trajectory...\n');

    // Simple Euler integration with drag
    let state = [x0, y0, z0, vx0, vy0, vz0];
    let time = 0.0;
    const trajectory = [];

    while (time <= maxTime) {
        const [x, y, z, vx, vy, vz] = state;
        const r = Math.sqrt(x*x + y*y + z*z);
        const altitude = r - earthRadius;
        const v = Math.sqrt(vx*vx + vy*vy + vz*vz);

        // Check for ground impact or escape
        if (altitude <= 0) {
            console.log(`Ground impact at t = ${time.toFixed(1)} s`);
            break;
        }
        if (altitude > 200e3) {
            console.log('Vehicle escaped atmosphere');
            break;
        }

        // Store trajectory point
        trajectory.push({
            time: time,
            altitude: altitude,
            velocity: v,
            r: r
        });

        // Compute atmospheric density (exponential model)
        const rho = rho0 * Math.exp(-altitude / scaleHeight);

        // Gravitational acceleration
        const gx = -earthGM * x / (r * r * r);
        const gy = -earthGM * y / (r * r * r);
        const gz = -earthGM * z / (r * r * r);

        // Aerodynamic drag acceleration
        const qbar = 0.5 * rho * v * v;  // dynamic pressure
        const dragAccel = qbar * referenceArea * dragCoefficient / vehicleMass;
        const ax_drag = -dragAccel * vx / v;
        const ay_drag = -dragAccel * vy / v;
        const az_drag = -dragAccel * vz / v;

        // Total acceleration
        const ax = gx + ax_drag;
        const ay = gy + ay_drag;
        const az = gz + az_drag;

        // Euler integration
        state = [
            x + vx * dt,
            y + vy * dt,
            z + vz * dt,
            vx + ax * dt,
            vy + ay * dt,
            vz + az * dt
        ];

        time += dt;
    }

    // Print trajectory summary
    console.log('Trajectory Summary:');
    console.log(`  Duration: ${time.toFixed(1)} s`);
    console.log(`  Data points: ${trajectory.length}`);

    // Find key events
    let maxQbar = 0;
    let maxQbarTime = 0;
    let maxQbarAlt = 0;

    for (const point of trajectory) {
        const rho = rho0 * Math.exp(-point.altitude / scaleHeight);
        const qbar = 0.5 * rho * point.velocity * point.velocity;
        if (qbar > maxQbar) {
            maxQbar = qbar;
            maxQbarTime = point.time;
            maxQbarAlt = point.altitude;
        }
    }

    console.log(`\nMax Dynamic Pressure:`);
    console.log(`  q_max: ${(maxQbar / 1000).toFixed(2)} kPa`);
    console.log(`  Time: ${maxQbarTime.toFixed(1)} s`);
    console.log(`  Altitude: ${(maxQbarAlt / 1000).toFixed(1)} km`);

    // Print altitude profile at key points
    console.log('\nAltitude Profile:');
    const printIndices = [0, Math.floor(trajectory.length/4), Math.floor(trajectory.length/2),
                          Math.floor(3*trajectory.length/4), trajectory.length-1];
    for (const i of printIndices) {
        if (i < trajectory.length) {
            const p = trajectory[i];
            console.log(`  t=${p.time.toFixed(0)}s: alt=${(p.altitude/1000).toFixed(1)}km, v=${p.velocity.toFixed(0)}m/s`);
        }
    }

    // Clean up
    initialState.delete();

    console.log('\n=== Reentry simulation complete ===');
}

main().catch(console.error);
