/**
 * Example 33: Satellite Separation using Differential Drag
 *
 * Ported from: examples/tudatpy/propagation/separation_satellites_diff_drag.py
 *
 * This example demonstrates separation of satellites in a formation
 * using differential drag (different drag areas) for orbit maneuvering
 * without propulsion.
 *
 * Key concepts:
 * - Differential drag control
 * - Formation flying
 * - Relative motion dynamics
 * - Atmospheric density effects
 * - Along-track separation
 *
 * Run with: node 33_separation_diff_drag.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Satellite Separation using Differential Drag ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const GM_EARTH = 3.986004418e14;  // Earth gravitational parameter [m^3/s^2]
    const R_EARTH = 6.378137e6;       // Earth radius [m]
    const J2 = 1.08263e-3;            // J2 coefficient

    // Atmospheric model (exponential)
    const rho0 = 1.225;               // Sea level density [kg/m^3]
    const H = 8500;                   // Scale height [m]

    // Satellite parameters
    const satellite1 = {
        name: 'Sat-A (High Drag)',
        mass: 100,              // kg
        area: 2.0,              // m^2 (larger area - high drag)
        Cd: 2.2                 // Drag coefficient
    };

    const satellite2 = {
        name: 'Sat-B (Low Drag)',
        mass: 100,              // kg
        area: 0.5,              // m^2 (smaller area - low drag)
        Cd: 2.2                 // Drag coefficient
    };

    // Orbit parameters (LEO with significant drag)
    const altitude = 350e3;     // 350 km - significant drag environment
    const inclination = 51.6 * Math.PI / 180;  // ISS-like inclination

    const r0 = R_EARTH + altitude;
    const v0 = Math.sqrt(GM_EARTH / r0);
    const orbitalPeriod = 2 * Math.PI * Math.sqrt(r0 * r0 * r0 / GM_EARTH);

    console.log('Orbit Configuration:');
    console.log(`  Altitude: ${altitude / 1000} km`);
    console.log(`  Inclination: ${(inclination * 180 / Math.PI).toFixed(1)} deg`);
    console.log(`  Orbital velocity: ${(v0 / 1000).toFixed(3)} km/s`);
    console.log(`  Orbital period: ${(orbitalPeriod / 60).toFixed(1)} minutes`);

    console.log('\nSatellite Configuration:');
    console.log(`  ${satellite1.name}: mass=${satellite1.mass} kg, area=${satellite1.area} m^2`);
    console.log(`  ${satellite2.name}: mass=${satellite2.mass} kg, area=${satellite2.area} m^2`);

    // Ballistic coefficients (higher = less drag effect)
    const BC1 = satellite1.mass / (satellite1.Cd * satellite1.area);
    const BC2 = satellite2.mass / (satellite2.Cd * satellite2.area);
    console.log(`\nBallistic Coefficients:`);
    console.log(`  ${satellite1.name}: ${BC1.toFixed(1)} kg/m^2`);
    console.log(`  ${satellite2.name}: ${BC2.toFixed(1)} kg/m^2`);

    // Initial states (both satellites start together)
    // Using orbital elements, then convert to Cartesian
    const RAAN = 0;
    const argPeriapsis = 0;
    const trueAnomaly = 0;

    function orbitalToCartesian(a, e, i, RAAN, omega, nu) {
        // Position in orbital plane
        const p = a * (1 - e * e);
        const r = p / (1 + e * Math.cos(nu));

        const xOrb = r * Math.cos(nu);
        const yOrb = r * Math.sin(nu);

        // Velocity in orbital plane
        const h = Math.sqrt(GM_EARTH * p);
        const vxOrb = -GM_EARTH / h * Math.sin(nu);
        const vyOrb = GM_EARTH / h * (e + Math.cos(nu));

        // Rotation matrices
        const cosRAAN = Math.cos(RAAN), sinRAAN = Math.sin(RAAN);
        const cosI = Math.cos(i), sinI = Math.sin(i);
        const cosOmega = Math.cos(omega), sinOmega = Math.sin(omega);

        // Combined rotation
        const Px = cosRAAN * cosOmega - sinRAAN * sinOmega * cosI;
        const Py = sinRAAN * cosOmega + cosRAAN * sinOmega * cosI;
        const Pz = sinOmega * sinI;
        const Qx = -cosRAAN * sinOmega - sinRAAN * cosOmega * cosI;
        const Qy = -sinRAAN * sinOmega + cosRAAN * cosOmega * cosI;
        const Qz = cosOmega * sinI;

        return [
            Px * xOrb + Qx * yOrb,
            Py * xOrb + Qy * yOrb,
            Pz * xOrb + Qz * yOrb,
            Px * vxOrb + Qx * vyOrb,
            Py * vxOrb + Qy * vyOrb,
            Pz * vxOrb + Qz * vyOrb
        ];
    }

    // Initial states (circular orbit)
    let state1 = orbitalToCartesian(r0, 0, inclination, RAAN, argPeriapsis, trueAnomaly);
    let state2 = orbitalToCartesian(r0, 0, inclination, RAAN, argPeriapsis, trueAnomaly);

    // Atmospheric density at altitude
    function atmosphericDensity(alt) {
        if (alt < 0) alt = 0;
        // More realistic density model for LEO
        if (alt < 150e3) {
            return rho0 * Math.exp(-alt / H);
        } else {
            // Adjusted for higher altitudes
            const rho150 = rho0 * Math.exp(-150e3 / H);
            const H_high = 50000;  // Larger scale height at high altitude
            return rho150 * Math.exp(-(alt - 150e3) / H_high);
        }
    }

    // Compute acceleration with J2 and drag
    function computeAcceleration(state, satellite) {
        const pos = [state[0], state[1], state[2]];
        const vel = [state[3], state[4], state[5]];

        const r = Math.sqrt(pos[0]*pos[0] + pos[1]*pos[1] + pos[2]*pos[2]);
        const alt = r - R_EARTH;

        // Two-body acceleration
        const r3 = r * r * r;
        const aGrav = [
            -GM_EARTH * pos[0] / r3,
            -GM_EARTH * pos[1] / r3,
            -GM_EARTH * pos[2] / r3
        ];

        // J2 perturbation
        const z2 = pos[2] * pos[2];
        const r2 = r * r;
        const factor = 1.5 * J2 * GM_EARTH * R_EARTH * R_EARTH / (r2 * r3);
        const aJ2 = [
            factor * pos[0] * (5 * z2 / r2 - 1),
            factor * pos[1] * (5 * z2 / r2 - 1),
            factor * pos[2] * (5 * z2 / r2 - 3)
        ];

        // Drag acceleration
        const rho = atmosphericDensity(alt);
        const vRel = vel;  // Assuming no atmospheric rotation for simplicity
        const v = Math.sqrt(vRel[0]*vRel[0] + vRel[1]*vRel[1] + vRel[2]*vRel[2]);

        const dragFactor = -0.5 * rho * satellite.Cd * satellite.area / satellite.mass * v;
        const aDrag = [
            dragFactor * vRel[0],
            dragFactor * vRel[1],
            dragFactor * vRel[2]
        ];

        return [
            aGrav[0] + aJ2[0] + aDrag[0],
            aGrav[1] + aJ2[1] + aDrag[1],
            aGrav[2] + aJ2[2] + aDrag[2]
        ];
    }

    // Derivatives for a satellite
    function derivatives(state, satellite) {
        const acc = computeAcceleration(state, satellite);
        return [state[3], state[4], state[5], acc[0], acc[1], acc[2]];
    }

    // RK4 integrator
    function rk4Step(state, satellite, dt) {
        const k1 = derivatives(state, satellite);
        const s1 = state.map((v, i) => v + 0.5 * dt * k1[i]);
        const k2 = derivatives(s1, satellite);
        const s2 = state.map((v, i) => v + 0.5 * dt * k2[i]);
        const k3 = derivatives(s2, satellite);
        const s3 = state.map((v, i) => v + dt * k3[i]);
        const k4 = derivatives(s3, satellite);
        return state.map((v, i) => v + (dt / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    }

    // Compute along-track separation
    function computeSeparation(state1, state2) {
        // Relative position
        const relPos = [
            state2[0] - state1[0],
            state2[1] - state1[1],
            state2[2] - state1[2]
        ];

        // Total separation distance
        const totalSep = Math.sqrt(relPos[0]*relPos[0] + relPos[1]*relPos[1] + relPos[2]*relPos[2]);

        // Along-track direction (velocity direction of sat1)
        const v1 = Math.sqrt(state1[3]*state1[3] + state1[4]*state1[4] + state1[5]*state1[5]);
        const alongTrack = [state1[3]/v1, state1[4]/v1, state1[5]/v1];

        // Along-track separation (dot product)
        const alongTrackSep = relPos[0]*alongTrack[0] + relPos[1]*alongTrack[1] + relPos[2]*alongTrack[2];

        // Radial direction
        const r1 = Math.sqrt(state1[0]*state1[0] + state1[1]*state1[1] + state1[2]*state1[2]);
        const radial = [state1[0]/r1, state1[1]/r1, state1[2]/r1];
        const radialSep = relPos[0]*radial[0] + relPos[1]*radial[1] + relPos[2]*radial[2];

        // Cross-track (perpendicular to both)
        const crossTrack = [
            alongTrack[1]*radial[2] - alongTrack[2]*radial[1],
            alongTrack[2]*radial[0] - alongTrack[0]*radial[2],
            alongTrack[0]*radial[1] - alongTrack[1]*radial[0]
        ];
        const crossTrackSep = relPos[0]*crossTrack[0] + relPos[1]*crossTrack[1] + relPos[2]*crossTrack[2];

        return { total: totalSep, alongTrack: alongTrackSep, radial: radialSep, crossTrack: crossTrackSep };
    }

    // Simulation
    console.log('\n--- Running Separation Simulation ---\n');

    const dt = 10;  // Time step [s]
    const numDays = 7;  // Simulation duration
    const totalTime = numDays * 86400;

    const separationHistory = [];
    const altitudeHistory = [];
    let t = 0;

    while (t <= totalTime) {
        // Compute separation
        const sep = computeSeparation(state1, state2);

        // Compute altitudes
        const r1 = Math.sqrt(state1[0]*state1[0] + state1[1]*state1[1] + state1[2]*state1[2]);
        const r2 = Math.sqrt(state2[0]*state2[0] + state2[1]*state2[1] + state2[2]*state2[2]);
        const alt1 = r1 - R_EARTH;
        const alt2 = r2 - R_EARTH;

        // Store periodically
        if (Math.floor(t) % 3600 === 0 || separationHistory.length === 0) {
            separationHistory.push({
                t: t / 86400,  // days
                total: sep.total / 1000,  // km
                alongTrack: sep.alongTrack / 1000,
                radial: sep.radial,  // m
                crossTrack: sep.crossTrack  // m
            });

            altitudeHistory.push({
                t: t / 86400,
                alt1: alt1 / 1000,
                alt2: alt2 / 1000
            });
        }

        // Integrate both satellites
        state1 = rk4Step(state1, satellite1, dt);
        state2 = rk4Step(state2, satellite2, dt);
        t += dt;
    }

    // Results
    console.log(`Simulation completed: ${numDays} days`);
    console.log(`Data points: ${separationHistory.length}`);

    // Final separation
    const finalSep = separationHistory[separationHistory.length - 1];
    console.log(`\n=== Final Separation (after ${numDays} days) ===`);
    console.log(`  Total: ${finalSep.total.toFixed(2)} km`);
    console.log(`  Along-track: ${finalSep.alongTrack.toFixed(2)} km`);
    console.log(`  Radial: ${finalSep.radial.toFixed(1)} m`);
    console.log(`  Cross-track: ${finalSep.crossTrack.toFixed(1)} m`);

    // Altitude decay
    const initialAlt = altitudeHistory[0];
    const finalAlt = altitudeHistory[altitudeHistory.length - 1];
    console.log(`\n=== Altitude Decay ===`);
    console.log(`  ${satellite1.name}:`);
    console.log(`    Initial: ${initialAlt.alt1.toFixed(2)} km`);
    console.log(`    Final: ${finalAlt.alt1.toFixed(2)} km`);
    console.log(`    Decay: ${(initialAlt.alt1 - finalAlt.alt1).toFixed(2)} km`);
    console.log(`  ${satellite2.name}:`);
    console.log(`    Initial: ${initialAlt.alt2.toFixed(2)} km`);
    console.log(`    Final: ${finalAlt.alt2.toFixed(2)} km`);
    console.log(`    Decay: ${(initialAlt.alt2 - finalAlt.alt2).toFixed(2)} km`);

    // Separation rate
    const sepRate = (finalSep.alongTrack - separationHistory[0].alongTrack) / numDays;
    console.log(`\n=== Separation Rate ===`);
    console.log(`  Along-track rate: ${sepRate.toFixed(2)} km/day`);
    console.log(`  Time to 100 km separation: ${(100 / Math.abs(sepRate)).toFixed(1)} days`);

    // Differential drag analysis
    console.log('\n--- Differential Drag Analysis ---');

    const rho_avg = atmosphericDensity(altitude);
    const dragAccel1 = 0.5 * rho_avg * v0 * v0 * satellite1.Cd * satellite1.area / satellite1.mass;
    const dragAccel2 = 0.5 * rho_avg * v0 * v0 * satellite2.Cd * satellite2.area / satellite2.mass;
    const diffDrag = dragAccel1 - dragAccel2;

    console.log(`At ${altitude / 1000} km altitude:`);
    console.log(`  Atmospheric density: ${rho_avg.toExponential(3)} kg/m^3`);
    console.log(`  Drag acceleration (${satellite1.name}): ${(dragAccel1 * 1e6).toFixed(4)} mm/s^2`);
    console.log(`  Drag acceleration (${satellite2.name}): ${(dragAccel2 * 1e6).toFixed(4)} mm/s^2`);
    console.log(`  Differential drag: ${(diffDrag * 1e6).toFixed(4)} mm/s^2`);

    // Theoretical separation estimate
    const theoreticalSep = 0.5 * diffDrag * totalTime * totalTime;
    console.log(`\nTheoretical separation (constant drag): ${(theoreticalSep / 1000).toFixed(2)} km`);

    // Sample trajectory
    console.log('\n--- Separation Over Time ---');
    console.log('Day | Along-track [km] | Altitude Diff [m]');
    console.log('----|------------------|------------------');
    for (let i = 0; i < separationHistory.length; i += Math.max(1, Math.floor(separationHistory.length / 8))) {
        const sep = separationHistory[i];
        const alt = altitudeHistory[i];
        const altDiff = (alt.alt2 - alt.alt1) * 1000;  // m
        console.log(`${sep.t.toFixed(1).padStart(3)} | ${sep.alongTrack.toFixed(2).padStart(16)} | ${altDiff.toFixed(1).padStart(17)}`);
    }

    // Mission design implications
    console.log('\n--- Mission Design Implications ---');
    console.log('Differential drag can be used for:');
    console.log('  - Formation initialization (spreading satellites)');
    console.log('  - Station-keeping adjustments');
    console.log('  - Fuel-free orbit maneuvering');
    console.log('  - Constellation phasing');
    console.log('\nLimitations:');
    console.log('  - Only works in atmospheres (< ~600 km for Earth)');
    console.log('  - Causes altitude decay');
    console.log('  - Slow compared to propulsive maneuvers');
    console.log('  - Requires attitude control for area changes');

    console.log('\n=== Differential Drag Separation Complete ===');
}

main().catch(console.error);
