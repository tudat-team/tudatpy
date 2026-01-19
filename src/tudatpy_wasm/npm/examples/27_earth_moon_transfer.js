/**
 * Example 27: Thrust-Based Earth-Moon Transfer
 *
 * Ported from: examples/tudatpy/propagation/thrust_between_Earth_Moon.py
 *
 * This example demonstrates continuous thrust propagation in the Earth-Moon
 * system with mass depletion and multi-body gravitational perturbations.
 *
 * Key concepts:
 * - Continuous low-thrust propulsion
 * - Multi-body gravitational environment (Earth, Moon, Sun)
 * - Mass propagation coupled with thrust
 * - Hybrid termination conditions
 * - Dependent variables (altitude, mass)
 *
 * Run with: node 27_earth_moon_transfer.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Earth-Moon Transfer with Continuous Thrust ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const GM_EARTH = 3.986004418e14;  // Earth gravitational parameter [m^3/s^2]
    const GM_MOON = 4.9028e12;        // Moon gravitational parameter [m^3/s^2]
    const GM_SUN = 1.32712440018e20;  // Sun gravitational parameter [m^3/s^2]
    const R_EARTH = 6.378137e6;       // Earth radius [m]
    const MOON_DISTANCE = 3.844e8;    // Average Earth-Moon distance [m]
    const JULIAN_DAY = 86400;         // Seconds per day

    // Vehicle parameters
    const initialMass = 5000;           // kg
    const thrustMagnitude = 10.0;       // N
    const specificImpulse = 5000;       // s
    const g0 = 9.80665;                 // Standard gravity [m/s^2]

    // Compute mass flow rate: mdot = T / (Isp * g0)
    const massFlowRate = thrustMagnitude / (specificImpulse * g0);

    console.log('Vehicle Configuration:');
    console.log(`  Initial mass: ${initialMass} kg`);
    console.log(`  Thrust: ${thrustMagnitude} N`);
    console.log(`  Specific impulse: ${specificImpulse} s`);
    console.log(`  Mass flow rate: ${(massFlowRate * 1000).toFixed(4)} g/s`);
    console.log(`  Initial acceleration: ${(thrustMagnitude / initialMass * 1000).toFixed(4)} mm/s^2`);

    // Initial state: 8000 km altitude, circular velocity
    const r0 = 8.0e6;  // Initial distance from Earth center [m]
    const v0 = 7.5e3;  // Initial velocity [m/s] (slightly above circular)

    // Circular velocity at this altitude
    const vCircular = Math.sqrt(GM_EARTH / r0);
    console.log(`\nInitial Orbit:`);
    console.log(`  Altitude: ${((r0 - R_EARTH) / 1000).toFixed(0)} km`);
    console.log(`  Velocity: ${(v0 / 1000).toFixed(3)} km/s`);
    console.log(`  Circular velocity: ${(vCircular / 1000).toFixed(3)} km/s`);
    console.log(`  Energy parameter: ${v0 > vCircular ? 'hyperbolic-ish' : 'elliptical'}`);

    // Simulation parameters
    const simulationDays = 30;
    const dt = 60;  // Time step [s]
    const totalSteps = Math.floor(simulationDays * JULIAN_DAY / dt);

    // State vector: [x, y, z, vx, vy, vz, mass]
    let state = [r0, 0, 0, 0, v0, 0, initialMass];

    // Termination conditions
    const maxAltitude = 100e6;    // 100,000 km
    const minMass = 4000;         // 4000 kg (burned 1000 kg)

    console.log('\nTermination Conditions:');
    console.log(`  Max altitude: ${(maxAltitude / 1e6).toFixed(0)} million km`);
    console.log(`  Min mass: ${minMass} kg (max burn: ${initialMass - minMass} kg)`);
    console.log(`  Max time: ${simulationDays} days`);

    // Simple Moon position model (circular orbit in XY plane)
    const moonPeriod = 27.3 * JULIAN_DAY;  // Sidereal month
    const moonOmega = 2 * Math.PI / moonPeriod;

    function getMoonPosition(t) {
        const angle = moonOmega * t;
        return [
            MOON_DISTANCE * Math.cos(angle),
            MOON_DISTANCE * Math.sin(angle),
            0
        ];
    }

    // Simplified Sun position (very distant, nearly fixed direction)
    const SUN_DISTANCE = 1.496e11;  // 1 AU
    function getSunPosition(t) {
        // Sun moves slowly, approximate as fixed for 30 days
        const angle = 2 * Math.PI * t / (365.25 * JULIAN_DAY);
        return [
            SUN_DISTANCE * Math.cos(angle),
            SUN_DISTANCE * Math.sin(angle),
            0
        ];
    }

    // Compute gravitational acceleration from a body
    function gravityAccel(pos, bodyPos, bodyGM) {
        const dx = bodyPos[0] - pos[0];
        const dy = bodyPos[1] - pos[1];
        const dz = bodyPos[2] - pos[2];
        const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
        const r3 = r * r * r;
        return [bodyGM * dx / r3, bodyGM * dy / r3, bodyGM * dz / r3];
    }

    // Compute thrust acceleration (aligned with velocity)
    function thrustAccel(vel, mass) {
        const v = Math.sqrt(vel[0]*vel[0] + vel[1]*vel[1] + vel[2]*vel[2]);
        if (v < 1e-10) return [0, 0, 0];
        const accel = thrustMagnitude / mass;
        return [
            accel * vel[0] / v,
            accel * vel[1] / v,
            accel * vel[2] / v
        ];
    }

    // RK4 integration for state + mass
    function derivatives(t, s) {
        const pos = [s[0], s[1], s[2]];
        const vel = [s[3], s[4], s[5]];
        const mass = s[6];

        // Gravitational accelerations
        const aEarth = gravityAccel(pos, [0, 0, 0], GM_EARTH);
        const moonPos = getMoonPosition(t);
        const aMoon = gravityAccel(pos, moonPos, GM_MOON);
        const sunPos = getSunPosition(t);
        const aSun = gravityAccel(pos, sunPos, GM_SUN);

        // Thrust acceleration
        const aThrust = thrustAccel(vel, mass);

        // Total acceleration
        const ax = aEarth[0] + aMoon[0] + aSun[0] + aThrust[0];
        const ay = aEarth[1] + aMoon[1] + aSun[1] + aThrust[1];
        const az = aEarth[2] + aMoon[2] + aSun[2] + aThrust[2];

        // Mass rate (negative because mass decreases)
        const mdot = -massFlowRate;

        return [vel[0], vel[1], vel[2], ax, ay, az, mdot];
    }

    function rk4Step(t, s, h) {
        const k1 = derivatives(t, s);
        const s1 = s.map((v, i) => v + 0.5 * h * k1[i]);
        const k2 = derivatives(t + 0.5 * h, s1);
        const s2 = s.map((v, i) => v + 0.5 * h * k2[i]);
        const k3 = derivatives(t + 0.5 * h, s2);
        const s3 = s.map((v, i) => v + h * k3[i]);
        const k4 = derivatives(t + h, s3);
        return s.map((v, i) => v + (h / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    }

    // Storage for trajectory data
    const trajectory = [];
    const altitudeHistory = [];
    const massHistory = [];

    console.log('\n--- Running Propagation ---\n');

    let t = 0;
    let terminationReason = 'max_time';

    for (let step = 0; step <= totalSteps; step++) {
        // Compute current altitude
        const r = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
        const altitude = r - R_EARTH;

        // Store data periodically
        if (step % 100 === 0) {
            trajectory.push({
                t: t,
                x: state[0],
                y: state[1],
                z: state[2],
                vx: state[3],
                vy: state[4],
                vz: state[5]
            });
            altitudeHistory.push({ t: t / JULIAN_DAY, altitude: altitude / 1000 });
            massHistory.push({ t: t / JULIAN_DAY, mass: state[6] });
        }

        // Check termination conditions
        if (altitude > maxAltitude) {
            terminationReason = 'max_altitude';
            console.log(`Termination: Reached altitude ${(altitude / 1e6).toFixed(2)} million km`);
            break;
        }
        if (state[6] < minMass) {
            terminationReason = 'min_mass';
            console.log(`Termination: Mass dropped to ${state[6].toFixed(1)} kg`);
            break;
        }

        // Integrate one step
        state = rk4Step(t, state, dt);
        t += dt;
    }

    if (terminationReason === 'max_time') {
        console.log(`Termination: Reached max simulation time (${simulationDays} days)`);
    }

    // Final state summary
    const finalR = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
    const finalV = Math.sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
    const finalAltitude = finalR - R_EARTH;
    const propellantUsed = initialMass - state[6];

    console.log('\n=== Final State ===');
    console.log(`Time elapsed: ${(t / JULIAN_DAY).toFixed(2)} days`);
    console.log(`Final altitude: ${(finalAltitude / 1000).toFixed(0)} km`);
    console.log(`Final velocity: ${(finalV / 1000).toFixed(3)} km/s`);
    console.log(`Final mass: ${state[6].toFixed(1)} kg`);
    console.log(`Propellant used: ${propellantUsed.toFixed(1)} kg`);

    // Orbital energy analysis
    const specificEnergy = 0.5 * finalV * finalV - GM_EARTH / finalR;
    if (specificEnergy < 0) {
        const sma = -GM_EARTH / (2 * specificEnergy);
        console.log(`Orbit type: Elliptical`);
        console.log(`Semi-major axis: ${(sma / 1000).toFixed(0)} km`);
    } else {
        console.log(`Orbit type: Hyperbolic (escaping Earth)`);
        console.log(`Excess velocity: ${(Math.sqrt(2 * specificEnergy) / 1000).toFixed(3)} km/s`);
    }

    // Show Moon distance at final time
    const moonPos = getMoonPosition(t);
    const moonDist = Math.sqrt(
        Math.pow(state[0] - moonPos[0], 2) +
        Math.pow(state[1] - moonPos[1], 2) +
        Math.pow(state[2] - moonPos[2], 2)
    );
    console.log(`Distance to Moon: ${(moonDist / 1000).toFixed(0)} km`);

    // Print trajectory summary
    console.log('\n--- Trajectory Summary ---');
    console.log(`Total data points: ${trajectory.length}`);

    // Find min and max altitudes
    const minAlt = Math.min(...altitudeHistory.map(h => h.altitude));
    const maxAlt = Math.max(...altitudeHistory.map(h => h.altitude));
    console.log(`Altitude range: ${minAlt.toFixed(0)} - ${maxAlt.toFixed(0)} km`);

    // Delta-V calculation
    const deltaV = specificImpulse * g0 * Math.log(initialMass / state[6]);
    console.log(`Total delta-V applied: ${(deltaV / 1000).toFixed(3)} km/s`);

    // Comparison with Hohmann transfer
    const hohmannDV1 = Math.sqrt(GM_EARTH / r0) * (Math.sqrt(2 * MOON_DISTANCE / (r0 + MOON_DISTANCE)) - 1);
    const hohmannDV2 = Math.sqrt(GM_EARTH / MOON_DISTANCE) * (1 - Math.sqrt(2 * r0 / (r0 + MOON_DISTANCE)));
    const hohmannTotal = Math.abs(hohmannDV1) + Math.abs(hohmannDV2);

    console.log('\n--- Comparison with Hohmann Transfer ---');
    console.log(`Hohmann transfer delta-V: ${(hohmannTotal / 1000).toFixed(3)} km/s`);
    console.log(`Low-thrust delta-V used: ${(deltaV / 1000).toFixed(3)} km/s`);
    console.log(`Difference: ${((deltaV - hohmannTotal) / 1000).toFixed(3)} km/s`);
    console.log('Note: Low-thrust requires more delta-V but can use efficient propulsion');

    console.log('\n=== Earth-Moon Transfer Complete ===');
}

main().catch(console.error);
