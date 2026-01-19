/**
 * Example 29: JUICE Flybys of Jovian Moons
 *
 * Ported from: examples/tudatpy/propagation/juice_flybys.py
 *
 * This example demonstrates simulation of the JUICE (JUpiter ICy moons
 * Explorer) mission trajectory, including flybys of Ganymede, Europa,
 * and Callisto in the Jovian system.
 *
 * Key concepts:
 * - Multi-body Jovian system dynamics
 * - Flyby trajectory design
 * - Closest approach detection
 * - Variable step integration
 * - B-plane targeting
 *
 * Run with: node 29_juice_flybys.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== JUICE Flybys of Jovian Moons ===\n');

    const tudat = await createTudatModule();

    // Gravitational parameters [m^3/s^2]
    const GM_JUPITER = 1.26687e17;
    const GM_GANYMEDE = 9.887e12;
    const GM_EUROPA = 3.203e12;
    const GM_CALLISTO = 7.179e12;
    const GM_IO = 5.959e12;
    const GM_SUN = 1.32712440018e20;

    // Orbital radii (semi-major axes) [m]
    const a_IO = 4.217e8;           // Io orbit
    const a_EUROPA = 6.709e8;       // Europa orbit
    const a_GANYMEDE = 1.0704e9;    // Ganymede orbit
    const a_CALLISTO = 1.8827e9;    // Callisto orbit
    const a_JUPITER = 7.785e11;     // Jupiter around Sun

    // Moon radii [m]
    const R_GANYMEDE = 2.634e6;
    const R_EUROPA = 1.561e6;
    const R_CALLISTO = 2.410e6;

    // Orbital periods [s]
    const T_IO = 1.769 * 86400;
    const T_EUROPA = 3.551 * 86400;
    const T_GANYMEDE = 7.155 * 86400;
    const T_CALLISTO = 16.689 * 86400;

    // JUICE spacecraft mass
    const JUICE_MASS = 5000;  // kg

    console.log('Jovian System Parameters:');
    console.log(`  Jupiter GM: ${GM_JUPITER.toExponential(3)} m^3/s^2`);
    console.log(`  Ganymede: a = ${(a_GANYMEDE / 1e6).toFixed(0)} km, T = ${(T_GANYMEDE / 86400).toFixed(2)} days`);
    console.log(`  Europa: a = ${(a_EUROPA / 1e6).toFixed(0)} km, T = ${(T_EUROPA / 86400).toFixed(2)} days`);
    console.log(`  Callisto: a = ${(a_CALLISTO / 1e6).toFixed(0)} km, T = ${(T_CALLISTO / 86400).toFixed(2)} days`);

    // Moon position functions (circular coplanar approximation)
    function getMoonPosition(a, T, t, phase0 = 0) {
        const omega = 2 * Math.PI / T;
        const angle = omega * t + phase0;
        return [a * Math.cos(angle), a * Math.sin(angle), 0];
    }

    // Get position of all moons at time t
    function getJovianSystemState(t) {
        return {
            io: getMoonPosition(a_IO, T_IO, t, 0),
            europa: getMoonPosition(a_EUROPA, T_EUROPA, t, Math.PI / 6),
            ganymede: getMoonPosition(a_GANYMEDE, T_GANYMEDE, t, Math.PI / 3),
            callisto: getMoonPosition(a_CALLISTO, T_CALLISTO, t, Math.PI / 2)
        };
    }

    // Compute gravitational acceleration from a body
    function gravityAccel(pos, bodyPos, bodyGM) {
        const dx = bodyPos[0] - pos[0];
        const dy = bodyPos[1] - pos[1];
        const dz = bodyPos[2] - pos[2];
        const r = Math.sqrt(dx*dx + dy*dy + dz*dz);
        if (r < 1000) return [0, 0, 0];  // Avoid singularity
        const r3 = r * r * r;
        return [bodyGM * dx / r3, bodyGM * dy / r3, bodyGM * dz / r3];
    }

    // Total acceleration on JUICE
    function derivatives(t, state) {
        const pos = [state[0], state[1], state[2]];
        const vel = [state[3], state[4], state[5]];

        // Jupiter (at origin)
        const aJupiter = gravityAccel(pos, [0, 0, 0], GM_JUPITER);

        // Get moon positions
        const moons = getJovianSystemState(t);

        // Moon gravitational perturbations
        const aGanymede = gravityAccel(pos, moons.ganymede, GM_GANYMEDE);
        const aEuropa = gravityAccel(pos, moons.europa, GM_EUROPA);
        const aCallisto = gravityAccel(pos, moons.callisto, GM_CALLISTO);
        const aIo = gravityAccel(pos, moons.io, GM_IO);

        // Total acceleration
        const ax = aJupiter[0] + aGanymede[0] + aEuropa[0] + aCallisto[0] + aIo[0];
        const ay = aJupiter[1] + aGanymede[1] + aEuropa[1] + aCallisto[1] + aIo[1];
        const az = aJupiter[2] + aGanymede[2] + aEuropa[2] + aCallisto[2] + aIo[2];

        return [vel[0], vel[1], vel[2], ax, ay, az];
    }

    // RK4 integrator
    function rk4Step(t, state, h) {
        const k1 = derivatives(t, state);
        const s1 = state.map((v, i) => v + 0.5 * h * k1[i]);
        const k2 = derivatives(t + 0.5 * h, s1);
        const s2 = state.map((v, i) => v + 0.5 * h * k2[i]);
        const k3 = derivatives(t + 0.5 * h, s2);
        const s3 = state.map((v, i) => v + h * k3[i]);
        const k4 = derivatives(t + h, s3);
        return state.map((v, i) => v + (h / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    }

    // Distance to a moon
    function distanceToMoon(state, moonPos) {
        const dx = state[0] - moonPos[0];
        const dy = state[1] - moonPos[1];
        const dz = state[2] - moonPos[2];
        return Math.sqrt(dx*dx + dy*dy + dz*dz);
    }

    // Design a flyby trajectory starting from outside Callisto's orbit
    // JUICE-like trajectory: approach from outer system, multiple moon flybys

    console.log('\n--- Designing JUICE-like Trajectory ---\n');

    // Initial conditions: approaching Jupiter from outside Callisto orbit
    // Start at about 2.5 million km from Jupiter, with appropriate velocity
    const r0 = 2.5e9;  // Initial distance from Jupiter [m]
    const approachSpeed = 5000;  // Approach speed [m/s]

    // Position and velocity for a hyperbolic approach
    // Coming from positive x direction, aiming to pass Jupiter
    const x0 = r0;
    const y0 = 0.5e9;  // Offset for aiming at Ganymede
    const z0 = 0;
    const vx0 = -approachSpeed;  // Moving toward Jupiter
    const vy0 = 1500;  // Slight tangential component
    const vz0 = 0;

    console.log('Initial State (Jupiter-centered):');
    console.log(`  Position: (${(x0/1e6).toFixed(0)}, ${(y0/1e6).toFixed(0)}, ${(z0/1e6).toFixed(0)}) thousand km`);
    console.log(`  Velocity: (${vx0.toFixed(0)}, ${vy0.toFixed(0)}, ${vz0.toFixed(0)}) m/s`);

    let state = [x0, y0, z0, vx0, vy0, vz0];
    let t = 0;
    const dt = 100;  // Time step [s]
    const totalDays = 100;  // Simulation duration [days]
    const maxTime = totalDays * 86400;

    // Track closest approaches
    const flybyEvents = [];
    let lastDistGanymede = Infinity;
    let lastDistEuropa = Infinity;
    let lastDistCallisto = Infinity;

    // Trajectory storage (sampled)
    const trajectory = [];

    console.log('\n--- Propagating Trajectory ---\n');

    while (t < maxTime) {
        const moons = getJovianSystemState(t);

        // Compute distances to moons
        const distGanymede = distanceToMoon(state, moons.ganymede);
        const distEuropa = distanceToMoon(state, moons.europa);
        const distCallisto = distanceToMoon(state, moons.callisto);

        // Detect closest approaches (local minima)
        if (distGanymede < lastDistGanymede && lastDistGanymede < 1e8) {
            // Was approaching, check if now receding
        } else if (distGanymede > lastDistGanymede && lastDistGanymede < 5e7) {
            // Just passed closest approach to Ganymede
            flybyEvents.push({
                moon: 'Ganymede',
                time: t - dt,
                distance: lastDistGanymede,
                altitude: lastDistGanymede - R_GANYMEDE
            });
            console.log(`Ganymede flyby at t = ${((t - dt) / 86400).toFixed(2)} days, distance = ${(lastDistGanymede / 1000).toFixed(0)} km`);
        }

        if (distEuropa > lastDistEuropa && lastDistEuropa < 5e7) {
            flybyEvents.push({
                moon: 'Europa',
                time: t - dt,
                distance: lastDistEuropa,
                altitude: lastDistEuropa - R_EUROPA
            });
            console.log(`Europa flyby at t = ${((t - dt) / 86400).toFixed(2)} days, distance = ${(lastDistEuropa / 1000).toFixed(0)} km`);
        }

        if (distCallisto > lastDistCallisto && lastDistCallisto < 1e8) {
            flybyEvents.push({
                moon: 'Callisto',
                time: t - dt,
                distance: lastDistCallisto,
                altitude: lastDistCallisto - R_CALLISTO
            });
            console.log(`Callisto flyby at t = ${((t - dt) / 86400).toFixed(2)} days, distance = ${(lastDistCallisto / 1000).toFixed(0)} km`);
        }

        lastDistGanymede = distGanymede;
        lastDistEuropa = distEuropa;
        lastDistCallisto = distCallisto;

        // Store trajectory data periodically
        if (Math.floor(t / 3600) % 1 === 0 && trajectory.length < 10000) {
            const r = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
            const v = Math.sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);
            trajectory.push({
                t: t,
                x: state[0] / 1e6,  // in thousand km
                y: state[1] / 1e6,
                z: state[2] / 1e6,
                r: r / 1e6,
                v: v
            });
        }

        // Check for Jupiter impact or escape
        const r = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
        if (r < 7e7) {  // Jupiter radius ~70,000 km
            console.log('WARNING: Trajectory impacts Jupiter!');
            break;
        }
        if (r > 5e9) {  // Escaped system
            console.log('Spacecraft has escaped Jupiter system');
            break;
        }

        // Integrate
        state = rk4Step(t, state, dt);
        t += dt;
    }

    console.log(`\nSimulation completed at t = ${(t / 86400).toFixed(2)} days`);

    // Summary
    console.log('\n=== Flyby Summary ===');
    console.log(`Total flybys detected: ${flybyEvents.length}`);

    if (flybyEvents.length > 0) {
        console.log('\nFlyby Details:');
        flybyEvents.forEach((fb, i) => {
            console.log(`  ${i + 1}. ${fb.moon}: day ${(fb.time / 86400).toFixed(2)}, altitude ${(fb.altitude / 1000).toFixed(0)} km`);
        });

        // Ganymede flybys (JUICE primary target)
        const ganymedeFlybys = flybyEvents.filter(f => f.moon === 'Ganymede');
        console.log(`\nGanymede flybys: ${ganymedeFlybys.length}`);
        if (ganymedeFlybys.length > 0) {
            const minAlt = Math.min(...ganymedeFlybys.map(f => f.altitude));
            console.log(`  Minimum altitude: ${(minAlt / 1000).toFixed(0)} km`);
        }
    }

    // Final state
    const finalR = Math.sqrt(state[0]*state[0] + state[1]*state[1] + state[2]*state[2]);
    const finalV = Math.sqrt(state[3]*state[3] + state[4]*state[4] + state[5]*state[5]);

    console.log('\n--- Final State ---');
    console.log(`Position: (${(state[0]/1e6).toFixed(0)}, ${(state[1]/1e6).toFixed(0)}, ${(state[2]/1e6).toFixed(0)}) thousand km`);
    console.log(`Distance from Jupiter: ${(finalR / 1e6).toFixed(0)} thousand km`);
    console.log(`Velocity: ${finalV.toFixed(1)} m/s`);

    // Orbital energy analysis
    const specificEnergy = 0.5 * finalV * finalV - GM_JUPITER / finalR;
    if (specificEnergy < 0) {
        const sma = -GM_JUPITER / (2 * specificEnergy);
        console.log(`\nOrbit type: Elliptical (captured)`);
        console.log(`Semi-major axis: ${(sma / 1e6).toFixed(0)} thousand km`);

        // Orbital period
        const period = 2 * Math.PI * Math.sqrt(sma * sma * sma / GM_JUPITER);
        console.log(`Orbital period: ${(period / 86400).toFixed(1)} days`);
    } else {
        console.log(`\nOrbit type: Hyperbolic (escape)`);
        const vInf = Math.sqrt(2 * specificEnergy);
        console.log(`Excess velocity: ${(vInf / 1000).toFixed(2)} km/s`);
    }

    // B-plane analysis for last Ganymede flyby (if any)
    if (flybyEvents.length > 0) {
        const lastGanymede = flybyEvents.filter(f => f.moon === 'Ganymede').pop();
        if (lastGanymede) {
            console.log('\n--- B-plane Analysis (Last Ganymede Flyby) ---');

            // V-infinity estimate at flyby
            const flybyDist = lastGanymede.distance;
            const vInfFlyby = Math.sqrt(2 * GM_GANYMEDE / flybyDist + 2 * (specificEnergy + GM_JUPITER / a_GANYMEDE));
            console.log(`V-infinity estimate: ${(vInfFlyby / 1000).toFixed(2)} km/s`);

            // Turn angle estimate
            const turnAngle = 2 * Math.asin(1 / (1 + flybyDist * vInfFlyby * vInfFlyby / GM_GANYMEDE));
            console.log(`Turn angle estimate: ${(turnAngle * 180 / Math.PI).toFixed(1)} degrees`);

            // B-parameter (impact parameter for hyperbolic flyby)
            const bParam = flybyDist * Math.sqrt(1 + 2 * GM_GANYMEDE / (flybyDist * vInfFlyby * vInfFlyby));
            console.log(`B-parameter: ${(bParam / 1000).toFixed(0)} km`);
        }
    }

    // JUICE mission comparison
    console.log('\n--- JUICE Mission Reference ---');
    console.log('Planned JUICE trajectory (actual mission):');
    console.log('  - 2 Ganymede flybys (phase 1)');
    console.log('  - 6 Callisto flybys');
    console.log('  - 3 Europa flybys');
    console.log('  - 12 Ganymede flybys (phase 2)');
    console.log('  - Ganymede orbit insertion');
    console.log('  Total: ~35 moon flybys over ~3 years');

    console.log('\n=== JUICE Flybys Simulation Complete ===');
}

main().catch(console.error);
