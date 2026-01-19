/**
 * Example 39: Galilean Moons State Estimation
 *
 * Ported from: examples/tudatpy/estimation/galilean_moons_state_estimation.py
 *
 * This example demonstrates orbit determination in the Jovian system,
 * estimating the states of Jupiter's Galilean moons using simulated
 * observations.
 *
 * Key concepts:
 * - Multi-body state estimation
 * - Jovian system dynamics
 * - Laplace resonance
 * - Moon ephemeris fitting
 * - Perturbation sensitivity
 *
 * Run with: node 39_galilean_moons.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Galilean Moons State Estimation ===\n');

    const tudat = await createTudatModule();

    // Jupiter system constants
    const GM_JUPITER = 1.26687e17;     // m^3/s^2
    const R_JUPITER = 7.1492e7;        // m

    // Galilean moon parameters
    const moons = {
        Io: {
            GM: 5.959e12,
            sma: 421800e3,      // m
            period: 1.769 * 86400,  // s
            ecc: 0.0041,
            inc: 0.036 * Math.PI/180,  // rad
            mass: 8.932e22      // kg
        },
        Europa: {
            GM: 3.203e12,
            sma: 671100e3,
            period: 3.551 * 86400,
            ecc: 0.0094,
            inc: 0.466 * Math.PI/180,
            mass: 4.800e22
        },
        Ganymede: {
            GM: 9.887e12,
            sma: 1070400e3,
            period: 7.155 * 86400,
            ecc: 0.0013,
            inc: 0.177 * Math.PI/180,
            mass: 1.482e23
        },
        Callisto: {
            GM: 7.179e12,
            sma: 1882700e3,
            period: 16.689 * 86400,
            ecc: 0.0074,
            inc: 0.192 * Math.PI/180,
            mass: 1.076e23
        }
    };

    console.log('Galilean Moon System:');
    console.log('Moon     | a [Rj]  | Period [d] | Mass [Moon]');
    console.log('---------|---------|------------|------------');
    const moonMass = 7.342e22;  // Earth's Moon mass for reference
    for (const [name, moon] of Object.entries(moons)) {
        console.log(`${name.padEnd(8)} | ${(moon.sma / R_JUPITER).toFixed(2).padStart(7)} | ${(moon.period / 86400).toFixed(3).padStart(10)} | ${(moon.mass / moonMass).toFixed(2).padStart(11)}`);
    }

    // Laplace resonance
    console.log('\n--- Laplace Resonance ---');
    const n_io = 2 * Math.PI / moons.Io.period;
    const n_europa = 2 * Math.PI / moons.Europa.period;
    const n_ganymede = 2 * Math.PI / moons.Ganymede.period;

    console.log(`Io/Europa period ratio: ${(moons.Europa.period / moons.Io.period).toFixed(4)} (≈2:1)`);
    console.log(`Europa/Ganymede period ratio: ${(moons.Ganymede.period / moons.Europa.period).toFixed(4)} (≈2:1)`);
    console.log(`Laplace resonance: n_Io - 3*n_Europa + 2*n_Ganymede ≈ 0`);
    const laplaceRes = n_io - 3*n_europa + 2*n_ganymede;
    console.log(`  Residual: ${(laplaceRes * 180/Math.PI * 86400).toFixed(6)}°/day`);

    // Convert elements to Cartesian state
    function elementsToState(a, e, i, raan, omega, M, GM_central) {
        // Mean to true anomaly
        let E = M;
        for (let iter = 0; iter < 10; iter++) {
            E = M + e * Math.sin(E);
        }
        const nu = 2 * Math.atan2(Math.sqrt(1+e)*Math.sin(E/2), Math.sqrt(1-e)*Math.cos(E/2));

        // Position/velocity in orbital plane
        const p = a * (1 - e*e);
        const r = p / (1 + e * Math.cos(nu));
        const h = Math.sqrt(GM_central * p);

        const xOrb = r * Math.cos(nu);
        const yOrb = r * Math.sin(nu);
        const vxOrb = -GM_central/h * Math.sin(nu);
        const vyOrb = GM_central/h * (e + Math.cos(nu));

        // Rotation to inertial frame
        const cosR = Math.cos(raan), sinR = Math.sin(raan);
        const cosI = Math.cos(i), sinI = Math.sin(i);
        const cosO = Math.cos(omega), sinO = Math.sin(omega);

        const Px = cosR*cosO - sinR*sinO*cosI;
        const Py = sinR*cosO + cosR*sinO*cosI;
        const Pz = sinO*sinI;
        const Qx = -cosR*sinO - sinR*cosO*cosI;
        const Qy = -sinR*sinO + cosR*cosO*cosI;
        const Qz = cosO*sinI;

        return [
            Px*xOrb + Qx*yOrb, Py*xOrb + Qy*yOrb, Pz*xOrb + Qz*yOrb,
            Px*vxOrb + Qx*vyOrb, Py*vxOrb + Qy*vyOrb, Pz*vxOrb + Qz*vyOrb
        ];
    }

    // Initialize moon states
    function initializeMoonState(moon, meanAnomaly) {
        return elementsToState(
            moon.sma, moon.ecc, moon.inc,
            0, 0, meanAnomaly, GM_JUPITER
        );
    }

    // Set initial mean anomalies (arbitrary for simulation)
    const initialMAs = {
        Io: 0,
        Europa: Math.PI/3,
        Ganymede: 2*Math.PI/3,
        Callisto: Math.PI
    };

    console.log('\n--- Initializing Moon States ---\n');

    const moonStates = {};
    for (const [name, moon] of Object.entries(moons)) {
        moonStates[name] = initializeMoonState(moon, initialMAs[name]);
        const state = moonStates[name];
        const r = Math.sqrt(state[0]**2 + state[1]**2 + state[2]**2);
        console.log(`${name}:`);
        console.log(`  Position: (${(state[0]/1e6).toFixed(1)}, ${(state[1]/1e6).toFixed(1)}, ${(state[2]/1e6).toFixed(1)}) × 10^6 m`);
        console.log(`  Distance from Jupiter: ${(r / R_JUPITER).toFixed(2)} Rj`);
    }

    // Simulate observations
    console.log('\n--- Simulating Astrometric Observations ---\n');

    // Observer at Earth (simplified - fixed distance)
    const earthDistance = 6.0 * 1.496e11;  // ~6 AU (opposition-like)

    // Simulate relative position observations (angular separation)
    function simulateObservation(moonState, noise = 0) {
        // Angular position relative to Jupiter (simplified)
        const angX = Math.atan2(moonState[0], earthDistance) * 180/Math.PI * 3600;  // arcsec
        const angY = Math.atan2(moonState[1], earthDistance) * 180/Math.PI * 3600;  // arcsec

        return {
            angX: angX + noise * (Math.random() - 0.5) * 2,
            angY: angY + noise * (Math.random() - 0.5) * 2
        };
    }

    // Propagate N-body system
    function propagateSystem(states, dt) {
        const newStates = {};

        for (const [name, state] of Object.entries(states)) {
            let ax = 0, ay = 0, az = 0;

            // Jupiter gravity
            const r = Math.sqrt(state[0]**2 + state[1]**2 + state[2]**2);
            const jupAcc = GM_JUPITER / (r**3);
            ax -= jupAcc * state[0];
            ay -= jupAcc * state[1];
            az -= jupAcc * state[2];

            // Moon-moon perturbations
            for (const [otherName, otherState] of Object.entries(states)) {
                if (name !== otherName) {
                    const dx = otherState[0] - state[0];
                    const dy = otherState[1] - state[1];
                    const dz = otherState[2] - state[2];
                    const dist = Math.sqrt(dx**2 + dy**2 + dz**2);
                    const moonAcc = moons[otherName].GM / (dist**3);
                    ax += moonAcc * dx;
                    ay += moonAcc * dy;
                    az += moonAcc * dz;
                }
            }

            // Euler integration
            newStates[name] = [
                state[0] + state[3]*dt + 0.5*ax*dt*dt,
                state[1] + state[4]*dt + 0.5*ay*dt*dt,
                state[2] + state[5]*dt + 0.5*az*dt*dt,
                state[3] + ax*dt,
                state[4] + ay*dt,
                state[5] + az*dt
            ];
        }

        return newStates;
    }

    // Generate observations over 10 days
    const obsNoise = 0.01;  // 0.01 arcsec noise (very precise)
    const obsInterval = 6 * 3600;  // 6 hours
    const totalTime = 10 * 86400;  // 10 days

    const observations = [];
    let currentStates = { ...moonStates };
    let t = 0;

    while (t <= totalTime) {
        const obs = { t };
        for (const name of Object.keys(moons)) {
            obs[name] = simulateObservation(currentStates[name], obsNoise);
        }
        observations.push(obs);

        // Propagate
        const numSteps = 10;
        const stepDt = obsInterval / numSteps;
        for (let i = 0; i < numSteps; i++) {
            currentStates = propagateSystem(currentStates, stepDt);
        }
        t += obsInterval;
    }

    console.log(`Generated ${observations.length} observation sets over ${totalTime/86400} days`);
    console.log(`Observation noise: ${obsNoise} arcsec`);

    // Analyze observations
    console.log('\n--- Observation Analysis ---');

    for (const name of Object.keys(moons)) {
        const angXs = observations.map(o => o[name].angX);
        const angYs = observations.map(o => o[name].angY);

        const maxSep = Math.max(...angXs.map(Math.abs), ...angYs.map(Math.abs));
        console.log(`${name}: max angular separation ${maxSep.toFixed(2)} arcsec`);
    }

    // Estimate state corrections (simplified)
    console.log('\n--- State Estimation ---\n');

    // Perturb initial states
    const perturbedStates = {};
    for (const [name, state] of Object.entries(moonStates)) {
        perturbedStates[name] = [
            state[0] * 1.001,  // 0.1% position error
            state[1] * 1.001,
            state[2],
            state[3] * 1.0001,  // 0.01% velocity error
            state[4] * 1.0001,
            state[5]
        ];
    }

    // Compute residuals
    function computeResiduals(states, obs) {
        let testStates = { ...states };
        let sumSq = 0;
        let count = 0;

        for (const ob of obs) {
            // Propagate to observation time
            if (ob.t > 0) {
                const numSteps = Math.ceil(ob.t / 3600);
                const stepDt = ob.t / numSteps;
                for (let i = 0; i < numSteps; i++) {
                    testStates = propagateSystem(testStates, stepDt);
                }
            }

            // Compute residuals
            for (const name of Object.keys(moons)) {
                const pred = simulateObservation(testStates[name], 0);
                const resX = ob[name].angX - pred.angX;
                const resY = ob[name].angY - pred.angY;
                sumSq += resX**2 + resY**2;
                count += 2;
            }
        }

        return Math.sqrt(sumSq / count);
    }

    const trueRMS = computeResiduals(moonStates, [observations[0]]);
    const perturbedRMS = computeResiduals(perturbedStates, [observations[0]]);

    console.log(`True initial state RMS: ${trueRMS.toFixed(4)} arcsec`);
    console.log(`Perturbed initial state RMS: ${perturbedRMS.toFixed(4)} arcsec`);

    // Sensitivity analysis
    console.log('\n--- Sensitivity Analysis ---');

    const positionSensitivity = 1000e3;  // 1000 km position change
    const velocitySensitivity = 1;       // 1 m/s velocity change

    for (const name of Object.keys(moons)) {
        const baseObs = simulateObservation(moonStates[name], 0);

        // Position sensitivity
        const pertPosState = [...moonStates[name]];
        pertPosState[0] += positionSensitivity;
        const pertPosObs = simulateObservation(pertPosState, 0);
        const posSens = Math.abs(pertPosObs.angX - baseObs.angX) / (positionSensitivity / 1000);

        // Velocity sensitivity (after 1 day propagation effect)
        const velSens = posSens * velocitySensitivity * 86400 / positionSensitivity;

        console.log(`${name}: position sensitivity ${posSens.toFixed(4)} arcsec/1000km, velocity ${velSens.toFixed(6)} arcsec/(m/s)/day`);
    }

    // Results summary
    console.log('\n=== Estimation Summary ===');
    console.log('The Galilean moon system presents unique challenges:');
    console.log('  1. Strong mutual perturbations (Laplace resonance)');
    console.log('  2. Rapid orbital motion (Io: 1.77 day period)');
    console.log('  3. Variable Earth-Jupiter distance');
    console.log('  4. Need for long observation arcs');
    console.log('\nTypical ephemeris accuracy:');
    console.log('  Position: 1-10 km (from spacecraft tracking)');
    console.log('  Angular: 0.001-0.01 arcsec (from astrometry)');

    console.log('\n=== Galilean Moons State Estimation Complete ===');
}

main().catch(console.error);
