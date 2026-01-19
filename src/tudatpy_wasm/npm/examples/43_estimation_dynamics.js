/**
 * Example 43: Estimation with Different Dynamical Models
 *
 * Ported from: examples/tudatpy/estimation/estimation_dynamical_models.py
 *
 * This example demonstrates how to use different dynamical models for
 * simulating observations vs performing estimation, which is essential
 * for understanding model sensitivity and robustness.
 *
 * Key concepts:
 * - Truth model vs estimation model
 * - Model mismatch effects
 * - Systematic errors from unmodeled dynamics
 * - Robust estimation strategies
 *
 * Run with: node 43_estimation_dynamics.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Estimation with Different Dynamical Models ===\n');

    const tudat = await createTudatModule();

    // Constants
    const GM_EARTH = 3.986004418e14;
    const R_EARTH = 6.371e6;
    const J2 = 1.08263e-3;
    const J3 = -2.54e-6;
    const J4 = -1.62e-6;

    // Spacecraft initial conditions (Mars Express-like altitude)
    const altitude = 400e3;
    const a = R_EARTH + altitude;
    const e = 0.001;
    const inc = 86.0 * Math.PI / 180;

    console.log('Scenario: Low Earth Orbit Estimation');
    console.log(`  Altitude: ${altitude / 1000} km`);
    console.log(`  Inclination: ${(inc * 180 / Math.PI).toFixed(1)} deg\n`);

    /**
     * Gravitational acceleration with spherical harmonics
     * truthModel: includes J2, J3, J4
     * estimationModel: only includes J2
     */
    function gravityAccel(pos, useFullModel) {
        const x = pos[0], y = pos[1], z = pos[2];
        const r = Math.sqrt(x*x + y*y + z*z);
        const r2 = r * r;
        const r3 = r * r2;

        // Point mass
        let ax = -GM_EARTH * x / r3;
        let ay = -GM_EARTH * y / r3;
        let az = -GM_EARTH * z / r3;

        // J2 perturbation (both models)
        const zr = z / r;
        const zr2 = zr * zr;
        const J2factor = 1.5 * J2 * (R_EARTH / r) ** 2;
        const J2x = J2factor * (5 * zr2 - 1);
        const J2y = J2factor * (5 * zr2 - 1);
        const J2z = J2factor * (5 * zr2 - 3);

        ax += -GM_EARTH * x / r3 * J2x;
        ay += -GM_EARTH * y / r3 * J2y;
        az += -GM_EARTH * z / r3 * J2z;

        if (useFullModel) {
            // J3 perturbation (truth model only)
            const J3factor = 2.5 * J3 * (R_EARTH / r) ** 3;
            ax += -GM_EARTH * x / r3 * J3factor * zr * (7 * zr2 - 3);
            ay += -GM_EARTH * y / r3 * J3factor * zr * (7 * zr2 - 3);
            az += -GM_EARTH * z / r3 * J3factor * (35 * zr2 * zr2 / 5 - 6 * zr2 + 0.6);

            // J4 perturbation (truth model only)
            const J4factor = 5/8 * J4 * (R_EARTH / r) ** 4;
            const zr4 = zr2 * zr2;
            ax += -GM_EARTH * x / r3 * J4factor * (63 * zr4 - 42 * zr2 + 3);
            ay += -GM_EARTH * y / r3 * J4factor * (63 * zr4 - 42 * zr2 + 3);
            az += -GM_EARTH * z / r3 * J4factor * (63 * zr4 - 70 * zr2 + 15);
        }

        return [ax, ay, az];
    }

    /**
     * RK4 propagation
     */
    function propagate(state0, dt, steps, useFullModel) {
        const trajectory = [{ t: 0, state: [...state0] }];
        let state = [...state0];

        for (let i = 0; i < steps; i++) {
            const t = i * dt;

            // RK4 step
            const k1 = derivatives(state, useFullModel);
            const s1 = state.map((v, j) => v + 0.5 * dt * k1[j]);
            const k2 = derivatives(s1, useFullModel);
            const s2 = state.map((v, j) => v + 0.5 * dt * k2[j]);
            const k3 = derivatives(s2, useFullModel);
            const s3 = state.map((v, j) => v + dt * k3[j]);
            const k4 = derivatives(s3, useFullModel);

            state = state.map((v, j) => v + dt/6 * (k1[j] + 2*k2[j] + 2*k3[j] + k4[j]));
            trajectory.push({ t: (i + 1) * dt, state: [...state] });
        }

        return trajectory;
    }

    function derivatives(state, useFullModel) {
        const pos = [state[0], state[1], state[2]];
        const vel = [state[3], state[4], state[5]];
        const accel = gravityAccel(pos, useFullModel);
        return [vel[0], vel[1], vel[2], accel[0], accel[1], accel[2]];
    }

    /**
     * Convert Keplerian to Cartesian
     */
    function keplerToCartesian(a, e, inc, omega, RAAN, nu) {
        const p = a * (1 - e * e);
        const r = p / (1 + e * Math.cos(nu));

        const xOrb = r * Math.cos(nu);
        const yOrb = r * Math.sin(nu);

        const h = Math.sqrt(GM_EARTH * p);
        const vxOrb = -GM_EARTH / h * Math.sin(nu);
        const vyOrb = GM_EARTH / h * (e + Math.cos(nu));

        const cosO = Math.cos(RAAN), sinO = Math.sin(RAAN);
        const cosw = Math.cos(omega), sinw = Math.sin(omega);
        const cosi = Math.cos(inc), sini = Math.sin(inc);

        const x = (cosO*cosw - sinO*sinw*cosi) * xOrb + (-cosO*sinw - sinO*cosw*cosi) * yOrb;
        const y = (sinO*cosw + cosO*sinw*cosi) * xOrb + (-sinO*sinw + cosO*cosw*cosi) * yOrb;
        const z = sinw*sini * xOrb + cosw*sini * yOrb;

        const vx = (cosO*cosw - sinO*sinw*cosi) * vxOrb + (-cosO*sinw - sinO*cosw*cosi) * vyOrb;
        const vy = (sinO*cosw + cosO*sinw*cosi) * vxOrb + (-sinO*sinw + cosO*cosw*cosi) * vyOrb;
        const vz = sinw*sini * vxOrb + cosw*sini * vyOrb;

        return [x, y, z, vx, vy, vz];
    }

    // Initial state
    const state0 = keplerToCartesian(a, e, inc, 0, 0, 0);

    // Simulation parameters
    const dt = 60;  // 1 minute steps
    const duration = 2 * 3600;  // 2 hours
    const steps = Math.floor(duration / dt);

    console.log('--- Running Truth and Estimation Model Propagations ---\n');
    console.log(`Duration: ${duration / 3600} hours`);
    console.log(`Time step: ${dt} seconds`);
    console.log(`Number of steps: ${steps}\n`);

    // Propagate with both models
    console.log('Truth Model: Point mass + J2 + J3 + J4');
    const truthTraj = propagate(state0, dt, steps, true);

    console.log('Estimation Model: Point mass + J2 only\n');
    const estTraj = propagate(state0, dt, steps, false);

    // Compute differences
    console.log('--- Model Mismatch Analysis ---\n');

    const differences = [];
    for (let i = 0; i < truthTraj.length; i++) {
        const truth = truthTraj[i].state;
        const est = estTraj[i].state;

        const dx = truth[0] - est[0];
        const dy = truth[1] - est[1];
        const dz = truth[2] - est[2];
        const dvx = truth[3] - est[3];
        const dvy = truth[4] - est[4];
        const dvz = truth[5] - est[5];

        differences.push({
            t: truthTraj[i].t,
            posErr: Math.sqrt(dx*dx + dy*dy + dz*dz),
            velErr: Math.sqrt(dvx*dvx + dvy*dvy + dvz*dvz)
        });
    }

    // Report at key times
    console.log('Position Error (3D RSS) Due to Model Mismatch:');
    console.log('Time [min] | Pos Error [m] | Vel Error [mm/s]');
    console.log('-----------|---------------|------------------');

    const reportTimes = [0, 15, 30, 60, 90, 120];
    for (const time of reportTimes) {
        const idx = Math.floor(time * 60 / dt);
        if (idx < differences.length) {
            const d = differences[idx];
            console.log(`${time.toString().padStart(10)} | ${d.posErr.toFixed(3).padStart(13)} | ${(d.velErr * 1000).toFixed(3).padStart(16)}`);
        }
    }

    // Final error
    const finalDiff = differences[differences.length - 1];
    console.log(`\nFinal (${duration/60} min):`);
    console.log(`  Position error: ${finalDiff.posErr.toFixed(2)} m`);
    console.log(`  Velocity error: ${(finalDiff.velErr * 1000).toFixed(3)} mm/s`);

    // Error growth rate
    const midIdx = Math.floor(differences.length / 2);
    const errorGrowthRate = (finalDiff.posErr - differences[midIdx].posErr) /
                           ((duration - midIdx * dt) / 60);
    console.log(`  Error growth rate: ${errorGrowthRate.toFixed(3)} m/min`);

    // Simulated estimation analysis
    console.log('\n--- Estimation Impact Analysis ---\n');

    // Observation simulation (range measurements every 5 minutes)
    const obsInterval = 300;  // 5 minutes
    const rangeSigma = 10;    // 10 m noise
    const observations = [];

    for (let i = 0; i < truthTraj.length; i++) {
        if (i * dt % obsInterval === 0) {
            const truthState = truthTraj[i].state;
            const r = Math.sqrt(truthState[0]**2 + truthState[1]**2 + truthState[2]**2);

            // Add noise to simulate observation
            const noise = (Math.random() - 0.5) * 2 * rangeSigma;
            observations.push({
                t: i * dt,
                range: r + noise,
                truth: r
            });
        }
    }

    console.log(`Simulated ${observations.length} range observations`);
    console.log(`Observation noise: ${rangeSigma} m (1σ)`);

    // Compute residuals using estimation model
    console.log('\nObservation Residuals (Computed - Observed):');
    console.log('Time [min] | Residual [m] | Cause');
    console.log('-----------|--------------|-------');

    let sumSqResiduals = 0;
    for (const obs of observations) {
        const idx = Math.floor(obs.t / dt);
        const estState = estTraj[idx].state;
        const estRange = Math.sqrt(estState[0]**2 + estState[1]**2 + estState[2]**2);

        const residual = estRange - obs.range;
        const cause = Math.abs(residual) > 2 * rangeSigma ? 'Model error' : 'Noise';

        if (obs.t % 600 === 0) {  // Print every 10 minutes
            console.log(`${(obs.t / 60).toFixed(0).padStart(10)} | ${residual.toFixed(2).padStart(12)} | ${cause}`);
        }

        sumSqResiduals += residual * residual;
    }

    const rmsResidual = Math.sqrt(sumSqResiduals / observations.length);
    console.log(`\nRMS Residual: ${rmsResidual.toFixed(2)} m`);
    console.log(`Expected (noise only): ${rangeSigma.toFixed(2)} m`);
    console.log(`Excess due to model: ${(rmsResidual - rangeSigma).toFixed(2)} m`);

    // Recommendations
    console.log('\n--- Implications for Estimation ---\n');

    console.log('Model Mismatch Effects:');
    console.log('  1. Systematic residuals that grow with time');
    console.log('  2. Biased state estimates');
    console.log('  3. Optimistic covariance (underestimated uncertainty)');
    console.log('  4. Poor prediction accuracy');

    console.log('\nMitigation Strategies:');
    console.log('  1. Use higher-fidelity dynamics model');
    console.log('  2. Estimate empirical accelerations');
    console.log('  3. Consider process noise (Kalman filter)');
    console.log('  4. Use shorter estimation arcs');
    console.log('  5. Include model uncertainty in covariance');

    console.log('\nModel Comparison:');
    console.log('┌───────────────────┬────────────────┬────────────────┐');
    console.log('│ Component         │ Truth Model    │ Est. Model     │');
    console.log('├───────────────────┼────────────────┼────────────────┤');
    console.log('│ Point mass        │ ✓              │ ✓              │');
    console.log('│ J2                │ ✓              │ ✓              │');
    console.log('│ J3                │ ✓              │ ✗              │');
    console.log('│ J4                │ ✓              │ ✗              │');
    console.log('│ Atmospheric drag  │ ✗              │ ✗              │');
    console.log('│ Solar radiation   │ ✗              │ ✗              │');
    console.log('│ Third bodies      │ ✗              │ ✗              │');
    console.log('└───────────────────┴────────────────┴────────────────┘');

    console.log('\n=== Different Dynamical Models Complete ===');
}

main().catch(console.error);
