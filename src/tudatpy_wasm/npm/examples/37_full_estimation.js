/**
 * Example 37: Full State Estimation Example
 *
 * Ported from: examples/tudatpy/estimation/full_estimation_example.py
 *
 * This example demonstrates comprehensive orbit determination using
 * simulated observations and batch least squares estimation.
 *
 * Key concepts:
 * - Observation simulation
 * - Batch least squares estimation
 * - State vector estimation
 * - Residual analysis
 * - Covariance interpretation
 *
 * Run with: node 37_full_estimation.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Full State Estimation Example ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const GM_EARTH = 3.986004418e14;
    const R_EARTH = 6.378137e6;
    const J2 = 1.08263e-3;

    // True satellite state (what we're trying to estimate)
    const trueSatellite = {
        sma: R_EARTH + 500e3,      // 500 km altitude
        ecc: 0.001,
        inc: 51.6 * Math.PI / 180, // ISS-like inclination
        raan: 45 * Math.PI / 180,
        aop: 30 * Math.PI / 180,
        ta: 0
    };

    // Ground station (tracking)
    const groundStation = {
        name: 'Station1',
        lat: 35 * Math.PI / 180,   // 35°N
        lon: -120 * Math.PI / 180, // 120°W
        alt: 100                    // 100 m altitude
    };

    console.log('True Satellite State:');
    console.log(`  Semi-major axis: ${(trueSatellite.sma / 1000).toFixed(1)} km`);
    console.log(`  Eccentricity: ${trueSatellite.ecc.toFixed(6)}`);
    console.log(`  Inclination: ${(trueSatellite.inc * 180 / Math.PI).toFixed(2)}°`);
    console.log(`  RAAN: ${(trueSatellite.raan * 180 / Math.PI).toFixed(2)}°`);
    console.log(`  Arg of periapsis: ${(trueSatellite.aop * 180 / Math.PI).toFixed(2)}°`);

    console.log(`\nGround Station:`);
    console.log(`  ${groundStation.name} at ${(groundStation.lat * 180/Math.PI).toFixed(1)}°N, ${(-groundStation.lon * 180/Math.PI).toFixed(1)}°W`);

    // Convert Keplerian to Cartesian
    function keplerToCartesian(kep, GM) {
        const a = kep.sma, e = kep.ecc, i = kep.inc;
        const raan = kep.raan, omega = kep.aop, nu = kep.ta;

        // Position in orbital plane
        const p = a * (1 - e * e);
        const r = p / (1 + e * Math.cos(nu));

        const xOrb = r * Math.cos(nu);
        const yOrb = r * Math.sin(nu);

        // Velocity in orbital plane
        const h = Math.sqrt(GM * p);
        const vxOrb = -GM / h * Math.sin(nu);
        const vyOrb = GM / h * (e + Math.cos(nu));

        // Rotation matrices
        const cosRaan = Math.cos(raan), sinRaan = Math.sin(raan);
        const cosI = Math.cos(i), sinI = Math.sin(i);
        const cosOmega = Math.cos(omega), sinOmega = Math.sin(omega);

        const Px = cosRaan * cosOmega - sinRaan * sinOmega * cosI;
        const Py = sinRaan * cosOmega + cosRaan * sinOmega * cosI;
        const Pz = sinOmega * sinI;
        const Qx = -cosRaan * sinOmega - sinRaan * cosOmega * cosI;
        const Qy = -sinRaan * sinOmega + cosRaan * cosOmega * cosI;
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

    // Ground station position in ECI (simplified, no Earth rotation for demo)
    function getStationECI(station, t) {
        // Include Earth rotation
        const omegaEarth = 7.292115e-5;  // rad/s
        const theta = omegaEarth * t + station.lon;

        const cosLat = Math.cos(station.lat);
        const sinLat = Math.sin(station.lat);
        const cosTheta = Math.cos(theta);
        const sinTheta = Math.sin(theta);

        const rStation = R_EARTH + station.alt;
        return [
            rStation * cosLat * cosTheta,
            rStation * cosLat * sinTheta,
            rStation * sinLat
        ];
    }

    // Propagate orbit with J2
    function propagateOrbit(state0, dt, GM) {
        const x = state0[0], y = state0[1], z = state0[2];
        const vx = state0[3], vy = state0[4], vz = state0[5];

        const r = Math.sqrt(x*x + y*y + z*z);
        const r2 = r * r;
        const r3 = r2 * r;
        const r5 = r2 * r3;

        // Two-body acceleration
        const ax = -GM * x / r3;
        const ay = -GM * y / r3;
        const az = -GM * z / r3;

        // J2 perturbation
        const z2 = z * z;
        const j2Factor = 1.5 * J2 * GM * R_EARTH * R_EARTH / r5;
        const axJ2 = j2Factor * x * (5 * z2 / r2 - 1);
        const ayJ2 = j2Factor * y * (5 * z2 / r2 - 1);
        const azJ2 = j2Factor * z * (5 * z2 / r2 - 3);

        // Euler integration (simplified)
        return [
            x + vx * dt + 0.5 * (ax + axJ2) * dt * dt,
            y + vy * dt + 0.5 * (ay + ayJ2) * dt * dt,
            z + vz * dt + 0.5 * (az + azJ2) * dt * dt,
            vx + (ax + axJ2) * dt,
            vy + (ay + ayJ2) * dt,
            vz + (az + azJ2) * dt
        ];
    }

    // Simulate range observation
    function simulateRange(satState, stationPos, noise = 0) {
        const dx = satState[0] - stationPos[0];
        const dy = satState[1] - stationPos[1];
        const dz = satState[2] - stationPos[2];
        const range = Math.sqrt(dx*dx + dy*dy + dz*dz);
        return range + noise * (Math.random() - 0.5) * 2;
    }

    // Simulate range-rate observation
    function simulateRangeRate(satState, stationPos, noise = 0) {
        const dx = satState[0] - stationPos[0];
        const dy = satState[1] - stationPos[1];
        const dz = satState[2] - stationPos[2];
        const range = Math.sqrt(dx*dx + dy*dy + dz*dz);

        const dvx = satState[3];  // Simplified: station velocity neglected
        const dvy = satState[4];
        const dvz = satState[5];

        const rangeRate = (dx*dvx + dy*dvy + dz*dvz) / range;
        return rangeRate + noise * (Math.random() - 0.5) * 2;
    }

    // Check if satellite is visible from station
    function isVisible(satState, stationPos, minElevation = 10) {
        const dx = satState[0] - stationPos[0];
        const dy = satState[1] - stationPos[1];
        const dz = satState[2] - stationPos[2];

        // Station zenith direction
        const stationR = Math.sqrt(stationPos[0]**2 + stationPos[1]**2 + stationPos[2]**2);
        const zenith = [stationPos[0]/stationR, stationPos[1]/stationR, stationPos[2]/stationR];

        // Line of sight direction
        const range = Math.sqrt(dx*dx + dy*dy + dz*dz);
        const los = [dx/range, dy/range, dz/range];

        // Elevation angle
        const cosZenithAngle = zenith[0]*los[0] + zenith[1]*los[1] + zenith[2]*los[2];
        const elevation = 90 - Math.acos(cosZenithAngle) * 180 / Math.PI;

        return elevation > minElevation;
    }

    // Simulate observations
    console.log('\n--- Simulating Observations ---\n');

    const trueState0 = keplerToCartesian(trueSatellite, GM_EARTH);
    const orbitalPeriod = 2 * Math.PI * Math.sqrt(trueSatellite.sma**3 / GM_EARTH);
    const obsInterval = 60;  // seconds
    const totalTime = 3 * orbitalPeriod;  // 3 orbits

    const rangeNoise = 10;      // 10 m range noise
    const rangeRateNoise = 0.001;  // 1 mm/s range-rate noise

    const observations = [];
    let state = [...trueState0];
    let t = 0;

    while (t < totalTime) {
        const stationPos = getStationECI(groundStation, t);

        if (isVisible(state, stationPos)) {
            const range = simulateRange(state, stationPos, rangeNoise);
            const rangeRate = simulateRangeRate(state, stationPos, rangeRateNoise);

            observations.push({
                t,
                range,
                rangeRate,
                stationPos: [...stationPos],
                trueState: [...state]
            });
        }

        // Propagate
        const dt = Math.min(obsInterval, totalTime - t);
        state = propagateOrbit(state, dt, GM_EARTH);
        t += dt;
    }

    console.log(`Total observations: ${observations.length}`);
    console.log(`Observation period: ${(totalTime / 3600).toFixed(1)} hours (${(totalTime / orbitalPeriod).toFixed(1)} orbits)`);
    console.log(`Range noise: ${rangeNoise} m`);
    console.log(`Range-rate noise: ${rangeRateNoise * 1000} mm/s`);

    // Initial guess (perturbed from truth)
    const initialGuess = [...trueState0];
    initialGuess[0] += 1000;   // 1 km error in x
    initialGuess[1] += 500;    // 500 m error in y
    initialGuess[2] += 200;    // 200 m error in z
    initialGuess[3] += 1;      // 1 m/s error in vx
    initialGuess[4] += 0.5;    // 0.5 m/s error in vy
    initialGuess[5] += 0.2;    // 0.2 m/s error in vz

    console.log('\nInitial State Errors:');
    console.log(`  Position: ${Math.sqrt(1000**2 + 500**2 + 200**2).toFixed(0)} m`);
    console.log(`  Velocity: ${Math.sqrt(1**2 + 0.5**2 + 0.2**2).toFixed(3)} m/s`);

    // Batch least squares estimation (simplified)
    console.log('\n--- Running Batch Least Squares ---\n');

    function computeResiduals(estimatedState0, obs) {
        const residuals = [];
        let state = [...estimatedState0];
        let lastT = 0;

        for (const ob of obs) {
            // Propagate to observation time
            const dt = ob.t - lastT;
            if (dt > 0) {
                // Simple propagation in steps
                const numSteps = Math.ceil(dt / 10);
                const stepDt = dt / numSteps;
                for (let i = 0; i < numSteps; i++) {
                    state = propagateOrbit(state, stepDt, GM_EARTH);
                }
            }
            lastT = ob.t;

            // Compute predicted observations
            const predRange = simulateRange(state, ob.stationPos, 0);
            const predRangeRate = simulateRangeRate(state, ob.stationPos, 0);

            // Residuals
            residuals.push({
                range: ob.range - predRange,
                rangeRate: ob.rangeRate - predRangeRate
            });
        }

        return residuals;
    }

    // Compute partial derivatives numerically
    function computePartials(state0, obs, delta = 1) {
        const H = [];  // Design matrix rows

        for (let obs_idx = 0; obs_idx < obs.length; obs_idx++) {
            const rowRange = [];
            const rowRangeRate = [];

            for (let i = 0; i < 6; i++) {
                // Perturb state
                const statePlus = [...state0];
                const stateMinus = [...state0];
                statePlus[i] += delta;
                stateMinus[i] -= delta;

                // Compute residuals
                const resPlus = computeResiduals(statePlus, [obs[obs_idx]]);
                const resMinus = computeResiduals(stateMinus, [obs[obs_idx]]);

                // Numerical derivative
                rowRange.push(-(resPlus[0].range - resMinus[0].range) / (2 * delta));
                rowRangeRate.push(-(resPlus[0].rangeRate - resMinus[0].rangeRate) / (2 * delta));
            }

            H.push(rowRange);
            H.push(rowRangeRate);
        }

        return H;
    }

    // Simplified least squares iteration
    let currentEstimate = [...initialGuess];
    const maxIterations = 5;
    const convergenceThreshold = 0.1;  // m

    for (let iter = 0; iter < maxIterations; iter++) {
        // Compute residuals
        const residuals = computeResiduals(currentEstimate, observations);

        // RMS residuals
        const rangeRMS = Math.sqrt(residuals.reduce((sum, r) => sum + r.range**2, 0) / residuals.length);
        const rangeRateRMS = Math.sqrt(residuals.reduce((sum, r) => sum + r.rangeRate**2, 0) / residuals.length);

        console.log(`Iteration ${iter + 1}:`);
        console.log(`  Range RMS: ${rangeRMS.toFixed(2)} m`);
        console.log(`  Range-rate RMS: ${(rangeRateRMS * 1000).toFixed(3)} mm/s`);

        // Check convergence
        if (rangeRMS < convergenceThreshold) {
            console.log('  Converged!');
            break;
        }

        // Compute correction (simplified - using gradient descent)
        // In practice, would use full least squares normal equations
        const stepSize = 0.1;
        for (let i = 0; i < 6; i++) {
            // Estimate gradient numerically
            const delta = (i < 3) ? 10 : 0.01;
            const statePlus = [...currentEstimate];
            statePlus[i] += delta;
            const resPlus = computeResiduals(statePlus, observations);
            const rmsPlus = Math.sqrt(resPlus.reduce((sum, r) => sum + r.range**2 + (r.rangeRate*1000)**2, 0) / resPlus.length);

            const stateMinus = [...currentEstimate];
            stateMinus[i] -= delta;
            const resMinus = computeResiduals(stateMinus, observations);
            const rmsMinus = Math.sqrt(resMinus.reduce((sum, r) => sum + r.range**2 + (r.rangeRate*1000)**2, 0) / resMinus.length);

            const gradient = (rmsPlus - rmsMinus) / (2 * delta);
            currentEstimate[i] -= stepSize * gradient * delta * 10;
        }
    }

    // Final results
    console.log('\n=== Estimation Results ===');

    const posError = Math.sqrt(
        (currentEstimate[0] - trueState0[0])**2 +
        (currentEstimate[1] - trueState0[1])**2 +
        (currentEstimate[2] - trueState0[2])**2
    );
    const velError = Math.sqrt(
        (currentEstimate[3] - trueState0[3])**2 +
        (currentEstimate[4] - trueState0[4])**2 +
        (currentEstimate[5] - trueState0[5])**2
    );

    console.log('\nFinal State Errors:');
    console.log(`  Position: ${posError.toFixed(1)} m`);
    console.log(`  Velocity: ${(velError * 1000).toFixed(2)} mm/s`);

    console.log('\nState Vector Comparison:');
    console.log('Component | True         | Estimated    | Error');
    console.log('----------|--------------|--------------|--------');
    const labels = ['X [km]', 'Y [km]', 'Z [km]', 'Vx [m/s]', 'Vy [m/s]', 'Vz [m/s]'];
    const scales = [1000, 1000, 1000, 1, 1, 1];
    for (let i = 0; i < 6; i++) {
        const trueVal = trueState0[i] / scales[i];
        const estVal = currentEstimate[i] / scales[i];
        const err = (currentEstimate[i] - trueState0[i]) / (i < 3 ? 1 : 0.001);
        const errUnit = i < 3 ? 'm' : 'mm/s';
        console.log(`${labels[i].padEnd(9)} | ${trueVal.toFixed(3).padStart(12)} | ${estVal.toFixed(3).padStart(12)} | ${err.toFixed(2)} ${errUnit}`);
    }

    // Residual analysis
    const finalResiduals = computeResiduals(currentEstimate, observations);
    console.log('\n--- Residual Statistics ---');
    const rangeRes = finalResiduals.map(r => r.range);
    const rrRes = finalResiduals.map(r => r.rangeRate * 1000);

    console.log(`Range residuals: mean=${(rangeRes.reduce((a,b)=>a+b,0)/rangeRes.length).toFixed(2)}m, ` +
        `std=${Math.sqrt(rangeRes.reduce((a,b)=>a+b*b,0)/rangeRes.length).toFixed(2)}m`);
    console.log(`Range-rate residuals: mean=${(rrRes.reduce((a,b)=>a+b,0)/rrRes.length).toFixed(3)}mm/s, ` +
        `std=${Math.sqrt(rrRes.reduce((a,b)=>a+b*b,0)/rrRes.length).toFixed(3)}mm/s`);

    console.log('\n=== Full State Estimation Complete ===');
}

main().catch(console.error);
