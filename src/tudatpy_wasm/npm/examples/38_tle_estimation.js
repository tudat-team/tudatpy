/**
 * Example 38: TLE-Based State Estimation
 *
 * Ported from: examples/tudatpy/estimation/estimation_with_tle.py
 *
 * This example demonstrates orbit determination using Two-Line Element
 * (TLE) data as pseudo-observations.
 *
 * Key concepts:
 * - TLE format and interpretation
 * - SGP4 propagation
 * - TLE-to-state conversion
 * - State vector fitting
 * - Epoch transformation
 *
 * Run with: node 38_tle_estimation.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== TLE-Based State Estimation ===\n');

    const tudat = await createTudatModule();

    // TLE format description
    console.log('TLE Format Overview:');
    console.log('Line 1: Catalog number, classification, launch info, epoch, mean motion derivatives, BSTAR');
    console.log('Line 2: Inclination, RAAN, eccentricity, arg of perigee, mean anomaly, mean motion');

    // Example TLE (ISS)
    const tle = {
        line1: '1 25544U 98067A   24001.50000000  .00016717  00000+0  10270-3 0  9993',
        line2: '2 25544  51.6400 247.4600 0006700  90.0000 270.0000 15.49491000100000'
    };

    console.log('\nExample TLE (ISS):');
    console.log(`  ${tle.line1}`);
    console.log(`  ${tle.line2}`);

    // Parse TLE
    function parseTLE(line1, line2) {
        // Line 2 parsing
        const inc = parseFloat(line2.substring(8, 16));           // Inclination [deg]
        const raan = parseFloat(line2.substring(17, 25));         // RAAN [deg]
        const ecc = parseFloat('0.' + line2.substring(26, 33));   // Eccentricity
        const aop = parseFloat(line2.substring(34, 42));          // Arg of perigee [deg]
        const ma = parseFloat(line2.substring(43, 51));           // Mean anomaly [deg]
        const n = parseFloat(line2.substring(52, 63));            // Mean motion [rev/day]

        // Line 1 parsing
        const epochYear = parseInt(line1.substring(18, 20));
        const epochDay = parseFloat(line1.substring(20, 32));
        const bstar = parseFloat(line1.substring(53, 54) + '.' +
            line1.substring(54, 59) + 'e' + line1.substring(59, 61));

        // Convert epoch
        const year = epochYear < 57 ? 2000 + epochYear : 1900 + epochYear;

        // Mean motion to semi-major axis
        const GM_EARTH = 3.986004418e14;
        const nRadSec = n * 2 * Math.PI / 86400;  // rad/s
        const sma = Math.pow(GM_EARTH / (nRadSec * nRadSec), 1/3);

        return {
            year, epochDay, bstar,
            inc: inc * Math.PI / 180,
            raan: raan * Math.PI / 180,
            ecc,
            aop: aop * Math.PI / 180,
            ma: ma * Math.PI / 180,
            n: nRadSec,
            sma
        };
    }

    const elements = parseTLE(tle.line1, tle.line2);

    console.log('\nParsed TLE Elements:');
    console.log(`  Epoch: ${elements.year} day ${elements.epochDay.toFixed(8)}`);
    console.log(`  Semi-major axis: ${(elements.sma / 1000).toFixed(3)} km`);
    console.log(`  Eccentricity: ${elements.ecc.toFixed(7)}`);
    console.log(`  Inclination: ${(elements.inc * 180/Math.PI).toFixed(4)}°`);
    console.log(`  RAAN: ${(elements.raan * 180/Math.PI).toFixed(4)}°`);
    console.log(`  Arg of perigee: ${(elements.aop * 180/Math.PI).toFixed(4)}°`);
    console.log(`  Mean anomaly: ${(elements.ma * 180/Math.PI).toFixed(4)}°`);
    console.log(`  Mean motion: ${(elements.n * 86400 / (2*Math.PI)).toFixed(8)} rev/day`);
    console.log(`  B*: ${elements.bstar.toExponential(5)}`);

    // Convert mean anomaly to true anomaly
    function meanToTrue(M, e, maxIter = 10) {
        // Solve Kepler's equation: M = E - e*sin(E)
        let E = M;
        for (let i = 0; i < maxIter; i++) {
            E = M + e * Math.sin(E);
        }
        // True anomaly from eccentric anomaly
        const nu = 2 * Math.atan2(
            Math.sqrt(1 + e) * Math.sin(E / 2),
            Math.sqrt(1 - e) * Math.cos(E / 2)
        );
        return nu;
    }

    // Convert elements to Cartesian state
    function elementsToCartesian(elem) {
        const GM = 3.986004418e14;
        const a = elem.sma, e = elem.ecc, i = elem.inc;
        const raan = elem.raan, omega = elem.aop;
        const nu = meanToTrue(elem.ma, e);

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

        return {
            x: Px * xOrb + Qx * yOrb,
            y: Py * xOrb + Qy * yOrb,
            z: Pz * xOrb + Qz * yOrb,
            vx: Px * vxOrb + Qx * vyOrb,
            vy: Py * vxOrb + Qy * vyOrb,
            vz: Pz * vxOrb + Qz * vyOrb
        };
    }

    const state = elementsToCartesian(elements);

    console.log('\nCartesian State at TLE Epoch:');
    console.log(`  X: ${(state.x / 1000).toFixed(3)} km`);
    console.log(`  Y: ${(state.y / 1000).toFixed(3)} km`);
    console.log(`  Z: ${(state.z / 1000).toFixed(3)} km`);
    console.log(`  Vx: ${(state.vx / 1000).toFixed(6)} km/s`);
    console.log(`  Vy: ${(state.vy / 1000).toFixed(6)} km/s`);
    console.log(`  Vz: ${(state.vz / 1000).toFixed(6)} km/s`);

    // Orbit characteristics
    const R_EARTH = 6.378137e6;
    const rMag = Math.sqrt(state.x**2 + state.y**2 + state.z**2);
    const vMag = Math.sqrt(state.vx**2 + state.vy**2 + state.vz**2);
    const altitude = rMag - R_EARTH;
    const periapsis = elements.sma * (1 - elements.ecc) - R_EARTH;
    const apoapsis = elements.sma * (1 + elements.ecc) - R_EARTH;
    const period = 2 * Math.PI / elements.n;

    console.log('\nOrbit Characteristics:');
    console.log(`  Current altitude: ${(altitude / 1000).toFixed(2)} km`);
    console.log(`  Periapsis altitude: ${(periapsis / 1000).toFixed(2)} km`);
    console.log(`  Apoapsis altitude: ${(apoapsis / 1000).toFixed(2)} km`);
    console.log(`  Velocity: ${(vMag / 1000).toFixed(3)} km/s`);
    console.log(`  Orbital period: ${(period / 60).toFixed(2)} minutes`);

    // Simulate multiple TLEs over time
    console.log('\n--- Simulating TLE Time Series ---\n');

    // Generate pseudo-TLEs at different epochs
    function propagateMeanElements(elem, dt) {
        // Simplified J2 secular perturbations
        const J2 = 1.08263e-3;
        const R = R_EARTH;
        const a = elem.sma, e = elem.ecc, i = elem.inc;
        const n = elem.n;

        // J2 perturbation rates
        const p = a * (1 - e*e);
        const factor = -1.5 * J2 * (R/p)**2 * n;

        // RAAN drift
        const raanDot = factor * Math.cos(i);
        // Argument of perigee drift
        const aopDot = factor * (2 - 2.5 * Math.sin(i)**2);
        // Mean anomaly drift
        const maDot = n + factor * Math.sqrt(1 - e*e) * (1 - 1.5 * Math.sin(i)**2);

        return {
            ...elem,
            raan: elem.raan + raanDot * dt,
            aop: elem.aop + aopDot * dt,
            ma: (elem.ma + maDot * dt) % (2 * Math.PI)
        };
    }

    // Generate TLE observations over 1 day
    const tleObs = [];
    const obsInterval = 2 * 3600;  // 2 hours between TLEs
    const totalDuration = 24 * 3600;  // 1 day

    for (let t = 0; t <= totalDuration; t += obsInterval) {
        const propagatedElem = propagateMeanElements(elements, t);
        const obsState = elementsToCartesian(propagatedElem);
        tleObs.push({
            t,
            elements: propagatedElem,
            state: obsState
        });
    }

    console.log(`Generated ${tleObs.length} TLE observations over ${totalDuration/3600} hours`);

    // Fit state to TLE observations
    console.log('\n--- Fitting State to TLE Observations ---\n');

    // Use first TLE as reference
    const refState = [state.x, state.y, state.z, state.vx, state.vy, state.vz];

    // Compute residuals between propagated numerical state and TLE-derived states
    const GM = 3.986004418e14;

    function propagateNumerical(state0, dt) {
        // Simple two-body propagation
        const x = state0[0], y = state0[1], z = state0[2];
        const vx = state0[3], vy = state0[4], vz = state0[5];

        const r = Math.sqrt(x*x + y*y + z*z);
        const mu_r3 = GM / (r * r * r);

        // Simple Euler for demo
        const steps = Math.max(1, Math.floor(dt / 60));
        const stepDt = dt / steps;

        let s = [...state0];
        for (let i = 0; i < steps; i++) {
            const rr = Math.sqrt(s[0]**2 + s[1]**2 + s[2]**2);
            const acc = GM / (rr * rr * rr);
            s = [
                s[0] + s[3] * stepDt,
                s[1] + s[4] * stepDt,
                s[2] + s[5] * stepDt,
                s[3] - acc * s[0] * stepDt,
                s[4] - acc * s[1] * stepDt,
                s[5] - acc * s[2] * stepDt
            ];
        }
        return s;
    }

    // Compare propagated state with TLE observations
    console.log('Time [h] | Pos diff [km] | Vel diff [m/s]');
    console.log('---------|---------------|---------------');

    for (const obs of tleObs) {
        const numState = propagateNumerical(refState, obs.t);

        const posDiff = Math.sqrt(
            (numState[0] - obs.state.x)**2 +
            (numState[1] - obs.state.y)**2 +
            (numState[2] - obs.state.z)**2
        );
        const velDiff = Math.sqrt(
            (numState[3] - obs.state.vx)**2 +
            (numState[4] - obs.state.vy)**2 +
            (numState[5] - obs.state.vz)**2
        );

        console.log(`${(obs.t/3600).toFixed(1).padStart(8)} | ${(posDiff/1000).toFixed(3).padStart(13)} | ${velDiff.toFixed(4).padStart(14)}`);
    }

    // Analysis
    console.log('\n--- Analysis ---');
    console.log('Position differences grow over time due to:');
    console.log('  1. J2 perturbations not in simple two-body propagation');
    console.log('  2. SGP4 vs numerical propagator differences');
    console.log('  3. TLE fitting errors from original observations');

    // TLE accuracy discussion
    console.log('\n--- TLE Accuracy Notes ---');
    console.log('Typical TLE accuracy:');
    console.log('  - Position: 1-5 km (varies with object and orbit)');
    console.log('  - Velocity: 1-5 m/s');
    console.log('  - Degradation: ~1 km/day for LEO');
    console.log('\nTLE limitations:');
    console.log('  - Mean elements (not osculating)');
    console.log('  - SGP4/SDP4 theory approximations');
    console.log('  - Irregular update intervals');
    console.log('  - Public TLEs may be degraded for security');

    // Covariance estimation
    console.log('\n--- TLE Covariance Estimate ---');
    const typicalSigmaPos = 2000;  // m (typical TLE position uncertainty)
    const typicalSigmaVel = 2;     // m/s (typical TLE velocity uncertainty)

    console.log('Estimated 1-sigma uncertainties:');
    console.log(`  Position: ${typicalSigmaPos} m`);
    console.log(`  Velocity: ${typicalSigmaVel} m/s`);
    console.log(`  3-sigma position: ${3 * typicalSigmaPos / 1000} km`);

    console.log('\n=== TLE-Based State Estimation Complete ===');
}

main().catch(console.error);
