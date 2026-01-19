/**
 * Example 34: Cassini-1 MGA Trajectory Optimization
 *
 * Ported from: examples/tudatpy/mission_design/cassini1_mga_optimization.py
 *
 * This example demonstrates optimization of a Cassini-like multiple
 * gravity assist trajectory using impulsive maneuvers.
 *
 * Key concepts:
 * - Multiple gravity assist (MGA) trajectory design
 * - Impulsive delta-V optimization
 * - EVVEJSA sequence (Earth-Venus-Venus-Earth-Jupiter-Saturn)
 * - Powered vs unpowered flybys
 * - Porkchop plot methodology
 *
 * Run with: node 34_cassini_mga.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Cassini-1 MGA Trajectory Optimization ===\n');

    const tudat = await createTudatModule();

    // Constants
    const AU = 1.496e11;                // Astronomical unit [m]
    const GM_SUN = 1.32712440018e20;    // Sun gravitational parameter [m^3/s^2]
    const JULIAN_DAY = 86400;

    // Planetary data
    const planets = {
        Earth: {
            a: 1.000 * AU, e: 0.0167, GM: 3.986e14, R: 6.378e6,
            SOI: 9.25e8  // Sphere of influence [m]
        },
        Venus: {
            a: 0.723 * AU, e: 0.0068, GM: 3.249e14, R: 6.052e6,
            SOI: 6.16e8
        },
        Jupiter: {
            a: 5.203 * AU, e: 0.0489, GM: 1.267e17, R: 7.149e7,
            SOI: 4.82e10
        },
        Saturn: {
            a: 9.537 * AU, e: 0.0565, GM: 3.793e16, R: 6.027e7,
            SOI: 5.46e10
        }
    };

    // Cassini-1 sequence: Earth -> Venus -> Venus -> Earth -> Jupiter -> Saturn
    const sequence = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn'];

    console.log('Cassini-1 MGA Sequence:');
    console.log(`  ${sequence.join(' -> ')}\n`);

    console.log('Planetary Parameters:');
    for (const [name, data] of Object.entries(planets)) {
        console.log(`  ${name}: a = ${(data.a / AU).toFixed(3)} AU, e = ${data.e}`);
    }

    // Get planet position at time t (simplified circular coplanar)
    function getPlanetPosition(planet, t) {
        const data = planets[planet];
        const n = Math.sqrt(GM_SUN / (data.a * data.a * data.a));  // Mean motion
        const theta = n * t;  // True anomaly (circular)
        return {
            x: data.a * Math.cos(theta),
            y: data.a * Math.sin(theta),
            z: 0,
            vx: -data.a * n * Math.sin(theta),
            vy: data.a * n * Math.cos(theta),
            vz: 0
        };
    }

    // Solve Lambert problem (simplified Battin method)
    function solveLambert(r1, r2, tof, GM, isPrograde = true) {
        const r1Mag = Math.sqrt(r1.x*r1.x + r1.y*r1.y + r1.z*r1.z);
        const r2Mag = Math.sqrt(r2.x*r2.x + r2.y*r2.y + r2.z*r2.z);

        // Angle between position vectors
        const cosTA = (r1.x*r2.x + r1.y*r2.y + r1.z*r2.z) / (r1Mag * r2Mag);
        let sinTA = Math.sqrt(1 - cosTA*cosTA);
        if (!isPrograde) sinTA = -sinTA;

        const A = Math.sqrt(r1Mag * r2Mag * (1 + cosTA));
        if (Math.abs(A) < 1e-10) return null;

        // Iterative solution for z parameter
        let z = 0;
        const maxIter = 100;
        const tol = 1e-8;

        for (let iter = 0; iter < maxIter; iter++) {
            let C, S;
            if (z > 0) {
                const sqrtZ = Math.sqrt(z);
                C = (1 - Math.cos(sqrtZ)) / z;
                S = (sqrtZ - Math.sin(sqrtZ)) / Math.pow(sqrtZ, 3);
            } else if (z < 0) {
                const sqrtNZ = Math.sqrt(-z);
                C = (1 - Math.cosh(sqrtNZ)) / z;
                S = (Math.sinh(sqrtNZ) - sqrtNZ) / Math.pow(sqrtNZ, 3);
            } else {
                C = 0.5;
                S = 1/6;
            }

            const y = r1Mag + r2Mag + A * (z*S - 1) / Math.sqrt(C);
            if (y < 0) {
                z = z * 0.5;
                continue;
            }

            const x = Math.sqrt(y / C);
            const tCalc = (x*x*x * S + A * Math.sqrt(y)) / Math.sqrt(GM);

            const error = tCalc - tof;
            if (Math.abs(error) < tol * tof) break;

            // Newton-Raphson update
            const dTdz = (x*x*x * (C - 1.5*S/C) / (2*C) + A/8 * (3*S*Math.sqrt(y)/C + A/x)) / Math.sqrt(GM);
            z = z - error / dTdz;
        }

        // Compute velocities
        let C, S;
        if (z > 0) {
            const sqrtZ = Math.sqrt(z);
            C = (1 - Math.cos(sqrtZ)) / z;
            S = (sqrtZ - Math.sin(sqrtZ)) / Math.pow(sqrtZ, 3);
        } else {
            C = 0.5;
            S = 1/6;
        }

        const y = r1Mag + r2Mag + A * (z*S - 1) / Math.sqrt(C);
        const f = 1 - y / r1Mag;
        const g = A * Math.sqrt(y / GM);
        const gDot = 1 - y / r2Mag;

        const v1 = {
            x: (r2.x - f * r1.x) / g,
            y: (r2.y - f * r1.y) / g,
            z: (r2.z - f * r1.z) / g
        };

        const v2 = {
            x: (gDot * r2.x - r1.x) / g,
            y: (gDot * r2.y - r1.y) / g,
            z: (gDot * r2.z - r1.z) / g
        };

        return { v1, v2 };
    }

    // Compute gravity assist turn angle and delta-V
    function gravityAssist(vInf_in, vInf_out, planet) {
        const data = planets[planet];

        // V-infinity magnitudes
        const vInfInMag = Math.sqrt(vInf_in.x**2 + vInf_in.y**2 + vInf_in.z**2);
        const vInfOutMag = Math.sqrt(vInf_out.x**2 + vInf_out.y**2 + vInf_out.z**2);

        // Turn angle
        const dot = vInf_in.x*vInf_out.x + vInf_in.y*vInf_out.y + vInf_in.z*vInf_out.z;
        const cosAngle = dot / (vInfInMag * vInfOutMag);
        const turnAngle = Math.acos(Math.max(-1, Math.min(1, cosAngle)));

        // Minimum flyby radius for this turn angle (unpowered)
        // delta = 2 * arcsin(1 / (1 + rp * v_inf^2 / GM))
        // Solving for rp: rp = GM / v_inf^2 * (1/sin(delta/2) - 1)
        const avgVInf = (vInfInMag + vInfOutMag) / 2;
        const sinHalfTurn = Math.sin(turnAngle / 2);
        const minRp = data.GM / (avgVInf * avgVInf) * (1 / sinHalfTurn - 1);

        // Check if flyby is possible (rp > planet radius)
        const isPossible = minRp > data.R * 1.1;  // 10% safety margin

        // If V-infinity magnitude changes, powered flyby needed
        const deltaVInf = Math.abs(vInfOutMag - vInfInMag);

        return {
            turnAngle: turnAngle * 180 / Math.PI,
            minRp: minRp,
            altitude: (minRp - data.R) / 1000,  // km
            isPossible,
            deltaVInf
        };
    }

    // Evaluate trajectory for given times of flight
    function evaluateTrajectory(launchEpoch, tofs) {
        let totalDeltaV = 0;
        const legs = [];

        let t = launchEpoch;

        for (let i = 0; i < sequence.length - 1; i++) {
            const departPlanet = sequence[i];
            const arrivePlanet = sequence[i + 1];
            const tof = tofs[i] * JULIAN_DAY;

            // Get planet positions
            const r1 = getPlanetPosition(departPlanet, t);
            const r2 = getPlanetPosition(arrivePlanet, t + tof);

            // Solve Lambert problem
            const lambert = solveLambert(r1, r2, tof, GM_SUN);
            if (!lambert) return { totalDeltaV: Infinity, legs: [] };

            const leg = {
                from: departPlanet,
                to: arrivePlanet,
                tof: tofs[i],
                vDepart: lambert.v1,
                vArrive: lambert.v2
            };

            // Compute delta-V for this leg
            if (i === 0) {
                // Launch from Earth: delta-V is v_infinity at departure
                const vPlanet = { x: r1.vx, y: r1.vy, z: r1.vz };
                const vInf = {
                    x: lambert.v1.x - vPlanet.x,
                    y: lambert.v1.y - vPlanet.y,
                    z: lambert.v1.z - vPlanet.z
                };
                const vInfMag = Math.sqrt(vInf.x**2 + vInf.y**2 + vInf.z**2);

                // C3 = v_infinity^2
                const C3 = vInfMag * vInfMag;
                leg.C3 = C3;
                leg.vInfLaunch = vInfMag;

                // Delta-V from 200 km parking orbit
                const parkingRadius = planets.Earth.R + 200e3;
                const vPark = Math.sqrt(planets.Earth.GM / parkingRadius);
                const vEscape = Math.sqrt(vPark*vPark + vInfMag*vInfMag);
                leg.deltaV = vEscape - vPark;
                totalDeltaV += leg.deltaV;
            } else {
                // Gravity assist at intermediate planet
                const prevLeg = legs[i - 1];
                const vPlanet = getPlanetPosition(departPlanet, t);

                // V-infinity in (from previous leg)
                const vInf_in = {
                    x: prevLeg.vArrive.x - vPlanet.vx,
                    y: prevLeg.vArrive.y - vPlanet.vy,
                    z: prevLeg.vArrive.z - vPlanet.vz
                };

                // V-infinity out (to next leg)
                const vInf_out = {
                    x: lambert.v1.x - vPlanet.vx,
                    y: lambert.v1.y - vPlanet.vy,
                    z: lambert.v1.z - vPlanet.vz
                };

                const assist = gravityAssist(vInf_in, vInf_out, departPlanet);
                leg.flyby = assist;

                // Add powered flyby delta-V if needed
                if (assist.deltaVInf > 100) {  // More than 100 m/s change
                    leg.deltaV = assist.deltaVInf;
                    totalDeltaV += assist.deltaVInf;
                } else {
                    leg.deltaV = 50;  // TCM
                    totalDeltaV += 50;
                }
            }

            // Saturn orbit insertion
            if (arrivePlanet === 'Saturn' && i === sequence.length - 2) {
                const vPlanet = { x: r2.vx, y: r2.vy, z: r2.vz };
                const vInf = {
                    x: lambert.v2.x - vPlanet.x,
                    y: lambert.v2.y - vPlanet.y,
                    z: lambert.v2.z - vPlanet.z
                };
                const vInfMag = Math.sqrt(vInf.x**2 + vInf.y**2 + vInf.z**2);
                leg.vInfArrive = vInfMag;

                // Saturn orbit insertion (into 1 million km orbit)
                const captureRadius = 1e9;
                const vCapture = Math.sqrt(planets.Saturn.GM / captureRadius);
                const vArrival = Math.sqrt(vCapture*vCapture + vInfMag*vInfMag);
                const soiDeltaV = vArrival - vCapture;
                leg.soiDeltaV = soiDeltaV;
                totalDeltaV += soiDeltaV;
            }

            legs.push(leg);
            t += tof;
        }

        return { totalDeltaV, legs };
    }

    // Grid search optimization
    console.log('\n--- Running Grid Search Optimization ---\n');

    // Search space (times of flight in days)
    const tofRanges = [
        [150, 250],   // Earth -> Venus
        [200, 400],   // Venus -> Venus (resonant)
        [100, 300],   // Venus -> Earth
        [600, 1200],  // Earth -> Jupiter
        [900, 1800]   // Jupiter -> Saturn
    ];

    const candidates = [
        { name: 'A', tofs: [200, 300, 200, 800, 1200] },
        { name: 'B', tofs: [180, 350, 180, 900, 1400] },
        { name: 'C', tofs: [220, 280, 220, 1000, 1100] },
        { name: 'D', tofs: [193, 424, 55, 500, 1280] },  // Cassini actual (approx)
    ];

    let bestResult = null;
    let bestDeltaV = Infinity;

    for (const candidate of candidates) {
        const result = evaluateTrajectory(0, candidate.tofs);
        const totalTof = candidate.tofs.reduce((a, b) => a + b, 0);

        console.log(`Candidate ${candidate.name}:`);
        console.log(`  TOFs: [${candidate.tofs.join(', ')}] days`);
        console.log(`  Total flight time: ${(totalTof / 365.25).toFixed(2)} years`);
        console.log(`  Total delta-V: ${(result.totalDeltaV / 1000).toFixed(3)} km/s`);

        if (result.totalDeltaV < bestDeltaV) {
            bestDeltaV = result.totalDeltaV;
            bestResult = { ...result, candidate };
        }

        // Print flyby details
        for (const leg of result.legs) {
            if (leg.flyby) {
                console.log(`    ${leg.from} flyby: turn angle = ${leg.flyby.turnAngle.toFixed(1)}°, ` +
                    `altitude = ${leg.flyby.altitude.toFixed(0)} km`);
            }
        }
        console.log('');
    }

    // Best result summary
    console.log('=== Optimization Result ===');
    console.log(`\nBest candidate: ${bestResult.candidate.name}`);
    console.log(`Total delta-V: ${(bestDeltaV / 1000).toFixed(3)} km/s`);

    console.log('\nLeg Details:');
    let cumTime = 0;
    for (let i = 0; i < bestResult.legs.length; i++) {
        const leg = bestResult.legs[i];
        cumTime += bestResult.candidate.tofs[i];
        console.log(`  ${leg.from} -> ${leg.to}:`);
        console.log(`    TOF: ${bestResult.candidate.tofs[i]} days`);
        console.log(`    Cumulative: ${(cumTime / 365.25).toFixed(2)} years`);
        console.log(`    Delta-V: ${(leg.deltaV / 1000).toFixed(3)} km/s`);
        if (leg.C3) console.log(`    C3: ${(leg.C3 / 1e6).toFixed(2)} km²/s²`);
        if (leg.soiDeltaV) console.log(`    SOI: ${(leg.soiDeltaV / 1000).toFixed(3)} km/s`);
    }

    // Actual Cassini comparison
    console.log('\n--- Actual Cassini-Huygens Mission ---');
    console.log('Launch: October 15, 1997');
    console.log('Trajectory: VVEJGA (Venus-Venus-Earth-Jupiter GA)');
    console.log('Saturn arrival: July 1, 2004');
    console.log('Total flight time: 6.7 years');
    console.log('Launch C3: 16.6 km²/s²');
    console.log('Total delta-V: ~2.0 km/s');

    console.log('\n=== Cassini-1 MGA Optimization Complete ===');
}

main().catch(console.error);
