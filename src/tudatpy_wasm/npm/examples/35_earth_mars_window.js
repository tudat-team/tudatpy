/**
 * Example 35: Earth-Mars Transfer Window Analysis
 *
 * Ported from: examples/tudatpy/mission_design/earth_mars_transfer_window.py
 *
 * This example demonstrates analysis of Earth-Mars transfer windows,
 * computing porkchop plots for launch/arrival date combinations.
 *
 * Key concepts:
 * - Transfer window analysis
 * - Porkchop plots
 * - Synodic period
 * - Lambert problem grid search
 * - C3 and arrival V-infinity
 *
 * Run with: node 35_earth_mars_window.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Earth-Mars Transfer Window Analysis ===\n');

    const tudat = await createTudatModule();

    // Constants
    const AU = 1.496e11;                // Astronomical unit [m]
    const GM_SUN = 1.32712440018e20;    // Sun gravitational parameter [m^3/s^2]
    const JULIAN_DAY = 86400;
    const JULIAN_YEAR = 365.25 * JULIAN_DAY;

    // Planetary orbital elements (simplified circular coplanar)
    const earth = {
        a: 1.000 * AU,
        T: 365.25 * JULIAN_DAY,  // Orbital period [s]
        GM: 3.986e14
    };

    const mars = {
        a: 1.524 * AU,
        T: 687 * JULIAN_DAY,  // Orbital period [s]
        GM: 4.283e13
    };

    // Synodic period (time between transfer windows)
    const synodicPeriod = 1 / Math.abs(1/earth.T - 1/mars.T);
    console.log('Orbital Parameters:');
    console.log(`  Earth: a = ${(earth.a / AU).toFixed(3)} AU, T = ${(earth.T / JULIAN_DAY).toFixed(1)} days`);
    console.log(`  Mars:  a = ${(mars.a / AU).toFixed(3)} AU, T = ${(mars.T / JULIAN_DAY).toFixed(1)} days`);
    console.log(`  Synodic period: ${(synodicPeriod / JULIAN_DAY).toFixed(1)} days (${(synodicPeriod / JULIAN_YEAR).toFixed(2)} years)`);

    // Get planet position at time t
    function getPlanetPosition(planet, t) {
        const n = 2 * Math.PI / planet.T;  // Mean motion
        const theta = n * t;
        const r = planet.a;
        return {
            x: r * Math.cos(theta),
            y: r * Math.sin(theta),
            z: 0,
            vx: -r * n * Math.sin(theta),
            vy: r * n * Math.cos(theta),
            vz: 0
        };
    }

    // Simplified Lambert solver
    function solveLambert(r1, r2, tof, GM) {
        const r1Mag = Math.sqrt(r1.x**2 + r1.y**2 + r1.z**2);
        const r2Mag = Math.sqrt(r2.x**2 + r2.y**2 + r2.z**2);

        // Transfer angle
        const cosTheta = (r1.x*r2.x + r1.y*r2.y) / (r1Mag * r2Mag);
        const sinTheta = (r1.x*r2.y - r1.y*r2.x) / (r1Mag * r2Mag);
        const theta = Math.atan2(sinTheta, cosTheta);

        // Chord length
        const c = Math.sqrt(r1Mag**2 + r2Mag**2 - 2*r1Mag*r2Mag*cosTheta);

        // Semi-perimeter
        const s = (r1Mag + r2Mag + c) / 2;

        // Minimum energy transfer
        const aMin = s / 2;
        const tMin = Math.sqrt(2) / 4 * Math.sqrt(s**3 / GM) * (Math.PI - Math.sqrt(1 - c/s) + Math.asin(Math.sqrt(c/s)));

        // Use iterative method for actual TOF
        // Simplified: estimate semi-major axis from TOF
        const n_transfer = Math.PI / tof * 2;  // Approximate
        const a_est = Math.pow(GM / (n_transfer**2), 1/3);
        const a = Math.max(aMin, a_est);

        // Solve for transfer orbit parameters
        const p = (4 * a * (s - r1Mag) * (s - r2Mag) / (c**2)) * Math.pow(Math.sin((theta)/2), 2);

        // Lagrange coefficients
        const f = 1 - r2Mag / p * (1 - cosTheta);
        const g = r1Mag * r2Mag * sinTheta / Math.sqrt(GM * p);
        const fDot = Math.sqrt(GM / p) * Math.tan(theta/2) * ((1 - cosTheta)/p - 1/r1Mag - 1/r2Mag);
        const gDot = 1 - r1Mag / p * (1 - cosTheta);

        // Velocities
        const v1 = {
            x: (r2.x - f * r1.x) / g,
            y: (r2.y - f * r1.y) / g,
            z: 0
        };

        const v2 = {
            x: fDot * r1.x + gDot * v1.x,
            y: fDot * r1.y + gDot * v1.y,
            z: 0
        };

        return { v1, v2, a };
    }

    // Compute transfer parameters for given launch and arrival epochs
    function computeTransfer(launchEpoch, arrivalEpoch) {
        const tof = arrivalEpoch - launchEpoch;
        if (tof <= 0) return null;

        const r1 = getPlanetPosition(earth, launchEpoch);
        const r2 = getPlanetPosition(mars, arrivalEpoch);

        try {
            const lambert = solveLambert(r1, r2, tof, GM_SUN);
            if (!lambert) return null;

            // Departure V-infinity (relative to Earth)
            const vInfDep = {
                x: lambert.v1.x - r1.vx,
                y: lambert.v1.y - r1.vy,
                z: 0
            };
            const vInfDepMag = Math.sqrt(vInfDep.x**2 + vInfDep.y**2);

            // C3 (characteristic energy)
            const C3 = vInfDepMag ** 2;

            // Arrival V-infinity (relative to Mars)
            const vMars = getPlanetPosition(mars, arrivalEpoch);
            const vInfArr = {
                x: lambert.v2.x - vMars.vx,
                y: lambert.v2.y - vMars.vy,
                z: 0
            };
            const vInfArrMag = Math.sqrt(vInfArr.x**2 + vInfArr.y**2);

            // Delta-V for Earth departure (from 200 km parking orbit)
            const parkingRadius = 6578e3;
            const vPark = Math.sqrt(earth.GM / parkingRadius);
            const vEscape = Math.sqrt(vPark**2 + vInfDepMag**2);
            const dvDepart = vEscape - vPark;

            // Delta-V for Mars orbit insertion (into 400 km orbit)
            const marsOrbitRadius = 3796e3;
            const vMarsOrbit = Math.sqrt(mars.GM / marsOrbitRadius);
            const vMarsArrive = Math.sqrt(vMarsOrbit**2 + vInfArrMag**2);
            const dvArrive = vMarsArrive - vMarsOrbit;

            return {
                tof: tof / JULIAN_DAY,
                C3: C3 / 1e6,  // km²/s²
                vInfDep: vInfDepMag / 1000,  // km/s
                vInfArr: vInfArrMag / 1000,  // km/s
                dvDepart: dvDepart / 1000,   // km/s
                dvArrive: dvArrive / 1000,   // km/s
                dvTotal: (dvDepart + dvArrive) / 1000,
                a: lambert.a / AU  // AU
            };
        } catch (e) {
            return null;
        }
    }

    // Generate porkchop plot data
    console.log('\n--- Generating Porkchop Plot Data ---\n');

    // Define search window (one synodic period)
    const startLaunch = 0;
    const endLaunch = synodicPeriod;
    const launchStep = 10 * JULIAN_DAY;

    // TOF range
    const minTOF = 150 * JULIAN_DAY;
    const maxTOF = 400 * JULIAN_DAY;
    const tofStep = 10 * JULIAN_DAY;

    const porkchopData = [];
    let minC3 = Infinity;
    let minC3Launch = 0, minC3Arrival = 0;
    let minTotal = Infinity;
    let minTotalLaunch = 0, minTotalArrival = 0;

    for (let launch = startLaunch; launch <= endLaunch; launch += launchStep) {
        for (let tof = minTOF; tof <= maxTOF; tof += tofStep) {
            const arrival = launch + tof;
            const result = computeTransfer(launch, arrival);

            if (result && result.C3 > 0 && result.C3 < 100) {
                porkchopData.push({
                    launch: launch / JULIAN_DAY,
                    arrival: arrival / JULIAN_DAY,
                    ...result
                });

                if (result.C3 < minC3) {
                    minC3 = result.C3;
                    minC3Launch = launch;
                    minC3Arrival = arrival;
                }

                if (result.dvTotal < minTotal) {
                    minTotal = result.dvTotal;
                    minTotalLaunch = launch;
                    minTotalArrival = arrival;
                }
            }
        }
    }

    console.log(`Computed ${porkchopData.length} transfer opportunities`);

    // Analyze optimal transfers
    console.log('\n=== Optimal Transfer Analysis ===');

    console.log('\nMinimum C3 (departure energy):');
    const minC3Transfer = computeTransfer(minC3Launch, minC3Arrival);
    console.log(`  Launch: day ${(minC3Launch / JULIAN_DAY).toFixed(0)}`);
    console.log(`  Arrival: day ${(minC3Arrival / JULIAN_DAY).toFixed(0)}`);
    console.log(`  TOF: ${minC3Transfer.tof.toFixed(0)} days`);
    console.log(`  C3: ${minC3Transfer.C3.toFixed(2)} km²/s²`);
    console.log(`  V∞ departure: ${minC3Transfer.vInfDep.toFixed(2)} km/s`);
    console.log(`  V∞ arrival: ${minC3Transfer.vInfArr.toFixed(2)} km/s`);
    console.log(`  Total ΔV: ${minC3Transfer.dvTotal.toFixed(2)} km/s`);

    console.log('\nMinimum total ΔV:');
    const minDvTransfer = computeTransfer(minTotalLaunch, minTotalArrival);
    console.log(`  Launch: day ${(minTotalLaunch / JULIAN_DAY).toFixed(0)}`);
    console.log(`  Arrival: day ${(minTotalArrival / JULIAN_DAY).toFixed(0)}`);
    console.log(`  TOF: ${minDvTransfer.tof.toFixed(0)} days`);
    console.log(`  C3: ${minDvTransfer.C3.toFixed(2)} km²/s²`);
    console.log(`  V∞ departure: ${minDvTransfer.vInfDep.toFixed(2)} km/s`);
    console.log(`  V∞ arrival: ${minDvTransfer.vInfArr.toFixed(2)} km/s`);
    console.log(`  Total ΔV: ${minDvTransfer.dvTotal.toFixed(2)} km/s`);

    // Hohmann transfer comparison
    console.log('\n--- Hohmann Transfer Reference ---');
    const hohmannA = (earth.a + mars.a) / 2;
    const hohmannTOF = Math.PI * Math.sqrt(hohmannA**3 / GM_SUN);
    const v1 = Math.sqrt(GM_SUN / earth.a);
    const vT1 = Math.sqrt(GM_SUN * (2/earth.a - 1/hohmannA));
    const v2 = Math.sqrt(GM_SUN / mars.a);
    const vT2 = Math.sqrt(GM_SUN * (2/mars.a - 1/hohmannA));
    const hohmannDV = Math.abs(vT1 - v1) + Math.abs(v2 - vT2);

    console.log(`Hohmann transfer:`);
    console.log(`  TOF: ${(hohmannTOF / JULIAN_DAY).toFixed(0)} days`);
    console.log(`  ΔV1: ${((vT1 - v1) / 1000).toFixed(2)} km/s`);
    console.log(`  ΔV2: ${((v2 - vT2) / 1000).toFixed(2)} km/s`);
    console.log(`  Total: ${(hohmannDV / 1000).toFixed(2)} km/s`);

    // Sample porkchop data
    console.log('\n--- Porkchop Plot Sample (C3 values) ---');
    console.log('Launch[day] | TOF[day] | C3[km²/s²] | V∞arr[km/s] | ΔVtot[km/s]');
    console.log('------------|----------|------------|-------------|------------');

    // Sample at regular intervals
    const sampleInterval = Math.floor(porkchopData.length / 15);
    for (let i = 0; i < porkchopData.length; i += sampleInterval) {
        const d = porkchopData[i];
        console.log(`${d.launch.toFixed(0).padStart(11)} | ${d.tof.toFixed(0).padStart(8)} | ${d.C3.toFixed(2).padStart(10)} | ${d.vInfArr.toFixed(2).padStart(11)} | ${d.dvTotal.toFixed(2).padStart(11)}`);
    }

    // Transfer window characteristics
    console.log('\n--- Transfer Window Characteristics ---');

    // Find range of good opportunities (C3 < 15 km²/s²)
    const goodTransfers = porkchopData.filter(d => d.C3 < 15);
    if (goodTransfers.length > 0) {
        const firstGood = goodTransfers[0];
        const lastGood = goodTransfers[goodTransfers.length - 1];
        console.log(`\nTransfers with C3 < 15 km²/s²:`);
        console.log(`  Count: ${goodTransfers.length}`);
        console.log(`  Launch window: day ${firstGood.launch.toFixed(0)} - ${lastGood.launch.toFixed(0)}`);
        console.log(`  Window duration: ~${((lastGood.launch - firstGood.launch)).toFixed(0)} days`);
    }

    // Type I vs Type II transfers
    const typeI = porkchopData.filter(d => d.tof < 250);
    const typeII = porkchopData.filter(d => d.tof >= 250);
    console.log(`\nTransfer types:`);
    console.log(`  Type I (short, < 250 days): ${typeI.length} opportunities`);
    console.log(`  Type II (long, >= 250 days): ${typeII.length} opportunities`);

    if (typeI.length > 0) {
        const bestTypeI = typeI.reduce((a, b) => a.dvTotal < b.dvTotal ? a : b);
        console.log(`  Best Type I: TOF=${bestTypeI.tof.toFixed(0)}d, ΔV=${bestTypeI.dvTotal.toFixed(2)} km/s`);
    }
    if (typeII.length > 0) {
        const bestTypeII = typeII.reduce((a, b) => a.dvTotal < b.dvTotal ? a : b);
        console.log(`  Best Type II: TOF=${bestTypeII.tof.toFixed(0)}d, ΔV=${bestTypeII.dvTotal.toFixed(2)} km/s`);
    }

    // Mission design notes
    console.log('\n--- Mission Design Notes ---');
    console.log('Earth-Mars transfer windows occur every ~26 months');
    console.log('Typical characteristics:');
    console.log('  - Type I: shorter, higher C3, less efficient but faster');
    console.log('  - Type II: longer, lower C3, more efficient but slower');
    console.log('  - Optimal window width: ~30-60 days for C3 < 15 km²/s²');
    console.log('  - Mars arrival velocity: 2-4 km/s typical');

    console.log('\n=== Earth-Mars Transfer Window Analysis Complete ===');
}

main().catch(console.error);
