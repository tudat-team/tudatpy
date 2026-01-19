/**
 * Example 07: Orbital Element Conversions
 *
 * This example demonstrates all the orbital element conversion functions
 * available in Tudat WASM, including Keplerian, Cartesian, and other
 * element sets.
 *
 * Run with: node 07_orbital_elements.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Orbital Element Conversions ===\n');

    const tudat = await createTudatModule();

    const earthGM = 3.986004418e14;  // m³/s²
    const earthRadius = 6378137.0;   // m

    // ========================================
    // Keplerian to Cartesian Conversion
    // ========================================
    console.log('1. Keplerian to Cartesian Conversion');
    console.log('=====================================\n');

    // Define various orbit types
    const orbits = [
        {
            name: 'ISS (LEO circular)',
            a: earthRadius + 400e3,
            e: 0.0001,
            i: 51.6,
            aop: 0,
            raan: 0,
            ta: 0
        },
        {
            name: 'GTO (elliptical)',
            a: (earthRadius + 200e3 + earthRadius + 35786e3) / 2,
            e: 0.73,
            i: 28.5,
            aop: 180,
            raan: 45,
            ta: 0
        },
        {
            name: 'Molniya (highly elliptical)',
            a: 26600e3,
            e: 0.74,
            i: 63.4,
            aop: 270,
            raan: 90,
            ta: 0
        },
        {
            name: 'Sun-synchronous (polar)',
            a: earthRadius + 800e3,
            e: 0.001,
            i: 98.6,
            aop: 0,
            raan: 0,
            ta: 45
        }
    ];

    for (const orbit of orbits) {
        console.log(`${orbit.name}:`);
        console.log(`  Keplerian: a=${(orbit.a/1000).toFixed(1)} km, e=${orbit.e}, i=${orbit.i}°`);

        const kepler = new tudat.Vector6d();
        kepler.set(0, orbit.a);
        kepler.set(1, orbit.e);
        kepler.set(2, orbit.i * Math.PI / 180);
        kepler.set(3, orbit.aop * Math.PI / 180);
        kepler.set(4, orbit.raan * Math.PI / 180);
        kepler.set(5, orbit.ta * Math.PI / 180);

        const cartesian = tudat.astro.element_conversion.keplerian_to_cartesian(kepler, earthGM);

        const pos = [cartesian.get(0), cartesian.get(1), cartesian.get(2)];
        const vel = [cartesian.get(3), cartesian.get(4), cartesian.get(5)];
        const r = Math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2);
        const v = Math.sqrt(vel[0]**2 + vel[1]**2 + vel[2]**2);

        console.log(`  Cartesian: r=${(r/1000).toFixed(1)} km, v=${(v/1000).toFixed(3)} km/s`);
        console.log(`  Position: [${(pos[0]/1000).toFixed(1)}, ${(pos[1]/1000).toFixed(1)}, ${(pos[2]/1000).toFixed(1)}] km`);
        console.log(`  Velocity: [${vel[0].toFixed(1)}, ${vel[1].toFixed(1)}, ${vel[2].toFixed(1)}] m/s`);
        console.log('');

        kepler.delete();
        cartesian.delete();
    }

    // ========================================
    // Cartesian to Keplerian Conversion
    // ========================================
    console.log('\n2. Cartesian to Keplerian Conversion');
    console.log('=====================================\n');

    // Create a state and convert back
    const testState = new tudat.Vector6d();
    testState.set(0, 6800000);   // x [m]
    testState.set(1, 1000000);   // y [m]
    testState.set(2, 500000);    // z [m]
    testState.set(3, -500);      // vx [m/s]
    testState.set(4, 7200);      // vy [m/s]
    testState.set(5, 1000);      // vz [m/s]

    console.log('Input Cartesian state:');
    console.log(`  Position: [${testState.get(0)/1000}, ${testState.get(1)/1000}, ${testState.get(2)/1000}] km`);
    console.log(`  Velocity: [${testState.get(3)}, ${testState.get(4)}, ${testState.get(5)}] m/s`);

    const keplerOut = tudat.astro.element_conversion.cartesian_to_keplerian(testState, earthGM);

    console.log('\nOutput Keplerian elements:');
    console.log(`  Semi-major axis: ${(keplerOut.get(0)/1000).toFixed(2)} km`);
    console.log(`  Eccentricity: ${keplerOut.get(1).toFixed(6)}`);
    console.log(`  Inclination: ${(keplerOut.get(2) * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  Arg of periapsis: ${(keplerOut.get(3) * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  RAAN: ${(keplerOut.get(4) * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  True anomaly: ${(keplerOut.get(5) * 180 / Math.PI).toFixed(4)}°`);

    testState.delete();
    keplerOut.delete();

    // ========================================
    // Anomaly Conversions
    // ========================================
    console.log('\n\n3. Anomaly Conversions');
    console.log('======================\n');

    const eccentricities = [0.0, 0.1, 0.5, 0.9];

    console.log('True Anomaly → Eccentric Anomaly → Mean Anomaly:');
    console.log('------------------------------------------------');

    for (const e of eccentricities) {
        console.log(`\nEccentricity = ${e}:`);

        const trueAnomalies = [0, 45, 90, 135, 180];  // degrees

        for (const taDeg of trueAnomalies) {
            const ta = taDeg * Math.PI / 180;

            // True to Eccentric anomaly
            const E = tudat.astro.element_conversion.true_to_eccentric_anomaly(ta, e);

            // Eccentric to Mean anomaly
            const M = tudat.astro.element_conversion.eccentric_to_mean_anomaly(E, e);

            console.log(`  ν=${taDeg.toString().padStart(3)}° → E=${(E * 180 / Math.PI).toFixed(2).padStart(7)}° → M=${(M * 180 / Math.PI).toFixed(2).padStart(7)}°`);
        }
    }

    console.log('\n\nMean Anomaly → Eccentric Anomaly (Kepler equation):');
    console.log('---------------------------------------------------');

    for (const e of [0.1, 0.5, 0.9]) {
        console.log(`\nEccentricity = ${e}:`);

        const meanAnomalies = [0, 45, 90, 135, 180];  // degrees

        for (const mDeg of meanAnomalies) {
            const M = mDeg * Math.PI / 180;

            // Mean to Eccentric anomaly (iterative solution)
            const E = tudat.astro.element_conversion.mean_to_eccentric_anomaly(e, M);

            // Back to true anomaly
            const ta = tudat.astro.element_conversion.eccentric_to_true_anomaly(E, e);

            console.log(`  M=${mDeg.toString().padStart(3)}° → E=${(E * 180 / Math.PI).toFixed(2).padStart(7)}° → ν=${(ta * 180 / Math.PI).toFixed(2).padStart(7)}°`);
        }
    }

    // ========================================
    // Spherical Coordinates
    // ========================================
    console.log('\n\n4. Coordinate System Conversions');
    console.log('=================================\n');

    console.log('Cartesian to Spherical:');

    const cartesianPoints = [
        [1, 0, 0],
        [0, 1, 0],
        [0, 0, 1],
        [1, 1, 1],
        [6378137, 0, 0]  // Point on equator
    ];

    for (const cart of cartesianPoints) {
        const cartVec = new tudat.Vector3d();
        cartVec.set(0, cart[0]);
        cartVec.set(1, cart[1]);
        cartVec.set(2, cart[2]);

        const spherical = tudat.astro.frame_conversion.cartesian_to_spherical(cartVec);

        const r = spherical.get(0);
        const lat = spherical.get(1) * 180 / Math.PI;
        const lon = spherical.get(2) * 180 / Math.PI;

        console.log(`  [${cart[0]}, ${cart[1]}, ${cart[2]}] → r=${r.toFixed(2)}, lat=${lat.toFixed(2)}°, lon=${lon.toFixed(2)}°`);

        cartVec.delete();
        spherical.delete();
    }

    // ========================================
    // Orbital Energy and Velocity
    // ========================================
    console.log('\n\n5. Orbital Energy and Velocity');
    console.log('===============================\n');

    const altitudes = [200, 400, 800, 35786, 384400];  // km
    const names = ['LEO-low', 'ISS', 'LEO-high', 'GEO', 'Moon'];

    console.log('Circular orbit velocities and energies:');
    console.log('---------------------------------------');

    for (let i = 0; i < altitudes.length; i++) {
        const r = earthRadius + altitudes[i] * 1000;
        const v_circ = Math.sqrt(earthGM / r);
        const v_esc = Math.sqrt(2 * earthGM / r);
        const energy = -earthGM / (2 * r);
        const period = 2 * Math.PI * Math.sqrt(r**3 / earthGM);

        console.log(`\n${names[i]} (${altitudes[i]} km altitude):`);
        console.log(`  Orbital velocity: ${(v_circ/1000).toFixed(3)} km/s`);
        console.log(`  Escape velocity: ${(v_esc/1000).toFixed(3)} km/s`);
        console.log(`  Specific energy: ${(energy/1e6).toFixed(3)} MJ/kg`);

        if (period < 86400) {
            console.log(`  Period: ${(period/60).toFixed(2)} min`);
        } else {
            console.log(`  Period: ${(period/86400).toFixed(2)} days`);
        }
    }

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
