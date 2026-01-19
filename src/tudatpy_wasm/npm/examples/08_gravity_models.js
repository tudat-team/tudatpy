/**
 * Example 08: Gravity Models and Spherical Harmonics
 *
 * This example demonstrates gravity field modeling including point mass,
 * spherical harmonics, and third-body perturbations.
 *
 * Run with: node 08_gravity_models.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Gravity Models and Spherical Harmonics ===\n');

    const tudat = await createTudatModule();

    const earthGM = 3.986004418e14;  // m³/s²
    const earthRadius = 6378137.0;   // m

    // Earth J2-J6 coefficients (unnormalized)
    const earthJ2 = 1.08262668e-3;
    const earthJ3 = -2.53265648e-6;
    const earthJ4 = -1.61962159e-6;
    const earthJ5 = -2.27296082e-7;
    const earthJ6 = 5.40681239e-7;

    // ========================================
    // Point Mass Gravity
    // ========================================
    console.log('1. Point Mass Gravity Model');
    console.log('===========================\n');

    const altitudes = [200, 400, 800, 2000, 35786];  // km

    console.log('Gravitational acceleration at different altitudes:');
    console.log('--------------------------------------------------');

    for (const alt of altitudes) {
        const r = earthRadius + alt * 1000;
        const g = earthGM / (r * r);

        console.log(`  ${alt.toString().padStart(5)} km altitude: g = ${g.toFixed(4)} m/s² (${(g/9.81*100).toFixed(2)}% of surface)`);
    }

    // ========================================
    // J2 Perturbation Effects
    // ========================================
    console.log('\n\n2. J2 Perturbation Effects');
    console.log('==========================\n');

    console.log('The J2 coefficient causes:');
    console.log('  - Nodal regression (RAAN drift)');
    console.log('  - Apsidal precession (argument of periapsis drift)');
    console.log('  - Mean motion variation\n');

    // Calculate J2 effects for different orbits
    const orbits = [
        { name: 'ISS', alt: 400, inc: 51.6 },
        { name: 'Sun-sync', alt: 800, inc: 98.6 },
        { name: 'GPS', alt: 20200, inc: 55.0 },
        { name: 'GEO', alt: 35786, inc: 0.0 }
    ];

    console.log('RAAN and Argument of Periapsis drift rates:');
    console.log('-------------------------------------------');

    for (const orbit of orbits) {
        const a = earthRadius + orbit.alt * 1000;
        const n = Math.sqrt(earthGM / Math.pow(a, 3));  // mean motion
        const i = orbit.inc * Math.PI / 180;

        // RAAN drift rate (rad/s)
        const raanDot = -1.5 * n * earthJ2 * Math.pow(earthRadius / a, 2) * Math.cos(i);

        // Argument of periapsis drift rate (rad/s)
        const aopDot = 0.75 * n * earthJ2 * Math.pow(earthRadius / a, 2) * (5 * Math.cos(i)**2 - 1);

        // Convert to deg/day
        const raanDotDegDay = raanDot * 180 / Math.PI * 86400;
        const aopDotDegDay = aopDot * 180 / Math.PI * 86400;

        console.log(`\n  ${orbit.name} (${orbit.alt} km, ${orbit.inc}° inc):`);
        console.log(`    RAAN drift: ${raanDotDegDay.toFixed(4)}°/day`);
        console.log(`    AoP drift: ${aopDotDegDay.toFixed(4)}°/day`);

        // Sun-synchronous condition check
        const requiredRaanDot = 360.0 / 365.25;  // deg/day for sun-sync
        if (Math.abs(raanDotDegDay + requiredRaanDot) < 0.1) {
            console.log(`    → Sun-synchronous orbit!`);
        }
    }

    // ========================================
    // Higher Order Harmonics
    // ========================================
    console.log('\n\n3. Higher Order Spherical Harmonics');
    console.log('====================================\n');

    console.log('Earth zonal harmonic coefficients:');
    console.log(`  J2 = ${earthJ2.toExponential(6)} (dominant oblateness)`);
    console.log(`  J3 = ${earthJ3.toExponential(6)} (pear-shaped)`);
    console.log(`  J4 = ${earthJ4.toExponential(6)}`);
    console.log(`  J5 = ${earthJ5.toExponential(6)}`);
    console.log(`  J6 = ${earthJ6.toExponential(6)}`);

    // Relative magnitudes
    console.log('\nRelative magnitudes (compared to J2):');
    console.log(`  |J3/J2| = ${Math.abs(earthJ3/earthJ2).toExponential(2)}`);
    console.log(`  |J4/J2| = ${Math.abs(earthJ4/earthJ2).toExponential(2)}`);
    console.log(`  |J5/J2| = ${Math.abs(earthJ5/earthJ2).toExponential(2)}`);
    console.log(`  |J6/J2| = ${Math.abs(earthJ6/earthJ2).toExponential(2)}`);

    // ========================================
    // Legendre Polynomials
    // ========================================
    console.log('\n\n4. Legendre Polynomials');
    console.log('=======================\n');

    console.log('Legendre polynomials P_n(x) for gravity field expansion:');

    const xValues = [-1.0, -0.5, 0.0, 0.5, 1.0];

    console.log('\n  x      P_0    P_1    P_2     P_3      P_4');
    console.log('  ----   ----   ----   -----   ------   ------');

    for (const x of xValues) {
        // Compute Legendre polynomials
        const P0 = 1;
        const P1 = x;
        const P2 = 0.5 * (3 * x * x - 1);
        const P3 = 0.5 * (5 * x * x * x - 3 * x);
        const P4 = 0.125 * (35 * x**4 - 30 * x * x + 3);

        console.log(`  ${x.toFixed(1).padStart(4)}   ${P0.toFixed(2).padStart(4)}   ${P1.toFixed(2).padStart(4)}   ${P2.toFixed(2).padStart(5)}   ${P3.toFixed(3).padStart(6)}   ${P4.toFixed(3).padStart(6)}`);
    }

    // ========================================
    // Third-Body Perturbations
    // ========================================
    console.log('\n\n5. Third-Body Perturbations');
    console.log('===========================\n');

    // Third body data
    const thirdBodies = [
        { name: 'Moon', gm: 4.9028e12, distance: 384400e3 },
        { name: 'Sun', gm: 1.327e20, distance: 1.496e11 }
    ];

    console.log('Third-body acceleration at different satellite altitudes:');
    console.log('---------------------------------------------------------\n');

    const satAltitudes = [400, 35786, 384400];  // LEO, GEO, Moon distance

    for (const alt of satAltitudes) {
        const satR = earthRadius + alt * 1000;
        const earthAccel = earthGM / (satR * satR);

        console.log(`Satellite at ${alt} km altitude (Earth accel: ${earthAccel.toExponential(3)} m/s²):`);

        for (const body of thirdBodies) {
            // Simplified third-body acceleration (tidal acceleration)
            // a ≈ 2 * GM_third * r_sat / r_third³
            const tidalAccel = 2 * body.gm * satR / Math.pow(body.distance, 3);
            const ratio = tidalAccel / earthAccel;

            console.log(`  ${body.name}: ${tidalAccel.toExponential(3)} m/s² (${(ratio * 100).toExponential(2)}% of Earth)`);
        }
        console.log('');
    }

    // ========================================
    // Gravity Field at Different Locations
    // ========================================
    console.log('\n6. Gravity Variation with Latitude');
    console.log('===================================\n');

    // Surface gravity varies with latitude due to:
    // 1. Oblateness (equatorial radius > polar radius)
    // 2. Centrifugal acceleration (equator spins faster)
    // 3. J2 effect

    const earthOmega = 7.2921159e-5;  // rad/s (Earth rotation rate)
    const latitudes = [0, 30, 45, 60, 90];

    console.log('Surface gravity variation with latitude:');
    console.log('----------------------------------------');

    for (const lat of latitudes) {
        const phi = lat * Math.PI / 180;

        // Approximate radius at latitude (oblate spheroid)
        const f = 1/298.257;  // flattening
        const r_lat = earthRadius * (1 - f * Math.sin(phi)**2);

        // Point mass gravity
        const g_point = earthGM / (r_lat * r_lat);

        // J2 correction (approximate)
        const g_J2 = g_point * (1 + 1.5 * earthJ2 * (3 * Math.sin(phi)**2 - 1));

        // Centrifugal reduction at surface
        const centrifugal = earthOmega * earthOmega * r_lat * Math.cos(phi)**2;

        const g_total = g_J2 - centrifugal;

        console.log(`  ${lat.toString().padStart(2)}° latitude: g = ${g_total.toFixed(5)} m/s² (centrifugal: -${(centrifugal*1000).toFixed(3)} mm/s²)`);
    }

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
