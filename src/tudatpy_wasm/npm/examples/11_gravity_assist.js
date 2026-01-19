/**
 * Example 11: Gravity Assist Maneuvers
 *
 * This example demonstrates gravity assist (flyby) trajectory design,
 * including unpowered and powered swing-by maneuvers.
 *
 * Run with: node 11_gravity_assist.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Gravity Assist Maneuvers ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Gravity Assist Fundamentals
    // ========================================
    console.log('1. Gravity Assist Fundamentals');
    console.log('==============================\n');

    console.log('A gravity assist uses a planet\'s gravity and motion to');
    console.log('change a spacecraft\'s velocity without using propellant.\n');

    console.log('Key equations:');
    console.log('  - Turn angle: δ = 2·arcsin(1/(1 + r_p·v_∞²/μ))');
    console.log('  - ΔV (heliocentric) ≈ 2·v_∞·sin(δ/2)');
    console.log('  - Max ΔV = 2·v_planet (theoretical limit)\n');

    // Planetary data
    const planets = {
        'Venus': { gm: 3.249e14, radius: 6052e3, v_orb: 35020, r_orb: 1.082e11 },
        'Earth': { gm: 3.986e14, radius: 6378e3, v_orb: 29780, r_orb: 1.496e11 },
        'Mars': { gm: 4.283e13, radius: 3396e3, v_orb: 24130, r_orb: 2.279e11 },
        'Jupiter': { gm: 1.267e17, radius: 71492e3, v_orb: 13060, r_orb: 7.785e11 },
        'Saturn': { gm: 3.794e16, radius: 60268e3, v_orb: 9680, r_orb: 14.27e11 }
    };

    // ========================================
    // Turn Angle Calculation
    // ========================================
    console.log('2. Turn Angle vs Flyby Distance');
    console.log('================================\n');

    // Calculate turn angle for different V_infinity and periapsis distances
    function calculateTurnAngle(v_inf, r_p, gm) {
        const e = 1 + r_p * v_inf * v_inf / gm;  // eccentricity
        return 2 * Math.asin(1 / e);  // turn angle
    }

    console.log('Earth flyby (v_∞ = 5 km/s):');
    console.log('---------------------------');

    const v_inf = 5000;  // m/s
    const earthGM = planets['Earth'].gm;
    const earthR = planets['Earth'].radius;

    const altitudes = [200, 500, 1000, 2000, 5000, 10000];  // km

    console.log('Altitude [km]    Turn Angle [°]    ΔV (helio) [km/s]');
    console.log('-------------    --------------    -----------------');

    for (const alt of altitudes) {
        const r_p = earthR + alt * 1000;
        const delta = calculateTurnAngle(v_inf, r_p, earthGM);
        const deltaV = 2 * v_inf * Math.sin(delta / 2);

        console.log(`${alt.toString().padStart(8)}         ${(delta * 180 / Math.PI).toFixed(2).padStart(10)}         ${(deltaV / 1000).toFixed(3).padStart(10)}`);
    }

    // ========================================
    // Planetary Comparison
    // ========================================
    console.log('\n\n3. Gravity Assist Capability by Planet');
    console.log('=======================================\n');

    console.log('Maximum ΔV for 500 km periapsis altitude, v_∞ = 10 km/s:');
    console.log('--------------------------------------------------------\n');

    const v_inf_test = 10000;  // m/s

    console.log('Planet      Radius [km]    GM [km³/s²]    Turn [°]    ΔV [km/s]');
    console.log('------      -----------    -----------    --------    ---------');

    for (const [name, data] of Object.entries(planets)) {
        const r_p = data.radius + 500e3;
        const delta = calculateTurnAngle(v_inf_test, r_p, data.gm);
        const deltaV = 2 * v_inf_test * Math.sin(delta / 2);

        console.log(`${name.padEnd(10)}  ${(data.radius/1000).toFixed(0).padStart(8)}      ${(data.gm/1e9).toExponential(2).padStart(10)}    ${(delta * 180 / Math.PI).toFixed(1).padStart(7)}     ${(deltaV / 1000).toFixed(2).padStart(8)}`);
    }

    // ========================================
    // Cassini Gravity Assists
    // ========================================
    console.log('\n\n4. Cassini Mission Gravity Assists (VVEJGA)');
    console.log('============================================\n');

    const cassiniFlybys = [
        { planet: 'Venus', date: '1998-04-26', v_inf: 6.3, alt: 284, dv: 7.0 },
        { planet: 'Venus', date: '1999-06-24', v_inf: 6.5, alt: 600, dv: 6.7 },
        { planet: 'Earth', date: '1999-08-18', v_inf: 8.9, alt: 1171, dv: 5.5 },
        { planet: 'Jupiter', date: '2000-12-30', v_inf: 11.6, alt: 9.72e6, dv: 2.2 }
    ];

    console.log('Flyby        Date          v_∞ [km/s]   Altitude [km]   ΔV [km/s]');
    console.log('------       ----------    ----------   -------------   ---------');

    for (const fb of cassiniFlybys) {
        console.log(`${fb.planet.padEnd(10)}   ${fb.date}      ${fb.v_inf.toFixed(1).padStart(6)}        ${fb.alt.toString().padStart(10)}       ${fb.dv.toFixed(1)}`);
    }

    console.log('\nTotal heliocentric velocity gain: ~21 km/s');
    console.log('(Equivalent to ~7 km/s ΔV if done by propulsion)');

    // ========================================
    // Powered Flyby
    // ========================================
    console.log('\n\n5. Powered Gravity Assist (Oberth Effect)');
    console.log('==========================================\n');

    console.log('Performing a burn at periapsis maximizes ΔV effectiveness');
    console.log('due to the Oberth effect: ΔE = v · Δv\n');

    // Compare powered vs unpowered flyby
    const v_inf_in = 8000;  // m/s incoming
    const r_p_powered = earthR + 300e3;  // 300 km periapsis
    const dv_burn = 500;  // m/s propulsive maneuver

    // Velocity at periapsis
    const v_p = Math.sqrt(v_inf_in**2 + 2 * earthGM / r_p_powered);

    // After burn
    const v_p_after = v_p + dv_burn;

    // Outgoing v_infinity
    const v_inf_out = Math.sqrt(v_p_after**2 - 2 * earthGM / r_p_powered);

    console.log('Powered Earth flyby comparison:');
    console.log('-------------------------------');
    console.log(`  Incoming v_∞: ${(v_inf_in/1000).toFixed(2)} km/s`);
    console.log(`  Periapsis altitude: 300 km`);
    console.log(`  Velocity at periapsis: ${(v_p/1000).toFixed(2)} km/s`);
    console.log(`  Propulsive burn: ${dv_burn} m/s`);
    console.log(`  Outgoing v_∞: ${(v_inf_out/1000).toFixed(2)} km/s`);
    console.log(`  v_∞ gain: ${((v_inf_out - v_inf_in)/1000).toFixed(2)} km/s`);
    console.log(`\n  Oberth multiplier: ${((v_inf_out - v_inf_in) / dv_burn).toFixed(2)}×`);

    // ========================================
    // B-Plane Targeting
    // ========================================
    console.log('\n\n6. B-Plane Targeting');
    console.log('====================\n');

    console.log('The B-plane is perpendicular to the incoming asymptote.');
    console.log('B-plane coordinates (B·R, B·T) define the flyby geometry:\n');
    console.log('  B = impact parameter = r_p × sqrt(1 + 2μ/(r_p × v_∞²))');
    console.log('  B·T = in-plane component (affects timing)');
    console.log('  B·R = out-of-plane component (affects inclination)\n');

    // Calculate B for different periapsis distances
    console.log('Impact parameter B for Earth flyby (v_∞ = 5 km/s):');
    console.log('--------------------------------------------------');

    for (const alt of [200, 500, 1000, 5000]) {
        const r_p = earthR + alt * 1000;
        const B = r_p * Math.sqrt(1 + 2 * earthGM / (r_p * v_inf * v_inf));

        console.log(`  Alt = ${alt} km: B = ${(B/1000).toFixed(0)} km`);
    }

    // ========================================
    // Escape and Capture
    // ========================================
    console.log('\n\n7. Escape and Capture ΔV');
    console.log('========================\n');

    console.log('For missions from/to planetary surfaces:\n');

    const moonGM = 4.9028e12;
    const moonRadius = 1737.4e3;

    // Calculate escape velocities
    const bodies = [
        { name: 'Moon', gm: moonGM, r: moonRadius },
        { name: 'Mars', gm: planets['Mars'].gm, r: planets['Mars'].radius },
        { name: 'Earth', gm: earthGM, r: earthR }
    ];

    console.log('Body        Escape velocity    Orbit velocity (100 km)');
    console.log('----        ---------------    -----------------------');

    for (const body of bodies) {
        const v_esc = Math.sqrt(2 * body.gm / body.r);
        const v_orb = Math.sqrt(body.gm / (body.r + 100e3));

        console.log(`${body.name.padEnd(10)}  ${(v_esc/1000).toFixed(3).padStart(10)} km/s      ${(v_orb/1000).toFixed(3).padStart(10)} km/s`);
    }

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
