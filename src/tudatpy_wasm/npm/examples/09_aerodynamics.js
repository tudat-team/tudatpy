/**
 * Example 09: Aerodynamics and Atmosphere Models
 *
 * This example demonstrates atmosphere models and aerodynamic force
 * calculations for spacecraft in planetary atmospheres.
 *
 * Run with: node 09_aerodynamics.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Aerodynamics and Atmosphere Models ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Exponential Atmosphere Model
    // ========================================
    console.log('1. Exponential Atmosphere Model');
    console.log('===============================\n');

    // Earth exponential atmosphere parameters
    const rho0 = 1.225;        // kg/m³ (sea level density)
    const H = 8500;            // m (scale height)
    const earthRadius = 6378137.0;

    console.log('Earth exponential atmosphere:');
    console.log(`  Sea level density: ${rho0} kg/m³`);
    console.log(`  Scale height: ${H/1000} km\n`);

    console.log('Density at different altitudes:');
    console.log('-------------------------------');

    const altitudes = [0, 50, 100, 150, 200, 300, 400, 500, 800, 1000];

    for (const alt of altitudes) {
        const h = alt * 1000;  // convert to meters
        const rho = rho0 * Math.exp(-h / H);

        console.log(`  ${alt.toString().padStart(4)} km: ρ = ${rho.toExponential(3)} kg/m³`);
    }

    // ========================================
    // NRLMSISE-00 Model (Approximation)
    // ========================================
    console.log('\n\n2. NRLMSISE-00 Atmosphere Model');
    console.log('================================\n');

    console.log('NRLMSISE-00 provides more accurate density including:');
    console.log('  - Diurnal variations');
    console.log('  - Solar activity (F10.7 index)');
    console.log('  - Geomagnetic activity (Ap index)');
    console.log('  - Latitude/longitude dependence\n');

    // Approximate NRLMSISE-00 densities at different conditions
    const conditions = [
        { name: 'Solar minimum, quiet', f107: 70, ap: 4, factor: 0.5 },
        { name: 'Solar average', f107: 140, ap: 10, factor: 1.0 },
        { name: 'Solar maximum, active', f107: 200, ap: 50, factor: 3.0 }
    ];

    console.log('Density variation at 400 km with solar activity:');
    console.log('------------------------------------------------');

    const baseRho400 = rho0 * Math.exp(-400000 / H);

    for (const cond of conditions) {
        const rho = baseRho400 * cond.factor;
        console.log(`  ${cond.name}:`);
        console.log(`    F10.7=${cond.f107}, Ap=${cond.ap}: ρ ≈ ${rho.toExponential(3)} kg/m³`);
    }

    // ========================================
    // Aerodynamic Forces
    // ========================================
    console.log('\n\n3. Aerodynamic Force Calculations');
    console.log('==================================\n');

    // Drag force: F_D = 0.5 * ρ * v² * C_D * A
    // Lift force: F_L = 0.5 * ρ * v² * C_L * A

    console.log('Drag force equation: F_D = ½ρv²C_D·A\n');

    // Typical spacecraft parameters
    const spacecraft = [
        { name: 'CubeSat (3U)', cd: 2.2, area: 0.03, mass: 4 },
        { name: 'ISS', cd: 2.0, area: 2500, mass: 420000 },
        { name: 'Starlink', cd: 2.2, area: 10, mass: 260 },
        { name: 'Space Shuttle', cd: 1.0, area: 250, mass: 100000 }
    ];

    const alt = 400;  // km
    const rho = rho0 * Math.exp(-alt * 1000 / H);
    const v = 7660;  // m/s (orbital velocity at 400 km)

    console.log(`At ${alt} km altitude (ρ = ${rho.toExponential(3)} kg/m³, v = ${v} m/s):`);
    console.log('------------------------------------------------------------------------\n');

    for (const sc of spacecraft) {
        const drag = 0.5 * rho * v * v * sc.cd * sc.area;
        const accel = drag / sc.mass;

        console.log(`${sc.name}:`);
        console.log(`  C_D = ${sc.cd}, A = ${sc.area} m², m = ${sc.mass} kg`);
        console.log(`  Drag force: ${drag.toExponential(3)} N`);
        console.log(`  Drag acceleration: ${accel.toExponential(3)} m/s²`);

        // Orbital decay rate approximation
        // dh/dt ≈ -ρ * v * C_D * A / (2m)
        const decayRate = rho * v * sc.cd * sc.area / (2 * sc.mass);
        console.log(`  Decay rate: ${(decayRate * 86400).toFixed(4)} m/day`);
        console.log('');
    }

    // ========================================
    // Ballistic Coefficient
    // ========================================
    console.log('\n4. Ballistic Coefficient');
    console.log('========================\n');

    console.log('Ballistic coefficient: β = m / (C_D · A)');
    console.log('Higher β → less affected by drag\n');

    console.log('Spacecraft ballistic coefficients:');
    console.log('----------------------------------');

    for (const sc of spacecraft) {
        const beta = sc.mass / (sc.cd * sc.area);
        console.log(`  ${sc.name}: β = ${beta.toFixed(2)} kg/m²`);
    }

    // ========================================
    // Dynamic Pressure
    // ========================================
    console.log('\n\n5. Dynamic Pressure During Reentry');
    console.log('===================================\n');

    console.log('Dynamic pressure: q = ½ρv²\n');

    // Reentry profile (simplified)
    const reentryPoints = [
        { alt: 120, v: 7800 },
        { alt: 100, v: 7600 },
        { alt: 80, v: 7200 },
        { alt: 60, v: 5000 },
        { alt: 40, v: 2000 },
        { alt: 20, v: 500 }
    ];

    console.log('Dynamic pressure during reentry:');
    console.log('--------------------------------');

    let maxQ = 0;
    let maxQAlt = 0;

    for (const point of reentryPoints) {
        const rho_point = rho0 * Math.exp(-point.alt * 1000 / H);
        const q = 0.5 * rho_point * point.v * point.v;

        if (q > maxQ) {
            maxQ = q;
            maxQAlt = point.alt;
        }

        console.log(`  ${point.alt.toString().padStart(3)} km, ${point.v} m/s: q = ${(q/1000).toFixed(2)} kPa`);
    }

    console.log(`\nMax Q occurs around ${maxQAlt} km altitude: ${(maxQ/1000).toFixed(2)} kPa`);

    // ========================================
    // Heating Rate
    // ========================================
    console.log('\n\n6. Aerodynamic Heating');
    console.log('======================\n');

    // Stagnation point heating rate (Chapman formula approximation)
    // q_dot ≈ k * sqrt(ρ/r_n) * v³

    console.log('Stagnation point heating rate ∝ √(ρ/r_n) · v³\n');

    const noseRadius = 0.3;  // m (nose radius)
    const k = 1.83e-8;  // Empirical constant (approximate)

    console.log(`Heating rate for nose radius = ${noseRadius} m:`);
    console.log('----------------------------------------------');

    for (const point of reentryPoints) {
        const rho_point = rho0 * Math.exp(-point.alt * 1000 / H);
        const qDot = k * Math.sqrt(rho_point / noseRadius) * Math.pow(point.v, 3);

        console.log(`  ${point.alt.toString().padStart(3)} km, ${point.v} m/s: q̇ = ${(qDot/1e6).toFixed(2)} MW/m²`);
    }

    // ========================================
    // Atmospheric Scale Heights
    // ========================================
    console.log('\n\n7. Scale Heights for Different Planets');
    console.log('======================================\n');

    const planets = [
        { name: 'Venus', H: 15900, rho0: 65, g: 8.87 },
        { name: 'Earth', H: 8500, rho0: 1.225, g: 9.81 },
        { name: 'Mars', H: 11100, rho0: 0.020, g: 3.71 },
        { name: 'Titan', H: 21000, rho0: 5.3, g: 1.35 }
    ];

    console.log('Planet      Scale Height    Surface Density    Surface g');
    console.log('------      ------------    ---------------    ---------');

    for (const planet of planets) {
        console.log(`${planet.name.padEnd(10)}  ${(planet.H/1000).toFixed(1).padStart(6)} km       ${planet.rho0.toFixed(3).padStart(8)} kg/m³     ${planet.g.toFixed(2)} m/s²`);
    }

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
