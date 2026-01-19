/**
 * Example 10: Solar Radiation Pressure
 *
 * This example demonstrates solar radiation pressure calculations,
 * including cannonball and flat plate models.
 *
 * Run with: node 10_radiation_pressure.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Solar Radiation Pressure ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const c = 299792458;  // Speed of light [m/s]
    const AU = 1.496e11;  // Astronomical Unit [m]
    const solarLuminosity = 3.828e26;  // [W]

    // Solar radiation pressure at 1 AU
    const solarFlux = solarLuminosity / (4 * Math.PI * AU * AU);  // W/m²
    const P_1AU = solarFlux / c;  // Radiation pressure [Pa]

    console.log('1. Solar Radiation Fundamentals');
    console.log('===============================\n');

    console.log(`Solar luminosity: ${solarLuminosity.toExponential(3)} W`);
    console.log(`Solar flux at 1 AU: ${solarFlux.toFixed(1)} W/m²`);
    console.log(`Radiation pressure at 1 AU: ${P_1AU.toExponential(3)} Pa`);
    console.log(`                           = ${(P_1AU * 1e6).toFixed(3)} μPa`);

    // ========================================
    // Radiation Pressure vs Distance
    // ========================================
    console.log('\n\n2. Radiation Pressure vs Distance from Sun');
    console.log('==========================================\n');

    const distances = [
        { name: 'Mercury', r: 0.387 },
        { name: 'Venus', r: 0.723 },
        { name: 'Earth', r: 1.0 },
        { name: 'Mars', r: 1.524 },
        { name: 'Jupiter', r: 5.203 },
        { name: 'Saturn', r: 9.537 }
    ];

    console.log('Radiation pressure scales as 1/r²:\n');
    console.log('Planet      Distance [AU]    Pressure [μPa]    Relative');
    console.log('------      -------------    --------------    --------');

    for (const planet of distances) {
        const P = P_1AU / (planet.r * planet.r);
        const relative = 1 / (planet.r * planet.r);
        console.log(`${planet.name.padEnd(10)}  ${planet.r.toFixed(3).padStart(8)}         ${(P * 1e6).toFixed(3).padStart(8)}         ${relative.toFixed(3)}`);
    }

    // ========================================
    // Cannonball Model
    // ========================================
    console.log('\n\n3. Cannonball SRP Model');
    console.log('=======================\n');

    console.log('For a spherical (cannonball) spacecraft:');
    console.log('  F_SRP = P · A · C_r');
    console.log('  where C_r is the radiation pressure coefficient (1 ≤ C_r ≤ 2)\n');

    // Typical spacecraft
    const spacecraft = [
        { name: 'GPS satellite', area: 20, mass: 2000, cr: 1.2 },
        { name: 'GEO comsat', area: 50, mass: 3000, cr: 1.5 },
        { name: 'Interplanetary probe', area: 10, mass: 500, cr: 1.3 },
        { name: 'Solar sail demo', area: 100, mass: 10, cr: 2.0 }
    ];

    console.log('SRP acceleration at 1 AU:');
    console.log('-------------------------\n');

    for (const sc of spacecraft) {
        const F_srp = P_1AU * sc.area * sc.cr;
        const a_srp = F_srp / sc.mass;

        console.log(`${sc.name}:`);
        console.log(`  Area = ${sc.area} m², Mass = ${sc.mass} kg, C_r = ${sc.cr}`);
        console.log(`  Force: ${(F_srp * 1e6).toFixed(3)} μN`);
        console.log(`  Acceleration: ${a_srp.toExponential(3)} m/s²`);
        console.log(`  Area-to-mass ratio: ${(sc.area/sc.mass * 1000).toFixed(2)} mm²/g`);
        console.log('');
    }

    // ========================================
    // Comparison with Gravity
    // ========================================
    console.log('\n4. SRP vs Gravitational Acceleration');
    console.log('=====================================\n');

    const sunGM = 1.327e20;  // m³/s²

    console.log('Ratio of SRP to solar gravity:');
    console.log('-------------------------------\n');

    // For a typical spacecraft
    const A_m = 0.01;  // area-to-mass ratio [m²/kg]
    const Cr = 1.5;

    for (const planet of distances) {
        const r = planet.r * AU;
        const a_grav = sunGM / (r * r);
        const P = P_1AU / (planet.r * planet.r);
        const a_srp = P * A_m * Cr;
        const ratio = a_srp / a_grav;

        console.log(`${planet.name}: a_SRP/a_grav = ${ratio.toExponential(2)} (A/m = ${A_m * 1000} mm²/g)`);
    }

    // ========================================
    // Shadow Effects
    // ========================================
    console.log('\n\n5. Eclipse/Shadow Effects');
    console.log('=========================\n');

    console.log('During eclipse, SRP force drops to zero.');
    console.log('Shadow types: Umbra (total), Penumbra (partial)\n');

    // Earth shadow calculation
    const earthRadius = 6378137;  // m
    const sunRadius = 6.96e8;     // m

    // Umbra cone half-angle
    const umbra_angle = Math.asin((sunRadius - earthRadius) / AU);
    const umbra_length = earthRadius / Math.sin(umbra_angle);

    console.log('Earth shadow geometry (at 1 AU from Sun):');
    console.log(`  Umbra cone half-angle: ${(umbra_angle * 180 / Math.PI).toFixed(4)}°`);
    console.log(`  Umbra length: ${(umbra_length / 1e6).toFixed(1)} × 10⁶ km`);

    // GEO satellite eclipse
    const geoAlt = 35786e3;  // m
    const geoRadius = earthRadius + geoAlt;

    console.log(`\nGEO satellite (${geoAlt/1e3} km altitude):`);
    console.log(`  Max shadow duration: ~70 minutes/day during equinox season`);
    console.log(`  Eclipse seasons: ~46 days around each equinox`);

    // ========================================
    // Solar Sail Dynamics
    // ========================================
    console.log('\n\n6. Solar Sail Performance');
    console.log('=========================\n');

    console.log('Solar sail characteristic acceleration: a_c = 2P/σ');
    console.log('where σ = m/A is the sail loading [kg/m²]\n');

    const sailLoadings = [1, 5, 10, 50, 100];  // g/m²

    console.log('Characteristic acceleration at 1 AU:');
    console.log('------------------------------------');

    for (const sigma_gm2 of sailLoadings) {
        const sigma = sigma_gm2 / 1000;  // Convert g/m² to kg/m²
        const a_c = 2 * P_1AU / sigma;

        console.log(`  σ = ${sigma_gm2.toString().padStart(3)} g/m²: a_c = ${(a_c * 1000).toFixed(4)} mm/s²`);
    }

    // Time to reach various destinations
    console.log('\nApproximate transfer times for optimized solar sail:');
    console.log('(Assuming characteristic acceleration of 1 mm/s² at 1 AU)\n');

    const destinations = [
        { name: 'Mars', time: '2-3 years' },
        { name: 'Jupiter', time: '5-7 years' },
        { name: 'Solar escape', time: '~10 years' },
        { name: 'Alpha Centauri', time: '~1000 years (0.1c laser sail: ~40 years)' }
    ];

    for (const dest of destinations) {
        console.log(`  ${dest.name}: ${dest.time}`);
    }

    // ========================================
    // Yarkovsky Effect
    // ========================================
    console.log('\n\n7. Yarkovsky Effect (Asteroids)');
    console.log('================================\n');

    console.log('Thermal re-radiation causes secular drift in asteroid orbits.');
    console.log('Direction depends on spin axis orientation:\n');
    console.log('  - Prograde spin → outward drift');
    console.log('  - Retrograde spin → inward drift\n');

    // Typical Yarkovsky drift
    const yarkovsky_drift = 1e-4;  // AU/Myr for ~1 km asteroid

    console.log(`Typical drift rate for 1 km asteroid: ~${yarkovsky_drift} AU/Myr`);
    console.log(`Important for asteroid orbit prediction over decades.`);

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
