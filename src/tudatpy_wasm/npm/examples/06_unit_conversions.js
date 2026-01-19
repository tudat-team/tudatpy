/**
 * Example 06: Unit Conversions and Physical Constants
 *
 * This example demonstrates the unit conversion utilities and
 * physical constants available in Tudat WASM.
 *
 * Run with: node 06_unit_conversions.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Unit Conversions and Physical Constants ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Physical Constants
    // ========================================
    console.log('Physical Constants:');
    console.log('-------------------');

    const constants = {
        'Speed of light': { value: tudat.constants.SPEED_OF_LIGHT, unit: 'm/s' },
        'Gravitational constant': { value: tudat.constants.GRAVITATIONAL_CONSTANT, unit: 'm³/(kg·s²)' },
        'Astronomical unit': { value: tudat.constants.ASTRONOMICAL_UNIT, unit: 'm' },
        'Earth GM': { value: tudat.constants.EARTH_GRAVITATIONAL_PARAMETER, unit: 'm³/s²' },
        'Sun GM': { value: tudat.constants.SUN_GRAVITATIONAL_PARAMETER, unit: 'm³/s²' },
        'Moon GM': { value: tudat.constants.MOON_GRAVITATIONAL_PARAMETER, unit: 'm³/s²' },
        'Earth radius': { value: tudat.constants.EARTH_EQUATORIAL_RADIUS, unit: 'm' },
        'Sea level gravity': { value: tudat.constants.SEA_LEVEL_GRAVITATIONAL_ACCELERATION, unit: 'm/s²' },
        'Stefan-Boltzmann': { value: tudat.constants.STEFAN_BOLTZMANN_CONSTANT, unit: 'W/(m²·K⁴)' },
        'Planck constant': { value: tudat.constants.PLANCK_CONSTANT, unit: 'J·s' }
    };

    for (const [name, data] of Object.entries(constants)) {
        if (data.value !== undefined) {
            console.log(`  ${name}: ${data.value.toExponential(6)} ${data.unit}`);
        }
    }

    // ========================================
    // Angle Conversions
    // ========================================
    console.log('\n\nAngle Conversions:');
    console.log('------------------');

    // Degrees to radians
    const angles_deg = [0, 30, 45, 60, 90, 180, 360];
    console.log('Degrees to Radians:');
    for (const deg of angles_deg) {
        const rad = tudat.math.unit_conversions.convert_degrees_to_radians(deg);
        console.log(`  ${deg}° = ${rad.toFixed(6)} rad`);
    }

    // Radians to degrees
    console.log('\nRadians to Degrees:');
    const angles_rad = [0, Math.PI/6, Math.PI/4, Math.PI/2, Math.PI, 2*Math.PI];
    for (const rad of angles_rad) {
        const deg = tudat.math.unit_conversions.convert_radians_to_degrees(rad);
        console.log(`  ${rad.toFixed(4)} rad = ${deg.toFixed(2)}°`);
    }

    // ========================================
    // Distance Conversions
    // ========================================
    console.log('\n\nDistance Conversions:');
    console.log('---------------------');

    // AU to meters
    const AU = tudat.constants.ASTRONOMICAL_UNIT;
    console.log('Astronomical Units to Meters:');
    const distances_au = [0.1, 0.5, 1.0, 1.52, 5.2, 30];  // Mercury-ish to Neptune-ish
    for (const au of distances_au) {
        const meters = au * AU;
        console.log(`  ${au} AU = ${(meters/1e9).toFixed(3)} × 10⁹ m`);
    }

    // Kilometers to meters (simple but common)
    console.log('\nKilometers to Meters:');
    const distances_km = [400, 35786, 384400];  // ISS, GEO, Moon
    const labels = ['ISS altitude', 'GEO altitude', 'Moon distance'];
    for (let i = 0; i < distances_km.length; i++) {
        console.log(`  ${labels[i]}: ${distances_km[i]} km = ${distances_km[i] * 1000} m`);
    }

    // ========================================
    // Time Conversions
    // ========================================
    console.log('\n\nTime Conversions:');
    console.log('-----------------');

    const SECONDS_PER_DAY = 86400;
    const SECONDS_PER_YEAR = 365.25 * SECONDS_PER_DAY;

    // Orbital periods in different units
    console.log('Orbital Periods:');
    const orbits = [
        { name: 'ISS', period_s: 92.68 * 60 },
        { name: 'GEO', period_s: 23.934 * 3600 },
        { name: 'Moon', period_s: 27.3 * SECONDS_PER_DAY },
        { name: 'Earth (year)', period_s: SECONDS_PER_YEAR },
        { name: 'Mars', period_s: 1.881 * SECONDS_PER_YEAR },
        { name: 'Jupiter', period_s: 11.86 * SECONDS_PER_YEAR }
    ];

    for (const orbit of orbits) {
        const days = orbit.period_s / SECONDS_PER_DAY;
        const years = orbit.period_s / SECONDS_PER_YEAR;
        if (days < 1) {
            console.log(`  ${orbit.name}: ${(orbit.period_s/60).toFixed(2)} min`);
        } else if (days < 365) {
            console.log(`  ${orbit.name}: ${days.toFixed(2)} days`);
        } else {
            console.log(`  ${orbit.name}: ${years.toFixed(2)} years`);
        }
    }

    // ========================================
    // Julian Date Conversions
    // ========================================
    console.log('\n\nJulian Date Conversions:');
    console.log('------------------------');

    // J2000 epoch
    const J2000_JD = 2451545.0;  // January 1, 2000, 12:00 TT
    console.log(`J2000 epoch: JD ${J2000_JD}`);

    // Convert seconds since J2000 to Julian Date
    const epochs = [
        { name: 'J2000', seconds: 0 },
        { name: '2024-01-01', seconds: 24 * 365.25 * SECONDS_PER_DAY },
        { name: '2030-01-01', seconds: 30 * 365.25 * SECONDS_PER_DAY }
    ];

    console.log('\nSeconds since J2000 to Julian Date:');
    for (const epoch of epochs) {
        const jd = J2000_JD + epoch.seconds / SECONDS_PER_DAY;
        console.log(`  ${epoch.name}: ${epoch.seconds.toExponential(4)} s → JD ${jd.toFixed(2)}`);
    }

    // ========================================
    // Velocity Conversions
    // ========================================
    console.log('\n\nVelocity Conversions:');
    console.log('---------------------');

    const velocities = [
        { name: 'Walking speed', ms: 1.4 },
        { name: 'Sound (sea level)', ms: 343 },
        { name: 'ISS orbital velocity', ms: 7660 },
        { name: 'Earth escape velocity', ms: 11186 },
        { name: 'Earth orbital velocity', ms: 29780 },
        { name: 'Solar escape (at Earth)', ms: 42100 }
    ];

    console.log('Velocities in m/s and km/s:');
    for (const v of velocities) {
        console.log(`  ${v.name}: ${v.ms.toFixed(1)} m/s = ${(v.ms/1000).toFixed(3)} km/s`);
    }

    // ========================================
    // Mass Conversions
    // ========================================
    console.log('\n\nPlanetary Masses (GM values):');
    console.log('-----------------------------');

    const G = tudat.constants.GRAVITATIONAL_CONSTANT;
    const bodies = [
        { name: 'Earth', gm: tudat.constants.EARTH_GRAVITATIONAL_PARAMETER },
        { name: 'Moon', gm: tudat.constants.MOON_GRAVITATIONAL_PARAMETER },
        { name: 'Sun', gm: tudat.constants.SUN_GRAVITATIONAL_PARAMETER }
    ];

    for (const body of bodies) {
        if (body.gm !== undefined) {
            const mass = body.gm / G;
            console.log(`  ${body.name}: GM = ${body.gm.toExponential(6)} m³/s² → M = ${mass.toExponential(3)} kg`);
        }
    }

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
