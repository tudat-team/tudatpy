/**
 * Example 19: Interplanetary Transfer Design
 *
 * This example demonstrates Hohmann transfers, porkchop plots,
 * and mission planning for interplanetary trajectories.
 *
 * Run with: node 19_interplanetary_transfer.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Interplanetary Transfer Design ===\n');

    const tudat = await createTudatModule();

    // Constants
    const sunGM = 1.32712440018e20;  // m³/s²
    const AU = 149597870700;          // m

    // Planetary data (simplified circular orbits)
    const planets = {
        Mercury: { sma: 0.387 * AU, period: 87.97, mu: 2.2032e13 },
        Venus:   { sma: 0.723 * AU, period: 224.7, mu: 3.2486e14 },
        Earth:   { sma: 1.000 * AU, period: 365.25, mu: 3.986e14 },
        Mars:    { sma: 1.524 * AU, period: 687.0, mu: 4.283e13 },
        Jupiter: { sma: 5.203 * AU, period: 4333, mu: 1.267e17 },
        Saturn:  { sma: 9.537 * AU, period: 10759, mu: 3.793e16 }
    };

    // ========================================
    // Hohmann Transfer
    // ========================================
    console.log('1. Hohmann Transfer (Minimum Energy)');
    console.log('=====================================\n');

    console.log('A Hohmann transfer is the most fuel-efficient two-impulse');
    console.log('transfer between coplanar circular orbits.\n');

    function hohmannTransfer(r1, r2, mu) {
        // Transfer orbit parameters
        const aTransfer = (r1 + r2) / 2;
        const eTransfer = Math.abs(r2 - r1) / (r1 + r2);

        // Velocities
        const v1Circular = Math.sqrt(mu / r1);
        const v2Circular = Math.sqrt(mu / r2);
        const v1Transfer = Math.sqrt(mu * (2/r1 - 1/aTransfer));
        const v2Transfer = Math.sqrt(mu * (2/r2 - 1/aTransfer));

        // Delta-V
        const dv1 = Math.abs(v1Transfer - v1Circular);
        const dv2 = Math.abs(v2Circular - v2Transfer);
        const totalDV = dv1 + dv2;

        // Transfer time (half orbital period)
        const tTransfer = Math.PI * Math.sqrt(Math.pow(aTransfer, 3) / mu);

        return {
            aTransfer, eTransfer,
            dv1, dv2, totalDV,
            tTransfer,
            v1Circular, v2Circular
        };
    }

    // Earth to Mars Hohmann transfer
    const earthMars = hohmannTransfer(planets.Earth.sma, planets.Mars.sma, sunGM);

    console.log('Earth to Mars Hohmann Transfer:');
    console.log(`  Departure orbit (Earth): ${(planets.Earth.sma / AU).toFixed(3)} AU`);
    console.log(`  Arrival orbit (Mars): ${(planets.Mars.sma / AU).toFixed(3)} AU`);
    console.log(`  Transfer semi-major axis: ${(earthMars.aTransfer / AU).toFixed(3)} AU`);
    console.log(`  Transfer eccentricity: ${earthMars.eTransfer.toFixed(4)}`);
    console.log(`  ΔV1 (departure): ${(earthMars.dv1 / 1000).toFixed(3)} km/s`);
    console.log(`  ΔV2 (arrival): ${(earthMars.dv2 / 1000).toFixed(3)} km/s`);
    console.log(`  Total ΔV: ${(earthMars.totalDV / 1000).toFixed(3)} km/s`);
    console.log(`  Transfer time: ${(earthMars.tTransfer / 86400).toFixed(1)} days`);

    // ========================================
    // Various Transfer Options
    // ========================================
    console.log('\n\n2. Transfer Options to Different Planets');
    console.log('=========================================\n');

    console.log('From Earth to:      ΔV (km/s)    Time (days)    Type');
    console.log('--------------      ---------    -----------    ----');

    const targets = ['Venus', 'Mars', 'Jupiter', 'Saturn'];

    for (const target of targets) {
        const transfer = hohmannTransfer(planets.Earth.sma, planets[target].sma, sunGM);
        const direction = planets[target].sma > planets.Earth.sma ? 'outer' : 'inner';
        console.log(`${target.padEnd(15)}     ${(transfer.totalDV / 1000).toFixed(3).padStart(6)}       ${(transfer.tTransfer / 86400).toFixed(0).padStart(7)}        ${direction}`);
    }

    // ========================================
    // Synodic Period and Launch Windows
    // ========================================
    console.log('\n\n3. Synodic Periods and Launch Windows');
    console.log('======================================\n');

    console.log('The synodic period is the time between successive alignments');
    console.log('of Earth and the target planet (optimal launch windows).\n');

    // Synodic period: 1/T_syn = |1/T1 - 1/T2|
    function synodicPeriod(T1, T2) {
        return 1 / Math.abs(1/T1 - 1/T2);
    }

    console.log('Planet     Orbital Period (days)    Synodic Period (days)');
    console.log('------     ---------------------    ---------------------');

    for (const [name, data] of Object.entries(planets)) {
        if (name !== 'Earth') {
            const syn = synodicPeriod(planets.Earth.period, data.period);
            console.log(`${name.padEnd(10)} ${data.period.toFixed(0).padStart(17)}            ${syn.toFixed(0).padStart(17)}`);
        }
    }

    // ========================================
    // Phase Angle
    // ========================================
    console.log('\n\n4. Phase Angle at Departure');
    console.log('===========================\n');

    console.log('The phase angle is the angular separation between Earth and');
    console.log('the target planet at departure.\n');

    // Phase angle for Hohmann transfer
    // During transfer, target moves: angle = n_target * t_transfer
    // Phase angle = 180° - angle_moved
    function phaseAngle(r1, r2, mu, targetPeriod) {
        const aTransfer = (r1 + r2) / 2;
        const tTransfer = Math.PI * Math.sqrt(Math.pow(aTransfer, 3) / mu);

        // Target angular motion during transfer
        const targetAngle = (tTransfer / (targetPeriod * 86400)) * 360;

        // For outer planet transfer, target should be ahead
        if (r2 > r1) {
            return 180 - targetAngle;
        } else {
            return 180 + targetAngle;
        }
    }

    console.log('Optimal phase angles (Hohmann transfer):');
    console.log('Target      Phase Angle    Position');
    console.log('------      -----------    --------');

    for (const target of targets) {
        const phase = phaseAngle(planets.Earth.sma, planets[target].sma, sunGM, planets[target].period);
        const position = planets[target].sma > planets.Earth.sma ? 'ahead' : 'behind';
        console.log(`${target.padEnd(11)} ${phase.toFixed(1).padStart(9)}°      ${position}`);
    }

    // ========================================
    // Porkchop Plot Concept
    // ========================================
    console.log('\n\n5. Porkchop Plot Concept');
    console.log('========================\n');

    console.log('A porkchop plot shows mission ΔV as a function of');
    console.log('launch date and arrival date (or time of flight).\n');

    // Simplified porkchop data for Earth-Mars
    console.log('Simplified Earth-Mars Porkchop (2026 launch window):');
    console.log('Launch Day    TOF (days)    C3 (km²/s²)    V∞ arr (km/s)');
    console.log('----------    ----------    -----------    -------------');

    // Example data points (illustrative)
    const porkchopData = [
        { launch: 'Day 0', tof: 180, c3: 15.2, vInf: 3.5 },
        { launch: 'Day 15', tof: 200, c3: 12.5, vInf: 3.2 },
        { launch: 'Day 30', tof: 220, c3: 10.8, vInf: 2.8 },
        { launch: 'Day 45', tof: 250, c3: 9.5, vInf: 2.5 },
        { launch: 'Day 60', tof: 280, c3: 11.2, vInf: 2.9 },
        { launch: 'Day 75', tof: 300, c3: 14.8, vInf: 3.4 }
    ];

    for (const point of porkchopData) {
        console.log(`${point.launch.padEnd(13)} ${point.tof.toString().padStart(10)}    ${point.c3.toFixed(1).padStart(11)}    ${point.vInf.toFixed(2).padStart(13)}`);
    }

    console.log('\nOptimal window: Day 45 (lowest total energy)');

    // ========================================
    // Type I and Type II Trajectories
    // ========================================
    console.log('\n\n6. Type I and Type II Trajectories');
    console.log('===================================\n');

    console.log('Transfer trajectories are classified by heliocentric angle:');
    console.log('  Type I: < 180° transfer angle (shorter)');
    console.log('  Type II: > 180° transfer angle (longer)\n');

    const transfer = hohmannTransfer(planets.Earth.sma, planets.Mars.sma, sunGM);

    console.log('Example: Earth-Mars');
    console.log(`  Hohmann (Type I): 180° transfer, ${(transfer.tTransfer / 86400).toFixed(0)} days`);
    console.log(`  Type II option: ~360° - 180° = longer path, ~500+ days`);
    console.log('\nType II transfers:');
    console.log('  + Sometimes lower C3 requirements');
    console.log('  + More flexible launch windows');
    console.log('  - Longer mission duration');
    console.log('  - More radiation exposure');

    // ========================================
    // Fast Transfers
    // ========================================
    console.log('\n\n7. Fast Transfers (Non-Hohmann)');
    console.log('===============================\n');

    console.log('Faster transfers require more ΔV but reduce mission time.\n');

    // Lambert problem concept for different TOFs
    console.log('Earth-Mars transfer options:');
    console.log('TOF (days)    ΔV departure    ΔV arrival    Total ΔV');
    console.log('----------    ------------    ----------    --------');

    // Illustrative data (actual values from Lambert solver would differ)
    const fastTransfers = [
        { tof: 120, dvDep: 5.2, dvArr: 4.8, total: 10.0 },
        { tof: 150, dvDep: 4.1, dvArr: 3.9, total: 8.0 },
        { tof: 180, dvDep: 3.5, dvArr: 3.2, total: 6.7 },
        { tof: 210, dvDep: 3.1, dvArr: 2.8, total: 5.9 },
        { tof: 259, dvDep: 2.9, dvArr: 2.6, total: 5.5 },  // Hohmann
        { tof: 300, dvDep: 3.2, dvArr: 2.5, total: 5.7 }
    ];

    for (const t of fastTransfers) {
        const note = t.tof === 259 ? ' (Hohmann)' : '';
        console.log(`${t.tof.toString().padStart(10)}    ${t.dvDep.toFixed(1).padStart(12)}    ${t.dvArr.toFixed(1).padStart(10)}    ${t.total.toFixed(1).padStart(8)}${note}`);
    }

    // ========================================
    // Gravity Assist Options
    // ========================================
    console.log('\n\n8. Gravity Assist Options');
    console.log('=========================\n');

    console.log('Gravity assists can reduce ΔV for outer planet missions.\n');

    console.log('Common gravity assist sequences:');
    console.log('  VEGA (Venus-Earth-Gravity Assist): Galileo, Cassini');
    console.log('  VEEGA (Venus-Earth-Earth): More energy gain');
    console.log('  MEGA (Mars-Earth): Ulysses solar polar mission');
    console.log('  JGA (Jupiter): Outer solar system missions\n');

    console.log('Example: Earth direct to Jupiter');
    const earthJupiter = hohmannTransfer(planets.Earth.sma, planets.Jupiter.sma, sunGM);
    console.log(`  Direct Hohmann ΔV: ${(earthJupiter.totalDV / 1000).toFixed(2)} km/s`);
    console.log(`  Transfer time: ${(earthJupiter.tTransfer / 86400 / 365.25).toFixed(2)} years\n`);

    console.log('With Venus-Venus-Earth-Jupiter (VVEJGA):');
    console.log('  Total ΔV: ~3.5 km/s (significantly less!)');
    console.log('  Transfer time: ~6 years (much longer)');
    console.log('  Used by: Galileo, Cassini');

    // ========================================
    // Mission Design Summary
    // ========================================
    console.log('\n\n9. Mars Mission Design Summary');
    console.log('==============================\n');

    console.log('Typical Mars mission profile:\n');

    console.log('Phase              Duration    ΔV Required');
    console.log('-----              --------    -----------');
    console.log('Earth parking orbit      -     LEO insertion');
    console.log('Trans-Mars injection     -     ~3.6 km/s');
    console.log('Cruise phase        180-300d         -');
    console.log('Mars orbit insertion     -     ~2.1 km/s (capture)');
    console.log('Aerobraking (if used)   months  minimal fuel');
    console.log('Science orbit            -     ~0.5 km/s (circularize)');

    console.log('\nTotal from LEO: ~6.2 km/s (without aerobraking)');
    console.log('Total from LEO: ~3.6 km/s (with aerobraking for capture)');

    // Next launch window
    console.log('\nMars launch windows occur every ~26 months');
    console.log('Recent/upcoming windows: 2022, 2024, 2026, 2028...');

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
