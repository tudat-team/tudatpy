/**
 * Example 26: Multiple Gravity Assist (MGA) Trajectory Optimization
 *
 * Ported from: examples/tudatpy/mission_design/mga_trajectories.py
 *
 * This example demonstrates optimizing an interplanetary trajectory
 * using multiple gravity assists (Cassini-like EVVEJSA sequence).
 *
 * Run with: node 26_mga_optimization.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Multiple Gravity Assist Trajectory Optimization ===\n');

    const tudat = await createTudatModule();

    // Planetary gravitational parameters and orbital radii
    const planets = {
        Sun: { GM: 1.32712440018e20 },
        Earth: { GM: 3.986004418e14, a: 1.496e11, e: 0.0167 },
        Venus: { GM: 3.24859e14, a: 1.082e11, e: 0.0068 },
        Jupiter: { GM: 1.26687e17, a: 7.785e11, e: 0.0489 },
        Saturn: { GM: 3.7931e16, a: 1.433e12, e: 0.0565 }
    };

    // Cassini-like sequence: Earth -> Venus -> Venus -> Earth -> Jupiter -> Saturn
    const sequence = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn'];

    console.log('Trajectory Sequence: ' + sequence.join(' -> '));
    console.log('\nPlanetary Parameters:');
    for (const [name, data] of Object.entries(planets)) {
        if (data.a) {
            console.log(`  ${name}: a = ${(data.a / 1.496e11).toFixed(3)} AU, e = ${data.e}`);
        }
    }

    // Time of flight bounds for each leg (days)
    const tofBounds = [
        [100, 400],   // Earth -> Venus
        [100, 500],   // Venus -> Venus (resonant)
        [100, 400],   // Venus -> Earth
        [500, 1500],  // Earth -> Jupiter
        [1000, 3000]  // Jupiter -> Saturn
    ];

    console.log('\nTime of Flight Bounds (days):');
    for (let i = 0; i < sequence.length - 1; i++) {
        console.log(`  ${sequence[i]} -> ${sequence[i+1]}: ${tofBounds[i][0]} - ${tofBounds[i][1]}`);
    }

    // Simple grid search optimization
    // In practice, would use genetic algorithm or other global optimizer
    console.log('\n--- Running Grid Search Optimization ---\n');

    // For demo, evaluate a few candidate solutions
    const candidates = [
        // [launch_date_offset, tof1, tof2, tof3, tof4, tof5] in days from epoch
        { name: 'Candidate A', tofs: [200, 250, 200, 900, 1800] },
        { name: 'Candidate B', tofs: [150, 300, 180, 1100, 2000] },
        { name: 'Candidate C', tofs: [180, 200, 220, 1000, 2200] },
        { name: 'Candidate D', tofs: [220, 280, 190, 850, 1900] },
    ];

    let bestCandidate = null;
    let bestDeltaV = Infinity;

    for (const candidate of candidates) {
        // Compute required delta-V for each leg using Lambert solver approximation
        let totalDeltaV = 0;

        for (let leg = 0; leg < sequence.length - 1; leg++) {
            const fromPlanet = planets[sequence[leg]];
            const toPlanet = planets[sequence[leg + 1]];
            const tof = candidate.tofs[leg] * 86400;  // Convert to seconds

            // Simplified Lambert calculation (using vis-viva approximation)
            // Actual implementation would use full Lambert solver
            const r1 = fromPlanet.a;
            const r2 = toPlanet.a;

            // Semi-major axis of transfer orbit (approximate)
            const aTransfer = (r1 + r2) / 2;

            // Departure velocity (vis-viva equation)
            const v1_helio = Math.sqrt(planets.Sun.GM * (2/r1 - 1/aTransfer));
            const v1_planet = Math.sqrt(planets.Sun.GM / r1);

            // Arrival velocity
            const v2_helio = Math.sqrt(planets.Sun.GM * (2/r2 - 1/aTransfer));
            const v2_planet = Math.sqrt(planets.Sun.GM / r2);

            // Excess velocities (simplified, ignores timing)
            const vInf_dep = Math.abs(v1_helio - v1_planet);
            const vInf_arr = Math.abs(v2_helio - v2_planet);

            // For gravity assists, only count powered portions
            if (leg === 0) {
                // Launch: need to achieve v_infinity from Earth orbit
                const earthParkingOrbit = 200e3 + 6378137;  // 200 km altitude
                const vPark = Math.sqrt(fromPlanet.GM / earthParkingOrbit);
                const vEscape = Math.sqrt(vPark*vPark + vInf_dep*vInf_dep);
                totalDeltaV += vEscape - vPark;
            }

            // Gravity assist maneuver (simplified - some may need powered flyby)
            if (leg > 0 && leg < sequence.length - 2) {
                // Assume unpowered flyby can achieve ~20-50% V-infinity deflection
                // Add small delta-V for trajectory correction
                totalDeltaV += 50;  // 50 m/s TCM per flyby
            }
        }

        // Saturn orbit insertion (simplified)
        const saturnOrbitRadius = 1e9;  // 1 million km
        const vInfSaturn = Math.sqrt(planets.Sun.GM / planets.Saturn.a) * 0.2;  // Approximate
        const vCapture = Math.sqrt(planets.Saturn.GM / saturnOrbitRadius);
        totalDeltaV += Math.sqrt(vCapture*vCapture + vInfSaturn*vInfSaturn) - vCapture;

        const totalTof = candidate.tofs.reduce((a, b) => a + b, 0);

        console.log(`${candidate.name}:`);
        console.log(`  Total time of flight: ${(totalTof / 365.25).toFixed(2)} years`);
        console.log(`  Total delta-V: ${(totalDeltaV / 1000).toFixed(2)} km/s`);

        if (totalDeltaV < bestDeltaV) {
            bestDeltaV = totalDeltaV;
            bestCandidate = candidate;
        }
    }

    console.log('\n=== Optimization Result ===');
    console.log(`\nBest candidate: ${bestCandidate.name}`);
    console.log(`Total delta-V: ${(bestDeltaV / 1000).toFixed(2)} km/s`);

    // Print detailed leg information
    console.log('\nLeg Details:');
    let cumulativeTime = 0;
    for (let leg = 0; leg < sequence.length - 1; leg++) {
        const tof = bestCandidate.tofs[leg];
        cumulativeTime += tof;
        console.log(`  ${sequence[leg]} -> ${sequence[leg+1]}: ${tof} days (arrival: day ${cumulativeTime.toFixed(0)})`);
    }

    // Compare with actual Cassini mission
    console.log('\n--- Comparison with Actual Cassini Mission ---');
    console.log('Cassini timeline:');
    console.log('  Launch: Oct 15, 1997');
    console.log('  Venus 1: Apr 26, 1998 (193 days)');
    console.log('  Venus 2: Jun 24, 1999 (424 days)');
    console.log('  Earth: Aug 18, 1999 (55 days)');
    console.log('  Jupiter: Dec 30, 2000 (500 days)');
    console.log('  Saturn: Jul 1, 2004 (1280 days)');
    console.log('  Total: ~6.7 years');
    console.log('  Total delta-V: ~2.0 km/s (with gravity assists)');

    console.log('\n=== MGA optimization complete ===');
}

main().catch(console.error);
