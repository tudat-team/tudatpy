/**
 * Example 41: Multiple Gravity Assist Trajectories
 *
 * Ported from: examples/tudatpy/mission_design/mga_trajectories.py
 *
 * This example demonstrates how Multiple Gravity Assist (MGA) transfer
 * trajectories can be simulated and analyzed. Three types of transfers
 * are shown:
 * - High-thrust transfer with unpowered legs
 * - High-thrust transfer with Deep Space Maneuvers (DSMs)
 * - Low-thrust transfer overview
 *
 * Key concepts:
 * - MGA trajectory design
 * - Unpowered vs powered legs
 * - Deep space maneuvers
 * - Delta-V budget computation
 * - Flyby altitude and turn angle
 *
 * Run with: node 41_mga_trajectories.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Multiple Gravity Assist Trajectories ===\n');

    const tudat = await createTudatModule();

    // Constants
    const AU = 1.496e11;                // Astronomical unit [m]
    const GM_SUN = 1.32712440018e20;    // Sun gravitational parameter [m^3/s^2]
    const JULIAN_DAY = 86400;           // Seconds per day

    // Planet data (semi-major axis, GM, radius)
    const planets = {
        Venus: { a: 0.723 * AU, GM: 3.249e14, R: 6.052e6, name: 'Venus' },
        Earth: { a: 1.000 * AU, GM: 3.986e14, R: 6.371e6, name: 'Earth' },
        Mars: { a: 1.524 * AU, GM: 4.283e13, R: 3.390e6, name: 'Mars' },
        Jupiter: { a: 5.203 * AU, GM: 1.267e17, R: 6.991e7, name: 'Jupiter' },
        Saturn: { a: 9.537 * AU, GM: 3.793e16, R: 5.823e7, name: 'Saturn' }
    };

    /**
     * Solve Lambert problem for transfer between two planets
     */
    function solveLambert(r1, r2, tof, GM, isRetrograde = false) {
        // Simple Lambert solver (single revolution)
        const r1Mag = Math.sqrt(r1[0]*r1[0] + r1[1]*r1[1] + r1[2]*r1[2]);
        const r2Mag = Math.sqrt(r2[0]*r2[0] + r2[1]*r2[1] + r2[2]*r2[2]);

        // Chord and semi-perimeter
        const c = Math.sqrt(
            (r2[0]-r1[0])**2 + (r2[1]-r1[1])**2 + (r2[2]-r1[2])**2
        );
        const s = (r1Mag + r2Mag + c) / 2;

        // Minimum energy ellipse
        const aMin = s / 2;
        const tofMin = Math.PI * Math.sqrt(aMin**3 / GM);

        // Estimate semi-major axis (iterate for accuracy)
        let a = aMin * Math.pow(tof / tofMin, 2/3);

        // Cross product for transfer plane normal
        const h = [
            r1[1]*r2[2] - r1[2]*r2[1],
            r1[2]*r2[0] - r1[0]*r2[2],
            r1[0]*r2[1] - r1[1]*r2[0]
        ];
        const hMag = Math.sqrt(h[0]**2 + h[1]**2 + h[2]**2);

        // Transfer angle
        const cosTheta = (r1[0]*r2[0] + r1[1]*r2[1] + r1[2]*r2[2]) / (r1Mag * r2Mag);
        let theta = Math.acos(Math.max(-1, Math.min(1, cosTheta)));
        if (isRetrograde ? h[2] > 0 : h[2] < 0) {
            theta = 2 * Math.PI - theta;
        }

        // Lagrange coefficients (simplified)
        const p = a * (1 - ((s - r1Mag - r2Mag) / c)**2);
        const f = 1 - r2Mag / p * (1 - cosTheta);
        const g = r1Mag * r2Mag * Math.sin(theta) / Math.sqrt(GM * p);

        // Departure velocity
        const v1 = [
            (r2[0] - f * r1[0]) / g,
            (r2[1] - f * r1[1]) / g,
            (r2[2] - f * r1[2]) / g
        ];

        // Arrival velocity (from energy)
        const v1Mag = Math.sqrt(v1[0]**2 + v1[1]**2 + v1[2]**2);
        const v2Mag = Math.sqrt(2 * GM / r2Mag - GM / a);

        // Direction at arrival
        const fDot = -Math.sqrt(GM / p) * Math.tan(theta / 2) *
                     ((1 - cosTheta) / p - 1/r1Mag - 1/r2Mag);
        const gDot = 1 - r1Mag / p * (1 - cosTheta);

        const v2 = [
            fDot * r1[0] + gDot * v1[0],
            fDot * r1[1] + gDot * v1[1],
            fDot * r1[2] + gDot * v1[2]
        ];

        return { v1, v2, a, theta };
    }

    /**
     * Get planet position at time t (circular coplanar approximation)
     */
    function getPlanetPosition(planet, t) {
        const T = 2 * Math.PI * Math.sqrt(planet.a**3 / GM_SUN);
        const theta = 2 * Math.PI * t / T;
        return [
            planet.a * Math.cos(theta),
            planet.a * Math.sin(theta),
            0
        ];
    }

    /**
     * Get planet velocity at time t
     */
    function getPlanetVelocity(planet, t) {
        const T = 2 * Math.PI * Math.sqrt(planet.a**3 / GM_SUN);
        const v = Math.sqrt(GM_SUN / planet.a);
        const theta = 2 * Math.PI * t / T;
        return [
            -v * Math.sin(theta),
            v * Math.cos(theta),
            0
        ];
    }

    /**
     * Compute flyby parameters
     */
    function computeFlyby(vInf_in, vInf_out, GM, R_min) {
        // V-infinity magnitudes
        const vInfInMag = Math.sqrt(vInf_in[0]**2 + vInf_in[1]**2 + vInf_in[2]**2);
        const vInfOutMag = Math.sqrt(vInf_out[0]**2 + vInf_out[1]**2 + vInf_out[2]**2);

        // Turn angle from dot product
        const dot = vInf_in[0]*vInf_out[0] + vInf_in[1]*vInf_out[1] + vInf_in[2]*vInf_out[2];
        const cosAlpha = dot / (vInfInMag * vInfOutMag);
        const alpha = Math.acos(Math.max(-1, Math.min(1, cosAlpha)));

        // Required periapsis for this turn angle (average v-infinity)
        const vInfAvg = (vInfInMag + vInfOutMag) / 2;
        const sinHalfAlpha = Math.sin(alpha / 2);
        const e = 1 / sinHalfAlpha;  // Eccentricity for hyperbolic flyby
        const rp = GM / (vInfAvg * vInfAvg) * (e - 1);

        // Delta-V from powered flyby (if needed to change v-infinity magnitude)
        const deltaVInf = Math.abs(vInfOutMag - vInfInMag);

        return {
            turnAngle: alpha * 180 / Math.PI,
            periapsis: rp,
            altitude: rp - R_min,
            vInfIn: vInfInMag,
            vInfOut: vInfOutMag,
            deltaV: deltaVInf
        };
    }

    // ========================================
    // PART 1: MGA with Unpowered Legs (EVVEJS)
    // ========================================
    console.log('--- Part 1: MGA Transfer with Unpowered Legs ---');
    console.log('Sequence: Earth -> Venus -> Venus -> Earth -> Jupiter -> Saturn\n');

    const sequence = ['Earth', 'Venus', 'Venus', 'Earth', 'Jupiter', 'Saturn'];

    // Times of flight for each leg [days]
    const tofs = [150, 300, 400, 800, 1500];

    // Departure epoch (arbitrary)
    let t = 0;
    let totalDeltaV = 0;
    const legs = [];

    console.log('Leg-by-Leg Analysis:\n');

    for (let i = 0; i < sequence.length - 1; i++) {
        const dep = planets[sequence[i]];
        const arr = planets[sequence[i + 1]];
        const tof = tofs[i] * JULIAN_DAY;

        // Get positions and velocities
        const r1 = getPlanetPosition(dep, t);
        const r2 = getPlanetPosition(arr, t + tof);
        const vPlanet1 = getPlanetVelocity(dep, t);
        const vPlanet2 = getPlanetVelocity(arr, t + tof);

        // Solve Lambert
        const lambert = solveLambert(r1, r2, tof, GM_SUN);

        // V-infinity at departure and arrival
        const vInfDep = [
            lambert.v1[0] - vPlanet1[0],
            lambert.v1[1] - vPlanet1[1],
            lambert.v1[2] - vPlanet1[2]
        ];
        const vInfArr = [
            lambert.v2[0] - vPlanet2[0],
            lambert.v2[1] - vPlanet2[1],
            lambert.v2[2] - vPlanet2[2]
        ];

        const vInfDepMag = Math.sqrt(vInfDep[0]**2 + vInfDep[1]**2 + vInfDep[2]**2);
        const vInfArrMag = Math.sqrt(vInfArr[0]**2 + vInfArr[1]**2 + vInfArr[2]**2);

        // Delta-V for this leg (departure maneuver for first leg only)
        let legDeltaV = 0;
        if (i === 0) {
            // Escape from LEO (400 km altitude)
            const r_LEO = 6371e3 + 400e3;
            const v_LEO = Math.sqrt(dep.GM / r_LEO);
            const v_escape = Math.sqrt(vInfDepMag**2 + 2 * dep.GM / r_LEO);
            legDeltaV = v_escape - v_LEO;
        }

        legs.push({
            from: sequence[i],
            to: sequence[i + 1],
            tof: tofs[i],
            vInfDep: vInfDepMag / 1000,
            vInfArr: vInfArrMag / 1000,
            deltaV: legDeltaV / 1000
        });

        console.log(`Leg ${i + 1}: ${sequence[i]} -> ${sequence[i + 1]}`);
        console.log(`  TOF: ${tofs[i]} days`);
        console.log(`  V∞ departure: ${(vInfDepMag / 1000).toFixed(2)} km/s`);
        console.log(`  V∞ arrival: ${(vInfArrMag / 1000).toFixed(2)} km/s`);
        if (legDeltaV > 0) {
            console.log(`  ΔV (escape): ${(legDeltaV / 1000).toFixed(2)} km/s`);
        }
        console.log('');

        totalDeltaV += legDeltaV;
        t += tof;
    }

    // Saturn orbit insertion
    const saturnInsertionDV = 1.5;  // km/s (typical)
    totalDeltaV += saturnInsertionDV * 1000;

    console.log(`Saturn Orbit Insertion ΔV: ${saturnInsertionDV.toFixed(2)} km/s`);
    console.log(`\nTotal Mission ΔV: ${(totalDeltaV / 1000).toFixed(2)} km/s`);
    console.log(`Total TOF: ${tofs.reduce((a, b) => a + b, 0)} days (${(tofs.reduce((a, b) => a + b, 0) / 365.25).toFixed(1)} years)`);

    // ========================================
    // PART 2: MGA with Deep Space Maneuvers
    // ========================================
    console.log('\n--- Part 2: MGA with Deep Space Maneuvers (DSMs) ---');
    console.log('Sequence: Earth -> Venus -> Earth (with DSM)\n');

    // Simplified DSM calculation
    // DSM occurs at a fraction of the leg
    const dsmFraction = 0.5;  // DSM at midpoint
    const dsmMagnitude = 0.3;  // km/s (typical DSM)

    console.log('Deep Space Maneuver Concept:');
    console.log('  - DSM location: 50% of leg');
    console.log(`  - DSM magnitude: ${dsmMagnitude} km/s`);
    console.log('  - Allows correction of trajectory without waiting for flyby');
    console.log('  - Useful when gravity assist alone cannot achieve desired change\n');

    // Example: Earth-Venus-Earth with DSM
    const eveDsm = {
        leg1: { from: 'Earth', to: 'Venus', tof: 150 },
        dsm: { deltav: dsmMagnitude, location: 0.5 },
        leg2: { from: 'Venus', to: 'Earth', tof: 250 }
    };

    console.log('EVE Transfer with DSM:');
    console.log(`  Leg 1: Earth -> Venus, TOF = ${eveDsm.leg1.tof} days`);
    console.log(`  DSM at ${eveDsm.dsm.location * 100}% of Leg 2: ΔV = ${eveDsm.dsm.deltav} km/s`);
    console.log(`  Leg 2: Venus -> Earth, TOF = ${eveDsm.leg2.tof} days`);
    console.log(`  Total TOF: ${eveDsm.leg1.tof + eveDsm.leg2.tof} days\n`);

    // Compare with no-DSM solution
    console.log('Comparison (no DSM vs with DSM):');
    console.log('  No DSM: Constrained to specific departure dates and TOFs');
    console.log('  With DSM: More flexibility in trajectory design');
    console.log('  Trade-off: Additional ΔV cost vs mission design flexibility');

    // ========================================
    // PART 3: Low-Thrust MGA Overview
    // ========================================
    console.log('\n--- Part 3: Low-Thrust MGA Overview ---');
    console.log('(Using hodographic shaping - see Example 30 for details)\n');

    console.log('Low-thrust MGA characteristics:');
    console.log('  - Continuous thrusting throughout legs');
    console.log('  - Uses electric propulsion (ion engines)');
    console.log('  - Higher specific impulse (Isp ~ 3000s vs 320s chemical)');
    console.log('  - Lower thrust (mN vs kN)');
    console.log('  - Longer transfer times but less propellant mass\n');

    // Compare propellant mass for same delta-V
    const totalDV = totalDeltaV / 1000;  // km/s
    const m0 = 5000;  // Initial mass [kg]
    const g0 = 9.81;  // m/s^2

    const ispChemical = 320;  // s
    const ispElectric = 3000; // s

    const mfChemical = m0 * Math.exp(-totalDV * 1000 / (ispChemical * g0));
    const mfElectric = m0 * Math.exp(-totalDV * 1000 / (ispElectric * g0));

    const propChemical = m0 - mfChemical;
    const propElectric = m0 - mfElectric;

    console.log(`For ΔV = ${totalDV.toFixed(2)} km/s and m0 = ${m0} kg:`);
    console.log(`  Chemical (Isp=${ispChemical}s): Propellant = ${propChemical.toFixed(0)} kg (${(propChemical/m0*100).toFixed(1)}%)`);
    console.log(`  Electric (Isp=${ispElectric}s): Propellant = ${propElectric.toFixed(0)} kg (${(propElectric/m0*100).toFixed(1)}%)`);
    console.log(`  Mass savings: ${(propChemical - propElectric).toFixed(0)} kg`);

    // ========================================
    // PART 4: Summary Tables
    // ========================================
    console.log('\n--- Summary ---\n');

    console.log('MGA Trajectory Types:');
    console.log('┌─────────────────────┬──────────────┬────────────┬─────────────┐');
    console.log('│ Type                │ Propulsion   │ ΔV Budget  │ TOF         │');
    console.log('├─────────────────────┼──────────────┼────────────┼─────────────┤');
    console.log('│ Unpowered legs      │ Chemical     │ Lowest     │ Constrained │');
    console.log('│ With DSMs           │ Chemical     │ Medium     │ Flexible    │');
    console.log('│ Low-thrust          │ Electric     │ Higher     │ Longer      │');
    console.log('└─────────────────────┴──────────────┴────────────┴─────────────┘');

    console.log('\nKey Design Parameters:');
    console.log('  - Departure epoch: Affects planetary alignment');
    console.log('  - TOF per leg: Determines transfer geometry');
    console.log('  - Flyby altitude: Controls turn angle');
    console.log('  - DSM location/magnitude: Fine-tune trajectory');
    console.log('  - Sequence: Planetary order affects total ΔV');

    console.log('\n=== MGA Trajectories Complete ===');
}

main().catch(console.error);
