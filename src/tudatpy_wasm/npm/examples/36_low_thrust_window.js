/**
 * Example 36: Low-Thrust Earth-Mars Transfer Window
 *
 * Ported from: examples/tudatpy/mission_design/low_thrust_earth_mars_transfer_window.py
 *
 * This example demonstrates analysis of low-thrust transfer windows
 * between Earth and Mars using electric propulsion.
 *
 * Key concepts:
 * - Low-thrust trajectory optimization
 * - Electric propulsion mission design
 * - Spiral trajectory approximation
 * - Thrust-coast-thrust profiles
 * - Specific impulse vs thrust tradeoffs
 *
 * Run with: node 36_low_thrust_window.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Low-Thrust Earth-Mars Transfer Window ===\n');

    const tudat = await createTudatModule();

    // Constants
    const AU = 1.496e11;                // Astronomical unit [m]
    const GM_SUN = 1.32712440018e20;    // Sun gravitational parameter [m^3/s^2]
    const JULIAN_DAY = 86400;
    const g0 = 9.80665;                 // Standard gravity [m/s²]

    // Spacecraft parameters
    const spacecraft = {
        initialMass: 1500,      // kg
        dryMass: 500,           // kg (payload + structure)
        thrust: 0.2,            // N (typical ion engine)
        Isp: 3000,              // s (specific impulse)
        power: 10000            // W (solar power at 1 AU)
    };

    // Planetary orbital elements
    const earth = { a: 1.000 * AU, T: 365.25 * JULIAN_DAY };
    const mars = { a: 1.524 * AU, T: 687 * JULIAN_DAY };

    console.log('Spacecraft Configuration:');
    console.log(`  Initial mass: ${spacecraft.initialMass} kg`);
    console.log(`  Dry mass: ${spacecraft.dryMass} kg`);
    console.log(`  Propellant: ${spacecraft.initialMass - spacecraft.dryMass} kg`);
    console.log(`  Thrust: ${spacecraft.thrust * 1000} mN`);
    console.log(`  Specific impulse: ${spacecraft.Isp} s`);
    console.log(`  Power: ${spacecraft.power / 1000} kW`);

    // Compute maximum delta-V from rocket equation
    const maxDeltaV = spacecraft.Isp * g0 * Math.log(spacecraft.initialMass / spacecraft.dryMass);
    console.log(`\nMaximum delta-V: ${(maxDeltaV / 1000).toFixed(2)} km/s`);

    // Initial acceleration
    const initialAccel = spacecraft.thrust / spacecraft.initialMass;
    console.log(`Initial acceleration: ${(initialAccel * 1000).toFixed(4)} mm/s²`);

    // Exhaust velocity
    const exhaustVel = spacecraft.Isp * g0;
    console.log(`Exhaust velocity: ${(exhaustVel / 1000).toFixed(1)} km/s`);

    // Mass flow rate
    const mdot = spacecraft.thrust / exhaustVel;
    console.log(`Mass flow rate: ${(mdot * 1e6).toFixed(3)} mg/s`);

    // Synodic period
    const synodicPeriod = 1 / Math.abs(1/earth.T - 1/mars.T);

    console.log('\n--- Low-Thrust Transfer Analysis ---\n');

    /**
     * For low-thrust transfers, we use a simplified spiral trajectory model.
     * The spacecraft spirals outward from Earth's orbit to Mars' orbit.
     *
     * Key equations:
     * - Δv ≈ |√(μ/r1) - √(μ/r2)| for tangential thrust
     * - TOF depends on thrust level and trajectory optimization
     */

    // Edelbaum approximation for low-thrust transfer
    function edelbaum_dv(r1, r2, delta_i = 0) {
        // Edelbaum's equation for low-thrust transfer between circular coplanar orbits
        const v1 = Math.sqrt(GM_SUN / r1);
        const v2 = Math.sqrt(GM_SUN / r2);

        // For coplanar case (delta_i = 0):
        // Δv = |v1 - v2| for pure tangential thrust
        // With plane change: Δv = √(v1² + v2² - 2v1v2cos(πδi/2))

        if (delta_i === 0) {
            return Math.abs(v1 - v2);
        } else {
            const cos_term = Math.cos(Math.PI * delta_i / 2);
            return Math.sqrt(v1*v1 + v2*v2 - 2*v1*v2*cos_term);
        }
    }

    // Minimum delta-V for Earth-Mars transfer
    const minDeltaV = edelbaum_dv(earth.a, mars.a);
    console.log(`Edelbaum minimum Δv: ${(minDeltaV / 1000).toFixed(2)} km/s`);

    // With typical plane change (Mars inclination ~1.85°)
    const marsInclination = 1.85 * Math.PI / 180;
    const deltaVWithPlaneChange = edelbaum_dv(earth.a, mars.a, marsInclination);
    console.log(`With plane change (${(marsInclination * 180/Math.PI).toFixed(2)}°): ${(deltaVWithPlaneChange / 1000).toFixed(2)} km/s`);

    // Estimate transfer time for continuous thrust
    function estimateTOF(deltaV, initialMass, thrust, Isp) {
        // For continuous tangential thrust:
        // TOF ≈ m0 * (1 - exp(-Δv/ve)) / mdot

        const ve = Isp * g0;
        const mdot = thrust / ve;
        const finalMass = initialMass * Math.exp(-deltaV / ve);
        const propellantMass = initialMass - finalMass;

        // Time = propellant mass / mass flow rate
        const burnTime = propellantMass / mdot;

        return {
            burnTime,
            finalMass,
            propellantUsed: propellantMass
        };
    }

    // Compute transfer parameters
    const transferParams = estimateTOF(minDeltaV, spacecraft.initialMass, spacecraft.thrust, spacecraft.Isp);

    console.log(`\nContinuous thrust transfer:`);
    console.log(`  Burn time: ${(transferParams.burnTime / JULIAN_DAY).toFixed(1)} days`);
    console.log(`  Final mass: ${transferParams.finalMass.toFixed(1)} kg`);
    console.log(`  Propellant used: ${transferParams.propellantUsed.toFixed(1)} kg`);

    // Add coast phases for realistic mission
    // Low-thrust missions often use thrust-coast-thrust profiles
    const coastFraction = 0.3;  // 30% coasting
    const totalTOF = transferParams.burnTime / (1 - coastFraction);
    console.log(`  With ${coastFraction * 100}% coast: ${(totalTOF / JULIAN_DAY).toFixed(1)} days total`);

    // Analyze transfer window sensitivity
    console.log('\n--- Transfer Window Sensitivity ---\n');

    // Phase angle at departure affects required delta-V
    function analyzePhaseAngle(phaseAngle) {
        // Additional delta-V needed if not optimal phase
        // Simplified: extra ~10% per 30° off optimal
        const optimalPhase = 44;  // degrees (Mars ahead of Earth)
        const phaseDiff = Math.abs(phaseAngle - optimalPhase);
        const penalty = 1 + 0.1 * Math.min(phaseDiff / 30, 2);
        return minDeltaV * penalty;
    }

    console.log('Phase angle sensitivity:');
    console.log('Phase[deg] | ΔV[km/s] | Propellant[kg] | Feasible');
    console.log('-----------|----------|----------------|----------');

    for (let phase = 0; phase <= 180; phase += 30) {
        const dv = analyzePhaseAngle(phase);
        const params = estimateTOF(dv, spacecraft.initialMass, spacecraft.thrust, spacecraft.Isp);
        const feasible = params.finalMass > spacecraft.dryMass;
        console.log(`${phase.toString().padStart(10)} | ${(dv / 1000).toFixed(2).padStart(8)} | ${params.propellantUsed.toFixed(1).padStart(14)} | ${feasible ? 'Yes' : 'No'}`);
    }

    // Optimal departure windows
    console.log('\n--- Optimal Departure Windows ---');

    // Window occurs when phase angle is optimal
    // Mars moves ~0.524°/day, Earth moves ~0.986°/day
    // Relative motion: ~0.462°/day
    const optimalPhaseAngle = 44;  // degrees
    const phaseRate = 360 / (synodicPeriod / JULIAN_DAY);  // deg/day
    const windowWidth = 60 / phaseRate;  // days for ±30° window

    console.log(`Optimal phase angle: ${optimalPhaseAngle}°`);
    console.log(`Phase rate: ${phaseRate.toFixed(3)}°/day`);
    console.log(`Window width (±30°): ~${windowWidth.toFixed(0)} days`);
    console.log(`Window occurs every: ${(synodicPeriod / JULIAN_DAY).toFixed(0)} days`);

    // Compare thrust levels
    console.log('\n--- Thrust Level Comparison ---');
    console.log('Thrust[mN] | TOF[days] | Prop[kg] | Final mass[kg]');
    console.log('-----------|-----------|----------|---------------');

    const thrustLevels = [50, 100, 200, 300, 500];
    for (const thrust of thrustLevels) {
        const params = estimateTOF(minDeltaV, spacecraft.initialMass, thrust / 1000, spacecraft.Isp);
        const tofDays = params.burnTime / JULIAN_DAY / (1 - coastFraction);
        console.log(`${thrust.toString().padStart(10)} | ${tofDays.toFixed(0).padStart(9)} | ${params.propellantUsed.toFixed(1).padStart(8)} | ${params.finalMass.toFixed(1).padStart(14)}`);
    }

    // Compare specific impulse
    console.log('\n--- Specific Impulse Comparison ---');
    console.log('Isp[s]  | Exhaust[km/s] | Prop[kg] | Final mass[kg]');
    console.log('--------|---------------|----------|---------------');

    const ispLevels = [1500, 2000, 3000, 4000, 5000];
    for (const isp of ispLevels) {
        const params = estimateTOF(minDeltaV, spacecraft.initialMass, spacecraft.thrust, isp);
        const ve = isp * g0 / 1000;
        console.log(`${isp.toString().padStart(7)} | ${ve.toFixed(1).padStart(13)} | ${params.propellantUsed.toFixed(1).padStart(8)} | ${params.finalMass.toFixed(1).padStart(14)}`);
    }

    // Mission comparison with chemical propulsion
    console.log('\n--- Comparison with Chemical Propulsion ---');

    // Chemical (bi-propellant)
    const chemIsp = 320;  // s
    const chemDeltaV = 4000;  // m/s (Hohmann + margin)
    const chemMassRatio = Math.exp(chemDeltaV / (chemIsp * g0));
    const chemPropellant = spacecraft.initialMass * (1 - 1/chemMassRatio);
    const chemFinalMass = spacecraft.initialMass - chemPropellant;

    console.log('\nChemical (Isp=320s):');
    console.log(`  Required Δv: ${chemDeltaV / 1000} km/s`);
    console.log(`  Propellant: ${chemPropellant.toFixed(0)} kg`);
    console.log(`  Final mass: ${chemFinalMass.toFixed(0)} kg`);
    console.log(`  TOF: ~260 days (Hohmann)`);

    console.log('\nElectric (Isp=3000s):');
    console.log(`  Required Δv: ${(minDeltaV / 1000).toFixed(2)} km/s`);
    console.log(`  Propellant: ${transferParams.propellantUsed.toFixed(0)} kg`);
    console.log(`  Final mass: ${transferParams.finalMass.toFixed(0)} kg`);
    console.log(`  TOF: ~${(totalTOF / JULIAN_DAY).toFixed(0)} days`);

    const massSaved = chemPropellant - transferParams.propellantUsed;
    console.log(`\nPropellant savings: ${massSaved.toFixed(0)} kg (${(massSaved/chemPropellant*100).toFixed(0)}%)`);
    console.log(`Time penalty: +${((totalTOF / JULIAN_DAY) - 260).toFixed(0)} days`);

    // Real mission examples
    console.log('\n--- Reference: Dawn Mission ---');
    console.log('Dawn spacecraft (Earth-Vesta-Ceres):');
    console.log('  Ion thrusters: 3 × NSTAR');
    console.log('  Thrust: 92 mN each');
    console.log('  Isp: 3100 s');
    console.log('  Xenon propellant: 425 kg');
    console.log('  Total Δv: 11.5 km/s');

    console.log('\n=== Low-Thrust Transfer Window Analysis Complete ===');
}

main().catch(console.error);
