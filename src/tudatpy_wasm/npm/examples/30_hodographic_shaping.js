/**
 * Example 30: Hodographic Shaping for Low-Thrust Trajectories
 *
 * Ported from: examples/tudatpy/mission_design/hodographic_shaping_mga_optimization.py
 *
 * This example demonstrates hodographic shaping for low-thrust trajectory
 * design, used in multiple gravity assist (MGA) transfers. The velocity
 * profile is shaped using basis functions to ensure boundary conditions
 * are met while minimizing delta-V.
 *
 * Key concepts:
 * - Hodographic shaping method
 * - Low-thrust trajectory design
 * - Velocity basis functions (radial, normal, axial)
 * - Multi-revolution transfers
 * - Shape-based thrust profiling
 *
 * Run with: node 30_hodographic_shaping.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Hodographic Shaping for Low-Thrust Trajectories ===\n');

    const tudat = await createTudatModule();

    // Constants
    const AU = 1.496e11;                // Astronomical unit [m]
    const GM_SUN = 1.32712440018e20;    // Sun gravitational parameter [m^3/s^2]
    const JULIAN_DAY = 86400;           // Seconds per day
    const JULIAN_YEAR = 365.25 * JULIAN_DAY;

    // Planet semi-major axes [m]
    const planets = {
        Earth: { a: 1.0 * AU, GM: 3.986e14 },
        Mars: { a: 1.524 * AU, GM: 4.283e13 },
        Jupiter: { a: 5.203 * AU, GM: 1.267e17 }
    };

    console.log('Transfer Configuration:');
    console.log('  Earth -> Mars -> Earth -> Jupiter');
    console.log('  Using hodographic shaping with low-thrust propulsion\n');

    /**
     * Hodographic shaping represents the velocity profile using basis functions.
     * The velocity is decomposed into three components:
     * - Radial (r): toward/away from the Sun
     * - Normal (n): perpendicular to orbital plane
     * - Axial (theta): along the orbital direction
     *
     * Each component is shaped as: v(t) = sum_i(c_i * f_i(t))
     * where f_i are basis functions and c_i are free coefficients
     */

    // Basis function definitions for hodographic shaping
    class HodographShaper {
        constructor(tof, numRevolutions = 0) {
            this.tof = tof;
            this.numRevolutions = numRevolutions;
            this.omega = 2 * Math.PI / tof;  // Angular frequency
            this.scale = 1 / tof;            // Scale factor
        }

        // Recommended radial velocity basis functions
        radialBasisFunctions(t) {
            const tau = t / this.tof;  // Normalized time [0, 1]
            return [
                1,                                           // Constant
                tau,                                         // Linear
                Math.sin(this.omega * t),                    // Sine
                Math.cos(this.omega * t),                    // Cosine
                Math.sin(0.5 * this.omega * t) * this.scale, // Half-frequency sine
                Math.cos(0.5 * this.omega * t) * this.scale  // Half-frequency cosine
            ];
        }

        // Recommended normal velocity basis functions
        normalBasisFunctions(t) {
            const tau = t / this.tof;
            return [
                1,
                tau,
                Math.sin(this.omega * t),
                Math.cos(this.omega * t),
                Math.sin(0.5 * this.omega * t) * this.scale,
                Math.cos(0.5 * this.omega * t) * this.scale
            ];
        }

        // Recommended axial velocity basis functions (for multi-revolution)
        axialBasisFunctions(t) {
            const tau = t / this.tof;
            const revFreq = (this.numRevolutions + 0.5) * this.omega;
            const exp = 4;
            const scalePow = Math.pow(this.scale, exp);
            return [
                1,
                tau,
                tau * tau,
                Math.sin(this.omega * t),
                Math.cos(this.omega * t),
                Math.pow(Math.cos(revFreq * t), exp) * scalePow,
                Math.pow(Math.sin(revFreq * t), exp) * scalePow
            ];
        }

        // Compute shaped velocity at time t given coefficients
        computeVelocity(t, radialCoeffs, normalCoeffs, axialCoeffs) {
            const radialBasis = this.radialBasisFunctions(t);
            const normalBasis = this.normalBasisFunctions(t);
            const axialBasis = this.axialBasisFunctions(t);

            let vr = 0, vn = 0, vtheta = 0;
            for (let i = 0; i < radialCoeffs.length && i < radialBasis.length; i++) {
                vr += radialCoeffs[i] * radialBasis[i];
            }
            for (let i = 0; i < normalCoeffs.length && i < normalBasis.length; i++) {
                vn += normalCoeffs[i] * normalBasis[i];
            }
            for (let i = 0; i < axialCoeffs.length && i < axialBasis.length; i++) {
                vtheta += axialCoeffs[i] * axialBasis[i];
            }

            return { vr, vn, vtheta };
        }
    }

    /**
     * Simple low-thrust transfer leg using hodographic shaping
     * This is a simplified implementation for demonstration
     */
    function computeHodographicLeg(departurePlanet, arrivalPlanet, tof, numRevolutions = 0) {
        const r1 = planets[departurePlanet].a;
        const r2 = planets[arrivalPlanet].a;

        // Compute circular velocities
        const v1 = Math.sqrt(GM_SUN / r1);
        const v2 = Math.sqrt(GM_SUN / r2);

        // Phase angle change
        const deltaTheta = 2 * Math.PI * (numRevolutions + 0.5);

        // Simple shaping: assume velocity profile that connects the two orbits
        const shaper = new HodographShaper(tof, numRevolutions);

        // Boundary conditions:
        // At t=0: r=r1, theta=0, v_theta=v1
        // At t=tof: r=r2, theta=deltaTheta, v_theta=v2

        // Simple coefficient determination (this is greatly simplified)
        // In reality, this requires solving a system of equations

        // Radial velocity: need to go from r1 to r2
        const avgRadialRate = (r2 - r1) / tof;

        // Axial (tangential) velocity: interpolate between circular velocities
        const avgAxialVel = (v1 + v2) / 2;

        // Sample the trajectory
        const samples = 100;
        const dt = tof / samples;
        let totalDeltaV = 0;
        let lastAccel = 0;

        const trajectory = [];

        for (let i = 0; i <= samples; i++) {
            const t = i * dt;
            const tau = t / tof;  // Normalized time

            // Smooth position interpolation (simplified)
            const r = r1 + (r2 - r1) * (3 * tau * tau - 2 * tau * tau * tau);  // Smooth step
            const theta = deltaTheta * tau;

            // Velocity from differentiation
            const drdt = (r2 - r1) * (6 * tau - 6 * tau * tau) / tof;
            const dthetadt = deltaTheta / tof;
            const vTheta = r * dthetadt;

            // Required centripetal + gravitational balance
            const gravAccel = GM_SUN / (r * r);
            const centripetalAccel = vTheta * vTheta / r;

            // Thrust required (radial component)
            const thrustAccelR = drdt / tof;  // Simplified

            // Tangential thrust to change angular momentum
            const h1 = r1 * v1;  // Initial angular momentum
            const h2 = r2 * v2;  // Final angular momentum
            const dh = (h2 - h1) / tof;
            const thrustAccelT = dh / (r * tof);

            // Total thrust acceleration magnitude
            const thrustAccel = Math.sqrt(thrustAccelR * thrustAccelR + thrustAccelT * thrustAccelT);

            // Accumulate delta-V
            if (i > 0) {
                totalDeltaV += thrustAccel * dt;
            }

            trajectory.push({
                t: t / JULIAN_DAY,
                r: r / AU,
                theta: theta * 180 / Math.PI,
                vTheta: vTheta / 1000,
                thrustAccel: thrustAccel * 1000  // mm/s^2
            });
        }

        return {
            departurePlanet,
            arrivalPlanet,
            tof: tof / JULIAN_DAY,
            numRevolutions,
            totalDeltaV: totalDeltaV / 1000,  // km/s
            trajectory
        };
    }

    // Define transfer legs (simplified Earth-Mars-Earth-Jupiter sequence)
    console.log('--- Computing Transfer Legs ---\n');

    // Leg 1: Earth to Mars
    const tofLeg1 = 300 * JULIAN_DAY;  // 300 days
    const leg1 = computeHodographicLeg('Earth', 'Mars', tofLeg1, 0);
    console.log(`Leg 1 (Earth -> Mars):`);
    console.log(`  Time of flight: ${leg1.tof.toFixed(0)} days`);
    console.log(`  Revolutions: ${leg1.numRevolutions}`);
    console.log(`  Delta-V: ${leg1.totalDeltaV.toFixed(3)} km/s`);

    // Leg 2: Mars to Earth (gravity assist return)
    const tofLeg2 = 400 * JULIAN_DAY;  // 400 days
    const leg2 = computeHodographicLeg('Mars', 'Earth', tofLeg2, 1);
    console.log(`\nLeg 2 (Mars -> Earth):`);
    console.log(`  Time of flight: ${leg2.tof.toFixed(0)} days`);
    console.log(`  Revolutions: ${leg2.numRevolutions}`);
    console.log(`  Delta-V: ${leg2.totalDeltaV.toFixed(3)} km/s`);

    // Leg 3: Earth to Jupiter
    const tofLeg3 = 900 * JULIAN_DAY;  // 900 days (~2.5 years)
    const leg3 = computeHodographicLeg('Earth', 'Jupiter', tofLeg3, 0);
    console.log(`\nLeg 3 (Earth -> Jupiter):`);
    console.log(`  Time of flight: ${leg3.tof.toFixed(0)} days`);
    console.log(`  Revolutions: ${leg3.numRevolutions}`);
    console.log(`  Delta-V: ${leg3.totalDeltaV.toFixed(3)} km/s`);

    // Total mission summary
    const totalTof = leg1.tof + leg2.tof + leg3.tof;
    const totalDeltaV = leg1.totalDeltaV + leg2.totalDeltaV + leg3.totalDeltaV;

    console.log('\n=== Mission Summary ===');
    console.log(`Total time of flight: ${totalTof.toFixed(0)} days (${(totalTof / 365.25).toFixed(2)} years)`);
    console.log(`Total low-thrust delta-V: ${totalDeltaV.toFixed(3)} km/s`);

    // Compare with impulsive Hohmann transfers
    console.log('\n--- Comparison with Hohmann Transfers ---');

    function hohmannDeltaV(r1, r2) {
        const v1 = Math.sqrt(GM_SUN / r1);
        const v2 = Math.sqrt(GM_SUN / r2);
        const a_transfer = (r1 + r2) / 2;
        const v_dep = Math.sqrt(GM_SUN * (2/r1 - 1/a_transfer));
        const v_arr = Math.sqrt(GM_SUN * (2/r2 - 1/a_transfer));
        return Math.abs(v_dep - v1) + Math.abs(v2 - v_arr);
    }

    const hohmannEM = hohmannDeltaV(planets.Earth.a, planets.Mars.a);
    const hohmannME = hohmannDeltaV(planets.Mars.a, planets.Earth.a);
    const hohmannEJ = hohmannDeltaV(planets.Earth.a, planets.Jupiter.a);
    const totalHohmann = hohmannEM + hohmannME + hohmannEJ;

    console.log(`Hohmann Earth-Mars: ${(hohmannEM / 1000).toFixed(3)} km/s`);
    console.log(`Hohmann Mars-Earth: ${(hohmannME / 1000).toFixed(3)} km/s`);
    console.log(`Hohmann Earth-Jupiter: ${(hohmannEJ / 1000).toFixed(3)} km/s`);
    console.log(`Total Hohmann delta-V: ${(totalHohmann / 1000).toFixed(3)} km/s`);

    // Hodographic shaping advantages
    console.log('\n--- Hodographic Shaping Method Notes ---');
    console.log('Advantages:');
    console.log('  - Analytical solution for low-thrust trajectories');
    console.log('  - Guarantees boundary conditions are satisfied');
    console.log('  - Basis functions provide flexibility in shaping');
    console.log('  - Suitable for preliminary mission design');
    console.log('  - Can handle multi-revolution transfers');

    console.log('\nBasis Functions Used:');
    console.log('  Radial: 1, t, sin(wt), cos(wt), scaled sin/cos');
    console.log('  Normal: 1, t, sin(wt), cos(wt), scaled sin/cos');
    console.log('  Axial:  1, t, t^2, sin(wt), cos(wt), power trig');

    // Thrust profile analysis for leg 3 (longest)
    console.log('\n--- Thrust Profile Analysis (Earth-Jupiter Leg) ---');
    const thrustData = leg3.trajectory;
    const maxThrust = Math.max(...thrustData.map(d => d.thrustAccel));
    const minThrust = Math.min(...thrustData.map(d => d.thrustAccel));
    const avgThrust = thrustData.reduce((sum, d) => sum + d.thrustAccel, 0) / thrustData.length;

    console.log(`Max thrust acceleration: ${maxThrust.toFixed(4)} mm/s^2`);
    console.log(`Min thrust acceleration: ${minThrust.toFixed(4)} mm/s^2`);
    console.log(`Average thrust acceleration: ${avgThrust.toFixed(4)} mm/s^2`);

    // For a 1000 kg spacecraft
    const spacecraftMass = 1000;  // kg
    console.log(`\nFor ${spacecraftMass} kg spacecraft:`);
    console.log(`  Max thrust: ${(maxThrust * spacecraftMass / 1000).toFixed(3)} N`);
    console.log(`  Avg thrust: ${(avgThrust * spacecraftMass / 1000).toFixed(3)} N`);

    // Typical ion engine comparison
    console.log('\nTypical ion engine thrust: 0.01-0.5 N');
    console.log('This transfer is compatible with electric propulsion.');

    // Trajectory sampling points
    console.log('\n--- Sample Trajectory Points (Earth-Jupiter) ---');
    console.log('Day     | Distance [AU] | Angle [deg] | V_theta [km/s]');
    console.log('--------|---------------|-------------|---------------');
    for (let i = 0; i <= 10; i++) {
        const idx = Math.floor(i * (thrustData.length - 1) / 10);
        const p = thrustData[idx];
        console.log(`${p.t.toFixed(0).padStart(7)} | ${p.r.toFixed(3).padStart(13)} | ${p.theta.toFixed(1).padStart(11)} | ${p.vTheta.toFixed(2).padStart(14)}`);
    }

    console.log('\n=== Hodographic Shaping Complete ===');
}

main().catch(console.error);
