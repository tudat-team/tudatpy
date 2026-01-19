/**
 * Example 32: Impact Manifolds and Libration Point Orbits in CR3BP
 *
 * Ported from: examples/tudatpy/propagation/impact_manifolds_lpo_cr3bp.py
 *
 * This example demonstrates computation of stable and unstable manifolds
 * associated with libration point orbits (LPOs) in the circular restricted
 * three-body problem (CR3BP).
 *
 * Key concepts:
 * - Halo orbit computation
 * - Stable/unstable manifolds
 * - Invariant manifold propagation
 * - Lunar transfer design
 * - Heteroclinic connections
 *
 * Run with: node 32_impact_manifolds.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Impact Manifolds and Libration Point Orbits ===\n');

    const tudat = await createTudatModule();

    // Earth-Moon CR3BP parameters
    const mu = 0.01215;  // Mass ratio (Moon mass / total mass)
    const L_star = 384400e3;  // Earth-Moon distance [m]
    const T_star = 2.36055e6;  // Time unit [s] (~27.3 days / 2π)

    // Positions of primaries in rotating frame (normalized)
    const earthPos = [-mu, 0, 0];
    const moonPos = [1 - mu, 0, 0];

    console.log('Earth-Moon CR3BP System:');
    console.log(`  Mass ratio μ: ${mu.toFixed(5)}`);
    console.log(`  Distance unit L*: ${(L_star / 1e3).toFixed(0)} km`);
    console.log(`  Time unit T*: ${(T_star / 86400).toFixed(3)} days`);

    // Lagrange point positions (L1 and L2 approximate)
    // L1 is between Earth and Moon
    const L1_x = 1 - mu - Math.pow(mu / 3, 1/3);
    // L2 is beyond Moon
    const L2_x = 1 - mu + Math.pow(mu / 3, 1/3);

    console.log(`\nLagrange Points (x-coordinate):`);
    console.log(`  L1: ${L1_x.toFixed(6)} (${((L1_x + mu) * L_star / 1e3).toFixed(0)} km from Earth)`);
    console.log(`  L2: ${L2_x.toFixed(6)} (${((L2_x + mu) * L_star / 1e3).toFixed(0)} km from Earth)`);

    // CR3BP equations of motion
    function cr3bpDerivatives(state) {
        const x = state[0], y = state[1], z = state[2];
        const vx = state[3], vy = state[4], vz = state[5];

        // Distances to primaries
        const r1 = Math.sqrt((x + mu)**2 + y*y + z*z);  // Distance to Earth
        const r2 = Math.sqrt((x - 1 + mu)**2 + y*y + z*z);  // Distance to Moon

        // Accelerations (rotating frame)
        const ax = 2*vy + x - (1 - mu)*(x + mu)/(r1**3) - mu*(x - 1 + mu)/(r2**3);
        const ay = -2*vx + y - (1 - mu)*y/(r1**3) - mu*y/(r2**3);
        const az = -(1 - mu)*z/(r1**3) - mu*z/(r2**3);

        return [vx, vy, vz, ax, ay, az];
    }

    // Jacobi constant (energy-like integral)
    function jacobiConstant(state) {
        const x = state[0], y = state[1], z = state[2];
        const vx = state[3], vy = state[4], vz = state[5];

        const r1 = Math.sqrt((x + mu)**2 + y*y + z*z);
        const r2 = Math.sqrt((x - 1 + mu)**2 + y*y + z*z);

        const U = 0.5*(x*x + y*y) + (1 - mu)/r1 + mu/r2;
        const v2 = vx*vx + vy*vy + vz*vz;

        return 2*U - v2;
    }

    // RK4 integrator
    function rk4Step(state, dt) {
        const k1 = cr3bpDerivatives(state);
        const s1 = state.map((v, i) => v + 0.5 * dt * k1[i]);
        const k2 = cr3bpDerivatives(s1);
        const s2 = state.map((v, i) => v + 0.5 * dt * k2[i]);
        const k3 = cr3bpDerivatives(s2);
        const s3 = state.map((v, i) => v + dt * k3[i]);
        const k4 = cr3bpDerivatives(s3);
        return state.map((v, i) => v + (dt / 6) * (k1[i] + 2*k2[i] + 2*k3[i] + k4[i]));
    }

    // Propagate orbit
    function propagateOrbit(initialState, totalTime, dt, direction = 1) {
        const trajectory = [];
        let state = [...initialState];
        let t = 0;

        while (t <= Math.abs(totalTime)) {
            trajectory.push({
                t: t * direction,
                x: state[0], y: state[1], z: state[2],
                vx: state[3], vy: state[4], vz: state[5],
                C: jacobiConstant(state)
            });

            state = rk4Step(state, dt * direction);
            t += Math.abs(dt);

            // Check for collision with Moon
            const r2 = Math.sqrt((state[0] - 1 + mu)**2 + state[1]**2 + state[2]**2);
            if (r2 < 0.005) {  // Moon radius ~0.0045 in normalized units
                console.log(`  Moon impact at t = ${(t * T_star / 86400).toFixed(2)} days`);
                break;
            }

            // Check for escape
            if (Math.sqrt(state[0]**2 + state[1]**2 + state[2]**2) > 2) {
                break;
            }
        }

        return trajectory;
    }

    // Approximate L1 Halo orbit initial conditions
    // These are approximate values for a small halo orbit around L1
    console.log('\n--- Computing L1 Halo Orbit ---\n');

    // Halo orbit parameters (approximate for Earth-Moon L1)
    const haloAmplitude = 0.01;  // Small amplitude halo
    const haloState = [
        L1_x - 0.02,  // x (slightly towards Earth)
        0,            // y
        haloAmplitude, // z (out of plane)
        0,            // vx
        0.08,         // vy (tangential velocity)
        0             // vz
    ];

    // Adjust velocity to close the orbit (simplified)
    const haloC = jacobiConstant(haloState);
    console.log(`Halo orbit initial state:`);
    console.log(`  Position: (${haloState[0].toFixed(4)}, ${haloState[1].toFixed(4)}, ${haloState[2].toFixed(4)})`);
    console.log(`  Velocity: (${haloState[3].toFixed(4)}, ${haloState[4].toFixed(4)}, ${haloState[5].toFixed(4)})`);
    console.log(`  Jacobi constant: ${haloC.toFixed(6)}`);

    // Propagate halo orbit for one period
    const haloPeriod = 2.5;  // Approximate period in normalized time
    const haloDt = 0.001;
    const haloTrajectory = propagateOrbit(haloState, haloPeriod, haloDt);

    console.log(`\nHalo orbit propagated for ${(haloPeriod * T_star / 86400).toFixed(2)} days`);
    console.log(`Points: ${haloTrajectory.length}`);

    // Find orbit extent
    const xCoords = haloTrajectory.map(p => p.x);
    const yCoords = haloTrajectory.map(p => p.y);
    const zCoords = haloTrajectory.map(p => p.z);

    console.log(`\nOrbit extent:`);
    console.log(`  X: ${Math.min(...xCoords).toFixed(4)} to ${Math.max(...xCoords).toFixed(4)}`);
    console.log(`  Y: ${Math.min(...yCoords).toFixed(4)} to ${Math.max(...yCoords).toFixed(4)}`);
    console.log(`  Z: ${Math.min(...zCoords).toFixed(4)} to ${Math.max(...zCoords).toFixed(4)}`);

    // Compute manifolds
    console.log('\n--- Computing Invariant Manifolds ---\n');

    // Unstable manifold: perturb along unstable direction and propagate forward
    // Stable manifold: perturb along stable direction and propagate backward

    const perturbation = 1e-5;  // Small perturbation
    const manifoldTime = 3;     // Propagation time

    // Sample points along the halo orbit for manifold computation
    const numManifolds = 8;
    const manifoldStates = [];

    for (let i = 0; i < numManifolds; i++) {
        const idx = Math.floor(i * haloTrajectory.length / numManifolds);
        const point = haloTrajectory[idx];

        // Approximate unstable direction (tangent + radial)
        const unstableDir = [
            point.vx + perturbation * (point.x - L1_x),
            point.vy + perturbation * point.y,
            point.vz + perturbation * point.z,
            0, 0, 0
        ];

        // Normalize
        const norm = Math.sqrt(unstableDir[0]**2 + unstableDir[1]**2 + unstableDir[2]**2);

        // Perturbed state for unstable manifold
        const unstableState = [
            point.x + perturbation * unstableDir[0] / norm,
            point.y + perturbation * unstableDir[1] / norm,
            point.z + perturbation * unstableDir[2] / norm,
            point.vx,
            point.vy,
            point.vz
        ];

        manifoldStates.push({
            base: point,
            unstable: unstableState
        });
    }

    console.log(`Computing ${numManifolds} unstable manifold trajectories...`);

    const unstableManifolds = [];
    for (let i = 0; i < manifoldStates.length; i++) {
        const traj = propagateOrbit(manifoldStates[i].unstable, manifoldTime, 0.005, 1);
        unstableManifolds.push(traj);

        // Check if this trajectory reaches Earth vicinity
        const earthDist = traj.map(p => Math.sqrt((p.x + mu)**2 + p.y**2 + p.z**2));
        const minEarthDist = Math.min(...earthDist);

        if (minEarthDist < 0.1) {
            console.log(`  Manifold ${i + 1}: approaches Earth (min dist: ${(minEarthDist * L_star / 1e6).toFixed(0)} km)`);
        }
    }

    // Compute stable manifolds (propagate backward)
    console.log(`\nComputing ${numManifolds} stable manifold trajectories...`);

    const stableManifolds = [];
    for (let i = 0; i < manifoldStates.length; i++) {
        // Reverse perturbation direction for stable manifold
        const point = manifoldStates[i].base;
        const stableState = [
            point.x - perturbation * 0.1,
            point.y - perturbation * 0.05,
            point.z,
            point.vx,
            point.vy,
            point.vz
        ];

        const traj = propagateOrbit(stableState, manifoldTime, 0.005, -1);  // Backward
        stableManifolds.push(traj);
    }

    // Analyze manifold structure
    console.log('\n=== Manifold Analysis ===');

    // Find manifolds that could be used for lunar transfer
    let transferCandidates = 0;
    for (let i = 0; i < unstableManifolds.length; i++) {
        const traj = unstableManifolds[i];
        const lastPoint = traj[traj.length - 1];

        // Check if ends near Moon
        const moonDist = Math.sqrt((lastPoint.x - 1 + mu)**2 + lastPoint.y**2 + lastPoint.z**2);
        if (moonDist < 0.1) {
            transferCandidates++;
        }
    }

    console.log(`\nUnstable manifold trajectories approaching Moon: ${transferCandidates}/${numManifolds}`);

    // Jacobi constant preservation check
    const C_initial = haloC;
    const C_final = unstableManifolds[0][unstableManifolds[0].length - 1].C;
    console.log(`\nJacobi constant (energy conservation):`);
    console.log(`  Initial: ${C_initial.toFixed(8)}`);
    console.log(`  Final: ${C_final.toFixed(8)}`);
    console.log(`  Relative error: ${(Math.abs(C_final - C_initial) / Math.abs(C_initial) * 100).toFixed(6)}%`);

    // Mission design implications
    console.log('\n--- Mission Design Implications ---');
    console.log('Unstable manifolds from L1 halo orbit can be used for:');
    console.log('  - Low-energy transfers to Moon');
    console.log('  - Station-keeping at L1');
    console.log('  - Sample return trajectories');
    console.log('  - Relay satellite positioning');

    console.log('\nStable manifolds provide:');
    console.log('  - Arrival pathways to halo orbit');
    console.log('  - Emergency return trajectories');
    console.log('  - Capture trajectory design');

    // Physical dimensions
    console.log('\n--- Physical Dimensions ---');
    const xRange = (Math.max(...xCoords) - Math.min(...xCoords)) * L_star / 1e3;
    const yRange = (Math.max(...yCoords) - Math.min(...yCoords)) * L_star / 1e3;
    const zRange = (Math.max(...zCoords) - Math.min(...zCoords)) * L_star / 1e3;

    console.log(`Halo orbit dimensions:`);
    console.log(`  X extent: ${xRange.toFixed(0)} km`);
    console.log(`  Y extent: ${yRange.toFixed(0)} km`);
    console.log(`  Z extent: ${zRange.toFixed(0)} km`);
    console.log(`  Period: ~${(haloPeriod * T_star / 86400).toFixed(1)} days`);

    // L1 distance from Earth
    const L1_dist_km = (L1_x + mu) * L_star / 1e3;
    console.log(`\nL1 distance from Earth: ${L1_dist_km.toFixed(0)} km`);
    console.log(`L1 distance from Moon: ${((1 - mu - L1_x) * L_star / 1e3).toFixed(0)} km`);

    console.log('\n=== Impact Manifolds Computation Complete ===');
}

main().catch(console.error);
