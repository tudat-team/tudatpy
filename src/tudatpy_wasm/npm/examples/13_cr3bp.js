/**
 * Example 13: Circular Restricted Three-Body Problem (CR3BP)
 *
 * This example demonstrates the CR3BP, including Lagrange points,
 * Jacobi constant, and halo orbits.
 *
 * Run with: node 13_cr3bp.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Circular Restricted Three-Body Problem ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // CR3BP Fundamentals
    // ========================================
    console.log('1. CR3BP Fundamentals');
    console.log('=====================\n');

    console.log('The CR3BP describes motion of a massless body under the');
    console.log('gravitational influence of two massive bodies in circular orbits.\n');

    console.log('Key parameters:');
    console.log('  - Mass ratio: μ = m₂/(m₁ + m₂)');
    console.log('  - Characteristic length: L = distance between primaries');
    console.log('  - Characteristic time: T = orbital period / 2π\n');

    // System parameters
    const systems = [
        { name: 'Earth-Moon', m1: 5.972e24, m2: 7.342e22, L: 384400e3 },
        { name: 'Sun-Earth', m1: 1.989e30, m2: 5.972e24, L: 1.496e11 },
        { name: 'Sun-Jupiter', m1: 1.989e30, m2: 1.898e27, L: 7.785e11 },
        { name: 'Mars-Phobos', m1: 6.417e23, m2: 1.06e16, L: 9.376e6 }
    ];

    console.log('System          Mass ratio μ         Characteristic length');
    console.log('------          ------------         ---------------------');

    for (const sys of systems) {
        const mu = sys.m2 / (sys.m1 + sys.m2);
        console.log(`${sys.name.padEnd(15)} ${mu.toExponential(4).padStart(12)}         ${(sys.L/1000).toExponential(4)} km`);
    }

    // ========================================
    // Lagrange Points
    // ========================================
    console.log('\n\n2. Lagrange Points');
    console.log('==================\n');

    console.log('Five equilibrium points exist in the rotating frame:\n');
    console.log('  L1: Between the primaries (unstable)');
    console.log('  L2: Beyond the smaller primary (unstable)');
    console.log('  L3: Beyond the larger primary (unstable)');
    console.log('  L4: Leading triangular point (stable for μ < 0.0385)');
    console.log('  L5: Trailing triangular point (stable for μ < 0.0385)\n');

    // Calculate L1, L2, L3 positions for Earth-Moon system
    const earthMoon = systems[0];
    const mu = earthMoon.m2 / (earthMoon.m1 + earthMoon.m2);

    console.log(`Earth-Moon system (μ = ${mu.toFixed(6)}):`);
    console.log('----------------------------------------');

    // L1 position (approximate)
    // r_L1 ≈ r₂ × (μ/3)^(1/3)
    const r_L1_approx = earthMoon.L * Math.pow(mu / 3, 1/3);
    const L1_from_Earth = earthMoon.L - r_L1_approx;

    // L2 position (approximate)
    const r_L2_approx = earthMoon.L * Math.pow(mu / 3, 1/3);
    const L2_from_Earth = earthMoon.L + r_L2_approx;

    // L3 position (approximate)
    const L3_from_Earth = -earthMoon.L * (1 + 7*mu/12);

    console.log(`  L1: ${(L1_from_Earth/1000).toFixed(0)} km from Earth (${(r_L1_approx/1000).toFixed(0)} km from Moon)`);
    console.log(`  L2: ${(L2_from_Earth/1000).toFixed(0)} km from Earth (${((L2_from_Earth - earthMoon.L)/1000).toFixed(0)} km beyond Moon)`);
    console.log(`  L3: ${(Math.abs(L3_from_Earth)/1000).toFixed(0)} km from Earth (opposite side of Moon)`);
    console.log(`  L4: Triangular point (60° ahead of Moon)`);
    console.log(`  L5: Triangular point (60° behind Moon)`);

    // ========================================
    // Jacobi Constant
    // ========================================
    console.log('\n\n3. Jacobi Constant (Energy Integral)');
    console.log('=====================================\n');

    console.log('The Jacobi constant is conserved in the rotating frame:');
    console.log('  C = -2U - (ẋ² + ẏ² + ż²)');
    console.log('where U is the effective potential.\n');

    // Calculate Jacobi constant at Lagrange points
    function effectivePotential(x, y, mu) {
        const r1 = Math.sqrt((x + mu)**2 + y**2);
        const r2 = Math.sqrt((x - 1 + mu)**2 + y**2);
        return -0.5 * (x**2 + y**2) - (1 - mu)/r1 - mu/r2;
    }

    console.log('Jacobi constant at Lagrange points (normalized units):');
    console.log('------------------------------------------------------');

    // L1 position (refined)
    const x_L1 = 1 - mu - Math.pow(mu/3, 1/3);
    const C_L1 = -2 * effectivePotential(x_L1, 0, mu);

    // L2 position
    const x_L2 = 1 - mu + Math.pow(mu/3, 1/3);
    const C_L2 = -2 * effectivePotential(x_L2, 0, mu);

    // L3 position
    const x_L3 = -1 - 5*mu/12;
    const C_L3 = -2 * effectivePotential(x_L3, 0, mu);

    // L4/L5 positions
    const x_L4 = 0.5 - mu;
    const y_L4 = Math.sqrt(3)/2;
    const C_L4 = -2 * effectivePotential(x_L4, y_L4, mu);

    console.log(`  C(L1) = ${C_L1.toFixed(6)}`);
    console.log(`  C(L2) = ${C_L2.toFixed(6)}`);
    console.log(`  C(L3) = ${C_L3.toFixed(6)}`);
    console.log(`  C(L4) = C(L5) = ${C_L4.toFixed(6)}`);

    console.log('\n  Note: L4 and L5 have the highest Jacobi constant,');
    console.log('        allowing the widest range of motion.\n');

    // ========================================
    // Zero-Velocity Curves
    // ========================================
    console.log('\n4. Zero-Velocity Curves');
    console.log('=======================\n');

    console.log('At zero velocity, the Jacobi constant equals -2U.');
    console.log('Zero-velocity curves define forbidden regions.\n');

    console.log('Energy ranges and accessible regions:');
    console.log('-------------------------------------');
    console.log(`  C > ${C_L1.toFixed(4)}: Motion confined near one primary`);
    console.log(`  C = ${C_L1.toFixed(4)}: Gateway opens at L1`);
    console.log(`  C = ${C_L2.toFixed(4)}: Gateway opens at L2 (escape possible)`);
    console.log(`  C = ${C_L3.toFixed(4)}: Gateway opens at L3`);
    console.log(`  C < ${C_L4.toFixed(4)}: Entire plane accessible`);

    // ========================================
    // Halo Orbits
    // ========================================
    console.log('\n\n5. Periodic Orbits Near Lagrange Points');
    console.log('=======================================\n');

    console.log('Families of periodic orbits exist near L1 and L2:\n');
    console.log('  - Lyapunov orbits: Planar, unstable');
    console.log('  - Halo orbits: 3D, out-of-plane oscillations');
    console.log('  - Lissajous orbits: Quasi-periodic, bounded motion\n');

    // Example halo orbit characteristics
    console.log('Typical halo orbit at Earth-Moon L2:');
    console.log('------------------------------------');
    console.log('  Period: ~14 days');
    console.log('  x-amplitude: ~10,000 km');
    console.log('  z-amplitude (out of plane): ~30,000 km');
    console.log('  Applications: Gateway station, deep space communication\n');

    console.log('Typical halo orbit at Sun-Earth L2:');
    console.log('-----------------------------------');
    console.log('  Period: ~180 days');
    console.log('  Amplitude: ~800,000 km');
    console.log('  Applications: JWST, WMAP, Gaia, Euclid');

    // ========================================
    // Stability and Station-Keeping
    // ========================================
    console.log('\n\n6. Stability and Station-Keeping');
    console.log('=================================\n');

    console.log('L1, L2, L3 are unstable (saddle points)');
    console.log('  → Instability timescale ~23 days for Earth-Moon L2');
    console.log('  → Requires ~1-10 m/s/year for station-keeping\n');

    console.log('L4, L5 are stable if μ < 0.0385 (Routh criterion)');
    console.log(`  → Earth-Moon: μ = ${mu.toFixed(4)} < 0.0385 ✓`);
    console.log(`  → Sun-Jupiter: μ = ${(systems[2].m2/(systems[2].m1+systems[2].m2)).toExponential(3)} < 0.0385 ✓`);
    console.log('  → Trojan asteroids naturally stable at L4/L5\n');

    // ========================================
    // Low-Energy Transfers
    // ========================================
    console.log('\n7. Low-Energy Transfers');
    console.log('=======================\n');

    console.log('The CR3BP enables fuel-efficient transfers via:');
    console.log('  - Weak stability boundary (ballistic capture)');
    console.log('  - Invariant manifolds (tubes connecting L-points)');
    console.log('  - Heteroclinic connections between periodic orbits\n');

    console.log('Example: Earth-Moon low-energy transfer');
    console.log('---------------------------------------');
    console.log('  Direct Hohmann transfer: ~3.1 km/s ΔV');
    console.log('  Low-energy transfer via L1: ~2.5 km/s ΔV');
    console.log('  Trade-off: Lower ΔV but ~3-4 month transfer time\n');

    console.log('Interplanetary Superhighway:');
    console.log('  Connecting manifolds of different systems');
    console.log('  (Sun-Earth L2 → Earth-Moon L1 → Moon)');

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
