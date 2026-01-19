/**
 * Example 12: Numerical Integrators
 *
 * This example demonstrates different numerical integration methods
 * available in Tudat WASM for solving ODEs.
 *
 * Run with: node 12_numerical_integrators.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Numerical Integrators ===\n');

    const tudat = await createTudatModule();

    // ========================================
    // Overview of Integrators
    // ========================================
    console.log('1. Available Integrators in Tudat');
    console.log('=================================\n');

    const integrators = [
        {
            name: 'Euler',
            order: 1,
            type: 'Fixed-step',
            stages: 1,
            use: 'Educational, simple problems'
        },
        {
            name: 'Runge-Kutta 4 (RK4)',
            order: 4,
            type: 'Fixed-step',
            stages: 4,
            use: 'General purpose, smooth problems'
        },
        {
            name: 'Runge-Kutta-Fehlberg 4(5)',
            order: '4(5)',
            type: 'Variable-step',
            stages: 6,
            use: 'Adaptive step size, moderate accuracy'
        },
        {
            name: 'Runge-Kutta-Fehlberg 7(8)',
            order: '7(8)',
            type: 'Variable-step',
            stages: 13,
            use: 'High accuracy orbit propagation'
        },
        {
            name: 'Dormand-Prince 8(7)',
            order: '8(7)',
            type: 'Variable-step',
            stages: 13,
            use: 'High accuracy, dense output'
        },
        {
            name: 'Adams-Bashforth-Moulton',
            order: 'Variable',
            type: 'Multi-step',
            stages: 'N/A',
            use: 'Long integrations, smooth problems'
        },
        {
            name: 'Bulirsch-Stoer',
            order: 'Adaptive',
            type: 'Extrapolation',
            stages: 'Variable',
            use: 'Very high accuracy requirements'
        }
    ];

    console.log('Name                        Order    Type           Use Case');
    console.log('----                        -----    ----           --------');

    for (const integ of integrators) {
        console.log(`${integ.name.padEnd(27)} ${integ.order.toString().padStart(5)}    ${integ.type.padEnd(13)}  ${integ.use}`);
    }

    // ========================================
    // Harmonic Oscillator Test
    // ========================================
    console.log('\n\n2. Test Problem: Harmonic Oscillator');
    console.log('=====================================\n');

    console.log('ODE: d²x/dt² = -ω²x');
    console.log('Exact solution: x(t) = A·cos(ωt + φ)\n');

    const omega = 1.0;  // angular frequency
    const A = 1.0;      // amplitude
    const phi = 0;      // phase

    // Initial conditions: x(0) = A, v(0) = 0
    let x = A;
    let v = 0;

    // Integration parameters
    const T = 2 * Math.PI / omega;  // One period
    const dt_values = [0.1, 0.01, 0.001];

    console.log('RK4 integration over one period:');
    console.log('--------------------------------');

    for (const dt of dt_values) {
        const nSteps = Math.round(T / dt);

        // Reset initial conditions
        x = A;
        v = 0;

        // RK4 integration
        for (let i = 0; i < nSteps; i++) {
            // RK4 stages for harmonic oscillator
            const k1x = v;
            const k1v = -omega * omega * x;

            const k2x = v + 0.5 * dt * k1v;
            const k2v = -omega * omega * (x + 0.5 * dt * k1x);

            const k3x = v + 0.5 * dt * k2v;
            const k3v = -omega * omega * (x + 0.5 * dt * k2x);

            const k4x = v + dt * k3v;
            const k4v = -omega * omega * (x + dt * k3x);

            x += dt * (k1x + 2*k2x + 2*k3x + k4x) / 6;
            v += dt * (k1v + 2*k2v + 2*k3v + k4v) / 6;
        }

        // Exact solution at t = T
        const x_exact = A * Math.cos(omega * T + phi);
        const error = Math.abs(x - x_exact);

        console.log(`  dt = ${dt.toFixed(3)}: x(T) = ${x.toFixed(8)}, error = ${error.toExponential(2)}`);
    }

    // ========================================
    // Two-Body Problem
    // ========================================
    console.log('\n\n3. Test Problem: Two-Body Orbit');
    console.log('================================\n');

    const GM = 3.986004418e14;
    const earthRadius = 6378137;

    // Circular orbit at 400 km altitude
    const r0 = earthRadius + 400e3;
    const v0 = Math.sqrt(GM / r0);

    console.log('Circular orbit propagation (RK4):');
    console.log('---------------------------------');
    console.log(`  Initial radius: ${(r0/1000).toFixed(1)} km`);
    console.log(`  Initial velocity: ${(v0/1000).toFixed(3)} km/s\n`);

    // Orbital period
    const T_orbit = 2 * Math.PI * Math.sqrt(r0**3 / GM);

    const stepsPerOrbit = [100, 1000, 10000];

    for (const nSteps of stepsPerOrbit) {
        const dt = T_orbit / nSteps;

        // Initial state [x, y, vx, vy]
        let state = [r0, 0, 0, v0];

        // Propagate one orbit with RK4
        for (let i = 0; i < nSteps; i++) {
            const [x, y, vx, vy] = state;
            const r = Math.sqrt(x*x + y*y);
            const r3 = r * r * r;

            // Accelerations
            const ax = -GM * x / r3;
            const ay = -GM * y / r3;

            // RK4 for 2D two-body problem
            const k1 = [vx, vy, ax, ay];

            const x2 = x + 0.5*dt*k1[0];
            const y2 = y + 0.5*dt*k1[1];
            const r2 = Math.sqrt(x2*x2 + y2*y2);
            const k2 = [vx + 0.5*dt*k1[2], vy + 0.5*dt*k1[3],
                       -GM*x2/(r2**3), -GM*y2/(r2**3)];

            const x3 = x + 0.5*dt*k2[0];
            const y3 = y + 0.5*dt*k2[1];
            const r3_ = Math.sqrt(x3*x3 + y3*y3);
            const k3 = [vx + 0.5*dt*k2[2], vy + 0.5*dt*k2[3],
                       -GM*x3/(r3_**3), -GM*y3/(r3_**3)];

            const x4 = x + dt*k3[0];
            const y4 = y + dt*k3[1];
            const r4 = Math.sqrt(x4*x4 + y4*y4);
            const k4 = [vx + dt*k3[2], vy + dt*k3[3],
                       -GM*x4/(r4**3), -GM*y4/(r4**3)];

            state = [
                x + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6,
                y + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6,
                vx + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6,
                vy + dt*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6
            ];
        }

        // Check orbit closure
        const r_final = Math.sqrt(state[0]**2 + state[1]**2);
        const error_r = Math.abs(r_final - r0);
        const error_pos = Math.sqrt((state[0] - r0)**2 + state[1]**2);

        console.log(`  ${nSteps} steps: Position error = ${error_pos.toExponential(2)} m, Radius error = ${error_r.toExponential(2)} m`);
    }

    // ========================================
    // Energy Conservation
    // ========================================
    console.log('\n\n4. Energy Conservation Test');
    console.log('===========================\n');

    console.log('Specific orbital energy: E = v²/2 - μ/r\n');

    const nSteps = 1000;
    const dt = T_orbit / nSteps;

    let state = [r0, 0, 0, v0];
    const E0 = v0*v0/2 - GM/r0;  // Initial energy

    let maxEnergyError = 0;

    for (let i = 0; i < nSteps; i++) {
        const [x, y, vx, vy] = state;
        const r = Math.sqrt(x*x + y*y);
        const v = Math.sqrt(vx*vx + vy*vy);
        const E = v*v/2 - GM/r;
        const relError = Math.abs((E - E0) / E0);
        maxEnergyError = Math.max(maxEnergyError, relError);

        // RK4 step (same as above)
        const r3 = r ** 3;
        const ax = -GM * x / r3;
        const ay = -GM * y / r3;

        const k1 = [vx, vy, ax, ay];
        const x2 = x + 0.5*dt*k1[0];
        const y2 = y + 0.5*dt*k1[1];
        const r2 = Math.sqrt(x2*x2 + y2*y2);
        const k2 = [vx + 0.5*dt*k1[2], vy + 0.5*dt*k1[3], -GM*x2/(r2**3), -GM*y2/(r2**3)];
        const x3 = x + 0.5*dt*k2[0];
        const y3 = y + 0.5*dt*k2[1];
        const r3_ = Math.sqrt(x3*x3 + y3*y3);
        const k3 = [vx + 0.5*dt*k2[2], vy + 0.5*dt*k2[3], -GM*x3/(r3_**3), -GM*y3/(r3_**3)];
        const x4 = x + dt*k3[0];
        const y4 = y + dt*k3[1];
        const r4 = Math.sqrt(x4*x4 + y4*y4);
        const k4 = [vx + dt*k3[2], vy + dt*k3[3], -GM*x4/(r4**3), -GM*y4/(r4**3)];

        state = [
            x + dt*(k1[0] + 2*k2[0] + 2*k3[0] + k4[0])/6,
            y + dt*(k1[1] + 2*k2[1] + 2*k3[1] + k4[1])/6,
            vx + dt*(k1[2] + 2*k2[2] + 2*k3[2] + k4[2])/6,
            vy + dt*(k1[3] + 2*k2[3] + 2*k3[3] + k4[3])/6
        ];
    }

    console.log(`Initial energy: ${E0.toExponential(6)} J/kg`);
    console.log(`Max relative energy error over 1 orbit: ${maxEnergyError.toExponential(2)}`);

    // ========================================
    // Integrator Selection Guide
    // ========================================
    console.log('\n\n5. Integrator Selection Guide');
    console.log('==============================\n');

    console.log('Recommended integrators by application:');
    console.log('---------------------------------------\n');

    const recommendations = [
        { app: 'LEO propagation (days)', integrator: 'RKF7(8) or DOPRI8(7)', step: '10-60 s' },
        { app: 'LEO with drag (days)', integrator: 'RKF4(5) variable step', step: 'adaptive' },
        { app: 'GEO/MEO (months)', integrator: 'DOPRI8(7)', step: '60-600 s' },
        { app: 'Interplanetary (years)', integrator: 'Adams-Bashforth-Moulton', step: 'adaptive' },
        { app: 'Close approaches', integrator: 'Bulirsch-Stoer', step: 'adaptive' },
        { app: 'Real-time simulation', integrator: 'RK4 fixed step', step: '1-10 s' },
        { app: 'Monte Carlo studies', integrator: 'RK4 or RKF4(5)', step: 'fixed/adaptive' }
    ];

    console.log('Application               Integrator                  Typical Step');
    console.log('-----------               ----------                  ------------');

    for (const rec of recommendations) {
        console.log(`${rec.app.padEnd(25)} ${rec.integrator.padEnd(27)} ${rec.step}`);
    }

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
