/**
 * Example 05: Solar System Propagation
 *
 * Ported from: examples/tudatpy/propagation/solar_system_propagation.py
 *
 * This example demonstrates multi-body propagation in the solar system,
 * computing planetary ephemerides using numerical integration.
 *
 * Run with: node 05_solar_system_propagation.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Solar System Propagation Example ===\n');

    const tudat = await createTudatModule();

    // Gravitational parameters (m^3/s^2)
    const sunGM = 1.32712440018e20;
    const earthGM = 3.986004418e14;
    const marsGM = 4.282837e13;
    const jupiterGM = 1.26686534e17;

    // Approximate orbital parameters (simplified circular orbits)
    const planets = {
        Earth: {
            sma: 1.496e11,       // 1 AU
            gm: earthGM,
            initialAngle: 0     // rad
        },
        Mars: {
            sma: 2.279e11,       // ~1.52 AU
            gm: marsGM,
            initialAngle: Math.PI / 4  // 45 degrees ahead
        },
        Jupiter: {
            sma: 7.785e11,       // ~5.2 AU
            gm: jupiterGM,
            initialAngle: Math.PI / 2  // 90 degrees ahead
        }
    };

    console.log('Solar System Configuration:');
    console.log('  Sun at origin');
    for (const [name, params] of Object.entries(planets)) {
        const period = 2 * Math.PI * Math.sqrt(Math.pow(params.sma, 3) / sunGM);
        console.log(`  ${name}: SMA = ${(params.sma/1e9).toFixed(3)} x10^9 m, Period = ${(period/86400/365.25).toFixed(2)} years`);
    }

    // Initialize planet states
    const states = {};
    for (const [name, params] of Object.entries(planets)) {
        const velocity = Math.sqrt(sunGM / params.sma);
        states[name] = {
            position: [
                params.sma * Math.cos(params.initialAngle),
                params.sma * Math.sin(params.initialAngle),
                0
            ],
            velocity: [
                -velocity * Math.sin(params.initialAngle),
                velocity * Math.cos(params.initialAngle),
                0
            ]
        };
    }

    // Simulation parameters
    const simulationDuration = 365.25 * 86400;  // 1 Earth year
    const numSteps = 1000;
    const dt = simulationDuration / numSteps;

    console.log(`\nSimulation: 1 Earth year (${numSteps} steps)`);
    console.log('Propagating planetary orbits...');

    // Store history
    const history = {
        time: [],
        Earth: { x: [], y: [] },
        Mars: { x: [], y: [] },
        Jupiter: { x: [], y: [] }
    };

    // Propagate using simple N-body integration
    for (let step = 0; step <= numSteps; step++) {
        const time = step * dt;
        history.time.push(time);

        // Store current positions
        for (const name of Object.keys(planets)) {
            history[name].x.push(states[name].position[0]);
            history[name].y.push(states[name].position[1]);
        }

        if (step < numSteps) {
            // Compute accelerations for each planet
            const accelerations = {};

            for (const name of Object.keys(planets)) {
                const pos = states[name].position;
                const r = Math.sqrt(pos[0]**2 + pos[1]**2 + pos[2]**2);

                // Sun's gravity
                accelerations[name] = [
                    -sunGM * pos[0] / Math.pow(r, 3),
                    -sunGM * pos[1] / Math.pow(r, 3),
                    -sunGM * pos[2] / Math.pow(r, 3)
                ];

                // Other planets' gravity (perturbations)
                for (const [otherName, otherParams] of Object.entries(planets)) {
                    if (otherName !== name) {
                        const otherPos = states[otherName].position;
                        const dx = otherPos[0] - pos[0];
                        const dy = otherPos[1] - pos[1];
                        const dz = otherPos[2] - pos[2];
                        const dist = Math.sqrt(dx**2 + dy**2 + dz**2);

                        accelerations[name][0] += otherParams.gm * dx / Math.pow(dist, 3);
                        accelerations[name][1] += otherParams.gm * dy / Math.pow(dist, 3);
                        accelerations[name][2] += otherParams.gm * dz / Math.pow(dist, 3);
                    }
                }
            }

            // Update states (velocity Verlet or simple Euler)
            for (const name of Object.keys(planets)) {
                // Update position
                states[name].position[0] += states[name].velocity[0] * dt + 0.5 * accelerations[name][0] * dt * dt;
                states[name].position[1] += states[name].velocity[1] * dt + 0.5 * accelerations[name][1] * dt * dt;
                states[name].position[2] += states[name].velocity[2] * dt + 0.5 * accelerations[name][2] * dt * dt;

                // Update velocity
                states[name].velocity[0] += accelerations[name][0] * dt;
                states[name].velocity[1] += accelerations[name][1] * dt;
                states[name].velocity[2] += accelerations[name][2] * dt;
            }
        }
    }

    // Analyze results
    console.log('\nResults after 1 Earth year:');

    for (const [name, params] of Object.entries(planets)) {
        const initialPos = { x: history[name].x[0], y: history[name].y[0] };
        const finalPos = { x: history[name].x[numSteps], y: history[name].y[numSteps] };

        const initialR = Math.sqrt(initialPos.x**2 + initialPos.y**2);
        const finalR = Math.sqrt(finalPos.x**2 + finalPos.y**2);

        const initialAngle = Math.atan2(initialPos.y, initialPos.x);
        const finalAngle = Math.atan2(finalPos.y, finalPos.x);

        // Angular motion
        let angularMotion = finalAngle - initialAngle;
        if (angularMotion < 0) angularMotion += 2 * Math.PI;

        const orbitalPeriod = 2 * Math.PI * Math.sqrt(Math.pow(params.sma, 3) / sunGM);
        const expectedAngularMotion = 2 * Math.PI * simulationDuration / orbitalPeriod;

        console.log(`\n  ${name}:`);
        console.log(`    Initial distance from Sun: ${(initialR/1e9).toFixed(3)} x10^9 m`);
        console.log(`    Final distance from Sun: ${(finalR/1e9).toFixed(3)} x10^9 m`);
        console.log(`    Distance change: ${((finalR - initialR)/1e6).toFixed(2)} x10^6 m`);
        console.log(`    Angular motion: ${(angularMotion * 180 / Math.PI).toFixed(1)} degrees`);
        console.log(`    Expected angular motion: ${(expectedAngularMotion * 180 / Math.PI % 360).toFixed(1)} degrees`);
    }

    // Compute Earth-Mars distance over time
    console.log('\nEarth-Mars distance analysis:');
    let minDist = Infinity;
    let maxDist = 0;
    let minTime = 0;
    let maxTime = 0;

    for (let i = 0; i < history.time.length; i++) {
        const dx = history.Mars.x[i] - history.Earth.x[i];
        const dy = history.Mars.y[i] - history.Earth.y[i];
        const dist = Math.sqrt(dx**2 + dy**2);

        if (dist < minDist) {
            minDist = dist;
            minTime = history.time[i];
        }
        if (dist > maxDist) {
            maxDist = dist;
            maxTime = history.time[i];
        }
    }

    console.log(`  Minimum distance: ${(minDist/1e9).toFixed(3)} x10^9 m at day ${(minTime/86400).toFixed(0)}`);
    console.log(`  Maximum distance: ${(maxDist/1e9).toFixed(3)} x10^9 m at day ${(maxTime/86400).toFixed(0)}`);

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
