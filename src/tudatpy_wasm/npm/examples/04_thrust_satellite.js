/**
 * Example 04: Thrust Satellite with Mass Propagation
 *
 * Ported from: examples/tudatpy/propagation/thrust_between_Earth_Moon.py
 *
 * This example demonstrates coupled translational and mass propagation
 * for a spacecraft with continuous thrust.
 *
 * Run with: node 04_thrust_satellite.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Thrust Satellite Example ===\n');

    const tudat = await createTudatModule();

    // Constants
    const earthGM = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;   // m
    const g0 = 9.80665;              // m/s^2 (standard gravity)

    // Spacecraft parameters
    const initialMass = 1000.0;      // kg
    const thrustMagnitude = 5.0;     // N (typical ion engine)
    const specificImpulse = 3000.0;  // s (typical ion engine)

    // Mass flow rate: dm/dt = -T / (g0 * Isp)
    const massFlowRate = thrustMagnitude / (g0 * specificImpulse);

    console.log('Spacecraft parameters:');
    console.log(`  Initial mass: ${initialMass} kg`);
    console.log(`  Thrust: ${thrustMagnitude} N`);
    console.log(`  Specific impulse: ${specificImpulse} s`);
    console.log(`  Mass flow rate: ${(massFlowRate * 1000).toFixed(4)} g/s`);

    // Initial orbit - GTO-like
    const periapsisAlt = 200.0e3;   // 200 km
    const apoapsisAlt = 35786.0e3;  // GEO altitude
    const periapsisRadius = earthRadius + periapsisAlt;
    const apoapsisRadius = earthRadius + apoapsisAlt;
    const semiMajorAxis = (periapsisRadius + apoapsisRadius) / 2;
    const eccentricity = (apoapsisRadius - periapsisRadius) / (apoapsisRadius + periapsisRadius);

    console.log('\nInitial orbit (GTO):');
    console.log(`  Periapsis altitude: ${periapsisAlt/1000} km`);
    console.log(`  Apoapsis altitude: ${apoapsisAlt/1000} km`);
    console.log(`  Semi-major axis: ${(semiMajorAxis/1000).toFixed(2)} km`);
    console.log(`  Eccentricity: ${eccentricity.toFixed(4)}`);

    // Start at periapsis
    const keplerianElements = new tudat.Vector6d();
    keplerianElements.set(0, semiMajorAxis);
    keplerianElements.set(1, eccentricity);
    keplerianElements.set(2, 0);  // equatorial
    keplerianElements.set(3, 0);  // arg of periapsis
    keplerianElements.set(4, 0);  // RAAN
    keplerianElements.set(5, 0);  // true anomaly (at periapsis)

    const initialState = tudat.astro.element_conversion.keplerian_to_cartesian(
        keplerianElements,
        earthGM
    );

    // Simulation parameters
    const orbitalPeriod = 2 * Math.PI * Math.sqrt(Math.pow(semiMajorAxis, 3) / earthGM);
    const simulationDuration = orbitalPeriod;  // One orbit
    const thrustDuration = orbitalPeriod / 2;  // Thrust for half orbit (around periapsis)

    console.log(`\nOrbital period: ${(orbitalPeriod/3600).toFixed(2)} hours`);
    console.log(`Thrust duration: ${(thrustDuration/3600).toFixed(2)} hours (half orbit)`);

    // Propagate with thrust
    console.log('\nPropagating with continuous thrust...');

    const numSteps = 500;
    const dt = simulationDuration / numSteps;

    let position = [initialState.get(0), initialState.get(1), initialState.get(2)];
    let velocity = [initialState.get(3), initialState.get(4), initialState.get(5)];
    let mass = initialMass;
    let time = 0;

    const history = [];
    history.push({
        time: 0,
        mass: mass,
        altitude: Math.sqrt(position[0]**2 + position[1]**2 + position[2]**2) - earthRadius,
        velocity: Math.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2)
    });

    // Simple Euler integration for demonstration
    // (In real applications, use proper numerical integrators)
    for (let i = 0; i < numSteps; i++) {
        const r = Math.sqrt(position[0]**2 + position[1]**2 + position[2]**2);
        const v = Math.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2);

        // Gravitational acceleration
        const gravAccel = [-earthGM * position[0] / Math.pow(r, 3),
                          -earthGM * position[1] / Math.pow(r, 3),
                          -earthGM * position[2] / Math.pow(r, 3)];

        // Thrust acceleration (in velocity direction, tangential thrust)
        let thrustAccel = [0, 0, 0];
        if (time < thrustDuration && mass > 100) {  // Thrust only during first half orbit
            const accelMag = thrustMagnitude / mass;
            thrustAccel = [accelMag * velocity[0] / v,
                          accelMag * velocity[1] / v,
                          accelMag * velocity[2] / v];
            mass -= massFlowRate * dt;
        }

        // Total acceleration
        const totalAccel = [gravAccel[0] + thrustAccel[0],
                           gravAccel[1] + thrustAccel[1],
                           gravAccel[2] + thrustAccel[2]];

        // Update state (Euler step)
        position = [position[0] + velocity[0] * dt,
                   position[1] + velocity[1] * dt,
                   position[2] + velocity[2] * dt];

        velocity = [velocity[0] + totalAccel[0] * dt,
                   velocity[1] + totalAccel[1] * dt,
                   velocity[2] + totalAccel[2] * dt];

        time += dt;

        // Store history
        history.push({
            time: time,
            mass: mass,
            altitude: Math.sqrt(position[0]**2 + position[1]**2 + position[2]**2) - earthRadius,
            velocity: Math.sqrt(velocity[0]**2 + velocity[1]**2 + velocity[2]**2)
        });
    }

    // Compute final orbital elements
    const finalState = new tudat.Vector6d();
    finalState.set(0, position[0]);
    finalState.set(1, position[1]);
    finalState.set(2, position[2]);
    finalState.set(3, velocity[0]);
    finalState.set(4, velocity[1]);
    finalState.set(5, velocity[2]);

    const finalKeplerian = tudat.astro.element_conversion.cartesian_to_keplerian(
        finalState,
        earthGM
    );

    console.log('\nResults after one orbit with thrust:');
    console.log(`  Final mass: ${mass.toFixed(2)} kg (consumed ${(initialMass - mass).toFixed(2)} kg)`);
    console.log(`  Initial SMA: ${(semiMajorAxis/1000).toFixed(2)} km`);
    console.log(`  Final SMA: ${(finalKeplerian.get(0)/1000).toFixed(2)} km`);
    console.log(`  SMA increase: ${((finalKeplerian.get(0) - semiMajorAxis)/1000).toFixed(2)} km`);
    console.log(`  Initial eccentricity: ${eccentricity.toFixed(4)}`);
    console.log(`  Final eccentricity: ${finalKeplerian.get(1).toFixed(4)}`);

    // Delta-V calculation
    const deltaV = specificImpulse * g0 * Math.log(initialMass / mass);
    console.log(`\nTotal delta-V provided: ${(deltaV/1000).toFixed(4)} km/s`);

    // Print some trajectory points
    console.log('\nSample trajectory points:');
    const sampleIndices = [0, 100, 200, 300, 400, 500];
    for (const i of sampleIndices) {
        if (i < history.length) {
            const h = history[i];
            console.log(`  t=${(h.time/3600).toFixed(2)} hr: alt=${(h.altitude/1000).toFixed(1)} km, v=${(h.velocity/1000).toFixed(3)} km/s, m=${h.mass.toFixed(1)} kg`);
        }
    }

    // Clean up
    keplerianElements.delete();
    initialState.delete();
    finalState.delete();
    finalKeplerian.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
