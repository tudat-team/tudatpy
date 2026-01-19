/**
 * Example 22: Covariance Propagation
 *
 * Ported from: examples/tudatpy/estimation/covariance_propagation_example.py
 *
 * This example demonstrates how to propagate state covariance (uncertainty)
 * along an orbit using the state transition matrix.
 *
 * Run with: node 22_covariance_propagation.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Covariance Propagation Example ===\n');

    const tudat = await createTudatModule();

    // Physical constants
    const earthGM = 3.986004418e14;  // m^3/s^2
    const earthRadius = 6378137.0;  // m

    // Define initial orbit (LEO circular)
    const altitude = 400.0e3;  // 400 km
    const semiMajorAxis = earthRadius + altitude;
    const orbitalPeriod = 2 * Math.PI * Math.sqrt(Math.pow(semiMajorAxis, 3) / earthGM);

    console.log('Orbit Parameters:');
    console.log(`  Altitude: ${(altitude / 1000).toFixed(1)} km`);
    console.log(`  Semi-major axis: ${(semiMajorAxis / 1000).toFixed(1)} km`);
    console.log(`  Orbital period: ${(orbitalPeriod / 60).toFixed(2)} minutes`);

    // Initial state covariance (6x6 matrix)
    // Diagonal: position uncertainty 100m, velocity uncertainty 0.1 m/s
    const positionSigma = 100.0;  // m
    const velocitySigma = 0.1;  // m/s

    // Create initial covariance matrix
    const initialCovariance = new tudat.MatrixXd(6, 6);
    for (let i = 0; i < 6; i++) {
        for (let j = 0; j < 6; j++) {
            initialCovariance.set(i, j, 0.0);
        }
    }
    // Set diagonal elements
    initialCovariance.set(0, 0, positionSigma * positionSigma);
    initialCovariance.set(1, 1, positionSigma * positionSigma);
    initialCovariance.set(2, 2, positionSigma * positionSigma);
    initialCovariance.set(3, 3, velocitySigma * velocitySigma);
    initialCovariance.set(4, 4, velocitySigma * velocitySigma);
    initialCovariance.set(5, 5, velocitySigma * velocitySigma);

    console.log('\nInitial Covariance (1-sigma):');
    console.log(`  Position: ${positionSigma.toFixed(1)} m`);
    console.log(`  Velocity: ${velocitySigma.toFixed(3)} m/s`);

    // For a simple two-body orbit, the state transition matrix can be
    // approximated using Clohessy-Wiltshire (CW) equations for circular orbits
    const n = Math.sqrt(earthGM / Math.pow(semiMajorAxis, 3));  // mean motion

    console.log(`  Mean motion: ${(n * 1000).toFixed(4)} mrad/s`);

    // Propagate covariance over one orbit
    const numSteps = 36;  // Every 10 degrees
    const dt = orbitalPeriod / numSteps;

    console.log('\nPropagating covariance...\n');
    console.log('Time (min) | Pos Sigma (m) | Vel Sigma (m/s)');
    console.log('-----------|---------------|----------------');

    // Propagate using linearized dynamics (CW approximation)
    // For circular orbit: x'' - 2n*y' - 3n²x = 0, y'' + 2n*x' = 0, z'' + n²z = 0
    for (let step = 0; step <= numSteps; step++) {
        const t = step * dt;
        const nt = n * t;

        // State transition matrix for CW equations
        // This is an approximation valid for circular orbits
        const phi = new Array(6).fill(null).map(() => new Array(6).fill(0));

        // In-plane (x-y) coupling
        phi[0][0] = 4 - 3 * Math.cos(nt);
        phi[0][1] = 0;
        phi[0][3] = Math.sin(nt) / n;
        phi[0][4] = (2 / n) * (1 - Math.cos(nt));

        phi[1][0] = 6 * (Math.sin(nt) - nt);
        phi[1][1] = 1;
        phi[1][3] = (2 / n) * (Math.cos(nt) - 1);
        phi[1][4] = (4 * Math.sin(nt) - 3 * nt) / n;

        phi[3][0] = 3 * n * Math.sin(nt);
        phi[3][3] = Math.cos(nt);
        phi[3][4] = 2 * Math.sin(nt);

        phi[4][0] = 6 * n * (Math.cos(nt) - 1);
        phi[4][3] = -2 * Math.sin(nt);
        phi[4][4] = 4 * Math.cos(nt) - 3;

        // Out-of-plane (z) - simple harmonic
        phi[2][2] = Math.cos(nt);
        phi[2][5] = Math.sin(nt) / n;
        phi[5][2] = -n * Math.sin(nt);
        phi[5][5] = Math.cos(nt);

        // Propagate covariance: P(t) = Phi * P0 * Phi^T
        // Compute Phi * P0
        const PhiP0 = new Array(6).fill(null).map(() => new Array(6).fill(0));
        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 6; j++) {
                for (let k = 0; k < 6; k++) {
                    PhiP0[i][j] += phi[i][k] * initialCovariance.get(k, j);
                }
            }
        }

        // Compute PhiP0 * Phi^T
        const Pt = new Array(6).fill(null).map(() => new Array(6).fill(0));
        for (let i = 0; i < 6; i++) {
            for (let j = 0; j < 6; j++) {
                for (let k = 0; k < 6; k++) {
                    Pt[i][j] += PhiP0[i][k] * phi[j][k];  // Phi^T
                }
            }
        }

        // Extract formal errors (sqrt of diagonal)
        const posSigma = Math.sqrt(Pt[0][0] + Pt[1][1] + Pt[2][2]) / Math.sqrt(3);
        const velSigma = Math.sqrt(Pt[3][3] + Pt[4][4] + Pt[5][5]) / Math.sqrt(3);

        if (step % 6 === 0) {  // Print every 60 degrees
            console.log(`${(t / 60).toFixed(1).padStart(10)} | ${posSigma.toFixed(1).padStart(13)} | ${velSigma.toFixed(4).padStart(14)}`);
        }
    }

    // Clean up
    initialCovariance.delete();

    console.log('\n=== Covariance propagation complete ===');
    console.log('\nNote: Covariance growth depends on:');
    console.log('  - Initial uncertainty');
    console.log('  - Orbital dynamics (secular growth in along-track direction)');
    console.log('  - Unmodeled perturbations (not included in this simple example)');
}

main().catch(console.error);
