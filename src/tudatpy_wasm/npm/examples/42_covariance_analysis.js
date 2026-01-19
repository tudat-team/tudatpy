/**
 * Example 42: Covariance Analysis
 *
 * Ported from: examples/tudatpy/estimation/covariance_estimated_parameters.py
 *
 * This example demonstrates covariance analysis for orbit estimation,
 * focusing on how to set up and analyze the covariance matrix of
 * estimated parameters without performing actual estimation.
 *
 * Key concepts:
 * - Covariance matrix computation
 * - Correlation coefficients
 * - Uncertainty ellipsoids
 * - Information content analysis
 * - Parameter observability
 *
 * Run with: node 42_covariance_analysis.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Covariance Analysis for Orbit Estimation ===\n');

    const tudat = await createTudatModule();

    // Constants
    const GM_EARTH = 3.986004418e14;
    const R_EARTH = 6.371e6;

    // DELFI-C3 satellite orbit parameters
    const altitude = 700e3;  // 700 km
    const a = R_EARTH + altitude;
    const e = 0.001;
    const i = 97.4 * Math.PI / 180;  // Sun-synchronous
    const omega = 0;
    const RAAN = 0;
    const nu = 0;

    // Observation parameters
    const observationInterval = 60;      // seconds
    const minElevation = 15 * Math.PI / 180;  // 15 degrees
    const measurementSigmaRangeRate = 0.001;  // 1 mm/s

    // Station coordinates (Delft, Netherlands)
    const stationLat = 52.0116 * Math.PI / 180;
    const stationLon = 4.3571 * Math.PI / 180;
    const stationAlt = 0;

    console.log('Scenario Setup:');
    console.log('  Satellite: DELFI-C3');
    console.log(`  Altitude: ${altitude / 1000} km`);
    console.log(`  Inclination: ${(i * 180 / Math.PI).toFixed(1)} deg (sun-synchronous)`);
    console.log(`  Observation: Range-rate, σ = ${measurementSigmaRangeRate * 1000} mm/s`);
    console.log(`  Min elevation: ${(minElevation * 180 / Math.PI).toFixed(0)} deg\n`);

    /**
     * Compute Cartesian state from Keplerian elements
     */
    function keplerToCartesian(a, e, i, omega, RAAN, nu, GM) {
        // Semi-latus rectum
        const p = a * (1 - e * e);
        const r = p / (1 + e * Math.cos(nu));

        // Position in orbital plane
        const xOrb = r * Math.cos(nu);
        const yOrb = r * Math.sin(nu);

        // Velocity in orbital plane
        const h = Math.sqrt(GM * p);
        const vxOrb = -GM / h * Math.sin(nu);
        const vyOrb = GM / h * (e + Math.cos(nu));

        // Rotation matrices
        const cosO = Math.cos(RAAN), sinO = Math.sin(RAAN);
        const cosw = Math.cos(omega), sinw = Math.sin(omega);
        const cosi = Math.cos(i), sini = Math.sin(i);

        // Transform to inertial frame
        const x = (cosO*cosw - sinO*sinw*cosi) * xOrb + (-cosO*sinw - sinO*cosw*cosi) * yOrb;
        const y = (sinO*cosw + cosO*sinw*cosi) * xOrb + (-sinO*sinw + cosO*cosw*cosi) * yOrb;
        const z = (sinw*sini) * xOrb + (cosw*sini) * yOrb;

        const vx = (cosO*cosw - sinO*sinw*cosi) * vxOrb + (-cosO*sinw - sinO*cosw*cosi) * vyOrb;
        const vy = (sinO*cosw + cosO*sinw*cosi) * vxOrb + (-sinO*sinw + cosO*cosw*cosi) * vyOrb;
        const vz = (sinw*sini) * vxOrb + (cosw*sini) * vyOrb;

        return [x, y, z, vx, vy, vz];
    }

    /**
     * Compute station position in ECEF
     */
    function stationECEF(lat, lon, alt) {
        const r = R_EARTH + alt;
        return [
            r * Math.cos(lat) * Math.cos(lon),
            r * Math.cos(lat) * Math.sin(lon),
            r * Math.sin(lat)
        ];
    }

    /**
     * Check visibility (simplified - Earth rotation not included for brevity)
     */
    function isVisible(satPos, stationPos, minEl) {
        // Vector from station to satellite
        const dx = satPos[0] - stationPos[0];
        const dy = satPos[1] - stationPos[1];
        const dz = satPos[2] - stationPos[2];
        const range = Math.sqrt(dx*dx + dy*dy + dz*dz);

        // Local up vector at station
        const stationR = Math.sqrt(stationPos[0]**2 + stationPos[1]**2 + stationPos[2]**2);
        const upX = stationPos[0] / stationR;
        const upY = stationPos[1] / stationR;
        const upZ = stationPos[2] / stationR;

        // Elevation
        const cosZenith = (dx*upX + dy*upY + dz*upZ) / range;
        const elevation = Math.PI/2 - Math.acos(cosZenith);

        return elevation > minEl;
    }

    /**
     * Compute range-rate observation partial derivatives
     * Simplified: partials w.r.t. position and velocity
     */
    function computePartials(satState, stationPos) {
        const dx = satState[0] - stationPos[0];
        const dy = satState[1] - stationPos[1];
        const dz = satState[2] - stationPos[2];
        const range = Math.sqrt(dx*dx + dy*dy + dz*dz);

        // Unit vector
        const ux = dx / range;
        const uy = dy / range;
        const uz = dz / range;

        // Range-rate partials (∂ρ̇/∂r and ∂ρ̇/∂v)
        // For range-rate, main contribution is from velocity projection
        const vx = satState[3];
        const vy = satState[4];
        const vz = satState[5];
        const rangeRate = ux*vx + uy*vy + uz*vz;

        // Partials w.r.t. position (through unit vector change)
        const drrdr = [
            (vx - rangeRate * ux) / range,
            (vy - rangeRate * uy) / range,
            (vz - rangeRate * uz) / range
        ];

        // Partials w.r.t. velocity (direct)
        const drrdv = [ux, uy, uz];

        return [...drrdr, ...drrdv];
    }

    /**
     * Build information matrix from observations
     */
    function buildInformationMatrix(states, stationPos, sigma) {
        const n = 6;  // State dimension
        const H = [];  // Design matrix rows

        for (const state of states) {
            if (isVisible([state[0], state[1], state[2]], stationPos, minElevation)) {
                H.push(computePartials(state, stationPos));
            }
        }

        if (H.length === 0) {
            return null;
        }

        // Information matrix: J = H^T * W * H, where W = 1/σ²
        const weight = 1 / (sigma * sigma);
        const J = Array(n).fill(null).map(() => Array(n).fill(0));

        for (const row of H) {
            for (let i = 0; i < n; i++) {
                for (let j = 0; j < n; j++) {
                    J[i][j] += weight * row[i] * row[j];
                }
            }
        }

        return { J, numObs: H.length };
    }

    /**
     * Invert matrix (simple 6x6)
     */
    function invertMatrix(A) {
        const n = A.length;
        const augmented = A.map((row, i) => [...row, ...Array(n).fill(0).map((_, j) => i === j ? 1 : 0)]);

        // Forward elimination
        for (let i = 0; i < n; i++) {
            let maxRow = i;
            for (let k = i + 1; k < n; k++) {
                if (Math.abs(augmented[k][i]) > Math.abs(augmented[maxRow][i])) {
                    maxRow = k;
                }
            }
            [augmented[i], augmented[maxRow]] = [augmented[maxRow], augmented[i]];

            for (let k = i + 1; k < n; k++) {
                const c = augmented[k][i] / augmented[i][i];
                for (let j = i; j < 2 * n; j++) {
                    augmented[k][j] -= c * augmented[i][j];
                }
            }
        }

        // Back substitution
        for (let i = n - 1; i >= 0; i--) {
            const c = augmented[i][i];
            for (let j = 0; j < 2 * n; j++) {
                augmented[i][j] /= c;
            }
            for (let k = i - 1; k >= 0; k--) {
                const c = augmented[k][i];
                for (let j = 0; j < 2 * n; j++) {
                    augmented[k][j] -= c * augmented[i][j];
                }
            }
        }

        return augmented.map(row => row.slice(n));
    }

    /**
     * Compute correlation matrix from covariance
     */
    function computeCorrelation(P) {
        const n = P.length;
        const C = Array(n).fill(null).map(() => Array(n).fill(0));

        for (let i = 0; i < n; i++) {
            for (let j = 0; j < n; j++) {
                C[i][j] = P[i][j] / Math.sqrt(P[i][i] * P[j][j]);
            }
        }

        return C;
    }

    // Simulate one orbit period
    console.log('--- Simulating Observations over One Orbit Period ---\n');

    const T = 2 * Math.PI * Math.sqrt(a * a * a / GM_EARTH);
    const dt = observationInterval;
    const numSteps = Math.floor(T / dt);

    const stationPos = stationECEF(stationLat, stationLon, stationAlt);
    const states = [];

    // Propagate orbit (simple 2-body)
    let state = keplerToCartesian(a, e, i, omega, RAAN, nu, GM_EARTH);

    for (let step = 0; step < numSteps; step++) {
        states.push([...state]);

        // Simple Kepler propagation for next step
        const t = (step + 1) * dt;
        const n = Math.sqrt(GM_EARTH / (a * a * a));
        const M = n * t;  // Mean anomaly

        // Solve Kepler's equation (Newton-Raphson)
        let E = M;
        for (let iter = 0; iter < 10; iter++) {
            E = E - (E - e * Math.sin(E) - M) / (1 - e * Math.cos(E));
        }

        // True anomaly
        const nuNew = 2 * Math.atan2(
            Math.sqrt(1 + e) * Math.sin(E / 2),
            Math.sqrt(1 - e) * Math.cos(E / 2)
        );

        state = keplerToCartesian(a, e, i, omega, RAAN, nuNew, GM_EARTH);
    }

    // Build information matrix
    const result = buildInformationMatrix(states, stationPos, measurementSigmaRangeRate);

    if (!result) {
        console.log('No observations available (satellite not visible)');
        return;
    }

    console.log(`Orbital period: ${(T / 60).toFixed(1)} minutes`);
    console.log(`Total observation epochs: ${numSteps}`);
    console.log(`Visible observations: ${result.numObs}`);
    console.log(`Visibility: ${(result.numObs / numSteps * 100).toFixed(1)}%\n`);

    // Compute covariance matrix
    const P = invertMatrix(result.J);

    // Extract uncertainties
    const sigmaX = Math.sqrt(P[0][0]);
    const sigmaY = Math.sqrt(P[1][1]);
    const sigmaZ = Math.sqrt(P[2][2]);
    const sigmaVx = Math.sqrt(P[3][3]);
    const sigmaVy = Math.sqrt(P[4][4]);
    const sigmaVz = Math.sqrt(P[5][5]);

    console.log('--- Covariance Analysis Results ---\n');

    console.log('Position Uncertainties (1σ):');
    console.log(`  σx: ${sigmaX.toFixed(1)} m`);
    console.log(`  σy: ${sigmaY.toFixed(1)} m`);
    console.log(`  σz: ${sigmaZ.toFixed(1)} m`);
    console.log(`  3D RSS: ${Math.sqrt(sigmaX**2 + sigmaY**2 + sigmaZ**2).toFixed(1)} m\n`);

    console.log('Velocity Uncertainties (1σ):');
    console.log(`  σvx: ${(sigmaVx * 1000).toFixed(2)} mm/s`);
    console.log(`  σvy: ${(sigmaVy * 1000).toFixed(2)} mm/s`);
    console.log(`  σvz: ${(sigmaVz * 1000).toFixed(2)} mm/s`);
    console.log(`  3D RSS: ${(Math.sqrt(sigmaVx**2 + sigmaVy**2 + sigmaVz**2) * 1000).toFixed(2)} mm/s\n`);

    // Correlation matrix
    const C = computeCorrelation(P);

    console.log('Correlation Matrix:');
    console.log('       x      y      z      vx     vy     vz');
    const labels = ['x ', 'y ', 'z ', 'vx', 'vy', 'vz'];
    for (let i = 0; i < 6; i++) {
        let row = labels[i] + ' ';
        for (let j = 0; j < 6; j++) {
            row += (C[i][j] >= 0 ? ' ' : '') + C[i][j].toFixed(2) + '  ';
        }
        console.log(row);
    }

    // Identify strong correlations
    console.log('\nStrong Correlations (|ρ| > 0.7):');
    let foundStrong = false;
    for (let i = 0; i < 6; i++) {
        for (let j = i + 1; j < 6; j++) {
            if (Math.abs(C[i][j]) > 0.7) {
                console.log(`  ${labels[i].trim()}-${labels[j].trim()}: ${C[i][j].toFixed(3)}`);
                foundStrong = true;
            }
        }
    }
    if (!foundStrong) {
        console.log('  None found');
    }

    // Uncertainty ellipsoid analysis
    console.log('\n--- Uncertainty Ellipsoid Analysis ---\n');

    // Position covariance eigenvalues (simplified)
    const posP = [[P[0][0], P[0][1], P[0][2]],
                  [P[1][0], P[1][1], P[1][2]],
                  [P[2][0], P[2][1], P[2][2]]];

    // Trace gives sum of eigenvalues
    const tracePos = posP[0][0] + posP[1][1] + posP[2][2];
    const avgSigmaPos = Math.sqrt(tracePos / 3);

    console.log(`Position ellipsoid:`);
    console.log(`  Average semi-axis: ${avgSigmaPos.toFixed(1)} m`);
    console.log(`  Volume (4/3 π σx σy σz): ${(4/3 * Math.PI * sigmaX * sigmaY * sigmaZ).toExponential(2)} m³`);

    // Velocity covariance
    const traceVel = P[3][3] + P[4][4] + P[5][5];
    const avgSigmaVel = Math.sqrt(traceVel / 3);

    console.log(`\nVelocity ellipsoid:`);
    console.log(`  Average semi-axis: ${(avgSigmaVel * 1000).toFixed(2)} mm/s`);

    // Observability analysis
    console.log('\n--- Observability Notes ---\n');
    console.log('Range-rate observations provide:');
    console.log('  - Good velocity determination (direct sensitivity)');
    console.log('  - Position from velocity integration over arc');
    console.log('  - Cross-track position weakly observable from single station');
    console.log('\nTo improve observability:');
    console.log('  - Add range observations');
    console.log('  - Use multiple ground stations');
    console.log('  - Longer observation arcs');
    console.log('  - Lower measurement noise');

    console.log('\n=== Covariance Analysis Complete ===');
}

main().catch(console.error);
