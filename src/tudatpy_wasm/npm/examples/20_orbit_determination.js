/**
 * Example 20: Orbit Determination Fundamentals
 *
 * This example demonstrates orbit determination concepts including
 * state estimation, covariance analysis, and observability.
 *
 * Run with: node 20_orbit_determination.js
 */

const createTudatModule = require('@tudat/tudatpy-wasm');

async function main() {
    console.log('=== Orbit Determination Fundamentals ===\n');

    const tudat = await createTudatModule();

    const earthGM = 3.986004418e14;
    const earthRadius = 6378137;

    // ========================================
    // Orbit Determination Overview
    // ========================================
    console.log('1. Orbit Determination Overview');
    console.log('===============================\n');

    console.log('Orbit determination estimates spacecraft state from observations.');
    console.log('Common observation types:');
    console.log('  - Range (distance to spacecraft)');
    console.log('  - Range-rate (Doppler velocity)');
    console.log('  - Angles (azimuth, elevation)');
    console.log('  - GPS positions');
    console.log('  - Optical tracking\n');

    console.log('Methods:');
    console.log('  - Least squares (batch processing)');
    console.log('  - Kalman filter (sequential processing)');
    console.log('  - Consider covariance analysis');
    console.log('  - Unscented Kalman filter (nonlinear)');

    // ========================================
    // State Vector and Covariance
    // ========================================
    console.log('\n\n2. State Vector and Covariance');
    console.log('==============================\n');

    console.log('The state vector typically includes:');
    console.log('  x = [x, y, z, vx, vy, vz]  (position and velocity)');
    console.log('  Extended: + drag coefficient, solar radiation pressure, etc.\n');

    // Example: ISS-like orbit
    const sma = 6778e3;
    const ecc = 0.0001;
    const inc = 51.6 * Math.PI / 180;

    const keplerElements = new tudat.Vector6d();
    keplerElements.set(0, sma);
    keplerElements.set(1, ecc);
    keplerElements.set(2, inc);
    keplerElements.set(3, 0);
    keplerElements.set(4, 0);
    keplerElements.set(5, 0);

    const state = tudat.astro.element_conversion.keplerian_to_cartesian(keplerElements, earthGM);

    console.log('Reference state (ISS-like):');
    console.log(`  Position: [${(state.get(0)/1000).toFixed(3)}, ${(state.get(1)/1000).toFixed(3)}, ${(state.get(2)/1000).toFixed(3)}] km`);
    console.log(`  Velocity: [${(state.get(3)/1000).toFixed(4)}, ${(state.get(4)/1000).toFixed(4)}, ${(state.get(5)/1000).toFixed(4)}] km/s\n`);

    // Example covariance matrix (diagonal for simplicity)
    console.log('Initial covariance (1-sigma uncertainties):');
    console.log('  Position: 100 m (each axis)');
    console.log('  Velocity: 0.1 m/s (each axis)\n');

    const sigmaPos = 100;     // m
    const sigmaVel = 0.1;     // m/s

    console.log('Covariance matrix P (diagonal):');
    console.log(`  P = diag([${sigmaPos**2}, ${sigmaPos**2}, ${sigmaPos**2}, ${sigmaVel**2}, ${sigmaVel**2}, ${sigmaVel**2}])`);
    console.log('  (Position variances in m², velocity variances in m²/s²)');

    // ========================================
    // Observation Modeling
    // ========================================
    console.log('\n\n3. Observation Modeling');
    console.log('=======================\n');

    // Ground station
    const stationLat = 28.5 * Math.PI / 180;   // Cape Canaveral
    const stationLon = -80.5 * Math.PI / 180;
    const stationAlt = 0;

    // Station position in ECEF (simplified)
    const stationX = earthRadius * Math.cos(stationLat) * Math.cos(stationLon);
    const stationY = earthRadius * Math.cos(stationLat) * Math.sin(stationLon);
    const stationZ = earthRadius * Math.sin(stationLat);

    console.log('Ground station (Cape Canaveral):');
    console.log(`  Latitude: ${(stationLat * 180 / Math.PI).toFixed(2)}°N`);
    console.log(`  Longitude: ${Math.abs(stationLon * 180 / Math.PI).toFixed(2)}°W`);
    console.log(`  ECEF: [${(stationX/1000).toFixed(2)}, ${(stationY/1000).toFixed(2)}, ${(stationZ/1000).toFixed(2)}] km\n`);

    // Range observation
    const dx = state.get(0) - stationX;
    const dy = state.get(1) - stationY;
    const dz = state.get(2) - stationZ;
    const range = Math.sqrt(dx*dx + dy*dy + dz*dz);

    // Range-rate observation
    const dvx = state.get(3);  // Simplified: station velocity = 0 in this frame
    const dvy = state.get(4);
    const dvz = state.get(5);
    const rangeRate = (dx*dvx + dy*dvy + dz*dvz) / range;

    console.log('Observation types and values:');
    console.log(`  Range: ${(range / 1000).toFixed(3)} km`);
    console.log(`  Range-rate: ${(rangeRate / 1000).toFixed(4)} km/s`);

    // Azimuth and elevation
    const az = Math.atan2(dy, dx);
    const el = Math.asin(dz / range);
    console.log(`  Azimuth: ${(az * 180 / Math.PI).toFixed(2)}°`);
    console.log(`  Elevation: ${(el * 180 / Math.PI).toFixed(2)}°`);

    // ========================================
    // Measurement Uncertainties
    // ========================================
    console.log('\n\n4. Measurement Uncertainties');
    console.log('============================\n');

    console.log('Typical measurement accuracies:\n');
    console.log('Measurement Type     Ground-based DSN    LEO Tracking');
    console.log('----------------     ----------------    ------------');
    console.log('Range                ~1 m                ~10 m');
    console.log('Range-rate           ~0.1 mm/s           ~1 mm/s');
    console.log('Azimuth              ~0.01°              ~0.1°');
    console.log('Elevation            ~0.01°              ~0.1°');
    console.log('GPS position         N/A                 ~1-10 m\n');

    console.log('Observation weights in least squares:');
    console.log('  W = R⁻¹ where R is measurement covariance');
    console.log('  Higher precision → larger weight');

    // ========================================
    // Partial Derivatives (H Matrix)
    // ========================================
    console.log('\n\n5. Partial Derivatives (H Matrix)');
    console.log('=================================\n');

    console.log('The observation matrix H relates observations to state:');
    console.log('  y = h(x) + noise');
    console.log('  H = ∂h/∂x (Jacobian)\n');

    // Range partials: ∂ρ/∂r = (r - r_station) / |r - r_station|
    const rangePartialX = dx / range;
    const rangePartialY = dy / range;
    const rangePartialZ = dz / range;

    console.log('Range partials (∂ρ/∂r):');
    console.log(`  ∂ρ/∂x = ${rangePartialX.toFixed(6)}`);
    console.log(`  ∂ρ/∂y = ${rangePartialY.toFixed(6)}`);
    console.log(`  ∂ρ/∂z = ${rangePartialZ.toFixed(6)}`);
    console.log('  ∂ρ/∂v = [0, 0, 0]  (range independent of velocity)\n');

    console.log('H matrix for range observation:');
    console.log(`  H_range = [${rangePartialX.toFixed(4)}, ${rangePartialY.toFixed(4)}, ${rangePartialZ.toFixed(4)}, 0, 0, 0]`);

    // ========================================
    // Dilution of Precision
    // ========================================
    console.log('\n\n6. Geometry and Dilution of Precision');
    console.log('======================================\n');

    console.log('Geometric Dilution of Precision (GDOP) measures');
    console.log('how observation geometry affects accuracy.\n');

    console.log('DOP types (GPS terminology):');
    console.log('  GDOP: Geometric (3D position + time)');
    console.log('  PDOP: Position (3D)');
    console.log('  HDOP: Horizontal (2D)');
    console.log('  VDOP: Vertical (1D)');
    console.log('  TDOP: Time\n');

    console.log('Good geometry: Low DOP (~1-2)');
    console.log('Poor geometry: High DOP (>6)');
    console.log('  σ_position = DOP × σ_range\n');

    console.log('Example: 4 GPS satellites');
    console.log('  If satellites clustered: PDOP ~ 10 → poor');
    console.log('  If well-distributed: PDOP ~ 1.5 → good');

    // ========================================
    // Batch Least Squares
    // ========================================
    console.log('\n\n7. Batch Least Squares Estimation');
    console.log('==================================\n');

    console.log('Batch least squares processes all observations together:\n');
    console.log('  1. Initial state guess: x₀');
    console.log('  2. Compute predicted observations: y_pred = h(x₀)');
    console.log('  3. Compute residuals: Δy = y_obs - y_pred');
    console.log('  4. Compute state correction:');
    console.log('     Δx = (HᵀWH)⁻¹ HᵀW Δy');
    console.log('  5. Update state: x₁ = x₀ + Δx');
    console.log('  6. Iterate until convergence\n');

    console.log('Final covariance:');
    console.log('  P = (HᵀWH)⁻¹');

    // Simplified example
    console.log('\nSimplified 1D example:');
    const truePos = 1000;  // m
    const measurements = [1005, 998, 1002, 997, 1001];  // 5 observations
    const measSigma = 5;   // m

    const meanObs = measurements.reduce((a, b) => a + b) / measurements.length;
    const estimatedSigma = measSigma / Math.sqrt(measurements.length);

    console.log(`  True position: ${truePos} m`);
    console.log(`  Observations: [${measurements.join(', ')}] m`);
    console.log(`  Measurement sigma: ${measSigma} m`);
    console.log(`  Estimated position: ${meanObs.toFixed(1)} m`);
    console.log(`  Estimated sigma: ${estimatedSigma.toFixed(2)} m (improved by √${measurements.length})`);

    // ========================================
    // Kalman Filter
    // ========================================
    console.log('\n\n8. Kalman Filter (Sequential Processing)');
    console.log('=========================================\n');

    console.log('The Kalman filter processes observations sequentially:\n');
    console.log('Prediction step (time update):');
    console.log('  x̂⁻ = Φ x̂⁺   (state propagation)');
    console.log('  P⁻ = Φ P⁺ Φᵀ + Q   (covariance propagation)\n');

    console.log('Update step (measurement update):');
    console.log('  K = P⁻ Hᵀ (H P⁻ Hᵀ + R)⁻¹   (Kalman gain)');
    console.log('  x̂⁺ = x̂⁻ + K (y - H x̂⁻)   (state update)');
    console.log('  P⁺ = (I - K H) P⁻   (covariance update)\n');

    // Simple Kalman filter example
    console.log('Example: Tracking position with range measurements\n');

    let x_est = 1050;  // Initial estimate (m)
    let P_est = 100 * 100;  // Initial variance (m²)
    const R = 25 * 25;  // Measurement variance (m²)
    const Q = 1;  // Process noise

    console.log('Step    Measurement    Estimate    Sigma (m)');
    console.log('----    -----------    --------    ---------');
    console.log(`Init    -              ${x_est.toFixed(1)}       ${Math.sqrt(P_est).toFixed(1)}`);

    for (let i = 0; i < measurements.length; i++) {
        // Prediction (no dynamics in this simple example)
        const x_pred = x_est;
        const P_pred = P_est + Q;

        // Update
        const K = P_pred / (P_pred + R);
        x_est = x_pred + K * (measurements[i] - x_pred);
        P_est = (1 - K) * P_pred;

        console.log(`${(i + 1).toString().padStart(4)}    ${measurements[i].toString().padStart(11)}    ${x_est.toFixed(1).padStart(8)}    ${Math.sqrt(P_est).toFixed(2).padStart(9)}`);
    }

    // ========================================
    // Observability
    // ========================================
    console.log('\n\n9. Observability Analysis');
    console.log('=========================\n');

    console.log('A state is observable if it can be determined from measurements.');
    console.log('The observability matrix O = [H; HΦ; HΦ²; ...]');
    console.log('System is observable if rank(O) = n (state dimension)\n');

    console.log('Common observability issues:');
    console.log('  - Range-only: Cannot determine position tangent to range');
    console.log('  - Single pass: Poor along-track determination');
    console.log('  - Coplanar stations: Poor out-of-plane determination\n');

    console.log('Solutions:');
    console.log('  - Multiple observation types (range + range-rate)');
    console.log('  - Multiple passes over time');
    console.log('  - Multiple ground stations');
    console.log('  - A priori information');

    // Clean up
    keplerElements.delete();
    state.delete();

    console.log('\n=== Example completed successfully ===');
}

main().catch(err => {
    console.error('Error:', err);
    process.exit(1);
});
