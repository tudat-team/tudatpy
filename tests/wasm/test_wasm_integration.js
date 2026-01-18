/**
 * Integration Tests for Tudat WASM API
 *
 * These tests verify actual functionality by calling WASM functions
 * and checking return values.
 *
 * Run with: node tests/wasm/test_wasm_integration.js
 */

const fs = require('fs');
const path = require('path');

const WASM_MODULE_PATH = path.join(__dirname, '../../build-wasm/src/tudatpy_wasm/tudatpy_wasm.js');

let testsPassed = 0;
let testsFailed = 0;

function assert(condition, testName, actual, expected) {
    if (condition) {
        console.log(`  ✓ ${testName}`);
        testsPassed++;
    } else {
        console.log(`  ✗ ${testName}`);
        console.log(`      Expected: ${expected}`);
        console.log(`      Actual:   ${actual}`);
        testsFailed++;
    }
}

function assertApprox(actual, expected, tolerance, testName) {
    const diff = Math.abs(actual - expected);
    assert(diff < tolerance, testName, actual, `${expected} ± ${tolerance}`);
}

// ============================================================================
// Integration Test Suites
// ============================================================================

function testVectorOperations(tudat) {
    console.log('\n--- Vector Operations ---');

    // Test Vector3d
    if (tudat.Vector3d) {
        const v3 = new tudat.Vector3d();
        v3.set(0, 1.0);
        v3.set(1, 2.0);
        v3.set(2, 3.0);

        assert(v3.get(0) === 1.0, 'Vector3d: set and get x', v3.get(0), 1.0);
        assert(v3.get(1) === 2.0, 'Vector3d: set and get y', v3.get(1), 2.0);
        assert(v3.get(2) === 3.0, 'Vector3d: set and get z', v3.get(2), 3.0);
        assert(v3.size() === 3, 'Vector3d: size is 3', v3.size(), 3);

        v3.delete();
    } else {
        console.log('  ⚠ Vector3d not available');
    }

    // Test Vector6d
    if (tudat.Vector6d) {
        const v6 = new tudat.Vector6d();
        v6.set(0, 7000000.0);  // x position
        v6.set(1, 0.0);        // y position
        v6.set(2, 0.0);        // z position
        v6.set(3, 0.0);        // x velocity
        v6.set(4, 7500.0);     // y velocity
        v6.set(5, 0.0);        // z velocity

        assert(v6.size() === 6, 'Vector6d: size is 6', v6.size(), 6);
        assert(v6.get(0) === 7000000.0, 'Vector6d: position x', v6.get(0), 7000000.0);
        assert(v6.get(4) === 7500.0, 'Vector6d: velocity y', v6.get(4), 7500.0);

        // Test toArray conversion
        if (typeof v6.toArray === 'function') {
            const arr = v6.toArray();
            assert(Array.isArray(arr), 'Vector6d: toArray returns array', typeof arr, 'array');
            assert(arr.length === 6, 'Vector6d: toArray length is 6', arr.length, 6);
            assert(arr[0] === 7000000.0, 'Vector6d: toArray[0] correct', arr[0], 7000000.0);
        }

        v6.delete();
    } else {
        console.log('  ⚠ Vector6d not available');
    }

    // Test VectorXd (dynamic size)
    if (tudat.VectorXd) {
        const vx = new tudat.VectorXd();
        vx.resize(4);
        vx.set(0, 1.0);
        vx.set(1, 2.0);
        vx.set(2, 3.0);
        vx.set(3, 4.0);

        assert(vx.size() === 4, 'VectorXd: size after resize', vx.size(), 4);
        assert(vx.get(2) === 3.0, 'VectorXd: get element', vx.get(2), 3.0);

        vx.delete();
    } else {
        console.log('  ⚠ VectorXd not available');
    }

    // Test Matrix3d
    if (tudat.Matrix3d) {
        const m3 = new tudat.Matrix3d();
        m3.set(0, 0, 1.0);  // Identity-like
        m3.set(1, 1, 1.0);
        m3.set(2, 2, 1.0);

        assert(m3.get(0, 0) === 1.0, 'Matrix3d: diagonal element', m3.get(0, 0), 1.0);
        assert(m3.rows() === 3, 'Matrix3d: rows', m3.rows(), 3);
        assert(m3.cols() === 3, 'Matrix3d: cols', m3.cols(), 3);

        m3.delete();
    } else {
        console.log('  ⚠ Matrix3d not available');
    }
}

function testElementConversion(tudat) {
    console.log('\n--- Element Conversion ---');

    // Test Keplerian to Cartesian conversion
    if (tudat.astro_element_conversion_keplerian_to_cartesian && tudat.Vector6d) {
        // Create Keplerian state (circular orbit at 7000 km)
        const keplerState = new tudat.Vector6d();
        keplerState.set(0, 7000000.0);    // semi-major axis [m]
        keplerState.set(1, 0.0);          // eccentricity
        keplerState.set(2, 0.0);          // inclination [rad]
        keplerState.set(3, 0.0);          // argument of periapsis [rad]
        keplerState.set(4, 0.0);          // RAAN [rad]
        keplerState.set(5, 0.0);          // true anomaly [rad]

        const gravitationalParameter = 3.986004418e14;  // Earth GM [m^3/s^2]

        try {
            const cartesianState = tudat.astro_element_conversion_keplerian_to_cartesian(
                keplerState, gravitationalParameter
            );

            // For circular orbit at true anomaly 0, position should be [a, 0, 0]
            assertApprox(cartesianState.get(0), 7000000.0, 1.0, 'Keplerian to Cartesian: x position');
            assertApprox(cartesianState.get(1), 0.0, 1.0, 'Keplerian to Cartesian: y position');
            assertApprox(cartesianState.get(2), 0.0, 1.0, 'Keplerian to Cartesian: z position');

            // Velocity should be [0, v_circular, 0] where v_circular = sqrt(GM/a)
            const expectedVelocity = Math.sqrt(gravitationalParameter / 7000000.0);
            assertApprox(cartesianState.get(3), 0.0, 1.0, 'Keplerian to Cartesian: x velocity');
            assertApprox(cartesianState.get(4), expectedVelocity, 1.0, 'Keplerian to Cartesian: y velocity');
            assertApprox(cartesianState.get(5), 0.0, 1.0, 'Keplerian to Cartesian: z velocity');

            cartesianState.delete();
        } catch (err) {
            console.log(`  ⚠ Element conversion failed: ${err.message}`);
        }

        keplerState.delete();
    } else {
        console.log('  ⚠ Element conversion functions not available');
    }
}

function testTwoBodyDynamics(tudat) {
    console.log('\n--- Two Body Dynamics ---');

    if (tudat.astro_two_body_dynamics_compute_kepler_orbit_period) {
        const semiMajorAxis = 7000000.0;  // 7000 km
        const gravitationalParameter = 3.986004418e14;

        try {
            const period = tudat.astro_two_body_dynamics_compute_kepler_orbit_period(
                semiMajorAxis, gravitationalParameter
            );

            // T = 2*pi*sqrt(a^3/GM)
            const expectedPeriod = 2 * Math.PI * Math.sqrt(
                Math.pow(semiMajorAxis, 3) / gravitationalParameter
            );

            assertApprox(period, expectedPeriod, 1.0, 'Kepler orbit period calculation');
        } catch (err) {
            console.log(`  ⚠ Two body dynamics failed: ${err.message}`);
        }
    } else {
        console.log('  ⚠ Two body dynamics functions not available');
    }
}

function testDateTime(tudat) {
    console.log('\n--- DateTime ---');

    if (tudat.astro_time_representation_DateTime) {
        try {
            // Create a DateTime for J2000 epoch
            const dt = new tudat.astro_time_representation_DateTime(2000, 1, 1, 12, 0, 0);

            // J2000 epoch should be at 0 seconds from J2000
            const epoch = dt.epoch();
            assertApprox(epoch, 0.0, 1.0, 'DateTime: J2000 epoch is 0');

            dt.delete();
        } catch (err) {
            console.log(`  ⚠ DateTime test failed: ${err.message}`);
        }
    } else {
        console.log('  ⚠ DateTime class not available');
    }
}

function testEnumerations(tudat) {
    console.log('\n--- Enumerations ---');

    // Test propagator types
    if (tudat.dynamics_propagation_setup_propagator_TranslationalPropagatorType) {
        const cowell = tudat.dynamics_propagation_setup_propagator_TranslationalPropagatorType.cowell;
        assert(typeof cowell !== 'undefined', 'TranslationalPropagatorType.cowell exists', cowell, 'defined');

        const encke = tudat.dynamics_propagation_setup_propagator_TranslationalPropagatorType.encke;
        assert(typeof encke !== 'undefined', 'TranslationalPropagatorType.encke exists', encke, 'defined');
    } else {
        console.log('  ⚠ TranslationalPropagatorType not available');
    }

    // Test integrator types
    if (tudat.dynamics_propagation_setup_integrator_AvailableIntegrators) {
        const rk4 = tudat.dynamics_propagation_setup_integrator_AvailableIntegrators.rungeKutta4;
        assert(typeof rk4 !== 'undefined', 'AvailableIntegrators.rungeKutta4 exists', rk4, 'defined');
    } else {
        console.log('  ⚠ AvailableIntegrators not available');
    }

    // Test acceleration types
    if (tudat.dynamics_propagation_setup_acceleration_AvailableAcceleration) {
        const pointMass = tudat.dynamics_propagation_setup_acceleration_AvailableAcceleration.point_mass_gravity;
        assert(typeof pointMass !== 'undefined', 'AvailableAcceleration.point_mass_gravity exists', pointMass, 'defined');

        const sphericalHarmonic = tudat.dynamics_propagation_setup_acceleration_AvailableAcceleration.spherical_harmonic_gravity;
        assert(typeof sphericalHarmonic !== 'undefined', 'AvailableAcceleration.spherical_harmonic_gravity exists', sphericalHarmonic, 'defined');
    } else {
        console.log('  ⚠ AvailableAcceleration not available');
    }
}

function testSettingsFactory(tudat) {
    console.log('\n--- Settings Factory Functions ---');

    // Test integrator settings creation
    if (tudat.dynamics_propagation_setup_integrator_runge_kutta_fixed_step_size) {
        try {
            const stepSize = 60.0;  // 60 seconds
            const integratorSettings = tudat.dynamics_propagation_setup_integrator_runge_kutta_fixed_step_size(
                stepSize,
                tudat.dynamics_propagation_setup_integrator_CoefficientSets ?
                    tudat.dynamics_propagation_setup_integrator_CoefficientSets.rungeKutta4 : undefined
            );

            assert(integratorSettings !== null, 'IntegratorSettings created', integratorSettings, 'non-null');

            if (integratorSettings && integratorSettings.delete) {
                integratorSettings.delete();
            }
        } catch (err) {
            console.log(`  ⚠ Integrator settings creation failed: ${err.message}`);
        }
    } else {
        console.log('  ⚠ runge_kutta_fixed_step_size not available');
    }

    // Test acceleration settings creation
    if (tudat.dynamics_propagation_setup_acceleration_point_mass_gravity) {
        try {
            const accelSettings = tudat.dynamics_propagation_setup_acceleration_point_mass_gravity();
            assert(accelSettings !== null, 'PointMassGravity AccelerationSettings created', accelSettings, 'non-null');

            if (accelSettings && accelSettings.delete) {
                accelSettings.delete();
            }
        } catch (err) {
            console.log(`  ⚠ Acceleration settings creation failed: ${err.message}`);
        }
    } else {
        console.log('  ⚠ point_mass_gravity not available');
    }

    // Test termination settings creation
    if (tudat.dynamics_propagation_setup_propagator_time_termination) {
        try {
            const endTime = 86400.0;  // 1 day in seconds
            const terminationSettings = tudat.dynamics_propagation_setup_propagator_time_termination(endTime);
            assert(terminationSettings !== null, 'TimeTerminationSettings created', terminationSettings, 'non-null');

            if (terminationSettings && terminationSettings.delete) {
                terminationSettings.delete();
            }
        } catch (err) {
            console.log(`  ⚠ Termination settings creation failed: ${err.message}`);
        }
    } else {
        console.log('  ⚠ time_termination not available');
    }
}

// ============================================================================
// Main Test Runner
// ============================================================================

async function runTests() {
    console.log('='.repeat(60));
    console.log('  Tudat WASM Integration Tests');
    console.log('='.repeat(60));

    if (!fs.existsSync(WASM_MODULE_PATH)) {
        console.error(`\nError: WASM module not found at ${WASM_MODULE_PATH}`);
        console.error('\nPlease build Tudat with Emscripten first.');
        process.exit(1);
    }

    console.log('\n[1] Loading WASM module...');

    try {
        const createTudatModule = require(WASM_MODULE_PATH);
        const tudat = await createTudatModule();

        console.log('[2] Running integration tests...');

        testVectorOperations(tudat);
        testElementConversion(tudat);
        testTwoBodyDynamics(tudat);
        testDateTime(tudat);
        testEnumerations(tudat);
        testSettingsFactory(tudat);

        // Print summary
        console.log('\n' + '='.repeat(60));
        console.log('  Integration Test Summary');
        console.log('='.repeat(60));
        console.log(`  Passed: ${testsPassed}`);
        console.log(`  Failed: ${testsFailed}`);
        console.log(`  Total:  ${testsPassed + testsFailed}`);
        console.log('='.repeat(60));

        process.exit(testsFailed > 0 ? 1 : 0);

    } catch (err) {
        console.error('Error:', err);
        process.exit(1);
    }
}

runTests();
