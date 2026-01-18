/**
 * Comprehensive Test Suite for Tudat WASM API
 *
 * This test suite verifies that all module bindings are working correctly.
 * Run with: node tests/wasm/test_wasm_api.js
 *
 * Prerequisites:
 *   Build Tudat with Emscripten first:
 *   mkdir build-wasm && cd build-wasm
 *   emcmake cmake ..
 *   make tudatpy_wasm -j$(nproc)
 */

const fs = require('fs');
const path = require('path');

// Path to the built WASM module
const WASM_MODULE_PATH = path.join(__dirname, '../../build-wasm/src/tudatpy_wasm/tudatpy_wasm.js');

// Test results tracking
let testsPassed = 0;
let testsFailed = 0;
const failedTests = [];

/**
 * Test assertion helper
 */
function assert(condition, testName, message = '') {
    if (condition) {
        console.log(`  ✓ ${testName}`);
        testsPassed++;
    } else {
        console.log(`  ✗ ${testName}${message ? ': ' + message : ''}`);
        testsFailed++;
        failedTests.push(testName);
    }
}

/**
 * Test module existence
 */
function testModuleExists(tudat, modulePath, testName) {
    const parts = modulePath.split('.');
    let obj = tudat;
    for (const part of parts) {
        if (obj && typeof obj[part] !== 'undefined') {
            obj = obj[part];
        } else {
            assert(false, testName, `Path ${modulePath} not found at ${part}`);
            return null;
        }
    }
    assert(true, testName);
    return obj;
}

/**
 * Test function existence
 */
function testFunctionExists(tudat, funcName, testName) {
    const func = tudat[funcName];
    assert(typeof func === 'function', testName, `${funcName} is not a function`);
    return func;
}

/**
 * Test enum existence and values
 */
function testEnumExists(tudat, enumName, expectedValues, testName) {
    const enumObj = tudat[enumName];
    if (!enumObj) {
        assert(false, testName, `Enum ${enumName} not found`);
        return;
    }

    let allFound = true;
    for (const val of expectedValues) {
        if (typeof enumObj[val] === 'undefined') {
            allFound = false;
            break;
        }
    }
    assert(allFound, testName);
}

/**
 * Test class existence
 */
function testClassExists(tudat, className, testName) {
    const cls = tudat[className];
    assert(typeof cls !== 'undefined', testName, `Class ${className} not found`);
    return cls;
}

// ============================================================================
// Test Suites
// ============================================================================

function testVectorTypes(tudat) {
    console.log('\n--- Vector and Matrix Types ---');

    testClassExists(tudat, 'Vector3d', 'Vector3d class');
    testClassExists(tudat, 'Vector6d', 'Vector6d class');
    testClassExists(tudat, 'Vector7d', 'Vector7d class');
    testClassExists(tudat, 'VectorXd', 'VectorXd class');
    testClassExists(tudat, 'Matrix3d', 'Matrix3d class');
    testClassExists(tudat, 'MatrixXd', 'MatrixXd class');
}

function testConstantsModule(tudat) {
    console.log('\n--- Constants Module ---');

    // Test physical constants
    testFunctionExists(tudat, 'constants_GRAVITATIONAL_CONSTANT',
        'constants: GRAVITATIONAL_CONSTANT exists');
    testFunctionExists(tudat, 'constants_SPEED_OF_LIGHT',
        'constants: SPEED_OF_LIGHT exists');
    testFunctionExists(tudat, 'constants_ASTRONOMICAL_UNIT',
        'constants: ASTRONOMICAL_UNIT exists');
}

function testMathModule(tudat) {
    console.log('\n--- Math Module ---');

    // Test interpolators
    testEnumExists(tudat, 'math_interpolators_AvailableLookupScheme',
        ['huntingAlgorithm', 'binarySearch'],
        'math.interpolators: AvailableLookupScheme enum');

    testEnumExists(tudat, 'math_interpolators_BoundaryInterpolationType',
        ['throw_exception_at_boundary', 'use_boundary_value'],
        'math.interpolators: BoundaryInterpolationType enum');

    // Test numerical integrators
    testEnumExists(tudat, 'math_numerical_integrators_AvailableIntegrators',
        ['rungeKutta4', 'euler', 'rungeKuttaVariableStepSize'],
        'math.numerical_integrators: AvailableIntegrators enum');

    // Test root finders
    testEnumExists(tudat, 'math_root_finders_AvailableRootFinders',
        ['bisection', 'newtonRaphson', 'secant'],
        'math.root_finders: AvailableRootFinders enum');
}

function testAstroModule(tudat) {
    console.log('\n--- Astro Module ---');

    // Test element conversion
    testEnumExists(tudat, 'astro_element_conversion_KeplerianElementIndices',
        ['semi_major_axis_index', 'eccentricity_index', 'inclination_index'],
        'astro.element_conversion: KeplerianElementIndices enum');

    testFunctionExists(tudat, 'astro_element_conversion_keplerian_to_cartesian',
        'astro.element_conversion: keplerian_to_cartesian function');

    testFunctionExists(tudat, 'astro_element_conversion_cartesian_to_keplerian',
        'astro.element_conversion: cartesian_to_keplerian function');

    // Test frame conversion
    testFunctionExists(tudat, 'astro_frame_conversion_inertial_to_rsw_rotation_matrix',
        'astro.frame_conversion: inertial_to_rsw_rotation_matrix function');

    // Test two body dynamics
    testFunctionExists(tudat, 'astro_two_body_dynamics_compute_kepler_orbit_period',
        'astro.two_body_dynamics: compute_kepler_orbit_period function');

    // Test time representation
    testClassExists(tudat, 'astro_time_representation_DateTime',
        'astro.time_representation: DateTime class');
}

function testDynamicsModule(tudat) {
    console.log('\n--- Dynamics Module ---');

    // Test environment
    testClassExists(tudat, 'dynamics_environment_Body',
        'dynamics.environment: Body class');

    testClassExists(tudat, 'dynamics_environment_SystemOfBodies',
        'dynamics.environment: SystemOfBodies class');

    // Test environment setup
    testClassExists(tudat, 'dynamics_environment_setup_BodySettings',
        'dynamics.environment_setup: BodySettings class');

    testFunctionExists(tudat, 'dynamics_environment_setup_get_default_body_settings',
        'dynamics.environment_setup: get_default_body_settings function');

    testFunctionExists(tudat, 'dynamics_environment_setup_create_system_of_bodies',
        'dynamics.environment_setup: create_system_of_bodies function');

    // Test atmosphere setup
    testEnumExists(tudat, 'dynamics_environment_setup_atmosphere_AtmosphereModelType',
        ['exponential_atmosphere', 'tabulated_atmosphere', 'nrlmsise00'],
        'dynamics.environment_setup.atmosphere: AtmosphereModelType enum');

    // Test gravity field setup
    testEnumExists(tudat, 'dynamics_environment_setup_gravity_field_GravityFieldType',
        ['central', 'central_spice', 'spherical_harmonic'],
        'dynamics.environment_setup.gravity_field: GravityFieldType enum');

    // Test ephemeris setup
    testEnumExists(tudat, 'dynamics_environment_setup_ephemeris_EphemerisType',
        ['approximate_jpl', 'direct_spice', 'tabulated', 'constant'],
        'dynamics.environment_setup.ephemeris: EphemerisType enum');

    // Test propagation setup - integrator
    testEnumExists(tudat, 'dynamics_propagation_setup_integrator_AvailableIntegrators',
        ['rungeKutta4', 'euler', 'rungeKuttaVariableStepSize'],
        'dynamics.propagation_setup.integrator: AvailableIntegrators enum');

    testFunctionExists(tudat, 'dynamics_propagation_setup_integrator_runge_kutta_fixed_step_size',
        'dynamics.propagation_setup.integrator: runge_kutta_fixed_step_size function');

    // Test propagation setup - acceleration
    testEnumExists(tudat, 'dynamics_propagation_setup_acceleration_AvailableAcceleration',
        ['undefined_acceleration', 'point_mass_gravity', 'spherical_harmonic_gravity',
         'aerodynamic', 'radiation_pressure'],
        'dynamics.propagation_setup.acceleration: AvailableAcceleration enum');

    testFunctionExists(tudat, 'dynamics_propagation_setup_acceleration_point_mass_gravity',
        'dynamics.propagation_setup.acceleration: point_mass_gravity function');

    // Test propagation setup - propagator
    testEnumExists(tudat, 'dynamics_propagation_setup_propagator_TranslationalPropagatorType',
        ['cowell', 'encke', 'gauss_keplerian'],
        'dynamics.propagation_setup.propagator: TranslationalPropagatorType enum');

    testFunctionExists(tudat, 'dynamics_propagation_setup_propagator_time_termination',
        'dynamics.propagation_setup.propagator: time_termination function');

    // Test simulator
    testClassExists(tudat, 'dynamics_simulator_SingleArcSimulator',
        'dynamics.simulator: SingleArcSimulator class');

    testFunctionExists(tudat, 'dynamics_simulator_create_dynamics_simulator',
        'dynamics.simulator: create_dynamics_simulator function');
}

function testEstimationModule(tudat) {
    console.log('\n--- Estimation Module ---');

    // Test observable models setup - links
    testEnumExists(tudat, 'estimation_observable_models_setup_links_LinkEndType',
        ['unidentified_link_end', 'transmitter', 'receiver'],
        'estimation.observable_models_setup.links: LinkEndType enum');

    // Test observable models setup - model settings
    testEnumExists(tudat, 'estimation_observable_models_setup_model_settings_ObservableType',
        ['one_way_range_type', 'angular_position_type', 'position_observable_type'],
        'estimation.observable_models_setup.model_settings: ObservableType enum');

    testFunctionExists(tudat, 'estimation_observable_models_setup_model_settings_one_way_range',
        'estimation.observable_models_setup.model_settings: one_way_range function');

    // Test observations setup - viability
    testEnumExists(tudat, 'estimation_observations_setup_viability_ObservationViabilityType',
        ['minimum_elevation_angle', 'body_avoidance_angle', 'body_occultation'],
        'estimation.observations_setup.viability: ObservationViabilityType enum');

    // Test estimation analysis
    testClassExists(tudat, 'estimation_estimation_analysis_Estimator',
        'estimation.estimation_analysis: Estimator class');
}

function testInterfaceModule(tudat) {
    console.log('\n--- Interface Module ---');

    // Test SPICE functions
    testFunctionExists(tudat, 'interface_spice_convert_julian_date_to_ephemeris_time',
        'interface.spice: convert_julian_date_to_ephemeris_time function');

    testFunctionExists(tudat, 'interface_spice_convert_ephemeris_time_to_julian_date',
        'interface.spice: convert_ephemeris_time_to_julian_date function');

    testFunctionExists(tudat, 'interface_spice_get_body_cartesian_state_at_epoch',
        'interface.spice: get_body_cartesian_state_at_epoch function');

    testFunctionExists(tudat, 'interface_spice_load_kernel',
        'interface.spice: load_kernel function');
}

function testDataModule(tudat) {
    console.log('\n--- Data Module ---');

    testFunctionExists(tudat, 'data_get_resource_path',
        'data: get_resource_path function');

    testFunctionExists(tudat, 'data_get_spice_kernel_path',
        'data: get_spice_kernel_path function');

    testFunctionExists(tudat, 'data_get_gravity_models_path',
        'data: get_gravity_models_path function');
}

function testTrajectoryDesignModule(tudat) {
    console.log('\n--- Trajectory Design Module ---');

    // Test transfer trajectory
    testEnumExists(tudat, 'trajectory_design_transfer_trajectory_TransferLegTypes',
        ['unpowered_unperturbed_leg_type', 'dsm_position_based_leg_type', 'dsm_velocity_based_leg_type'],
        'trajectory_design.transfer_trajectory: TransferLegTypes enum');

    testClassExists(tudat, 'trajectory_design_transfer_trajectory_TransferLeg',
        'trajectory_design.transfer_trajectory: TransferLeg class');

    testClassExists(tudat, 'trajectory_design_transfer_trajectory_TransferTrajectory',
        'trajectory_design.transfer_trajectory: TransferTrajectory class');
}

// ============================================================================
// Main Test Runner
// ============================================================================

async function runTests() {
    console.log('='.repeat(60));
    console.log('  Tudat WASM API Comprehensive Test Suite');
    console.log('='.repeat(60));

    // Check if the WASM module exists
    if (!fs.existsSync(WASM_MODULE_PATH)) {
        console.error(`\nError: WASM module not found at ${WASM_MODULE_PATH}`);
        console.error('\nPlease build Tudat with Emscripten first:');
        console.error('  mkdir build-wasm && cd build-wasm');
        console.error('  emcmake cmake ..');
        console.error('  make tudatpy_wasm -j$(nproc)');
        process.exit(1);
    }

    console.log('\n[1] Loading WASM module...');

    try {
        const createTudatModule = require(WASM_MODULE_PATH);
        const tudat = await createTudatModule();

        console.log('[2] Running tests...\n');

        // Run all test suites
        testVectorTypes(tudat);
        testConstantsModule(tudat);
        testMathModule(tudat);
        testAstroModule(tudat);
        testDynamicsModule(tudat);
        testEstimationModule(tudat);
        testInterfaceModule(tudat);
        testDataModule(tudat);
        testTrajectoryDesignModule(tudat);

        // Print summary
        console.log('\n' + '='.repeat(60));
        console.log('  Test Summary');
        console.log('='.repeat(60));
        console.log(`  Passed: ${testsPassed}`);
        console.log(`  Failed: ${testsFailed}`);
        console.log(`  Total:  ${testsPassed + testsFailed}`);

        if (failedTests.length > 0) {
            console.log('\n  Failed tests:');
            failedTests.forEach(t => console.log(`    - ${t}`));
        }

        console.log('='.repeat(60));

        process.exit(testsFailed > 0 ? 1 : 0);

    } catch (err) {
        console.error('Error loading WASM module:', err);
        process.exit(1);
    }
}

runTests();
