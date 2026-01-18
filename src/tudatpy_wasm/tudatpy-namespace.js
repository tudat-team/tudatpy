/**
 * Tudatpy WASM Namespace Builder
 *
 * This module transforms the flat Embind bindings into a nested namespace
 * structure that mirrors the Python API.
 *
 * Usage:
 *   import createTudatModule from './tudatpy_wasm.js';
 *
 *   async function main() {
 *       const Module = await createTudatModule();
 *       const tudat = buildTudatNamespaces(Module);
 *
 *       // Now use the familiar Python-like API:
 *       const keplerian = tudat.astro.element_conversion.cartesian_to_keplerian(
 *           state, gravitationalParameter
 *       );
 *   }
 */

/**
 * Builds nested namespace structure from flat Embind bindings.
 *
 * The Embind bindings use prefixed names like:
 *   - astro_element_conversion_cartesian_to_keplerian
 *   - dynamics_propagation_setup_integrator_runge_kutta_4
 *
 * This function creates:
 *   - tudat.astro.element_conversion.cartesian_to_keplerian
 *   - tudat.dynamics.propagation_setup.integrator.runge_kutta_4
 *
 * @param {Object} Module - The Emscripten module
 * @returns {Object} The tudat namespace object
 */
function buildTudatNamespaces(Module) {
    // Define the complete namespace hierarchy
    const tudat = {
        // Math module
        math: {
            geometry: {},
            interpolators: {},
            numerical_integrators: {},
            root_finders: {},
            statistics: {}
        },

        // Astrodynamics module
        astro: {
            element_conversion: {},
            frame_conversion: {},
            fundamentals: {},
            gravitation: {},
            polyhedron_utilities: {},
            time_representation: {},
            two_body_dynamics: {}
        },

        // Constants module
        constants: {},

        // Dynamics module (largest)
        dynamics: {
            environment: {},
            environment_setup: {
                aerodynamic_coefficients: {},
                atmosphere: {},
                ephemeris: {},
                gravity_field: {},
                gravity_field_variation: {},
                ground_station: {},
                radiation_pressure: {},
                rigid_body: {},
                rotation_model: {},
                shape: {},
                shape_deformation: {},
                vehicle_systems: {}
            },
            parameters: {},
            parameters_setup: {},
            propagation: {},
            propagation_setup: {
                acceleration: {},
                dependent_variable: {},
                integrator: {},
                mass_rate: {},
                propagator: {},
                thrust: {},
                torque: {}
            },
            simulator: {}
        },

        // Estimation module
        estimation: {
            estimation_analysis: {},
            observable_models: {
                observables_simulation: {}
            },
            observable_models_setup: {
                biases: {},
                light_time_corrections: {},
                links: {},
                model_settings: {}
            },
            observations: {
                observations_geometry: {},
                observations_processing: {}
            },
            observations_setup: {
                ancillary_settings: {},
                observations_dependent_variables: {},
                observations_simulation_settings: {},
                observations_wrapper: {},
                random_noise: {},
                viability: {}
            }
        },

        // Interface module
        interface: {
            spice: {}
        },

        // Data module
        data: {},

        // Trajectory design module
        trajectory_design: {
            shape_based_thrust: {},
            transfer_trajectory: {}
        },

        // Exceptions module
        exceptions: {
            spice_exceptions: {}
        }
    };

    /**
     * Navigate to a namespace by path (e.g., "astro_element_conversion")
     */
    function getNamespace(path) {
        const parts = path.split('_');
        let current = tudat;

        for (const part of parts) {
            if (current && typeof current === 'object' && part in current) {
                current = current[part];
            } else {
                return null;
            }
        }

        return current;
    }

    /**
     * Find the longest matching prefix that corresponds to a namespace
     */
    function findNamespaceAndName(key) {
        const parts = key.split('_');

        // Try progressively shorter prefixes to find the namespace
        for (let i = parts.length - 1; i > 0; i--) {
            const prefix = parts.slice(0, i).join('_');
            const ns = getNamespace(prefix);

            if (ns !== null) {
                const name = parts.slice(i).join('_');
                return { namespace: ns, name: name };
            }
        }

        return null;
    }

    // Map all Module properties to their namespaces
    for (const key of Object.keys(Module)) {
        // Skip internal properties
        if (key.startsWith('_') || key === 'tudat') continue;

        const value = Module[key];

        // Skip non-function/non-object properties (except for constants which are numbers)
        if (typeof value !== 'function' && typeof value !== 'object' && typeof value !== 'number') {
            continue;
        }

        // Find the namespace and name
        const result = findNamespaceAndName(key);
        if (result && result.name) {
            result.namespace[result.name] = value;
        }
    }

    // Attach to Module for easy access
    Module.tudat = tudat;

    return tudat;
}

/**
 * Helper to create Vector6d from array
 */
function createVector6d(Module, array) {
    if (array.length !== 6) {
        throw new Error('Vector6d requires exactly 6 elements');
    }
    const vec = new Module.Vector6d();
    for (let i = 0; i < 6; i++) {
        vec.set(i, array[i]);
    }
    return vec;
}

/**
 * Helper to create Vector3d from array
 */
function createVector3d(Module, array) {
    if (array.length !== 3) {
        throw new Error('Vector3d requires exactly 3 elements');
    }
    return new Module.Vector3d(array[0], array[1], array[2]);
}

/**
 * Helper to convert Vector6d to JavaScript array
 */
function vector6dToArray(vec) {
    return vec.toArray();
}

/**
 * Helper to convert Vector3d to JavaScript array
 */
function vector3dToArray(vec) {
    return vec.toArray();
}

// Export for both ES modules and CommonJS
if (typeof module !== 'undefined' && module.exports) {
    module.exports = {
        buildTudatNamespaces,
        createVector6d,
        createVector3d,
        vector6dToArray,
        vector3dToArray
    };
}

if (typeof window !== 'undefined') {
    window.buildTudatNamespaces = buildTudatNamespaces;
    window.createVector6d = createVector6d;
    window.createVector3d = createVector3d;
    window.vector6dToArray = vector6dToArray;
    window.vector3dToArray = vector3dToArray;
}
