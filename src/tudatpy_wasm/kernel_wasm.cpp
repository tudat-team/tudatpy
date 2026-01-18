/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Main entry point for Tudatpy WebAssembly bindings.
 *    This file mirrors the structure of kernel.cpp for PyBind11.
 *
 *    The Embind bindings are organized hierarchically:
 *    - tudat.math (interpolators, integrators, root_finders, geometry, statistics)
 *    - tudat.astro (element_conversion, frame_conversion, fundamentals, gravitation, etc.)
 *    - tudat.constants
 *    - tudat.dynamics (environment, environment_setup, propagation_setup, simulator, etc.)
 *    - tudat.estimation
 *    - tudat.interface (spice)
 *    - tudat.data
 *    - tudat.trajectory_design
 *    - tudat.exceptions
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <emscripten/val.h>
#include <set>
#include <string>
#include <vector>

// Core utility headers
#include "wasm_module.h"
#include "eigen_wasm.h"
#include "stl_wasm.h"
#include "shared_ptr_wasm.h"

// Include all module binding headers (will be added as implemented)
// #include "constants/expose_constants_wasm.h"
// #include "math/expose_math_wasm.h"
// #include "astro/expose_astro_wasm.h"
// #include "dynamics/expose_dynamics_wasm.h"
// #include "estimation/expose_estimation_wasm.h"
// #include "interface/expose_interface_wasm.h"
// #include "data/expose_data_wasm.h"
// #include "trajectory_design/expose_trajectory_design_wasm.h"
// #include "exceptions/expose_exceptions_wasm.h"

using namespace emscripten;

// ============================================================================
// Module Registry Access for JavaScript
// ============================================================================

/**
 * Returns all registered module paths as a JavaScript array.
 * This allows the JavaScript namespace builder to construct the
 * tudat.math.*, tudat.astro.*, etc. hierarchy.
 */
val getRegisteredModulePaths() {
    val arr = val::array();
    const std::set<std::string>& paths = tudatpy_wasm::ModuleRegistry::instance().getPaths();
    for (const auto& path : paths) {
        arr.call<void>("push", path);
    }
    return arr;
}

/**
 * Returns version information about the Tudat WASM build.
 */
val getTudatWasmInfo() {
    val info = val::object();
    info.set("version", "0.1.0");
#ifdef TUDAT_VERSION
    info.set("tudat_version", TUDAT_VERSION);
#else
    info.set("tudat_version", "unknown");
#endif
    info.set("build_type",
#ifdef NDEBUG
        "Release"
#else
        "Debug"
#endif
    );
    info.set("eigen_version", std::to_string(EIGEN_WORLD_VERSION) + "." +
                              std::to_string(EIGEN_MAJOR_VERSION) + "." +
                              std::to_string(EIGEN_MINOR_VERSION));
    return info;
}

// ============================================================================
// Namespace Builder JavaScript Code
// ============================================================================

/**
 * Returns JavaScript code that builds nested namespaces from flat bindings.
 * Call Module.buildNamespaces() after loading the WASM module.
 */
val getNamespaceBuilderCode() {
    return val(R"JS(
function buildNamespaces() {
    const Module = this;
    const tudat = {
        math: {
            geometry: {},
            interpolators: {},
            numerical_integrators: {},
            root_finders: {},
            statistics: {}
        },
        astro: {
            element_conversion: {},
            frame_conversion: {},
            fundamentals: {},
            gravitation: {},
            polyhedron_utilities: {},
            time_representation: {},
            two_body_dynamics: {}
        },
        constants: {},
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
        interface: {
            spice: {}
        },
        data: {},
        trajectory_design: {
            shape_based_thrust: {},
            transfer_trajectory: {}
        },
        exceptions: {
            spice_exceptions: {}
        }
    };

    // Build path-to-namespace mapping
    function getNamespace(path) {
        const parts = path.split('_');
        let current = tudat;
        for (const part of parts) {
            if (current && current[part]) {
                current = current[part];
            } else {
                return null;
            }
        }
        return current;
    }

    // Map all registered functions to their namespaces
    for (const key of Object.keys(Module)) {
        if (key.startsWith('_') || key === 'tudat') continue;
        if (typeof Module[key] !== 'function' && typeof Module[key] !== 'object') continue;

        // Find the last underscore that separates prefix from name
        const lastUnderscore = key.lastIndexOf('_');
        if (lastUnderscore <= 0) continue;

        const prefix = key.substring(0, lastUnderscore);
        const name = key.substring(lastUnderscore + 1);

        const ns = getNamespace(prefix);
        if (ns) {
            ns[name] = Module[key];
        }
    }

    Module.tudat = tudat;
    return tudat;
}
)JS");
}

// ============================================================================
// Main Embind Registration
// ============================================================================

EMSCRIPTEN_BINDINGS(tudatpy_wasm_kernel) {
    // Module info and utilities
    function("getRegisteredModulePaths", &getRegisteredModulePaths);
    function("getTudatWasmInfo", &getTudatWasmInfo);
    function("getNamespaceBuilderCode", &getNamespaceBuilderCode);

    // Constants for version info
    constant("TUDATPY_WASM_VERSION", std::string("0.1.0"));
}

#endif // __EMSCRIPTEN__
