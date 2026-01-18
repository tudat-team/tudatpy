/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Module hierarchy utilities for Emscripten Embind bindings.
 *    Provides namespace emulation since Embind has no submodule concept.
 */

#ifndef TUDATPY_WASM_MODULE_H
#define TUDATPY_WASM_MODULE_H

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <string>
#include <set>

namespace tudatpy_wasm {

/**
 * Registry for tracking all module paths for JavaScript namespace construction.
 * This singleton collects all registered module paths that can be queried
 * from JavaScript to build the nested namespace structure.
 */
class ModuleRegistry {
public:
    static ModuleRegistry& instance() {
        static ModuleRegistry reg;
        return reg;
    }

    void registerPath(const std::string& fullPath) {
        paths_.insert(fullPath);
    }

    const std::set<std::string>& getPaths() const { return paths_; }

private:
    ModuleRegistry() = default;
    std::set<std::string> paths_;
};

/**
 * Helper to create a prefixed binding name from module path and symbol name.
 * Example: namespaced("astro_element_conversion", "cartesian_to_keplerian")
 *          returns "astro_element_conversion_cartesian_to_keplerian"
 */
inline std::string namespaced(const std::string& modulePath, const std::string& name) {
    if (modulePath.empty()) {
        return name;
    }
    return modulePath + "_" + name;
}

} // namespace tudatpy_wasm

// ============================================================================
// Convenience Macros for Module Registration
// ============================================================================

/**
 * Declares the current module path for a binding file.
 * Must be used at file scope before using WASM_* macros.
 *
 * Example:
 *   WASM_MODULE_PATH("astro_element_conversion")
 */
#define WASM_MODULE_PATH(path) \
    namespace { \
        static const std::string __wasm_module_path__ = path; \
        struct __ModuleRegistrar { \
            __ModuleRegistrar() { tudatpy_wasm::ModuleRegistry::instance().registerPath(path); } \
        }; \
        static __ModuleRegistrar __module_registrar__; \
    }

/**
 * Registers a function with the current module prefix.
 *
 * Example:
 *   WASM_FUNCTION("cartesian_to_keplerian", &convertCartesianToKeplerian)
 *   // Registers as "astro_element_conversion_cartesian_to_keplerian" if
 *   // WASM_MODULE_PATH("astro_element_conversion") was used
 */
#define WASM_FUNCTION(name, func) \
    emscripten::function(tudatpy_wasm::namespaced(__wasm_module_path__, name).c_str(), func)

/**
 * Starts a class binding with the current module prefix.
 * Returns emscripten::class_<T> for chaining.
 *
 * Example:
 *   WASM_CLASS("IntegratorSettings", tni::IntegratorSettings<double>)
 *       .constructor<>()
 *       .function("getTimeStep", &tni::IntegratorSettings<double>::getTimeStep);
 */
#define WASM_CLASS(name, type) \
    emscripten::class_<type>(tudatpy_wasm::namespaced(__wasm_module_path__, name).c_str())

/**
 * Starts an enum binding with the current module prefix.
 *
 * Example:
 *   WASM_ENUM("CoefficientSets", tni::CoefficientSets)
 *       .value("rk_4", tni::rungeKutta4Classic)
 *       .value("rkf_45", tni::rungeKuttaFehlberg45);
 */
#define WASM_ENUM(name, type) \
    emscripten::enum_<type>(tudatpy_wasm::namespaced(__wasm_module_path__, name).c_str())

/**
 * Registers a constant/attribute with the current module prefix.
 *
 * Example:
 *   WASM_CONSTANT("SPEED_OF_LIGHT", physical_constants::SPEED_OF_LIGHT)
 */
#define WASM_CONSTANT(name, value) \
    emscripten::constant(tudatpy_wasm::namespaced(__wasm_module_path__, name).c_str(), value)

// ============================================================================
// JavaScript Namespace Builder (to be included in output JS)
// ============================================================================

/**
 * The following JavaScript code should be called after Module initialization
 * to build the nested namespace structure from flat prefixed names:
 *
 * function buildTudatNamespaces(Module) {
 *     const tudat = {};
 *     const paths = Module.getRegisteredPaths();
 *
 *     // Create namespace hierarchy
 *     for (const path of paths) {
 *         const parts = path.split('_');
 *         let current = tudat;
 *         for (const part of parts) {
 *             if (!current[part]) current[part] = {};
 *             current = current[part];
 *         }
 *     }
 *
 *     // Map prefixed functions to namespaces
 *     for (const key of Object.keys(Module)) {
 *         if (key.startsWith('_') || typeof Module[key] !== 'function') continue;
 *
 *         const lastUnderscore = key.lastIndexOf('_');
 *         if (lastUnderscore > 0) {
 *             const prefix = key.substring(0, lastUnderscore);
 *             const funcName = key.substring(lastUnderscore + 1);
 *
 *             // Navigate to namespace
 *             const parts = prefix.split('_');
 *             let target = tudat;
 *             let valid = true;
 *             for (const part of parts) {
 *                 if (target[part]) target = target[part];
 *                 else { valid = false; break; }
 *             }
 *             if (valid) target[funcName] = Module[key];
 *         }
 *     }
 *
 *     Module.tudat = tudat;
 *     return tudat;
 * }
 */

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_MODULE_H
