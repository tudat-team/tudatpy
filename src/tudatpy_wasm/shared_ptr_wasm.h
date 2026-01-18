/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Shared pointer handling utilities for Emscripten Embind bindings.
 *    Provides macros and helpers for registering classes with shared_ptr ownership.
 */

#ifndef TUDATPY_WASM_SHARED_PTR_H
#define TUDATPY_WASM_SHARED_PTR_H

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <memory>
#include <string>

namespace tudatpy_wasm {

// ============================================================================
// Shared Pointer Registration Helpers
// ============================================================================

/**
 * Helper struct for generating shared_ptr type names consistently
 */
template<typename T>
struct SharedPtrTypeName {
    static std::string get(const std::string& baseName) {
        return "shared_ptr_" + baseName;
    }
};

} // namespace tudatpy_wasm

// ============================================================================
// Convenience Macros for Class Registration with shared_ptr
// ============================================================================

/**
 * Begin a class binding with shared_ptr support.
 * This is the most common pattern for Tudat classes.
 *
 * Usage:
 *   WASM_CLASS_SHARED_BEGIN(Ephemeris, te::Ephemeris)
 *       .function("getCartesianState", &te::Ephemeris::getCartesianState)
 *   WASM_CLASS_END
 *
 * Expands to:
 *   class_<te::Ephemeris>("Ephemeris")
 *       .smart_ptr<std::shared_ptr<te::Ephemeris>>("shared_ptr_Ephemeris")
 *       .function(...)
 */
#define WASM_CLASS_SHARED_BEGIN(name, type) \
    emscripten::class_<type>(#name) \
        .smart_ptr<std::shared_ptr<type>>("shared_ptr_" #name)

/**
 * Begin a derived class binding with shared_ptr support.
 *
 * Usage:
 *   WASM_CLASS_SHARED_DERIVED_BEGIN(TleEphemeris, te::TleEphemeris, te::Ephemeris)
 *       .constructor<const std::string&, const std::string&>()
 *   WASM_CLASS_END
 */
#define WASM_CLASS_SHARED_DERIVED_BEGIN(name, type, base) \
    emscripten::class_<type, emscripten::base<base>>(#name) \
        .smart_ptr<std::shared_ptr<type>>("shared_ptr_" #name)

/**
 * End a class binding block.
 * Just a semicolon for readability.
 */
#define WASM_CLASS_END ;

/**
 * Add a smart_ptr constructor that creates shared_ptr.
 * Use when the class needs to be constructed via factory returning shared_ptr.
 *
 * Usage within a class binding:
 *   WASM_CLASS_SHARED_BEGIN(MyClass, MyType)
 *       WASM_SHARED_CONSTRUCTOR(MyClass, MyType, double, int)
 *   WASM_CLASS_END
 */
#define WASM_SHARED_CONSTRUCTOR(name, type, ...) \
    .smart_ptr_constructor("shared_ptr_" #name, \
        &std::make_shared<type, __VA_ARGS__>)

/**
 * Register a factory function that returns shared_ptr.
 *
 * Usage:
 *   WASM_FACTORY_FUNCTION("createBodySettings", createBodySettings)
 *
 * Where createBodySettings is a function returning std::shared_ptr<BodySettings>
 */
#define WASM_FACTORY_FUNCTION(name, func) \
    emscripten::function(name, func)

/**
 * Register a factory function with module prefix.
 * Assumes WASM_MODULE_PATH has been defined.
 *
 * Usage:
 *   WASM_MODULE_PATH("dynamics_environment_setup")
 *   ...
 *   WASM_PREFIXED_FACTORY("create_body_settings", createBodySettings)
 *   // Registers as "dynamics_environment_setup_create_body_settings"
 */
#define WASM_PREFIXED_FACTORY(name, func) \
    emscripten::function(tudatpy_wasm::namespaced(__wasm_module_path__, name).c_str(), func)

// ============================================================================
// Allow Subclass Registration
// ============================================================================

/**
 * Macro to register that a derived type can be passed as base type.
 * This is needed for polymorphic function parameters.
 *
 * Usage:
 *   WASM_ALLOW_SUBCLASS(TleEphemeris, te::TleEphemeris, te::Ephemeris)
 */
#define WASM_ALLOW_SUBCLASS(name, derived, base) \
    emscripten::class_<derived, emscripten::base<base>>(#name) \
        .smart_ptr<std::shared_ptr<derived>>("shared_ptr_" #name)

// ============================================================================
// Convenience for Settings Classes Pattern
// ============================================================================

/**
 * Many Tudat "Settings" classes follow a pattern:
 * - Base abstract class (e.g., IntegratorSettings)
 * - Multiple derived concrete classes
 * - Factory functions that return shared_ptr<Base>
 *
 * This macro helps register the common pattern.
 */
#define WASM_SETTINGS_BASE(name, type) \
    emscripten::class_<type>(#name) \
        .smart_ptr<std::shared_ptr<type>>("shared_ptr_" #name)

#define WASM_SETTINGS_DERIVED(name, type, base) \
    emscripten::class_<type, emscripten::base<base>>(#name) \
        .smart_ptr<std::shared_ptr<type>>("shared_ptr_" #name) \
        .allow_subclass<type, base>()

// ============================================================================
// Vector Registration for shared_ptr Types
// ============================================================================

/**
 * Register a vector of shared_ptr<T>.
 *
 * Usage:
 *   WASM_REGISTER_SHARED_VECTOR(AccelerationSettings, AccelerationSettings)
 *
 * Creates VectorSharedAccelerationSettings type
 */
#define WASM_REGISTER_SHARED_VECTOR(name, type) \
    emscripten::register_vector<std::shared_ptr<type>>("VectorShared" #name)

/**
 * Register a map from string to shared_ptr<T>.
 *
 * Usage:
 *   WASM_REGISTER_SHARED_MAP_STRING(Body, Body)
 *
 * Creates MapStringSharedBody type
 */
#define WASM_REGISTER_SHARED_MAP_STRING(name, type) \
    emscripten::register_map<std::string, std::shared_ptr<type>>("MapStringShared" #name)

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_SHARED_PTR_H
