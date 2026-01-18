/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    STL container registrations for Emscripten Embind bindings.
 *    Registers common vector and map types used in Tudat.
 */

#ifndef TUDATPY_WASM_STL_H
#define TUDATPY_WASM_STL_H

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include <vector>
#include <map>
#include <string>
#include <memory>
#include <Eigen/Core>

#include "eigen_wasm.h"

namespace tudatpy_wasm {

// Forward declarations of common Tudat types used in STL containers
// These will be properly defined in their respective binding files

} // namespace tudatpy_wasm

// Embind registrations are in stl_wasm.cpp (not in header to avoid duplicate symbols)

// ============================================================================
// Helper Macros for Registering Additional Container Types
// ============================================================================

/**
 * Register a std::vector<T> with a custom name.
 * Use when you need to register vectors of custom types.
 *
 * Example:
 *   REGISTER_WASM_VECTOR(AccelerationSettings, std::shared_ptr<AccelerationSettings>)
 */
#define REGISTER_WASM_VECTOR(name, type) \
    emscripten::register_vector<type>("Vector" #name)

/**
 * Register a std::map<K, V> with a custom name.
 *
 * Example:
 *   REGISTER_WASM_MAP(StringToBody, std::string, std::shared_ptr<Body>)
 */
#define REGISTER_WASM_MAP(name, key_type, value_type) \
    emscripten::register_map<key_type, value_type>("Map" #name)

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_STL_H
