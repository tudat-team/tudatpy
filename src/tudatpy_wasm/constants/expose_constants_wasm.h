/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Constants module bindings for WASM.
 *    Mirrors: src/tudatpy/constants/expose_constants.cpp
 */

#ifndef TUDATPY_WASM_CONSTANTS_H
#define TUDATPY_WASM_CONSTANTS_H

#ifdef __EMSCRIPTEN__

namespace tudatpy_wasm {
namespace constants {

// Note: For Embind, constants are registered in EMSCRIPTEN_BINDINGS blocks
// and don't need explicit expose functions like PyBind11

} // namespace constants
} // namespace tudatpy_wasm

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_CONSTANTS_H
