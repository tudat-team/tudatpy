/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Element conversion module bindings for WASM.
 *    Mirrors: src/tudatpy/astro/element_conversion/expose_element_conversion.cpp
 */

#ifndef TUDATPY_WASM_ELEMENT_CONVERSION_H
#define TUDATPY_WASM_ELEMENT_CONVERSION_H

#ifdef __EMSCRIPTEN__

namespace tudatpy_wasm {
namespace astro {
namespace element_conversion {

// Embind bindings are registered via EMSCRIPTEN_BINDINGS macro

} // namespace element_conversion
} // namespace astro
} // namespace tudatpy_wasm

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_ELEMENT_CONVERSION_H
