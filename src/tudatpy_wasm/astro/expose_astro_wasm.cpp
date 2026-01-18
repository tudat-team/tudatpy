/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../wasm_module.h"

// This file serves as the main entry point for the astro module bindings.
// All submodules are compiled separately and linked together.
// Submodules: element_conversion, frame_conversion, fundamentals, gravitation,
//             polyhedron_utilities, time_representation, two_body_dynamics

WASM_MODULE_PATH("astro")

EMSCRIPTEN_BINDINGS(tudatpy_astro) {
    using namespace emscripten;
    // All submodule bindings are registered through their individual
    // EMSCRIPTEN_BINDINGS blocks with prefixed names.
}

#endif
