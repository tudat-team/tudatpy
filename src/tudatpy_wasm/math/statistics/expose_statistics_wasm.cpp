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
#include "../../wasm_module.h"

#include <tudat/math/statistics.h>

namespace ts = tudat::statistics;

WASM_MODULE_PATH("math_statistics")

EMSCRIPTEN_BINDINGS(tudatpy_math_statistics) {
    using namespace emscripten;

    // Statistics functions would go here
    // Most statistics functions require complex return types
    // Keeping minimal for now
}

#endif
