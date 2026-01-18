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

// This file serves as the main entry point for the exceptions module bindings.
// Exception handling in WASM is different from native C++, so we provide
// basic exception type registration.

WASM_MODULE_PATH("exceptions")

EMSCRIPTEN_BINDINGS(tudatpy_exceptions) {
    using namespace emscripten;

    // Note: Emscripten has limited exception handling support.
    // Custom exception types can be registered here if needed.
    // Most exceptions will be converted to JavaScript errors automatically.
}

#endif
