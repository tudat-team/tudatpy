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
#include "../stl_wasm.h"

WASM_MODULE_PATH("data")

EMSCRIPTEN_BINDINGS(tudatpy_data) {
    using namespace emscripten;

    // Data module placeholder
    // Note: Most data reading functionality requires file system access
    // which is limited in browser environments.
    // Add bindings here as needed for browser-compatible data operations.
}

#endif
