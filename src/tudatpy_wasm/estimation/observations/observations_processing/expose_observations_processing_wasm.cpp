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
#include "../../../wasm_module.h"
#include "../../../stl_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/estimation_setup/observationsProcessing.h>

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observations_observations_processing")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_observations_processing) {
    using namespace emscripten;

    // ObservationCollection processing utilities
    // The main processing functions work on ObservationCollection objects
    // defined in the parent observations module
}

#endif
