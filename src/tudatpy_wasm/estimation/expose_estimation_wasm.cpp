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

// This file serves as the main entry point for the estimation module bindings.
// All submodules are compiled separately and linked together.
// The module hierarchy is managed through prefixed function/class names.

WASM_MODULE_PATH("estimation")

EMSCRIPTEN_BINDINGS(tudatpy_estimation) {
    using namespace emscripten;

    // The estimation module contains the following submodules:
    // - observable_models_setup: links, model_settings, light_time_corrections, biases
    // - observable_models: observables_simulation
    // - observations: observations_geometry, observations_processing
    // - observations_setup: ancillary_settings, observations_simulation_settings, random_noise,
    //                       viability, observations_dependent_variables, observations_wrapper
    // - estimation_analysis: Estimator, EstimationInput/Output, CovarianceAnalysis

    // All submodule bindings are registered through their individual
    // EMSCRIPTEN_BINDINGS blocks with prefixed names.
}

#endif
