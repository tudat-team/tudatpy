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
#include "../../eigen_wasm.h"
#include "../../stl_wasm.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/estimation_setup/fitOrbitToEphemeris.h>

namespace tss = tudat::simulation_setup;
namespace tep = tudat::estimatable_parameters;

WASM_MODULE_PATH("estimation_estimation_analysis_ephemeris_fit")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_estimation_analysis_ephemeris_fit) {
    using namespace emscripten;

    // createBestFitToCurrentEphemeris function
    // This fits an orbit to an existing ephemeris using estimation
    function("estimation_estimation_analysis_create_best_fit_to_ephemeris",
        &tss::createBestFitToCurrentEphemeris<double, double>);
}

#endif
