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
#include "../../../eigen_wasm.h"
#include "../../../stl_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/astro/observation_models/observationModel.h>

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observations_observations_geometry")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_observations_geometry) {
    using namespace emscripten;

    // ObservationAncilliarySimulationSettings class
    class_<tom::ObservationAncilliarySimulationSettings>(
        "estimation_observations_observations_geometry_ObservationAncilliarySimulationSettings")
        .smart_ptr<std::shared_ptr<tom::ObservationAncilliarySimulationSettings>>(
            "shared_ptr_ObservationAncilliarySimulationSettings")
        .constructor<>();
}

#endif
