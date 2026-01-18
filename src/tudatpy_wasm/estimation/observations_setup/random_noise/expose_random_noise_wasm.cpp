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

#include <tudat/simulation/estimation_setup/observationSimulationSettings.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("estimation_observations_setup_random_noise")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_setup_random_noise) {
    using namespace emscripten;

    // ObservationNoiseModel base class
    class_<tss::ObservationNoiseModel>("estimation_observations_setup_random_noise_ObservationNoiseModel")
        .smart_ptr<std::shared_ptr<tss::ObservationNoiseModel>>("shared_ptr_ObservationNoiseModel");

    // UnivariateGaussianObservationNoiseModel class
    class_<tss::UnivariateGaussianObservationNoiseModel, base<tss::ObservationNoiseModel>>(
        "estimation_observations_setup_random_noise_UnivariateGaussianObservationNoiseModel")
        .smart_ptr<std::shared_ptr<tss::UnivariateGaussianObservationNoiseModel>>(
            "shared_ptr_UnivariateGaussianObservationNoiseModel")
        .constructor<const double, const double, const int>()
        .function("getObservationNoise", &tss::UnivariateGaussianObservationNoiseModel::getObservationNoise);
}

#endif
