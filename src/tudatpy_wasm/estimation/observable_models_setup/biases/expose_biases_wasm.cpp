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

#include <tudat/simulation/estimation_setup/createObservationModel.h>

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observable_models_setup_biases")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observable_models_setup_biases) {
    using namespace emscripten;

    // ObservationBiasTypes enum
    enum_<tom::ObservationBiasTypes>("estimation_observable_models_setup_biases_ObservationBiasTypes")
        .value("multiple_observation_biases", tom::multiple_observation_biases)
        .value("constant_absolute_bias", tom::constant_absolute_bias)
        .value("constant_relative_bias", tom::constant_relative_bias)
        .value("arc_wise_constant_absolute_bias", tom::arc_wise_constant_absolute_bias)
        .value("arc_wise_constant_relative_bias", tom::arc_wise_constant_relative_bias)
        .value("constant_time_drift_bias", tom::constant_time_drift_bias)
        .value("arc_wise_time_drift_bias", tom::arc_wise_time_drift_bias)
        .value("constant_time_bias", tom::constant_time_bias)
        .value("arc_wise_time_bias", tom::arc_wise_time_bias);

    // ObservationBiasSettings base class
    class_<tom::ObservationBiasSettings>("estimation_observable_models_setup_biases_ObservationBiasSettings")
        .smart_ptr<std::shared_ptr<tom::ObservationBiasSettings>>(
            "shared_ptr_ObservationBiasSettings");

    // ConstantObservationBiasSettings
    class_<tom::ConstantObservationBiasSettings, base<tom::ObservationBiasSettings>>(
        "estimation_observable_models_setup_biases_ConstantObservationBiasSettings")
        .smart_ptr<std::shared_ptr<tom::ConstantObservationBiasSettings>>(
            "shared_ptr_ConstantObservationBiasSettings");

    // ArcWiseConstantObservationBiasSettings
    class_<tom::ArcWiseConstantObservationBiasSettings, base<tom::ObservationBiasSettings>>(
        "estimation_observable_models_setup_biases_ArcWiseConstantObservationBiasSettings")
        .smart_ptr<std::shared_ptr<tom::ArcWiseConstantObservationBiasSettings>>(
            "shared_ptr_ArcWiseConstantObservationBiasSettings");

    // MultipleObservationBiasSettings
    class_<tom::MultipleObservationBiasSettings, base<tom::ObservationBiasSettings>>(
        "estimation_observable_models_setup_biases_MultipleObservationBiasSettings")
        .smart_ptr<std::shared_ptr<tom::MultipleObservationBiasSettings>>(
            "shared_ptr_MultipleObservationBiasSettings");

    // Factory functions for bias settings
    function("estimation_observable_models_setup_biases_constant_absolute_bias",
        &tom::constantAbsoluteBias);

    function("estimation_observable_models_setup_biases_constant_relative_bias",
        &tom::constantRelativeBias);

    function("estimation_observable_models_setup_biases_multiple_biases",
        &tom::multipleObservationBiasSettings);
}

#endif
