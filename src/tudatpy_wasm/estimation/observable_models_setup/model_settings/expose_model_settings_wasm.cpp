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

#include <tudat/astro/observation_models/observableTypes.h>
#include <tudat/simulation/estimation_setup/createObservationModel.h>

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observable_models_setup_model_settings")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observable_models_setup_model_settings) {
    using namespace emscripten;

    // ObservableType enum
    enum_<tom::ObservableType>("estimation_observable_models_setup_model_settings_ObservableType")
        .value("undefined_observation_model", tom::undefined_observation_model)
        .value("one_way_range", tom::one_way_range)
        .value("angular_position", tom::angular_position)
        .value("position_observable", tom::position_observable)
        .value("one_way_doppler", tom::one_way_doppler)
        .value("one_way_differenced_range", tom::one_way_differenced_range)
        .value("n_way_range", tom::n_way_range)
        .value("two_way_doppler", tom::two_way_doppler)
        .value("euler_angle_313_observable", tom::euler_angle_313_observable)
        .value("velocity_observable", tom::velocity_observable)
        .value("relative_angular_position", tom::relative_angular_position)
        .value("n_way_differenced_range", tom::n_way_differenced_range)
        .value("relative_position_observable", tom::relative_position_observable)
        .value("dsn_one_way_averaged_doppler", tom::dsn_one_way_averaged_doppler)
        .value("dsn_n_way_averaged_doppler", tom::dsn_n_way_averaged_doppler)
        .value("doppler_measured_frequency", tom::doppler_measured_frequency)
        .value("dsn_n_way_range", tom::dsn_n_way_range);

    // ObservationModelSettings base class
    class_<tom::ObservationModelSettings>("estimation_observable_models_setup_model_settings_ObservationModelSettings")
        .smart_ptr<std::shared_ptr<tom::ObservationModelSettings>>(
            "shared_ptr_ObservationModelSettings");

    // Utility functions for observable types
    function("estimation_observable_models_setup_model_settings_get_observable_name",
        select_overload<std::string(const tom::ObservableType, const int)>(&tom::getObservableName));

    function("estimation_observable_models_setup_model_settings_get_observable_type",
        &tom::getObservableType);

    function("estimation_observable_models_setup_model_settings_get_observable_size",
        &tom::getObservableSize);

    // Factory functions for observation model settings
    function("estimation_observable_models_setup_model_settings_one_way_range",
        &tom::oneWayRangeSettings);

    function("estimation_observable_models_setup_model_settings_angular_position",
        &tom::angularPositionSettings);

    function("estimation_observable_models_setup_model_settings_relative_angular_position",
        &tom::relativeAngularPositionSettings);

    function("estimation_observable_models_setup_model_settings_position_observable",
        &tom::positionObservableSettings);

    function("estimation_observable_models_setup_model_settings_velocity_observable",
        &tom::velocityObservableSettings);

    function("estimation_observable_models_setup_model_settings_one_way_open_loop_doppler",
        &tom::oneWayOpenLoopDoppler);
}

#endif
