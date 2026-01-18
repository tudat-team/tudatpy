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

#include <tudat/simulation/estimation_setup/createLightTimeCorrection.h>

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observable_models_setup_light_time_corrections")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observable_models_setup_light_time_corrections) {
    using namespace emscripten;

    // LightTimeCorrectionType enum
    enum_<tom::LightTimeCorrectionType>(
        "estimation_observable_models_setup_light_time_corrections_LightTimeCorrectionType")
        .value("first_order_relativistic", tom::first_order_relativistic)
        .value("tabulated_tropospheric", tom::tabulated_tropospheric)
        .value("saastamoinen_tropospheric", tom::saastamoinen_tropospheric)
        .value("tabulated_ionospheric", tom::tabulated_ionospheric);

    // LightTimeCorrectionSettings base class
    class_<tom::LightTimeCorrectionSettings>(
        "estimation_observable_models_setup_light_time_corrections_LightTimeCorrectionSettings")
        .smart_ptr<std::shared_ptr<tom::LightTimeCorrectionSettings>>(
            "shared_ptr_LightTimeCorrectionSettings");

    // FirstOrderRelativisticLightTimeCorrectionSettings
    class_<tom::FirstOrderRelativisticLightTimeCorrectionSettings,
           base<tom::LightTimeCorrectionSettings>>(
        "estimation_observable_models_setup_light_time_corrections_FirstOrderRelativisticLightTimeCorrectionSettings")
        .smart_ptr<std::shared_ptr<tom::FirstOrderRelativisticLightTimeCorrectionSettings>>(
            "shared_ptr_FirstOrderRelativisticLightTimeCorrectionSettings");

    // Factory functions
    function("estimation_observable_models_setup_light_time_corrections_first_order_relativistic",
        &tom::firstOrderRelativisticLightTimeCorrectionSettings);
}

#endif
