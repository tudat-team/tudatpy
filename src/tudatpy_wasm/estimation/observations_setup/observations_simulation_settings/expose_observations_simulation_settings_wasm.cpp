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
namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observations_setup_observations_simulation_settings")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_setup_observations_simulation_settings) {
    using namespace emscripten;

    // ObservationSimulationSettings base class
    class_<tss::ObservationSimulationSettings<double>>(
        "estimation_observations_setup_observations_simulation_settings_ObservationSimulationSettings")
        .smart_ptr<std::shared_ptr<tss::ObservationSimulationSettings<double>>>(
            "shared_ptr_ObservationSimulationSettings");

    // TabulatedObservationSimulationSettings derived class
    class_<tss::TabulatedObservationSimulationSettings<double>,
           base<tss::ObservationSimulationSettings<double>>>(
        "estimation_observations_setup_observations_simulation_settings_TabulatedObservationSimulationSettings")
        .smart_ptr<std::shared_ptr<tss::TabulatedObservationSimulationSettings<double>>>(
            "shared_ptr_TabulatedObservationSimulationSettings");

    // Factory functions
    function("estimation_observations_setup_observations_simulation_settings_tabulated_settings",
        &tss::tabulatedObservationSimulationSettings<double>);

    function("estimation_observations_setup_observations_simulation_settings_continuous_arc_simulation_settings",
        &tss::perArcObservationSimulationSettingsList<double>);
}

#endif
