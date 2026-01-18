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
#include "../../stl_wasm.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/propagation_setup/createAccelerationModels.h>
#include <tudat/simulation/propagation_setup/createTorqueModel.h>
#include <tudat/simulation/propagation_setup/createMassRateModels.h>

namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("dynamics_propagation_setup")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup) {
    using namespace emscripten;

    // Factory function to create acceleration models
    function("dynamics_propagation_setup_create_acceleration_models",
        select_overload<tba::AccelerationMap(
            const tss::SystemOfBodies&,
            const tss::SelectedAccelerationMap&,
            const std::vector<std::string>&,
            const std::vector<std::string>&)>(
            &tss::createAccelerationModelsMap));

    // Factory function to create torque models
    function("dynamics_propagation_setup_create_torque_models",
        &tss::createTorqueModelsMap);

    // Factory function to create mass rate models
    function("dynamics_propagation_setup_create_mass_rate_models",
        &tss::createMassRateModelsMap);
}

#endif
