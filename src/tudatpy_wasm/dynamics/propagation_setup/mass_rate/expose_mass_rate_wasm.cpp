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

#include <tudat/simulation/propagation_setup/createMassRateModels.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_propagation_setup_mass_rate")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_mass_rate) {
    using namespace emscripten;

    // MassRateModelSettings base class
    class_<tss::MassRateModelSettings>("dynamics_propagation_setup_mass_rate_MassRateModelSettings")
        .smart_ptr<std::shared_ptr<tss::MassRateModelSettings>>("shared_ptr_MassRateModelSettings");

    // FromThrustMassRateSettings
    class_<tss::FromThrustMassRateSettings, base<tss::MassRateModelSettings>>(
        "dynamics_propagation_setup_mass_rate_FromThrustMassRateSettings")
        .smart_ptr<std::shared_ptr<tss::FromThrustMassRateSettings>>("shared_ptr_FromThrustMassRateSettings");

    // Factory functions
    function("dynamics_propagation_setup_mass_rate_from_thrust",
        &tss::fromThrustMassRate);
}

#endif
