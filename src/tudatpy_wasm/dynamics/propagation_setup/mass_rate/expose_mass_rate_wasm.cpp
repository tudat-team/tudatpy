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
#include <tudat/astro/basic_astro/massRateModel.h>

namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("dynamics_propagation_setup_mass_rate")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_mass_rate) {
    using namespace emscripten;

    // ========================================================================
    // AvailableMassRateModels enum
    // ========================================================================
    enum_<tba::AvailableMassRateModels>("dynamics_propagation_setup_mass_rate_AvailableMassRateModels")
        .value("undefined_mass_rate_type", tba::undefined_mass_rate_model)
        .value("custom_mass_rate_type", tba::custom_mass_rate_model)
        .value("from_thrust_mass_rate_type", tba::from_thrust_mass_rate_model);

    // ========================================================================
    // MassRateModelSettings base class
    // ========================================================================
    class_<tss::MassRateModelSettings>("dynamics_propagation_setup_mass_rate_MassRateModelSettings")
        .smart_ptr<std::shared_ptr<tss::MassRateModelSettings>>("shared_ptr_MassRateModelSettings");

    // ========================================================================
    // FromThrustMassRateSettings derived class
    // ========================================================================
    class_<tss::FromThrustMassRateSettings, base<tss::MassRateModelSettings>>(
        "dynamics_propagation_setup_mass_rate_FromThrustMassRateSettings")
        .smart_ptr<std::shared_ptr<tss::FromThrustMassRateSettings>>("shared_ptr_FromThrustMassRateSettings");

    // ========================================================================
    // CustomMassRateSettings derived class
    // ========================================================================
    class_<tss::CustomMassRateSettings, base<tss::MassRateModelSettings>>(
        "dynamics_propagation_setup_mass_rate_CustomMassRateSettings")
        .smart_ptr<std::shared_ptr<tss::CustomMassRateSettings>>("shared_ptr_CustomMassRateSettings");

    // ========================================================================
    // Factory functions
    // ========================================================================

    // From thrust mass rate
    function("dynamics_propagation_setup_mass_rate_from_thrust",
        &tss::fromThrustMassRate);

    // Custom mass rate (function of time)
    function("dynamics_propagation_setup_mass_rate_custom_mass_rate",
        &tss::customMassRate);
}

#endif
