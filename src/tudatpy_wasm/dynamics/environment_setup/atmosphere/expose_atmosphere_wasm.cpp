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

#include <tudat/simulation/environment_setup/createAtmosphereModel.h>

namespace tss = tudat::simulation_setup;
namespace ta = tudat::aerodynamics;

WASM_MODULE_PATH("dynamics_environment_setup_atmosphere")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_atmosphere) {
    using namespace emscripten;

    // AtmosphereTypes enum
    enum_<tss::AtmosphereTypes>("dynamics_environment_setup_atmosphere_AtmosphereTypes")
        .value("exponential_atmosphere", tss::exponential_atmosphere)
        .value("custom_constant_temperature_atmosphere", tss::custom_constant_temperature_atmosphere)
        .value("tabulated_atmosphere", tss::tabulated_atmosphere)
        .value("nrlmsise00", tss::nrlmsise00)
        .value("scaled_atmosphere", tss::scaled_atmosphere);

    // AtmosphereSettings base class
    class_<tss::AtmosphereSettings>("dynamics_environment_setup_atmosphere_AtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::AtmosphereSettings>>("shared_ptr_AtmosphereSettings")
        .function("getAtmosphereType", &tss::AtmosphereSettings::getAtmosphereType);

    // ExponentialAtmosphereSettings
    class_<tss::ExponentialAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_ExponentialAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::ExponentialAtmosphereSettings>>("shared_ptr_ExponentialAtmosphereSettings");

    // TabulatedAtmosphereSettings
    class_<tss::TabulatedAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_TabulatedAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::TabulatedAtmosphereSettings>>("shared_ptr_TabulatedAtmosphereSettings");

    // ScaledAtmosphereSettings
    class_<tss::ScaledAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_ScaledAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::ScaledAtmosphereSettings>>("shared_ptr_ScaledAtmosphereSettings");

    // Factory functions
    function("dynamics_environment_setup_atmosphere_exponential_atmosphere",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(const std::string&)>(
            &tss::exponentialAtmosphereSettings));

    function("dynamics_environment_setup_atmosphere_exponential_atmosphere_custom",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(const double, const double)>(
            &tss::exponentialAtmosphereSettings));

#if TUDAT_BUILD_WITH_NRLMSISE
    function("dynamics_environment_setup_atmosphere_nrlmsise00",
        &tss::nrlmsise00AtmosphereSettings);
#endif

    function("dynamics_environment_setup_atmosphere_scaled_by_constant",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(
            const std::shared_ptr<tss::AtmosphereSettings>,
            const double,
            const bool)>(&tss::scaledAtmosphereSettings));
}

#endif
