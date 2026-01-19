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

#include <tudat/simulation/environment_setup/createAtmosphereModel.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>
#include <tudat/basics/basicTypedefs.h>

namespace tss = tudat::simulation_setup;
namespace ta = tudat::aerodynamics;
namespace trf = tudat::reference_frames;
namespace tp = tudat::physical_constants;

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

    // AtmosphereDependentVariables enum
    enum_<tss::AtmosphereDependentVariables>("dynamics_environment_setup_atmosphere_AtmosphereDependentVariables")
        .value("tabulated_density", tss::density_dependent_atmosphere)
        .value("tabulated_pressure", tss::pressure_dependent_atmosphere)
        .value("tabulated_temperature", tss::temperature_dependent_atmosphere)
        .value("tabulated_gas_constant", tss::gas_constant_dependent_atmosphere)
        .value("tabulated_specific_heat_ratio", tss::specific_heat_ratio_dependent_atmosphere)
        .value("tabulated_molar_mass", tss::molar_mass_dependent_atmosphere);

    // ========================================================================
    // Wind model settings
    // ========================================================================
    class_<tss::WindModelSettings>(
        "dynamics_environment_setup_atmosphere_WindModelSettings")
        .smart_ptr<std::shared_ptr<tss::WindModelSettings>>(
            "shared_ptr_WindModelSettings");

    class_<tss::ConstantWindModelSettings, base<tss::WindModelSettings>>(
        "dynamics_environment_setup_atmosphere_ConstantWindModelSettings")
        .smart_ptr<std::shared_ptr<tss::ConstantWindModelSettings>>(
            "shared_ptr_ConstantWindModelSettings");

    class_<tss::CustomWindModelSettings, base<tss::WindModelSettings>>(
        "dynamics_environment_setup_atmosphere_CustomWindModelSettings")
        .smart_ptr<std::shared_ptr<tss::CustomWindModelSettings>>(
            "shared_ptr_CustomWindModelSettings");

    // ========================================================================
    // Atmosphere settings classes
    // ========================================================================
    class_<tss::AtmosphereSettings>("dynamics_environment_setup_atmosphere_AtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::AtmosphereSettings>>("shared_ptr_AtmosphereSettings")
        .function("getAtmosphereType", &tss::AtmosphereSettings::getAtmosphereType)
        .function("getWindSettings", &tss::AtmosphereSettings::getWindSettings)
        .function("setWindSettings", &tss::AtmosphereSettings::setWindSettings);

    class_<tss::ExponentialAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_ExponentialAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::ExponentialAtmosphereSettings>>(
            "shared_ptr_ExponentialAtmosphereSettings");

    class_<tss::CustomConstantTemperatureAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_CustomConstantTemperatureAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::CustomConstantTemperatureAtmosphereSettings>>(
            "shared_ptr_CustomConstantTemperatureAtmosphereSettings");

    class_<tss::TabulatedAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_TabulatedAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::TabulatedAtmosphereSettings>>(
            "shared_ptr_TabulatedAtmosphereSettings");

    class_<tss::ScaledAtmosphereSettings, base<tss::AtmosphereSettings>>(
        "dynamics_environment_setup_atmosphere_ScaledAtmosphereSettings")
        .smart_ptr<std::shared_ptr<tss::ScaledAtmosphereSettings>>(
            "shared_ptr_ScaledAtmosphereSettings");

    // ========================================================================
    // Wind model factory functions
    // ========================================================================
    function("dynamics_environment_setup_atmosphere_constant_wind_model",
        &tss::constantWindModelSettings);

    function("dynamics_environment_setup_atmosphere_custom_wind_model",
        &tss::customWindModelSettings);

    // ========================================================================
    // Atmosphere factory functions
    // ========================================================================

    // Exponential - predefined body name
    function("dynamics_environment_setup_atmosphere_exponential_predefined",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(const std::string&)>(
            &tss::exponentialAtmosphereSettings));

    // Exponential - full parameterization
    function("dynamics_environment_setup_atmosphere_exponential",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(
            const double, const double, const double, const double, const double)>(
            &tss::exponentialAtmosphereSettings));

    // Custom constant temperature atmosphere
    function("dynamics_environment_setup_atmosphere_custom_constant_temperature",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(
            const std::function<double(const double)>,
            const double, const double, const double)>(
            &tss::customConstantTemperatureAtmosphereSettings));

    // Tabulated atmosphere (requires file access - may have limitations in browser)
    function("dynamics_environment_setup_atmosphere_tabulated",
        &tss::tabulatedAtmosphereSettings);

#if TUDAT_BUILD_WITH_NRLMSISE
    // NRLMSISE-00 atmosphere model
    function("dynamics_environment_setup_atmosphere_nrlmsise00",
        &tss::nrlmsise00AtmosphereSettings);
#endif

    // Scaled atmosphere models
    function("dynamics_environment_setup_atmosphere_scaled_by_constant",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(
            const std::shared_ptr<tss::AtmosphereSettings>,
            const double,
            const bool)>(&tss::scaledAtmosphereSettings));

    function("dynamics_environment_setup_atmosphere_scaled_by_function",
        select_overload<std::shared_ptr<tss::AtmosphereSettings>(
            const std::shared_ptr<tss::AtmosphereSettings>,
            const std::function<double(const double)>,
            const bool)>(&tss::scaledAtmosphereSettings));
}

#endif
