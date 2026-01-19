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

#include <tudat/simulation/environment_setup/createEphemeris.h>

namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;

WASM_MODULE_PATH("dynamics_environment_setup_ephemeris")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_ephemeris) {
    using namespace emscripten;

    // EphemerisType enum
    enum_<tss::EphemerisType>("dynamics_environment_setup_ephemeris_EphemerisType")
        .value("approximate_planet_positions", tss::approximate_planet_positions)
        .value("direct_spice_ephemeris", tss::direct_spice_ephemeris)
        .value("tabulated_ephemeris", tss::tabulated_ephemeris)
        .value("interpolated_spice", tss::interpolated_spice)
        .value("constant_ephemeris", tss::constant_ephemeris)
        .value("kepler_ephemeris", tss::kepler_ephemeris)
        .value("custom_ephemeris", tss::custom_ephemeris)
        .value("direct_tle_ephemeris", tss::direct_tle_ephemeris)
        .value("scaled_ephemeris", tss::scaled_ephemeris);

    // EphemerisSettings base class
    class_<tss::EphemerisSettings>("dynamics_environment_setup_ephemeris_EphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::EphemerisSettings>>("shared_ptr_EphemerisSettings")
        .function("getEphemerisType", &tss::EphemerisSettings::getEphemerisType)
        .function("getFrameOrigin", &tss::EphemerisSettings::getFrameOrigin)
        .function("getFrameOrientation", &tss::EphemerisSettings::getFrameOrientation);

    // DirectSpiceEphemerisSettings
    class_<tss::DirectSpiceEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_DirectSpiceEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::DirectSpiceEphemerisSettings>>("shared_ptr_DirectSpiceEphemerisSettings");

    // ConstantEphemerisSettings
    class_<tss::ConstantEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_ConstantEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::ConstantEphemerisSettings>>("shared_ptr_ConstantEphemerisSettings");

    // KeplerEphemerisSettings
    class_<tss::KeplerEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_KeplerEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::KeplerEphemerisSettings>>("shared_ptr_KeplerEphemerisSettings");

    // TabulatedEphemerisSettings
    class_<tss::TabulatedEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_TabulatedEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::TabulatedEphemerisSettings>>("shared_ptr_TabulatedEphemerisSettings");

    // ScaledEphemerisSettings
    class_<tss::ScaledEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_ScaledEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::ScaledEphemerisSettings>>("shared_ptr_ScaledEphemerisSettings");

    // ApproximateJplEphemerisSettings
    class_<tss::ApproximateJplEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_ApproximateJplEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::ApproximateJplEphemerisSettings>>("shared_ptr_ApproximateJplEphemerisSettings");

    // InterpolatedSpiceEphemerisSettings
    class_<tss::InterpolatedSpiceEphemerisSettings, base<tss::DirectSpiceEphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_InterpolatedSpiceEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::InterpolatedSpiceEphemerisSettings>>(
            "shared_ptr_InterpolatedSpiceEphemerisSettings");

    // CustomEphemerisSettings
    class_<tss::CustomEphemerisSettings, base<tss::EphemerisSettings>>(
        "dynamics_environment_setup_ephemeris_CustomEphemerisSettings")
        .smart_ptr<std::shared_ptr<tss::CustomEphemerisSettings>>(
            "shared_ptr_CustomEphemerisSettings");

    // Factory functions
    function("dynamics_environment_setup_ephemeris_direct_spice",
        select_overload<std::shared_ptr<tss::EphemerisSettings>(
            const std::string, const std::string, const std::string)>(
            &tss::directSpiceEphemerisSettings));

    function("dynamics_environment_setup_ephemeris_approximate_jpl",
        select_overload<std::shared_ptr<tss::EphemerisSettings>(const std::string)>(
            &tss::approximateJplEphemerisSettings));

    function("dynamics_environment_setup_ephemeris_constant",
        &tss::constantEphemerisSettings);

    function("dynamics_environment_setup_ephemeris_kepler",
        &tss::keplerEphemerisSettings);

    function("dynamics_environment_setup_ephemeris_kepler_from_spice",
        &tss::keplerEphemerisFromSpiceSettings);

    function("dynamics_environment_setup_ephemeris_tabulated",
        select_overload<std::shared_ptr<tss::EphemerisSettings>(
            const std::map<double, Eigen::Vector6d>&, std::string, std::string)>(
            &tss::tabulatedEphemerisSettings));

    function("dynamics_environment_setup_ephemeris_interpolated_spice",
        &tss::interpolatedSpiceEphemerisSettings);

    function("dynamics_environment_setup_ephemeris_scaled_by_constant",
        select_overload<std::shared_ptr<tss::EphemerisSettings>(
            const std::shared_ptr<tss::EphemerisSettings>, const double, const bool)>(
            &tss::scaledEphemerisSettings));

    function("dynamics_environment_setup_ephemeris_scaled_by_vector",
        select_overload<std::shared_ptr<tss::EphemerisSettings>(
            const std::shared_ptr<tss::EphemerisSettings>, const Eigen::Vector6d, const bool)>(
            &tss::scaledEphemerisSettings));

    function("dynamics_environment_setup_ephemeris_custom_ephemeris",
        &tss::customEphemerisSettings);

    // SGP4/TLE ephemeris
    function("dynamics_environment_setup_ephemeris_sgp4",
        &tss::directTleEphemerisSettingsFromTleLines);
}

#endif
