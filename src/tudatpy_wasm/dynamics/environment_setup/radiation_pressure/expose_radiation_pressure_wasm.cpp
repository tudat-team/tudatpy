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

#include <tudat/simulation/environment_setup/createRadiationPressureInterface.h>
#include <tudat/simulation/environment_setup/createRadiationSourceModel.h>
#include <tudat/simulation/environment_setup/createRadiationPressureTargetModel.h>
#include <tudat/simulation/environment_setup/createSystemModel.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup_radiation_pressure")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_radiation_pressure) {
    using namespace emscripten;

    // RadiationPressureType enum (deprecated but still needed)
    enum_<tss::RadiationPressureType>("dynamics_environment_setup_radiation_pressure_RadiationPressureType")
        .value("cannon_ball_radiation_pressure_interface", tss::cannon_ball_radiation_pressure_interface);

    // RadiationPressureTargetModelType enum
    enum_<tss::RadiationPressureTargetModelType>("dynamics_environment_setup_radiation_pressure_RadiationPressureTargetModelType")
        .value("cannonball_target", tss::cannonball_target)
        .value("paneled_target", tss::paneled_target)
        .value("multi_type_target", tss::multi_type_target)
        .value("undefined_target", tss::undefined_target);

    // RadiationPressureInterfaceSettings base class (deprecated)
    class_<tss::RadiationPressureInterfaceSettings>(
        "dynamics_environment_setup_radiation_pressure_RadiationPressureInterfaceSettings")
        .smart_ptr<std::shared_ptr<tss::RadiationPressureInterfaceSettings>>(
            "shared_ptr_RadiationPressureInterfaceSettings")
        .function("getRadiationPressureType", &tss::RadiationPressureInterfaceSettings::getRadiationPressureType)
        .function("getSourceBody", &tss::RadiationPressureInterfaceSettings::getSourceBody)
        .function("getOccultingBodies", &tss::RadiationPressureInterfaceSettings::getOccultingBodies);

    // CannonBallRadiationPressureInterfaceSettings (deprecated)
    class_<tss::CannonBallRadiationPressureInterfaceSettings, base<tss::RadiationPressureInterfaceSettings>>(
        "dynamics_environment_setup_radiation_pressure_CannonBallRadiationPressureInterfaceSettings")
        .smart_ptr<std::shared_ptr<tss::CannonBallRadiationPressureInterfaceSettings>>(
            "shared_ptr_CannonBallRadiationPressureInterfaceSettings")
        .function("getArea", &tss::CannonBallRadiationPressureInterfaceSettings::getArea)
        .function("getRadiationPressureCoefficient",
            &tss::CannonBallRadiationPressureInterfaceSettings::getRadiationPressureCoefficient);

    // ========================================================================
    // Luminosity model settings
    // ========================================================================
    class_<tss::LuminosityModelSettings>(
        "dynamics_environment_setup_radiation_pressure_LuminosityModelSettings")
        .smart_ptr<std::shared_ptr<tss::LuminosityModelSettings>>(
            "shared_ptr_LuminosityModelSettings");

    function("dynamics_environment_setup_radiation_pressure_constant_luminosity",
        &tss::constantLuminosityModelSettings);

    function("dynamics_environment_setup_radiation_pressure_irradiance_based_constant_luminosity",
        &tss::irradianceBasedLuminosityModelSettings);

    // ========================================================================
    // Radiation source model settings
    // ========================================================================
    class_<tss::RadiationSourceModelSettings>(
        "dynamics_environment_setup_radiation_pressure_RadiationSourceModelSettings")
        .smart_ptr<std::shared_ptr<tss::RadiationSourceModelSettings>>(
            "shared_ptr_RadiationSourceModelSettings");

    function("dynamics_environment_setup_radiation_pressure_isotropic_radiation_source",
        &tss::isotropicPointRadiationSourceModelSettings);

    // ========================================================================
    // Radiation target model settings
    // ========================================================================
    class_<tss::RadiationPressureTargetModelSettings>(
        "dynamics_environment_setup_radiation_pressure_RadiationPressureTargetModelSettings")
        .smart_ptr<std::shared_ptr<tss::RadiationPressureTargetModelSettings>>(
            "shared_ptr_RadiationPressureTargetModelSettings");

    function("dynamics_environment_setup_radiation_pressure_cannonball_radiation_target",
        &tss::cannonballRadiationPressureTargetModelSettings);

    function("dynamics_environment_setup_radiation_pressure_paneled_radiation_target",
        select_overload<std::shared_ptr<tss::RadiationPressureTargetModelSettings>(
            const std::vector<std::string>&)>(
            &tss::paneledRadiationPressureTargetModelSettings));

    // ========================================================================
    // Panel reflection law settings
    // ========================================================================
    class_<tss::BodyPanelReflectionLawSettings>(
        "dynamics_environment_setup_radiation_pressure_BodyPanelReflectionLawSettings")
        .smart_ptr<std::shared_ptr<tss::BodyPanelReflectionLawSettings>>(
            "shared_ptr_BodyPanelReflectionLawSettings");

    function("dynamics_environment_setup_radiation_pressure_specular_diffuse_body_panel_reflection",
        &tss::specularDiffuseBodyPanelReflectionLawSettings);

    function("dynamics_environment_setup_radiation_pressure_lambertian_body_panel_reflection",
        &tss::lambertainBodyPanelReflectionLawSettings);

    // Deprecated cannonball function
    function("dynamics_environment_setup_radiation_pressure_cannonball",
        select_overload<std::shared_ptr<tss::RadiationPressureInterfaceSettings>(
            const std::string&, const double, const double, const std::vector<std::string>&)>(
            &tss::cannonBallRadiationPressureSettings));
}

#endif
