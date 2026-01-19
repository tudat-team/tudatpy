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
#include "../../eigen_wasm.h"
#include "../../stl_wasm.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup")

// Helper wrapper functions to avoid issues with default arguments
tss::BodyListSettings getDefaultBodySettingsWrapper(
    const std::vector<std::string>& bodies,
    const std::string& baseFrameOrigin,
    const std::string& baseFrameOrientation)
{
    return tss::getDefaultBodySettings(bodies, baseFrameOrigin, baseFrameOrientation);
}

tss::BodyListSettings getDefaultBodySettingsTimeLimitedWrapper(
    const std::vector<std::string>& bodies,
    const double initialTime,
    const double finalTime,
    const std::string& baseFrameOrigin,
    const std::string& baseFrameOrientation,
    const double timeStep)
{
    return tss::getDefaultBodySettings(bodies, initialTime, finalTime, baseFrameOrigin, baseFrameOrientation, timeStep);
}

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup) {
    using namespace emscripten;

    // BodySettings class
    class_<tss::BodySettings>("dynamics_environment_setup_BodySettings")
        .smart_ptr<std::shared_ptr<tss::BodySettings>>("shared_ptr_BodySettings")
        .property("constantMass", &tss::BodySettings::constantMass)
        .property("ephemerisSettings", &tss::BodySettings::ephemerisSettings)
        .property("gravityFieldSettings", &tss::BodySettings::gravityFieldSettings)
        .property("rotationModelSettings", &tss::BodySettings::rotationModelSettings)
        .property("atmosphereSettings", &tss::BodySettings::atmosphereSettings)
        .property("shapeModelSettings", &tss::BodySettings::shapeModelSettings)
        .property("radiationPressureSettings", &tss::BodySettings::radiationPressureSettings)
        .property("aerodynamicCoefficientSettings", &tss::BodySettings::aerodynamicCoefficientSettings)
        .property("gravityFieldVariationSettings", &tss::BodySettings::gravityFieldVariationSettings)
        .property("groundStationSettings", &tss::BodySettings::groundStationSettings);

    // BodyListSettings typedef
    class_<tss::BodyListSettings>("dynamics_environment_setup_BodyListSettings")
        .smart_ptr<std::shared_ptr<tss::BodyListSettings>>("shared_ptr_BodyListSettings")
        .function("get", &tss::BodyListSettings::get)
        .function("getFrameOrigin", &tss::BodyListSettings::getFrameOrigin)
        .function("getFrameOrientation", &tss::BodyListSettings::getFrameOrientation);

    // Factory function for default body settings
    function("dynamics_environment_setup_get_default_body_settings",
        &getDefaultBodySettingsWrapper);

    function("dynamics_environment_setup_get_default_body_settings_time_limited",
        &getDefaultBodySettingsTimeLimitedWrapper);

    // SystemOfBodies creation
    function("dynamics_environment_setup_create_system_of_bodies",
        &tss::createSystemOfBodies<double, double>);

    // Get default single body settings
    function("dynamics_environment_setup_get_default_single_body_settings",
        select_overload<std::shared_ptr<tss::BodySettings>(
            const std::string&, const std::string&)>(
            &tss::getDefaultSingleBodySettings));

    // Get default single body settings time-limited
    function("dynamics_environment_setup_get_default_single_body_settings_time_limited",
        select_overload<std::shared_ptr<tss::BodySettings>(
            const std::string&, const double, const double,
            const std::string&, const double)>(
            &tss::getDefaultSingleBodySettings));

    // Add aerodynamic coefficient interface
    function("dynamics_environment_setup_add_aerodynamic_coefficient_interface",
        &tss::addAerodynamicCoefficientInterface);

    // Add radiation pressure target model
    function("dynamics_environment_setup_add_radiation_pressure_target_model",
        &tss::addRadiationPressureTargetModel);

    // Add rotation model
    function("dynamics_environment_setup_add_rotation_model",
        &tss::addRotationModel);

    // Add gravity field model
    function("dynamics_environment_setup_add_gravity_field_model",
        &tss::addGravityFieldModel);

    // Add rigid body properties
    function("dynamics_environment_setup_add_rigid_body_properties",
        &tss::addRigidBodyProperties);

    // Add engine model
    function("dynamics_environment_setup_add_engine_model",
        &tss::addEngineModel);

    // Add flight conditions
    function("dynamics_environment_setup_add_flight_conditions",
        &tss::addFlightConditions);
}

#endif
