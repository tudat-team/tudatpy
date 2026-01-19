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
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/createRotationModel.h>
#include <tudat/simulation/environment_setup/defaultBodies.h>

namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("dynamics_environment_setup_rotation_model")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_rotation_model) {
    using namespace emscripten;

    // RotationModelType enum
    enum_<tss::RotationModelType>("dynamics_environment_setup_rotation_model_RotationModelType")
        .value("simple_rotational_model", tss::simple_rotation_model)
        .value("spice_rotation_model", tss::spice_rotation_model)
        .value("gcrs_to_itrs_rotation_model", tss::gcrs_to_itrs_rotation_model)
        .value("synchronous_rotation_model", tss::synchronous_rotation_model)
        .value("planetary_rotation_model", tss::planetary_rotation_model);

    // IAUConventions enum
    enum_<tba::IAUConventions>("dynamics_environment_setup_rotation_model_IAUConventions")
        .value("iau_2000_a", tba::iau_2000_a)
        .value("iau_2000_b", tba::iau_2000_b)
        .value("iau_2006", tba::iau_2006);

    // RotationModelSettings base class
    class_<tss::RotationModelSettings>("dynamics_environment_setup_rotation_model_RotationModelSettings")
        .smart_ptr<std::shared_ptr<tss::RotationModelSettings>>("shared_ptr_RotationModelSettings")
        .function("getRotationType", &tss::RotationModelSettings::getRotationType)
        .function("getOriginalFrame", &tss::RotationModelSettings::getOriginalFrame)
        .function("getTargetFrame", &tss::RotationModelSettings::getTargetFrame);

    // SimpleRotationModelSettings derived class
    class_<tss::SimpleRotationModelSettings, base<tss::RotationModelSettings>>(
        "dynamics_environment_setup_rotation_model_SimpleRotationModelSettings")
        .smart_ptr<std::shared_ptr<tss::SimpleRotationModelSettings>>(
            "shared_ptr_SimpleRotationModelSettings");

    // GcrsToItrsRotationModelSettings derived class
    class_<tss::GcrsToItrsRotationModelSettings, base<tss::RotationModelSettings>>(
        "dynamics_environment_setup_rotation_model_GcrsToItrsRotationModelSettings")
        .smart_ptr<std::shared_ptr<tss::GcrsToItrsRotationModelSettings>>(
            "shared_ptr_GcrsToItrsRotationModelSettings");

    // PlanetaryRotationModelSettings derived class
    class_<tss::PlanetaryRotationModelSettings, base<tss::RotationModelSettings>>(
        "dynamics_environment_setup_rotation_model_PlanetaryRotationModelSettings")
        .smart_ptr<std::shared_ptr<tss::PlanetaryRotationModelSettings>>(
            "shared_ptr_PlanetaryRotationModelSettings");

    // Factory functions
    function("dynamics_environment_setup_rotation_model_simple",
        select_overload<std::shared_ptr<tss::RotationModelSettings>(
            const std::string&, const std::string&, const Eigen::Matrix3d&,
            const double, const double)>(&tss::simpleRotationModelSettings));

    function("dynamics_environment_setup_rotation_model_simple_from_spice",
        &tss::simpleRotationModelFromSpiceSettings);

    function("dynamics_environment_setup_rotation_model_synchronous",
        &tss::synchronousRotationModelSettings);

    function("dynamics_environment_setup_rotation_model_spice",
        &tss::spiceRotationModelSettings);

    function("dynamics_environment_setup_rotation_model_constant",
        select_overload<std::shared_ptr<tss::RotationModelSettings>(
            const std::string&, const std::string&, const Eigen::Matrix3d&)>(
                &tss::constantRotationModelSettings));

    // GCRS to ITRS (high-accuracy Earth rotation)
    function("dynamics_environment_setup_rotation_model_gcrs_to_itrs",
        &tss::gcrsToItrsRotationModelSettings);

    // Mars high-accuracy rotation model
    function("dynamics_environment_setup_rotation_model_mars_high_accuracy",
        select_overload<std::shared_ptr<tss::RotationModelSettings>(
            const std::string&, const std::string&)>(
                &tss::getHighAccuracyMarsRotationModel));

    // Mars high-accuracy with custom angles
    function("dynamics_environment_setup_rotation_model_mars_high_accuracy_custom_angles",
        select_overload<std::shared_ptr<tss::RotationModelSettings>(
            const std::string&, const std::string&, const double, const double,
            const double, const double)>(
                &tss::getHighAccuracyMarsRotationModel));

    // Aerodynamic angle-based rotation
    function("dynamics_environment_setup_rotation_model_aerodynamic_angle_based",
        &tss::aerodynamicAngleRotationSettings);

    // Zero pitch moment aerodynamic angle-based
    function("dynamics_environment_setup_rotation_model_zero_pitch_moment_aerodynamic_angle_based",
        &tss::pitchTrimRotationSettings);

    // Custom inertial direction-based rotation
    function("dynamics_environment_setup_rotation_model_custom_inertial_direction_based",
        &tss::bodyFixedDirectionBasedRotationSettings);

    // Orbital state direction-based rotation
    function("dynamics_environment_setup_rotation_model_orbital_state_direction_based",
        &tss::orbitalStateBasedRotationSettings);
}

#endif
