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

#include <tudat/simulation/environment_setup/thrustSettings.h>

namespace tss = tudat::simulation_setup;

// Local ThrustFrames enum for backwards compatibility (matching Python bindings)
namespace tudat {
namespace simulation_setup {
enum ThrustFrames {
    unspecified_thrust_frame = -1,
    inertial_thrust_frame = 0,
    tnw_thrust_frame = 1
};
}
}

WASM_MODULE_PATH("dynamics_propagation_setup_thrust")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_thrust) {
    using namespace emscripten;

    // ========================================================================
    // ThrustMagnitudeTypes enum
    // ========================================================================
    enum_<tss::ThrustMagnitudeTypes>("dynamics_propagation_setup_thrust_ThrustMagnitudeTypes")
        .value("constant_thrust_magnitude", tss::constant_thrust_magnitude)
        .value("thrust_magnitude_from_time_function", tss::thrust_magnitude_from_time_function)
        .value("thrust_magnitude_from_dependent_variables", tss::thrust_magnitude_from_dependent_variables);

    // ========================================================================
    // ThrustDirectionTypes enum (legacy, for backwards compatibility)
    // ========================================================================
    enum_<tss::ThrustDirectionTypes>("dynamics_propagation_setup_thrust_ThrustDirectionGuidanceTypes")
        .value("colinear_with_state_segment_thrust_direction_type", tss::colinear_with_state_segment_thrust_direction)
        .value("thrust_direction_from_existing_body_orientation_type", tss::thrust_direction_from_existing_body_orientation)
        .value("custom_thrust_direction_type", tss::custom_thrust_direction)
        .value("custom_thrust_orientation_type", tss::custom_thrust_orientation)
        .value("mee_costate_based_thrust_direction_type", tss::mee_costate_based_thrust_direction);

    // ========================================================================
    // ThrustFrames enum (legacy, for backwards compatibility)
    // ========================================================================
    enum_<tss::ThrustFrames>("dynamics_propagation_setup_thrust_ThrustFrames")
        .value("unspecified_thrust_frame_type", tss::unspecified_thrust_frame)
        .value("inertial_thrust_frame_type", tss::inertial_thrust_frame)
        .value("tnw_thrust_frame_type", tss::tnw_thrust_frame);

    // ========================================================================
    // ThrustMagnitudeSettings base class and derived classes
    // ========================================================================

    // ThrustMagnitudeSettings base class
    class_<tss::ThrustMagnitudeSettings>("dynamics_propagation_setup_thrust_ThrustMagnitudeSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustMagnitudeSettings>>("shared_ptr_ThrustMagnitudeSettings")
        .property("thrustMagnitudeType", &tss::ThrustMagnitudeSettings::thrustMagnitudeType_)
        .property("thrustOriginId", &tss::ThrustMagnitudeSettings::thrustOriginId_);

    // ConstantThrustMagnitudeSettings derived class
    class_<tss::ConstantThrustMagnitudeSettings, base<tss::ThrustMagnitudeSettings>>(
        "dynamics_propagation_setup_thrust_ConstantThrustMagnitudeSettings")
        .smart_ptr<std::shared_ptr<tss::ConstantThrustMagnitudeSettings>>(
            "shared_ptr_ConstantThrustMagnitudeSettings")
        .property("thrustMagnitude", &tss::ConstantThrustMagnitudeSettings::thrustMagnitude_)
        .property("specificImpulse", &tss::ConstantThrustMagnitudeSettings::specificImpulse_);

    // CustomThrustMagnitudeSettings derived class
    class_<tss::CustomThrustMagnitudeSettings, base<tss::ThrustMagnitudeSettings>>(
        "dynamics_propagation_setup_thrust_CustomThrustMagnitudeSettings")
        .smart_ptr<std::shared_ptr<tss::CustomThrustMagnitudeSettings>>(
            "shared_ptr_CustomThrustMagnitudeSettings")
        .property("specificImpulseIsConstant", &tss::CustomThrustMagnitudeSettings::specificImpulseIsConstant_)
        .property("inputIsForce", &tss::CustomThrustMagnitudeSettings::inputIsForce_);

    // ========================================================================
    // ThrustDirectionSettings base class and derived classes (legacy)
    // ========================================================================

    // ThrustDirectionSettings base class
    class_<tss::ThrustDirectionSettings>("dynamics_propagation_setup_thrust_ThrustDirectionSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustDirectionSettings>>("shared_ptr_ThrustDirectionSettings")
        .property("thrustDirectionType", &tss::ThrustDirectionSettings::thrustDirectionType_)
        .property("relativeBody", &tss::ThrustDirectionSettings::relativeBody_);

    // ThrustDirectionFromStateGuidanceSettings derived class
    class_<tss::ThrustDirectionFromStateGuidanceSettings, base<tss::ThrustDirectionSettings>>(
        "dynamics_propagation_setup_thrust_ThrustDirectionFromStateGuidanceSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustDirectionFromStateGuidanceSettings>>(
            "shared_ptr_ThrustDirectionFromStateGuidanceSettings")
        .property("isColinearWithVelocity", &tss::ThrustDirectionFromStateGuidanceSettings::isColinearWithVelocity_)
        .property("directionIsOppositeToVector", &tss::ThrustDirectionFromStateGuidanceSettings::directionIsOppositeToVector_);

    // CustomThrustDirectionSettings derived class
    class_<tss::CustomThrustDirectionSettings, base<tss::ThrustDirectionSettings>>(
        "dynamics_propagation_setup_thrust_CustomThrustDirectionSettings")
        .smart_ptr<std::shared_ptr<tss::CustomThrustDirectionSettings>>(
            "shared_ptr_CustomThrustDirectionSettings");

    // CustomThrustOrientationSettings derived class
    class_<tss::CustomThrustOrientationSettings, base<tss::ThrustDirectionSettings>>(
        "dynamics_propagation_setup_thrust_CustomThrustOrientationSettings")
        .smart_ptr<std::shared_ptr<tss::CustomThrustOrientationSettings>>(
            "shared_ptr_CustomThrustOrientationSettings");

    // ========================================================================
    // Factory functions - Thrust Magnitude
    // ========================================================================

    // Constant thrust magnitude
    function("dynamics_propagation_setup_thrust_constant_thrust_magnitude",
        &tss::constantThrustMagnitudeSettings);

    // Custom thrust magnitude (force function + Isp function)
    function("dynamics_propagation_setup_thrust_custom_thrust_magnitude",
        &tss::fromFunctionThrustMagnitudeSettings);

    // Custom thrust magnitude with fixed Isp (force function + constant Isp)
    function("dynamics_propagation_setup_thrust_custom_thrust_magnitude_fixed_isp",
        &tss::fromFunctionThrustMagnitudeFixedIspSettings);

    // Custom thrust acceleration magnitude (acceleration function + Isp function)
    function("dynamics_propagation_setup_thrust_custom_thrust_acceleration_magnitude",
        &tss::customThrustAccelerationMagnitudeSettings);

    // Custom thrust acceleration magnitude with fixed Isp
    function("dynamics_propagation_setup_thrust_custom_thrust_acceleration_magnitude_fixed_isp",
        &tss::customThrustAccelerationMagnitudeFixedIspSettings);

    // ========================================================================
    // Factory functions - Thrust Direction (legacy, deprecated)
    // ========================================================================

    // Thrust direction from state guidance
    function("dynamics_propagation_setup_thrust_thrust_direction_from_state_guidance",
        &tss::thrustDirectionFromStateGuidanceSettings);

    // Thrust from existing body orientation (deprecated)
    function("dynamics_propagation_setup_thrust_thrust_from_existing_body_orientation",
        &tss::thrustFromExistingBodyOrientation);

    // Custom thrust orientation (deprecated)
    function("dynamics_propagation_setup_thrust_custom_thrust_orientation",
        select_overload<std::shared_ptr<tss::ThrustDirectionSettings>(
            const std::function<Eigen::Matrix3d(const double)>)>(
            &tss::customThrustOrientationSettings));

    // Custom thrust direction (deprecated)
    function("dynamics_propagation_setup_thrust_custom_thrust_direction",
        &tss::customThrustDirectionSettings);
}

#endif
