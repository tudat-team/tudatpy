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

WASM_MODULE_PATH("dynamics_propagation_setup_thrust")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_thrust) {
    using namespace emscripten;

    // ThrustDirectionTypes enum
    enum_<tss::ThrustDirectionTypes>("dynamics_propagation_setup_thrust_ThrustDirectionTypes")
        .value("colinear_with_state_segment_thrust_direction", tss::colinear_with_state_segment_thrust_direction)
        .value("thrust_direction_from_existing_body_orientation", tss::thrust_direction_from_existing_body_orientation)
        .value("custom_thrust_direction", tss::custom_thrust_direction)
        .value("custom_thrust_orientation", tss::custom_thrust_orientation)
        .value("mee_costate_based_thrust_direction", tss::mee_costate_based_thrust_direction);

    // ThrustMagnitudeTypes enum
    enum_<tss::ThrustMagnitudeTypes>("dynamics_propagation_setup_thrust_ThrustMagnitudeTypes")
        .value("constant_thrust_magnitude", tss::constant_thrust_magnitude)
        .value("thrust_magnitude_from_time_function", tss::thrust_magnitude_from_time_function)
        .value("thrust_magnitude_from_dependent_variables", tss::thrust_magnitude_from_dependent_variables);

    // ThrustDirectionSettings base class
    class_<tss::ThrustDirectionSettings>("dynamics_propagation_setup_thrust_ThrustDirectionSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustDirectionSettings>>("shared_ptr_ThrustDirectionSettings")
        .property("thrustDirectionType", &tss::ThrustDirectionSettings::thrustDirectionType_);

    // ThrustDirectionFromStateGuidanceSettings
    class_<tss::ThrustDirectionFromStateGuidanceSettings, base<tss::ThrustDirectionSettings>>(
        "dynamics_propagation_setup_thrust_ThrustDirectionFromStateGuidanceSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustDirectionFromStateGuidanceSettings>>(
            "shared_ptr_ThrustDirectionFromStateGuidanceSettings");

    // ThrustMagnitudeSettings base class
    class_<tss::ThrustMagnitudeSettings>("dynamics_propagation_setup_thrust_ThrustMagnitudeSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustMagnitudeSettings>>("shared_ptr_ThrustMagnitudeSettings")
        .property("thrustMagnitudeType", &tss::ThrustMagnitudeSettings::thrustMagnitudeType_);

    // ConstantThrustMagnitudeSettings
    class_<tss::ConstantThrustMagnitudeSettings, base<tss::ThrustMagnitudeSettings>>(
        "dynamics_propagation_setup_thrust_ConstantThrustMagnitudeSettings")
        .smart_ptr<std::shared_ptr<tss::ConstantThrustMagnitudeSettings>>(
            "shared_ptr_ConstantThrustMagnitudeSettings")
        .property("thrustMagnitude", &tss::ConstantThrustMagnitudeSettings::thrustMagnitude_)
        .property("specificImpulse", &tss::ConstantThrustMagnitudeSettings::specificImpulse_);

    // Factory functions
    function("dynamics_propagation_setup_thrust_thrust_direction_from_state_guidance",
        &tss::thrustDirectionFromStateGuidanceSettings);

    function("dynamics_propagation_setup_thrust_constant_thrust_magnitude",
        &tss::constantThrustMagnitudeSettings);

    function("dynamics_propagation_setup_thrust_from_function_thrust_magnitude",
        &tss::fromFunctionThrustMagnitudeSettings);
}

#endif
