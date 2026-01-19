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

#include <tudat/simulation/environment_setup.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup_vehicle_systems")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_vehicle_systems) {
    using namespace emscripten;

    // BodyPanelGeometrySettings base class
    class_<tss::BodyPanelGeometrySettings>(
        "dynamics_environment_setup_vehicle_systems_BodyPanelGeometrySettings")
        .smart_ptr<std::shared_ptr<tss::BodyPanelGeometrySettings>>(
            "shared_ptr_BodyPanelGeometrySettings");

    // FrameFixedBodyPanelGeometrySettings
    class_<tss::FrameFixedBodyPanelGeometrySettings, base<tss::BodyPanelGeometrySettings>>(
        "dynamics_environment_setup_vehicle_systems_FrameFixedBodyPanelGeometrySettings")
        .smart_ptr<std::shared_ptr<tss::FrameFixedBodyPanelGeometrySettings>>(
            "shared_ptr_FrameFixedBodyPanelGeometrySettings");

    // FrameVariableBodyPanelGeometrySettings
    class_<tss::FrameVariableBodyPanelGeometrySettings, base<tss::BodyPanelGeometrySettings>>(
        "dynamics_environment_setup_vehicle_systems_FrameVariableBodyPanelGeometrySettings")
        .smart_ptr<std::shared_ptr<tss::FrameVariableBodyPanelGeometrySettings>>(
            "shared_ptr_FrameVariableBodyPanelGeometrySettings");

    // BodyPanelSettings class
    class_<tss::BodyPanelSettings>(
        "dynamics_environment_setup_vehicle_systems_BodyPanelSettings")
        .smart_ptr<std::shared_ptr<tss::BodyPanelSettings>>(
            "shared_ptr_BodyPanelSettings");

    // FullPanelledBodySettings class
    class_<tss::FullPanelledBodySettings>(
        "dynamics_environment_setup_vehicle_systems_FullPanelledBodySettings")
        .smart_ptr<std::shared_ptr<tss::FullPanelledBodySettings>>(
            "shared_ptr_FullPanelledBodySettings");

    // MaterialProperties class
    class_<tss::MaterialProperties>(
        "dynamics_environment_setup_vehicle_systems_MaterialProperties")
        .smart_ptr<std::shared_ptr<tss::MaterialProperties>>(
            "shared_ptr_MaterialProperties");

    // Factory functions
    function("dynamics_environment_setup_vehicle_systems_frame_fixed_panel_geometry",
        select_overload<std::shared_ptr<tss::BodyPanelGeometrySettings>(
            const Eigen::Vector3d&, const double, const std::string&)>(
                &tss::frameFixedPanelGeometry));

    function("dynamics_environment_setup_vehicle_systems_body_tracking_panel_geometry",
        select_overload<std::shared_ptr<tss::BodyPanelGeometrySettings>(
            const std::string&, const bool, const double, const std::string&)>(
                &tss::bodyTrackingPanelGeometry));

    function("dynamics_environment_setup_vehicle_systems_body_panel_settings",
        select_overload<std::shared_ptr<tss::BodyPanelSettings>(
            std::shared_ptr<tss::BodyPanelGeometrySettings>,
            std::shared_ptr<tss::BodyPanelReflectionLawSettings>,
            std::string, std::shared_ptr<tss::MaterialProperties>)>(
                &tss::bodyPanelSettings));

    function("dynamics_environment_setup_vehicle_systems_full_panelled_body_settings",
        &tss::fullPanelledBodySettings);

    function("dynamics_environment_setup_vehicle_systems_box_wing_panelled_body_settings",
        &tss::bodyWingPanelledGeometry);

    function("dynamics_environment_setup_vehicle_systems_material_properties",
        &tss::materialProperties);
}

#endif
