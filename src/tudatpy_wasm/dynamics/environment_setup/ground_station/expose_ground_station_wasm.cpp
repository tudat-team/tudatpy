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

#include <tudat/simulation/environment_setup/createGroundStations.h>

namespace tss = tudat::simulation_setup;
namespace tcc = tudat::coordinate_conversions;

WASM_MODULE_PATH("dynamics_environment_setup_ground_station")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_ground_station) {
    using namespace emscripten;

    // PositionElementTypes enum
    enum_<tcc::PositionElementTypes>("dynamics_environment_setup_ground_station_PositionElementTypes")
        .value("cartesian_position", tcc::cartesian_position)
        .value("spherical_position", tcc::spherical_position)
        .value("geodetic_position", tcc::geodetic_position);

    // StationMotionModelTypes enum
    enum_<tss::StationMotionModelTypes>("dynamics_environment_setup_ground_station_StationMotionModelTypes")
        .value("linear_station_motion", tss::linear_station_motion)
        .value("piecewise_constant_station_motion", tss::piecewise_constant_station_motion)
        .value("custom_station_motion", tss::custom_station_motion)
        .value("body_deformation_station_motion", tss::body_deformation_station_motion);

    // GroundStationMotionSettings base class
    class_<tss::GroundStationMotionSettings>("dynamics_environment_setup_ground_station_GroundStationMotionSettings")
        .smart_ptr<std::shared_ptr<tss::GroundStationMotionSettings>>("shared_ptr_GroundStationMotionSettings")
        .function("getModelType", &tss::GroundStationMotionSettings::getModelType);

    // LinearGroundStationMotionSettings
    class_<tss::LinearGroundStationMotionSettings, base<tss::GroundStationMotionSettings>>(
        "dynamics_environment_setup_ground_station_LinearGroundStationMotionSettings")
        .smart_ptr<std::shared_ptr<tss::LinearGroundStationMotionSettings>>(
            "shared_ptr_LinearGroundStationMotionSettings");

    // GroundStationSettings class
    class_<tss::GroundStationSettings>("dynamics_environment_setup_ground_station_GroundStationSettings")
        .smart_ptr<std::shared_ptr<tss::GroundStationSettings>>("shared_ptr_GroundStationSettings")
        .function("getStationName", &tss::GroundStationSettings::getStationName)
        .function("getGroundStationPosition", &tss::GroundStationSettings::getGroundStationPosition)
        .function("getPositionElementType", &tss::GroundStationSettings::getPositionElementType);

    // PiecewiseConstantGroundStationMotionSettings
    class_<tss::PiecewiseConstantGroundStationMotionSettings, base<tss::GroundStationMotionSettings>>(
        "dynamics_environment_setup_ground_station_PiecewiseConstantGroundStationMotionSettings")
        .smart_ptr<std::shared_ptr<tss::PiecewiseConstantGroundStationMotionSettings>>(
            "shared_ptr_PiecewiseConstantGroundStationMotionSettings");

    // CustomGroundStationMotionSettings
    class_<tss::CustomGroundStationMotionSettings, base<tss::GroundStationMotionSettings>>(
        "dynamics_environment_setup_ground_station_CustomGroundStationMotionSettings")
        .smart_ptr<std::shared_ptr<tss::CustomGroundStationMotionSettings>>(
            "shared_ptr_CustomGroundStationMotionSettings");

    // ========================================================================
    // Factory functions
    // ========================================================================

    function("dynamics_environment_setup_ground_station_basic_station",
        &tss::groundStationSettings);

    function("dynamics_environment_setup_ground_station_linear_motion",
        &tss::linearGroundStationMotionSettings);

    function("dynamics_environment_setup_ground_station_piecewise_constant_station_motion",
        &tss::piecewiseConstantGroundStationMotionSettings);

    function("dynamics_environment_setup_ground_station_custom_station_motion",
        &tss::customGroundStationMotionSettings);

    function("dynamics_environment_setup_ground_station_body_deformation_motion",
        &tss::bodyDeformationStationMotionSettings);

}

#endif
