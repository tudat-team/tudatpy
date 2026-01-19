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

#include <tudat/astro/mission_segments/transferLeg.h>
#include <tudat/astro/mission_segments/transferNode.h>
#include <tudat/astro/mission_segments/transferTrajectory.h>
#include <tudat/astro/mission_segments/createTransferTrajectory.h>
#include <tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h>
#include <tudat/astro/low_thrust/shape_based/sphericalShapingLeg.h>
#include <tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h>
#include <tudat/math/root_finders.h>

namespace tms = tudat::mission_segments;
namespace tss = tudat::simulation_setup;
namespace tsbm = tudat::shape_based_methods;
namespace trf = tudat::root_finders;

WASM_MODULE_PATH("trajectory_design_transfer_trajectory")

EMSCRIPTEN_BINDINGS(tudatpy_trajectory_design_transfer_trajectory) {
    using namespace emscripten;

    // Default minimum pericenters constant
    constant("trajectory_design_transfer_trajectory_DEFAULT_MINIMUM_PERICENTERS",
        tms::DEFAULT_MINIMUM_PERICENTERS);

    // TransferLegTypes enum
    enum_<tms::TransferLegTypes>("trajectory_design_transfer_trajectory_TransferLegTypes")
        .value("unpowered_unperturbed_leg_type", tms::unpowered_unperturbed_leg)
        .value("dsm_position_based_leg_type", tms::dsm_position_based_leg)
        .value("dsm_velocity_based_leg_type", tms::dsm_velocity_based_leg)
        .value("hodographic_low_thrust_leg", tms::hodographic_low_thrust_leg)
        .value("spherical_shaping_low_thrust_leg", tms::spherical_shaping_low_thrust_leg);

    // TransferNodeTypes enum
    enum_<tms::TransferNodeTypes>("trajectory_design_transfer_trajectory_TransferNodeTypes")
        .value("swingby", tms::swingby)
        .value("escape_and_departure", tms::escape_and_departure)
        .value("capture_and_insertion", tms::capture_and_insertion);

    // TransferLeg base class
    class_<tms::TransferLeg>("trajectory_design_transfer_trajectory_TransferLeg")
        .smart_ptr<std::shared_ptr<tms::TransferLeg>>("shared_ptr_TransferLeg")
        .function("updateLegParameters", &tms::TransferLeg::updateLegParameters)
        .function("getLegDeltaV", &tms::TransferLeg::getLegDeltaV)
        .function("getTransferLegType", &tms::TransferLeg::getTransferLegType)
        .function("getDepartureVelocity", &tms::TransferLeg::getDepartureVelocity)
        .function("getArrivalVelocity", &tms::TransferLeg::getArrivalVelocity)
        .function("getLegTimeOfFlight", &tms::TransferLeg::getLegTimeOfFlight)
        .function("getLegDepartureTime", &tms::TransferLeg::getLegDepartureTime)
        .function("getLegArrivalTime", &tms::TransferLeg::getLegArrivalTime)
        .function("state_along_trajectory",
            select_overload<Eigen::Vector6d(const double)>(&tms::TransferLeg::getStateAlongTrajectory));

    // SphericalShapingLeg derived class
    class_<tsbm::SphericalShapingLeg, base<tms::TransferLeg>>(
        "trajectory_design_transfer_trajectory_SphericalShapingLeg")
        .smart_ptr<std::shared_ptr<tsbm::SphericalShapingLeg>>("shared_ptr_SphericalShapingLeg");

    // HodographicShapingLeg derived class
    class_<tsbm::HodographicShapingLeg, base<tms::TransferLeg>>(
        "trajectory_design_transfer_trajectory_HodographicShapingLeg")
        .smart_ptr<std::shared_ptr<tsbm::HodographicShapingLeg>>("shared_ptr_HodographicShapingLeg");

    // TransferNode base class
    class_<tms::TransferNode>("trajectory_design_transfer_trajectory_TransferNode")
        .smart_ptr<std::shared_ptr<tms::TransferNode>>("shared_ptr_TransferNode")
        .function("updateNodeParameters", &tms::TransferNode::updateNodeParameters)
        .function("getNodeDeltaV", &tms::TransferNode::getNodeDeltaV)
        .function("getTransferNodeType", &tms::TransferNode::getTransferNodeType)
        .function("nodeComputesOutgoingVelocity", &tms::TransferNode::nodeComputesOutgoingVelocity)
        .function("nodeComputesIncomingVelocity", &tms::TransferNode::nodeComputesIncomingVelocity)
        .function("getIncomingVelocity", &tms::TransferNode::getIncomingVelocity)
        .function("getOutgoingVelocity", &tms::TransferNode::getOutgoingVelocity);

    // TransferTrajectory class
    class_<tms::TransferTrajectory>("trajectory_design_transfer_trajectory_TransferTrajectory")
        .smart_ptr<std::shared_ptr<tms::TransferTrajectory>>("shared_ptr_TransferTrajectory")
        .function("evaluateTrajectory", &tms::TransferTrajectory::evaluateTrajectory)
        .function("getTotalDeltaV", &tms::TransferTrajectory::getTotalDeltaV)
        .function("getNodeDeltaV", &tms::TransferTrajectory::getNodeDeltaV)
        .function("getLegDeltaV", &tms::TransferTrajectory::getLegDeltaV)
        .function("getDeltaVPerNode", &tms::TransferTrajectory::getDeltaVPerNode)
        .function("getDeltaVPerLeg", &tms::TransferTrajectory::getDeltaVPerLeg)
        .function("getTotalTimeOfFlight", &tms::TransferTrajectory::getTotalTimeOfFlight)
        .function("getNumberOfNodes", &tms::TransferTrajectory::getNumberOfNodes)
        .function("getNumberOfLegs", &tms::TransferTrajectory::getNumberOfLegs)
        .function("getLegs", &tms::TransferTrajectory::getLegs)
        .function("getNodes", &tms::TransferTrajectory::getNodes)
        .function("getStatesAlongTrajectory",
            select_overload<std::map<double, Eigen::Vector6d>(const int)>(
                &tms::TransferTrajectory::getStatesAlongTrajectory))
        .function("getInertialThrustAccelerationsAlongTrajectory",
            select_overload<std::map<double, Eigen::Vector3d>(const int)>(
                &tms::TransferTrajectory::getInertialThrustAccelerationsAlongTrajectory))
        .function("getRswThrustAccelerationsAlongTrajectory",
            &tms::TransferTrajectory::getRswThrustAccelerationsAlongTrajectory)
        .function("getTnwThrustAccelerationsAlongTrajectory",
            &tms::TransferTrajectory::getTnwThrustAccelerationsAlongTrajectory);

    // TransferLegSettings base class
    class_<tms::TransferLegSettings>("trajectory_design_transfer_trajectory_TransferLegSettings")
        .smart_ptr<std::shared_ptr<tms::TransferLegSettings>>("shared_ptr_TransferLegSettings");

    // TransferNodeSettings base class
    class_<tms::TransferNodeSettings>("trajectory_design_transfer_trajectory_TransferNodeSettings")
        .smart_ptr<std::shared_ptr<tms::TransferNodeSettings>>("shared_ptr_TransferNodeSettings");

    // SwingbyNodeSettings
    class_<tms::SwingbyNodeSettings, base<tms::TransferNodeSettings>>(
        "trajectory_design_transfer_trajectory_SwingbyNodeSettings")
        .smart_ptr<std::shared_ptr<tms::SwingbyNodeSettings>>("shared_ptr_SwingbyNodeSettings");

    // EscapeAndDepartureNodeSettings
    class_<tms::EscapeAndDepartureNodeSettings, base<tms::TransferNodeSettings>>(
        "trajectory_design_transfer_trajectory_EscapeAndDepartureNodeSettings")
        .smart_ptr<std::shared_ptr<tms::EscapeAndDepartureNodeSettings>>(
            "shared_ptr_EscapeAndDepartureNodeSettings");

    // CaptureAndInsertionNodeSettings
    class_<tms::CaptureAndInsertionNodeSettings, base<tms::TransferNodeSettings>>(
        "trajectory_design_transfer_trajectory_CaptureAndInsertionNodeSettings")
        .smart_ptr<std::shared_ptr<tms::CaptureAndInsertionNodeSettings>>(
            "shared_ptr_CaptureAndInsertionNodeSettings");

    // Factory functions for leg settings
    function("trajectory_design_transfer_trajectory_unpowered_leg",
        &tms::unpoweredLeg);

    function("trajectory_design_transfer_trajectory_dsm_position_based_leg",
        &tms::dsmPositionBasedLeg);

    function("trajectory_design_transfer_trajectory_dsm_velocity_based_leg",
        &tms::dsmVelocityBasedLeg);

    // Spherical shaping leg
    function("trajectory_design_transfer_trajectory_spherical_shaping_leg",
        &tms::sphericalShapingLeg);

    // Hodographic shaping leg
    function("trajectory_design_transfer_trajectory_hodographic_shaping_leg",
        &tms::hodographicShapingLeg);

    // Factory functions for node settings
    function("trajectory_design_transfer_trajectory_swingby_node",
        &tms::swingbyNode);

    function("trajectory_design_transfer_trajectory_escape_and_departure_node",
        &tms::escapeAndDepartureNode);

    function("trajectory_design_transfer_trajectory_capture_and_insertion_node",
        &tms::captureAndInsertionNode);

    // Main trajectory creation function
    function("trajectory_design_transfer_trajectory_create_transfer_trajectory",
        &tms::createTransferTrajectory);

    // MGA trajectory settings functions
    function("trajectory_design_transfer_trajectory_mga_settings_unpowered_unperturbed_legs",
        select_overload<std::pair<std::vector<std::shared_ptr<tms::TransferLegSettings>>,
                                  std::vector<std::shared_ptr<tms::TransferNodeSettings>>>(
            const std::vector<std::string>&,
            const std::pair<double, double>,
            const std::pair<double, double>,
            const std::map<std::string, double>)>(&tms::getMgaTransferTrajectorySettingsWithoutDsm));

    function("trajectory_design_transfer_trajectory_mga_settings_dsm_position_based_legs",
        select_overload<std::pair<std::vector<std::shared_ptr<tms::TransferLegSettings>>,
                                  std::vector<std::shared_ptr<tms::TransferNodeSettings>>>(
            const std::vector<std::string>&,
            const std::pair<double, double>,
            const std::pair<double, double>,
            const std::map<std::string, double>)>(&tms::getMgaTransferTrajectorySettingsWithPositionBasedDsm));

    function("trajectory_design_transfer_trajectory_mga_settings_dsm_velocity_based_legs",
        select_overload<std::pair<std::vector<std::shared_ptr<tms::TransferLegSettings>>,
                                  std::vector<std::shared_ptr<tms::TransferNodeSettings>>>(
            const std::vector<std::string>&,
            const std::pair<double, double>,
            const std::pair<double, double>,
            const std::map<std::string, double>)>(&tms::getMgaTransferTrajectorySettingsWithVelocityBasedDsm));

    // MGA settings with spherical shaping legs
    function("trajectory_design_transfer_trajectory_mga_settings_spherical_shaping_legs",
        select_overload<std::pair<std::vector<std::shared_ptr<tms::TransferLegSettings>>,
                                  std::vector<std::shared_ptr<tms::TransferNodeSettings>>>(
            const std::vector<std::string>&,
            const std::shared_ptr<trf::RootFinderSettings>,
            const std::pair<double, double>,
            const std::pair<double, double>,
            const double,
            const double,
            const double,
            const std::map<std::string, double>)>(&tms::getMgaTransferTrajectorySettingsWithSphericalShapingThrust));

    // MGA settings with hodographic shaping legs (user-provided functions)
    function("trajectory_design_transfer_trajectory_mga_settings_hodographic_shaping_legs",
        select_overload<std::pair<std::vector<std::shared_ptr<tms::TransferLegSettings>>,
                                  std::vector<std::shared_ptr<tms::TransferNodeSettings>>>(
            const std::vector<std::string>&,
            const std::vector<tsbm::HodographicBasisFunctionList>&,
            const std::vector<tsbm::HodographicBasisFunctionList>&,
            const std::vector<tsbm::HodographicBasisFunctionList>&,
            const std::pair<double, double>,
            const std::pair<double, double>,
            const std::map<std::string, double>)>(&tms::getMgaTransferTrajectorySettingsWithHodographicShapingThrust));

    // MGA settings with hodographic shaping legs (recommended functions)
    function("trajectory_design_transfer_trajectory_mga_settings_hodographic_shaping_legs_with_recommended_functions",
        select_overload<std::pair<std::vector<std::shared_ptr<tms::TransferLegSettings>>,
                                  std::vector<std::shared_ptr<tms::TransferNodeSettings>>>(
            const std::vector<std::string>&,
            const std::vector<double>&,
            const std::vector<double>&,
            const std::pair<double, double>,
            const std::pair<double, double>,
            const std::map<std::string, double>)>(&tms::getMgaTransferTrajectorySettingsWithHodographicShapingThrust));

    // Utility function
    function("trajectory_design_transfer_trajectory_print_parameter_definitions",
        &tms::printTransferParameterDefinition);

    // Set low thrust acceleration
    function("trajectory_design_transfer_trajectory_set_low_thrust_acceleration",
        &tms::setLowThrustAcceleration);
}

#endif
