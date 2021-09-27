/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/mission_segments/createTransferTrajectory.h"

#include <pybind11/eigen.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <tudat/astro/mission_segments/createTransferTrajectory.h>
#include <tudat/simulation/propagation_setup/accelerationSettings.h>

namespace py = pybind11;
namespace tms = tudat::mission_segments;
namespace tss = tudat::simulation_setup;

namespace tudatpy {

void expose_transfer_trajectory(py::module &m) {

    py::enum_<tms::TransferLegTypes>(m, "TransferLegTypes")
            .value("unpowered_unperturbed_leg_type", tms::TransferLegTypes::unpowered_unperturbed_leg)
            .value("dsm_position_based_leg_type", tms::TransferLegTypes::dsm_position_based_leg)
            .value("dsm_velocity_based_leg_type", tms::TransferLegTypes::dsm_velocity_based_leg)
            .export_values();

    py::class_<
            tms::TransferNodeSettings,
            std::shared_ptr<tms::TransferNodeSettings> >(m, "TransferNodeSettings");

    py::class_<
            tms::SwingbyNodeSettings,
            std::shared_ptr<tms::SwingbyNodeSettings>,
            tms::TransferNodeSettings >(m, "SwingbyNodeSettings");

    py::class_<
            tms::EscapeAndDepartureNodeSettings,
            std::shared_ptr<tms::EscapeAndDepartureNodeSettings>,
            tms::TransferNodeSettings >(m, "EscapeAndDepartureNodeSettings");

    py::class_<
            tms::CaptureAndInsertionNodeSettings,
            std::shared_ptr<tms::CaptureAndInsertionNodeSettings>,
            tms::TransferNodeSettings >(m, "CaptureAndInsertionNodeSettings");

    py::class_<
            tms::TransferLegSettings,
            std::shared_ptr<tms::TransferLegSettings> >(m, "TransferLegSettings");

    m.def("mga_transfer_settings",
          py::overload_cast<
          const std::vector< std::string >&,
          const tms::TransferLegTypes,
          const std::pair< double, double >,
          const std::pair< double, double >,
          const std::map< std::string, double >
          >( &tms::getMgaTransferTrajectorySettings ),
          py::arg( "body_order" ),
          py::arg( "leg_type" ),
          py::arg( "departure_orbit" ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
          py::arg( "arrival_orbit" ) = std::make_pair( TUDAT_NAN, TUDAT_NAN ),
          py::arg( "minimum_pericenters" ) = tms::DEFAULT_MINIMUM_PERICENTERS
            );

    py::class_<
            tms::TransferTrajectory,
            std::shared_ptr<tms::TransferTrajectory> >(m, "TransferTrajectory")
            .def_property_readonly("delta_v", &tms::TransferTrajectory::getTotalDeltaV )
            .def("evaluate", &tms::TransferTrajectory::evaluateTrajectory,
                 py::arg( "times" ),
                 py::arg( "leg_parameters" ),
                 py::arg( "node_parameters" ) )
            .def("single_node_delta_v", &tms::TransferTrajectory::getNodeDeltaV,
                 py::arg( "node_index" ) )
            .def("single_leg_delta_v", &tms::TransferTrajectory::getLegDeltaV,
                 py::arg( "leg_index" ) )
            .def_property_readonly("per_node_delta_v", &tms::TransferTrajectory::getDeltaVPerNode )
            .def_property_readonly("per_leg_delta_v", &tms::TransferTrajectory::getDeltaVPerLeg )
            .def_property_readonly( "number_of_nodes", &tms::TransferTrajectory::getNumberOfNodes )
            .def_property_readonly( "number_of_legs", &tms::TransferTrajectory::getNumberOfLegs );


    m.def("unpowered_leg",
          &tms::unpoweredLeg );

    m.def("dsm_position_based_leg",
          &tms::dsmPositionBasedLeg );

    m.def("dsm_velocity_based_leg",
          &tms::dsmVelocityBasedLeg );

    m.def("swingby_node",
          &tms::swingbyNode,
          py::arg( "minimum_periapsis" ) = TUDAT_NAN);

    m.def("departure_node",
          &tms::escapeAndDepartureNode,
          py::arg( "departure_semi_major_axis" ),
          py::arg( "departure_eccentricity" ) );

    m.def("capture_node",
          &tms::captureAndInsertionNode,
          py::arg( "capture_semi_major_axis" ),
          py::arg( "capture_eccentricity" ) );


    m.def("print_parameter_definitions",
          &tms::printTransferParameterDefinition,
          py::arg( "leg_settings" ),
          py::arg( "node_settings" ) );

    m.def("create_transfer_trajectory",
          &tms::createTransferTrajectory,
          py::arg( "bodies" ),
          py::arg( "leg_settings" ),
          py::arg( "node_settings" ),
          py::arg( "node_names" ),
          py::arg( "central_body" ) );


    m.def("get_low_thrust_acceleration_settings",
          &tss::getLowThrustLegAccelerationSettings,
          py::arg("low_thrust_leg"),
          py::arg("bodies"),
          py::arg("body_to_propagate"),
          py::arg("specific_impulse_function"),
          py::arg("low_thrust_leg_initial_time") );
};

}// namespace tudatpy
