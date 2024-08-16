/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/mission_segments/createTransferTrajectory.h>
#include <tudat/math/root_finders.h>
#include <tudat/simulation/propagation_setup/accelerationSettings.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tms = tudat::mission_segments;
namespace tss = tudat::simulation_setup;
namespace tpc = tudat::physical_constants;
namespace trf = tudat::root_finders;
namespace tsbm = tudat::shape_based_methods;

namespace tudatpy {
    namespace trajectory_design {
        namespace transfer_trajectory {

            PYBIND11_MODULE(expose_transfer_trajectory, m) {
                m.attr("DEFAULT_MINIMUM_PERICENTERS") =
                    tms::DEFAULT_MINIMUM_PERICENTERS;

                py::enum_<tms::TransferLegTypes>(
                    m, "TransferLegTypes",
                    R"doc(Enumeration of available leg types.


	:member unpowered_unperturbed_leg_type:
	:member dsm_position_based_leg_type:
	:member dsm_velocity_based_leg_type:
	:member spherical_shaping_low_thrust_leg:
	:member hodographic_low_thrust_leg:
)doc")
                    .value(
                        "unpowered_unperturbed_leg_type",
                        tms::TransferLegTypes::unpowered_unperturbed_leg,
                        get_docstring(
                            "TransferLegTypes.unpowered_unperturbed_leg_type")
                            .c_str())
                    .value("dsm_position_based_leg_type",
                           tms::TransferLegTypes::dsm_position_based_leg,
                           get_docstring(
                               "TransferLegTypes.dsm_position_based_leg_type")
                               .c_str())
                    .value("dsm_velocity_based_leg_type",
                           tms::TransferLegTypes::dsm_velocity_based_leg,
                           get_docstring(
                               "TransferLegTypes.dsm_velocity_based_leg_type")
                               .c_str())
                    .value(
                        "spherical_shaping_low_thrust_leg",
                        tms::TransferLegTypes::spherical_shaping_low_thrust_leg,
                        get_docstring(
                            "TransferLegTypes.spherical_shaping_low_thrust_leg")
                            .c_str())
                    .value("hodographic_low_thrust_leg",
                           tms::TransferLegTypes::hodographic_low_thrust_leg,
                           get_docstring(
                               "TransferLegTypes.hodographic_low_thrust_leg")
                               .c_str())
                    .export_values();

                py::class_<tms::TransferLeg,
                           std::shared_ptr<tms::TransferLeg> >(
                    m, "TransferLeg",
                    R"doc(Base class for defining a transfer leg.

	Functional (base) class for transfer legs, requiring the leg type, departure body ephemeris and arrival body ephemeris.
	Transfer node classes requiring additional information must be created using an object derived from this class.

)doc")
                    .def("state_along_trajectory",
                         py::overload_cast<const double>(
                             &tms::TransferLeg::getStateAlongTrajectory),
                         py::arg("time_since_leg_beginning"),
                         get_docstring("TransferLeg.time_since_leg_beginning")
                             .c_str());


                py::class_<tsbm::SphericalShapingLeg,
                           std::shared_ptr<tsbm::SphericalShapingLeg>,
                           tms::TransferLeg>(
                    m, "SphericalShapingLeg",
                    R"doc(Class for defining low-thrust spherical-shaping leg.

	`TransferLeg` derived class for defining a low-thrust leg described using spherical shaping [3]_.

)doc");

                py::class_<tsbm::HodographicShapingLeg,
                           std::shared_ptr<tsbm::HodographicShapingLeg>,
                           tms::TransferLeg>(
                    m, "HodographicShapingLeg",
                    R"doc(Class for defining low-thrust hodographic-shaping leg.

	`TransferLeg` derived class for defining a low-thrust leg described using hodographic shaping [4]_.

)doc");

                py::class_<tms::TransferNodeSettings,
                           std::shared_ptr<tms::TransferNodeSettings> >(
                    m, "TransferNodeSettings",
                    R"doc(Base class for providing settings for transfer nodes.

	Functional (base) class for settings of transfer nodes that require no information in addition to their type.
	Transfer node classes requiring additional information must be created using an object derived from this class.

)doc");

                py::class_<tms::SwingbyNodeSettings,
                           std::shared_ptr<tms::SwingbyNodeSettings>,
                           tms::TransferNodeSettings>(
                    m, "SwingbyNodeSettings",
                    R"doc(Class for defining settings of swingby node.

	`TransferNodeSettings` derived class for providing settings for swingby nodes, which consist of the minimum periapsis
	radius.

)doc");

                py::class_<tms::EscapeAndDepartureNodeSettings,
                           std::shared_ptr<tms::EscapeAndDepartureNodeSettings>,
                           tms::TransferNodeSettings>(
                    m, "EscapeAndDepartureNodeSettings",
                    R"doc(Class for defining settings of escape and departure node.

	`TransferNodeSettings` derived class for providing settings for escape and departure nodes, which consist of the
	departure semi-major axis and eccentricity.

)doc");

                py::class_<
                    tms::CaptureAndInsertionNodeSettings,
                    std::shared_ptr<tms::CaptureAndInsertionNodeSettings>,
                    tms::TransferNodeSettings>(
                    m, "CaptureAndInsertionNodeSettings",
                    R"doc(Class for defining settings of capture and insertion node.

	`TransferNodeSettings` derived class for providing settings for capture and insertion nodes, which consist of the
	capture semi-major axis and eccentricity.

)doc");

                py::class_<tms::TransferLegSettings,
                           std::shared_ptr<tms::TransferLegSettings> >(
                    m, "TransferLegSettings",
                    R"doc(Base class for providing settings for transfer legs.

	Functional (base) class for settings of transfer legs that require no information in addition to their type.

)doc");

                m.def(
                    "mga_settings_unpowered_unperturbed_legs",
                    py::overload_cast<const std::vector<std::string>&,
                                      const std::pair<double, double>,
                                      const std::pair<double, double>,
                                      const std::map<std::string, double> >(
                        &tms::getMgaTransferTrajectorySettingsWithoutDsm),
                    py::arg("body_order"),
                    py::arg("departure_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("arrival_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("minimum_pericenters") =
                        tms::DEFAULT_MINIMUM_PERICENTERS,
                    R"doc(Function to get the legs and nodes settings of a transfer with just upowered legs.


	Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
	one initial node (departure or swingby), unpowered transfer leg(s) connected by swingby nodes, one final
	(capture or swingby) node.
	If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
	final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
	node.


	:param body_order:
		List of bodies to visit, including departure body, swingby bodies and arrival body.
	:param departure_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
		node as a swingby node (instead of a departure node).

	:param arrival_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
		node as a swingby node (instead of a capture node).

	:param minimum_pericenters:
		Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
		value. Default values from Izzo [1]_.

	:return:
		Tuple specifying the settings of each transfer leg and node.
)doc");

                m.def(
                    "mga_settings_dsm_position_based_legs",
                    py::overload_cast<const std::vector<std::string>&,
                                      const std::pair<double, double>,
                                      const std::pair<double, double>,
                                      const std::map<std::string, double> >(
                        &tms::
                            getMgaTransferTrajectorySettingsWithPositionBasedDsm),
                    py::arg("body_order"),
                    py::arg("departure_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("arrival_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("minimum_pericenters") =
                        tms::DEFAULT_MINIMUM_PERICENTERS,
                    R"doc(Function to get the legs and nodes settings of a transfer constituted by legs with 1 impulsive deep space maneuver (DSM)
described using the position formulation.


	Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
	one initial node (departure or swingby), position-based DSM transfer leg(s) connected by swingby nodes, one final
	(capture or swingby) node.
	If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
	final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
	node.


	:param body_order:
		List of bodies to visit, including departure body, swingby bodies and arrival body.
	:param departure_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
		node as a swingby node (instead of a departure node).

	:param arrival_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
		node as a swingby node (instead of a capture node).

	:param minimum_pericenters:
		Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
		value. Default values from Izzo [1]_.

	:return:
		Tuple specifying the settings of each transfer leg and node.
)doc");

                m.def(
                    "mga_settings_dsm_velocity_based_legs",
                    py::overload_cast<const std::vector<std::string>&,
                                      const std::pair<double, double>,
                                      const std::pair<double, double>,
                                      const std::map<std::string, double> >(
                        &tms::
                            getMgaTransferTrajectorySettingsWithVelocityBasedDsm),
                    py::arg("body_order"),
                    py::arg("departure_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("arrival_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("minimum_pericenters") =
                        tms::DEFAULT_MINIMUM_PERICENTERS,
                    R"doc(Function to get the legs and nodes settings of a transfer constituted by legs with 1 impulsive deep space maneuver (DSM)
described using the velocity formulation.


	Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
	one initial node (departure or swingby), velocity-based DSM transfer leg(s) connected by swingby nodes, one final
	(capture or swingby) node.
	If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
	final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
	node.


	:param body_order:
		List of bodies to visit, including departure body, swingby bodies and arrival body.
	:param departure_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
		node as a swingby node (instead of a departure node).

	:param arrival_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
		node as a swingby node (instead of a capture node).

	:param minimum_pericenters:
		Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
		value. Default values from Izzo [1]_.

	:return:
		Tuple specifying the settings of each transfer leg and node.
)doc");

                m.def(
                    "mga_settings_spherical_shaping_legs",
                    py::overload_cast<
                        const std::vector<std::string>&,
                        const std::shared_ptr<trf::RootFinderSettings>,
                        const std::pair<double, double>,
                        const std::pair<double, double>, const double,
                        const double, const double,
                        const std::map<std::string, double> >(
                        &tms::
                            getMgaTransferTrajectorySettingsWithSphericalShapingThrust),
                    py::arg("body_order"), py::arg("root_finder_settings"),
                    py::arg("departure_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("arrival_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("lower_bound_free_coefficient") = TUDAT_NAN,
                    py::arg("upper_bound_free_coefficient") = TUDAT_NAN,
                    py::arg("initial_value_free_coefficient") = TUDAT_NAN,
                    py::arg("minimum_pericenters") =
                        tms::DEFAULT_MINIMUM_PERICENTERS,
                    R"doc(Function to get the legs and nodes settings of a transfer constituted by low-thrust spherical shaping legs.


	Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
	one initial node (departure or swingby), spherical shaping leg(s) connected by swingby nodes, one final
	(capture or swingby) node.
	If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
	final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
	node.


	:param body_order:
		List of bodies to visit, including departure body, swingby bodies and arrival body.
	:param root_finder_settings:
		Settings of the root finder used by the spherical shaping leg when computing the value of the free coefficient
		that allows meeting the desired time of flight.

	:param departure_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
		node as a swingby node (instead of a departure node).

	:param arrival_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
		node as a swingby node (instead of a capture node).

	:param lower_bound_free_coefficient:
		Lower bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
		it must be specified if the selected root finder requires the definition of a lower bound.

	:param upper_bound_free_coefficient:
		Upper bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
		it must be specified if the selected root finder requires the definition of an upper bound.

	:param initial_value_free_coefficient:
		Initial guess for the free coeficient. Parameter is potentially used by the root finder:
		it must be specified if the selected root finder requires the definition of an initial guess.

	:param minimum_pericenters:
		Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
		value. Default values from Izzo [1]_.

	:return:
		Tuple specifying the settings of each transfer leg and node.
)doc");

                m.def(
                    "mga_settings_hodographic_shaping_legs",
                    py::overload_cast<
                        const std::vector<std::string>&,
                        const std::vector<tsbm::HodographicBasisFunctionList>&,
                        const std::vector<tsbm::HodographicBasisFunctionList>&,
                        const std::vector<tsbm::HodographicBasisFunctionList>&,
                        const std::pair<double, double>,
                        const std::pair<double, double>,
                        const std::map<std::string, double> >(
                        &tms::
                            getMgaTransferTrajectorySettingsWithHodographicShapingThrust),
                    py::arg("body_order"),
                    py::arg("radial_velocity_function_components_per_leg"),
                    py::arg("normal_velocity_function_components_per_leg"),
                    py::arg("axial_velocity_function_components_per_leg"),
                    py::arg("departure_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("arrival_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("minimum_pericenters") =
                        tms::DEFAULT_MINIMUM_PERICENTERS,
                    R"doc(Function to get the legs and nodes settings of a transfer constituted by low-thrust hodographic shaping legs,
with user-provided velocity shaping functions.


	Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
	one initial node (departure or swingby), hodographic shaping leg(s) connected by swingby nodes, one final
	(capture or swingby) node.
	If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
	final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
	node.


	:param body_order:
		List of bodies to visit, including departure body, swingby bodies and arrival body.
	:param radial_velocity_function_components_per_leg:
		List with the lists of radial velocity function components used in each leg.

	:param normal_velocity_function_components_per_leg:
		List with the lists of normal velocity function components used in each leg.

	:param axial_velocity_function_components_per_leg:
		List with the lists of axial velocity function components used in each leg.

	:param departure_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
		node as a swingby node (instead of a departure node).

	:param arrival_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
		node as a swingby node (instead of a capture node).

	:param minimum_pericenters:
		Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
		value. Default values from Izzo [1]_.

	:return:
		Tuple specifying the settings of each transfer leg and node.
)doc");

                m.def(
                    "mga_settings_hodographic_shaping_legs_with_recommended_"
                    "functions",
                    py::overload_cast<const std::vector<std::string>&,
                                      const std::vector<double>&,
                                      const std::vector<double>&,
                                      const std::pair<double, double>,
                                      const std::pair<double, double>,
                                      const std::map<std::string, double> >(
                        &tms::
                            getMgaTransferTrajectorySettingsWithHodographicShapingThrust),
                    py::arg("body_order"), py::arg("time_of_flight_per_leg"),
                    py::arg("number_of_revolutions_per_leg"),
                    py::arg("departure_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("arrival_orbit") =
                        std::make_pair(TUDAT_NAN, TUDAT_NAN),
                    py::arg("minimum_pericenters") =
                        tms::DEFAULT_MINIMUM_PERICENTERS,
                    R"doc(Function to get the legs and nodes settings of a transfer constituted by low-thrust hodographic shaping legs,
with user-provided velocity shaping functions.


	Function determines the legs and nodes settings of a multi-gravity assist transfer trajectory consisting of:
	one initial node (departure or swingby), hodographic shaping leg(s) connected by swingby nodes, one final
	(capture or swingby) node.
	If the departure orbit and arrival orbit are provided as arguments, the initial node is a departure node, and the
	final node a capture node. If the departure/arrival orbit is not specified, the initial/final node is a swingby
	node.


	:param body_order:
		List of bodies to visit, including departure body, swingby bodies and arrival body.
	:param radial_velocity_function_components_per_leg:
		List with the lists of radial velocity function components used in each leg.

	:param normal_velocity_function_components_per_leg:
		List with the lists of normal velocity function components used in each leg.

	:param axial_velocity_function_components_per_leg:
		List with the lists of axial velocity function components used in each leg.

	:param departure_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the departure orbit. The default values define the first
		node as a swingby node (instead of a departure node).

	:param arrival_orbit:
		Tuple of (semi-major axis, eccentricity) specifying the arrival orbit. The default values define the last
		node as a swingby node (instead of a capture node).

	:param minimum_pericenters:
		Minimum pericenter radii, where each body is specified as key and the respective minimum pericenter radius as
		value. Default values from Izzo [1]_.

	:return:
		Tuple specifying the settings of each transfer leg and node.
)doc");

                py::class_<tms::TransferTrajectory,
                           std::shared_ptr<tms::TransferTrajectory> >(
                    m, "TransferTrajectory",
                    R"doc(Class defining a transfer trajectory constituted by transfer legs and nodes.

	Class defining a transfer trajectory constituted by transfer legs and nodes. The object is tipically created using the `create_transfer_trajectory` function.
)doc")
                    .def_property_readonly(
                        "delta_v", &tms::TransferTrajectory::getTotalDeltaV,
                        R"doc(Total Delta V used in the transfer trajectory.
	)doc")
                    .def_property_readonly(
                        "time_of_flight",
                        &tms::TransferTrajectory::getTotalTimeOfFlight,
                        R"doc(Total time of flight of the transfer trajectory.
	)doc")
                    .def("evaluate",
                         &tms::TransferTrajectory::evaluateTrajectory,
                         py::arg("node_times"), py::arg("leg_parameters"),
                         py::arg("node_parameters"),
                         R"doc(Evaluate transfer trajectory.

	Function to evaluate the transfer trajectory, which consists of computing the transfer Delta V and time of
	flight, for the specified set of parameters.


	:param node_times:
		List of the time at each node.
	:param leg_parameters:
		List of lists with the parameters characterizing each leg. Each inner list corresponds to the
		parameters of one leg; if a leg does not require any parameter, its list can contain any value(s),
		therefore it is recommended to leave it empty.

	:param node_parameters:
		List of lists with the parameters characterizing each node. Each inner list corresponds to the
		parameters of one node; if a node does not require any parameter, its list can contain any value(s),
		therefore it is recommended to leave it empty.

	:return:
)doc")
                    .def(
                        "single_node_delta_v",
                        &tms::TransferTrajectory::getNodeDeltaV,
                        py::arg("node_index"),
                        R"doc(Retrieves the Delta V applied in the specified node.

	:param node_index:
		Index of the node for which the Delta V to be is retrieved.
	:return:
		Delta V for the specified node.
)doc")
                    .def(
                        "single_leg_delta_v",
                        &tms::TransferTrajectory::getLegDeltaV,
                        py::arg("leg_index"),
                        R"doc(Retrieves the Delta V applied in the specified leg.

	:param leg_index:
		Index of the leg for which the Delta V is to be retrieved.
	:return:
		Delta V for the specified leg.
)doc")
                    .def(
                        "states_along_trajectory",
                        py::overload_cast<const int>(
                            &tms::TransferTrajectory::getStatesAlongTrajectory),
                        py::arg("number_of_data_points_per_leg"),
                        R"doc(Returns the state history throughout the trajectory.

	Function that returns the state history throughout the trajectory, using the same number of data points in
	each leg. For each leg, the retrieved states are equally spaced in time.


	:param number_of_data_points_per_leg:
		Number of data points used to describe each leg.
	:return:
		Tuple of (state history, time history).
)doc")
                    .def(
                        "inertial_thrust_accelerations_along_trajectory",
                        py::overload_cast<const int>(
                            &tms::TransferTrajectory::
                                getInertialThrustAccelerationsAlongTrajectory),
                        py::arg("number_of_data_points_per_leg"),
                        R"doc(Returns the inertial thrust acceleration history throughout the trajectory.

	Function that returns the inertial thrust acceleration history throughout the trajectory, using the same number of data points in
	each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
	For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.


	:param number_of_data_points_per_leg:
		Number of data points used to describe each leg.
	:return:
		Tuple of (state history, time history).
)doc")
                    .def(
                        "rsw_thrust_accelerations_along_trajectory",
                        py::overload_cast<const int>(
                            &tms::TransferTrajectory::
                                getRswThrustAccelerationsAlongTrajectory),
                        py::arg("number_of_data_points_per_leg"),
                        R"doc(Returns the thrust acceleration history in the RSW frame throughout the trajectory.

	Function that returns the thrust acceleration history in the RSW frame throughout the trajectory, using the same number of data points in
	each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
	For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.


	:param number_of_data_points_per_leg:
		Number of data points used to describe each leg.
	:return:
		Tuple of (state history, time history).
)doc")
                    .def(
                        "tnw_thrust_accelerations_along_trajectory",
                        py::overload_cast<const int>(
                            &tms::TransferTrajectory::
                                getTnwThrustAccelerationsAlongTrajectory),
                        py::arg("number_of_data_points_per_leg"),
                        R"doc(Returns the thrust acceleration history in the TNW frame throughout the trajectory.

	Function that returns the thrust acceleration history in the TNW frame throughout the trajectory, using the same number of data points in
	each leg. For each leg, the retrieved thrust accelerations are equally spaced in time.
	For high-thrust legs (where only impulsive Delta Vs are applied) the thrust acceleration is always zero.


	:param number_of_data_points_per_leg:
		Number of data points used to describe each leg.
	:return:
		Tuple of (state history, time history).
)doc")
                    .def_property_readonly(
                        "delta_v_per_node",
                        &tms::TransferTrajectory::getDeltaVPerNode,
                        R"doc(List of the Delta V applied in each node.
	)doc")
                    .def_property_readonly(
                        "delta_v_per_leg",
                        &tms::TransferTrajectory::getDeltaVPerLeg,
                        R"doc(List of the Delta V applied in each leg.
	)doc")
                    .def_property_readonly(
                        "number_of_nodes",
                        &tms::TransferTrajectory::getNumberOfNodes,
                        R"doc(Number of nodes in the transfer trajectory.
	)doc")
                    .def_property_readonly(
                        "number_of_legs",
                        &tms::TransferTrajectory::getNumberOfLegs,
                        R"doc(Number of legs in the transfer trajectory.
	)doc")
                    .def_property_readonly(
                        "legs", &tms::TransferTrajectory::getLegs,
                        get_docstring("TransferTrajectory.legs").c_str());


                m.def(
                    "unpowered_leg", &tms::unpoweredLeg,
                    R"doc(Factory function for creating the settings of an unpowered leg.


	Factory function for creating the settings of an unpowered leg; the settings consist of just the leg type.
	Given the departure and arrival position, this leg computes the departure and arrival velocity using a Lambert
	targeter.
	The calculations performed in this leg do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Musegaas, 2012 [2]_.

	:return:
		Transfer leg settings.
)doc");

                m.def(
                    "dsm_position_based_leg", &tms::dsmPositionBasedLeg,
                    R"doc(Factory function for creating the settings of a transfer leg with 1 impulsive deep space maneuver (DSM) described using
the position formulation.


	Factory function for creating the settings of a transfer leg with 1 position-based DSM; the settings consist of just the leg type.
	Given the departure position and the DSM position this leg uses a Lambert targeter to compute the departure
	velocity and the velocity before the DSM. Given the DSM position and the arrival position, the leg also uses
	a Lambert targeter to compute the velocity after the DSM and the arrival velocity. The Delta V applied in the
	DSM is computed using the velocity before and after the DSM.
	The calculations performed in this leg do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Musegaas, 2012 [2]_.

	:return:
		Transfer leg settings.
)doc");

                m.def(
                    "dsm_velocity_based_leg", &tms::dsmVelocityBasedLeg,
                    R"doc(Factory function for creating the settings of a transfer leg with 1 impulsive deep space maneuver (DSM) described using
the velocity formulation.


	Factory function for creating the settings of a transfer leg with 1 velocity-based DSM; the settings consist of just the leg type.
	Given the departure position and velocity this leg "propagates" the Kepler elements until the instant of application
	of the DSM (giving the position at the DSM and the velocity before the DSM). Given the position of the DSM
	and the arrival position, it computes the velocity after the DSM
	(which is used to compute the Delta V applied in the DSM) and the arrival velocity using a Lambert targeter.
	The calculations performed in this leg do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Musegaas, 2012 [2]_.

	:return:
		Transfer leg settings.
)doc");

                m.def(
                    "spherical_shaping_leg", &tms::sphericalShapingLeg,
                    py::arg("root_finder_settings"),
                    py::arg("lower_bound_free_coefficient") = TUDAT_NAN,
                    py::arg("upper_bound_free_coefficient") = TUDAT_NAN,
                    py::arg("initial_value_free_coefficient") = TUDAT_NAN,
                    py::arg("time_to_azimuth_interpolator_step_size") =
                        tpc::JULIAN_DAY,
                    R"doc(Factory function for creating the settings of a low-thrust spherical shaping leg.


	Factory function for creating the settings of a low-thrust spherical shaping leg; the settings consist of
	variables necessary for setting up the root finder and variable to set up an interpolator.
	The trajectory is determined via spherical shaping, which shapes the position and time history throughout the
	transfer. The trajectory depends on a single parameter, which is selected using the root finder in order to meet
	a user-specified time of flight.
	The calculations performed in this leg do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Roegiers, 2014 [3]_.


	:param root_finder_settings:
		Settings of the root finder used by the spherical shaping leg when computing the value of the free coefficient
		that allows meeting the desired time of flight.

	:param lower_bound_free_coefficient:
		Lower bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
		it must be specified if the selected root finder requires the definition of a lower bound.

	:param upper_bound_free_coefficient:
		Upper bound of the possible values for the free coeficient. Parameter is potentially used by the root finder:
		it must be specified if the selected root finder requires the definition of an upper bound.

	:param initial_value_free_coefficient:
		Initial guess for the free coeficient. Parameter is potentially used by the root finder:
		it must be specified if the selected root finder requires the definition of an initial guess.

	:param time_to_azimuth_interpolator_step_size:
		Time step size used as reference to define the azimuth values at which the epoch is computed, when defining an
		interpolator to convert between epoch and azimuth.

	:return:
		Transfer leg settings.
)doc");

                m.def(
                    "hodographic_shaping_leg", &tms::hodographicShapingLeg,
                    py::arg("radial_velocity_function_components"),
                    py::arg("normal_velocity_function_components"),
                    py::arg("axial_velocity_function_components"),
                    R"doc(Factory function for creating the settings of a low-thrust hodographic shaping leg.


	Factory function for creating the settings of a low-thrust hodographic shaping leg; the settings consist of
	the functions used to shape the velocity throughout the transfer.
	Note that shape functions with at least 3 terms (i.e. 3 degrees of freedom) must be provided for each velocity
	component (radial, normal and axial); this is required to obtain a trajectory that satisfies the boundary
	conditions.
	The calculations performed in this leg do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Gondelach, 2012 [4]_.


	:param radial_velocity_function_components:
		List with components of the radial velocity shaping function, which determine the radial velocity profile
		throughout the transfer.

	:param normal_velocity_function_components:
		List with components of the normal velocity shaping function, which determine the normal velocity profile
		throughout the transfer.

	:param axial_velocity_function_components:
		List with components of the axial velocity shaping function, which determine the axial velocity profile
		throughout the transfer.

	:return:
		Transfer leg settings.
)doc");

                m.def(
                    "swingby_node", &tms::swingbyNode,
                    py::arg("minimum_periapsis") = TUDAT_NAN,
                    R"doc(Factory function for creating the settings of a swingby node.

	Factory function for creating the settings of a swingby node. The settings consist consist of the minimum
	allowed periapsis radius.
	The minimum periapsis radius can be set to infinity. In that case, the swingby does not affect the velocity of the
	spacecraft (e.g. might be relevant for swingbys of small bodies).
	The exact behavior of this node depends on the types of legs that precede and follow it. Given a known incoming
	and unknown outgoing velocity, the node forward propagates the gravity assist, possibly with a Delta V at
	the a periapsis. Given an unknown incoming and known outgoing velocity, the node backward propagates the
	gravity assist, possibly with a Delta V at the a periapsis. Given known incoming and outgoing velocities,
	the node computes the Delta V required to meet those.
	The calculations performed in this node do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Musegaas, 2012 [2]_.


	:param minimum_periapsis:
		Minimum periapsis radius. The minimum periapsis only needs to be specified if the types of swignby nodes that
		requires it is used. If that is the case and no minimum periapsis was selected an error is thrown.

	:return:
		Swingby node settings.
)doc");

                m.def(
                    "departure_node", &tms::escapeAndDepartureNode,
                    py::arg("departure_semi_major_axi    s"),
                    py::arg("departure_eccentricity"),
                    R"doc(Factory function for creating the settings of an escape or departure node.

	Factory function for creating the settings of an escape or departure node. The settings consist of the
	departure orbit eccentricity and semi-major axis.
	Given the initial orbit and the departure velocity, the node computes the Delta V that needs to be applied
	at the periapsis of the initial orbit to enter the escape trajectory.
	The calculations performed in this node do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Musegaas, 2012 [2]_.


	:param departure_semi_major_axis:
		Departure orbit semi-major axis.
	:param departure_eccentricity:
		Departure orbit eccentricity.
	:return:
		Escape or departure node settings.
)doc");

                m.def(
                    "capture_node", &tms::captureAndInsertionNode,
                    py::arg("capture_semi_major_axis"),
                    py::arg("capture_eccentricity"),
                    R"doc(Factory function for creating the settings of a capture or insertion node.

	Factory function for creating the settings of a capture or insertion node. The settings consist of the
	capture orbit eccentricity and semi-major axis.
	Given the the arrival velocity and the final orbit, the node computes the Delta V that needs to be applied
	at the periapsis of the final orbit to exit the capture trajectory.
	The calculations performed in this node do not involve any numerical integration, and are solved by
	(semi-)analytical models. For details see Musegaas, 2012 [2]_.


	:param capture_semi_major_axis:
		Capture orbit semi-major axis.
	:param capture_eccentricity:
		Capture orbit eccentricity.
	:return:
		Capture or insertion node settings.
)doc");

                m.def(
                    "print_parameter_definitions",
                    &tms::printTransferParameterDefinition,
                    py::arg("leg_settings"), py::arg("node_settings"),
                    R"doc(Prints the list of parameters required to define the transfer trajectory, according to the
specified node and leg settings.


	:param leg_settings:
		List of transfer leg settings.
	:param node_settings:
		List of transfer node settings.
	:return:
)doc");

                m.def(
                    "create_transfer_trajectory",
                    &tms::createTransferTrajectory, py::arg("bodies"),
                    py::arg("leg_settings"), py::arg("node_settings"),
                    py::arg("node_names"), py::arg("central_body"),
                    R"doc(Factory function for creating a transfer trajectory consisting of the specified sequence of transfer nodes and
transfer legs.


	The function creates a transfer trajectory based on the provided transfer nodes settings and transfer legs
	settings. The number of nodes should be equal to the number of legs plus 1.
	This function creates an instance of the `TransferTrajectory` class.


	:param bodies:
		System of bodies to be used in the transfer trajectory.
	:param leg_settings:
		List of transfer leg settings.
	:param node_settings:
		List of transfer node settings.
	:param node_names:
		Sequence of bodies used as transfer nodes.
	:param central_body:
		Central body with respect to which the two-body trajectory of the spacecraft
		is calculated.

	:return:
		Transfer trajectory object.
)doc");

                m.def("set_low_thrust_acceleration",
                      &tms::setLowThrustAcceleration, py::arg("transfer_leg"),
                      py::arg("bodies"), py::arg("body_name"),
                      py::arg("engine_name"),
                      get_docstring("set_low_thrust_acceleration").c_str());
            };

        }  // namespace transfer_trajectory
    }  // namespace trajectory_design
}  // namespace tudatpy
