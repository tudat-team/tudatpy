/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/aerodynamics.h>
#include <tudat/astro/ephemerides.h>
#include <tudat/astro/gravitation.h>
#include <tudat/astro/ground_stations/groundStation.h>
#include <tudat/basics/deprecationWarnings.h>
#include <tudat/simulation/environment_setup/body.h>
#include <tudat/simulation/environment_setup/createGroundStations.h>

#include "tudatpy/scalarTypes.h"

// namespace py = pybind11;
// namespace tba = tudat::basic_astrodynamics;
// namespace tss = tudat::simulation_setup;
// namespace tp = tudat::propagators;
// namespace tinterp = tudat::interpolators;
// namespace te = tudat::ephemerides;
// namespace tni = tudat::numerical_integrators;
// namespace trf = tudat::reference_frames;
// namespace tmrf = tudat::root_finders;


namespace py = pybind11;

namespace tba = tudat::basic_astrodynamics;
namespace ta = tudat::aerodynamics;
namespace tr = tudat::reference_frames;
namespace te = tudat::ephemerides;
namespace tgs = tudat::ground_stations;
namespace tr = tudat::reference_frames;
namespace tg = tudat::gravitation;
namespace trf = tudat::reference_frames;
namespace tss = tudat::simulation_setup;
namespace ti = tudat::interpolators;
namespace tsm = tudat::system_models;
namespace tom = tudat::observation_models;


namespace tudat {

    namespace aerodynamics {

        double getTotalSurfaceArea(
            const std::shared_ptr<HypersonicLocalInclinationAnalysis>
                coefficientGenerator) {
            double totalSurfaceArea = 0.0;
            for(int i = 0; i < coefficientGenerator->getNumberOfVehicleParts();
                i++) {
                totalSurfaceArea += std::fabs(
                    coefficientGenerator->getVehiclePart(i)->getTotalArea());
            }
            return totalSurfaceArea;
        }


        //! Function that saves the vehicle mesh data used for a
        //! HypersonicLocalInclinationAnalysis to a file
        std::pair<std::vector<Eigen::Vector3d>, std::vector<Eigen::Vector3d>>
        getVehicleMesh(const std::shared_ptr<HypersonicLocalInclinationAnalysis>
                           localInclinationAnalysis) {
            std::vector<boost::multi_array<Eigen::Vector3d, 2>> meshPoints =
                localInclinationAnalysis->getMeshPoints();
            std::vector<boost::multi_array<Eigen::Vector3d, 2>>
                meshSurfaceNormals =
                    localInclinationAnalysis->getPanelSurfaceNormals();


            //    boost::array< int, 3 > independentVariables;
            //    independentVariables[ 0 ] = 0;
            //    independentVariables[ 1 ] = 6;
            //    independentVariables[ 2 ] = 0;

            //    std::vector< std::vector< std::vector< double > > >
            //    pressureCoefficients =
            //            localInclinationAnalysis->getPressureCoefficientList(
            //            independentVariables );

            int counter = 0;
            std::vector<Eigen::Vector3d> meshPointsList;
            std::vector<Eigen::Vector3d> surfaceNormalsList;
            //    std::map< int, Eigen::Vector1d > pressureCoefficientsList;

            for(unsigned int i = 0; i < meshPoints.size(); i++) {
                for(unsigned int j = 0; j < meshPoints.at(i).shape()[0] - 1;
                    j++) {
                    for(unsigned int k = 0; k < meshPoints.at(i).shape()[1] - 1;
                        k++) {
                        meshPointsList.push_back(meshPoints[i][j][k]);
                        surfaceNormalsList.push_back(
                            meshSurfaceNormals[i][j][k]);
                        //                pressureCoefficientsList[ counter ] =
                        //                ( Eigen::Vector1d( ) <<
                        //                pressureCoefficients[ i ][ j ][ k ]
                        //                ).finished( );
                        counter++;
                    }
                }
            }

            return std::make_pair(meshPointsList, surfaceNormalsList);
        }

    }  // namespace aerodynamics

}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace environment {

            PYBIND11_MODULE(expose_environment, m) {
                py::enum_<ta::AerodynamicCoefficientsIndependentVariables>(
                    m, "AerodynamicCoefficientsIndependentVariables",
                    R"doc(Enumeration of the independent variables that can be used to compute aerodynamic coefficients.



	:member mach_number_dependent:
	:member angle_of_attack_dependent:
	:member sideslip_angle_dependent:
	:member altitude_dependent:
	:member time_dependent:
	:member control_surface_deflection_dependent:
	:member undefined_independent_variable:
)doc")
                    .value("mach_number_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               mach_number_dependent,
                           "")
                    .value("angle_of_attack_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               angle_of_attack_dependent,
                           "")
                    .value("sideslip_angle_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               angle_of_sideslip_dependent,
                           "")
                    .value("altitude_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               altitude_dependent,
                           "")
                    .value("time_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               time_dependent,
                           "")
                    .value("temperature_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               temperature_dependent,
                           "")
                    .value("velocity_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               velocity_dependent,
                           "")
                    .value("he_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               he_number_density_dependent,
                           "")
                    .value("o_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               o_number_density_dependent,
                           "")
                    .value("n2_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               n2_number_density_dependent,
                           "")
                    .value("o2_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               o2_number_density_dependent,
                           "")
                    .value("ar_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               ar_number_density_dependent,
                           "")
                    .value("h_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               h_number_density_dependent,
                           "")
                    .value("n_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               n_number_density_dependent,
                           "")
                    .value("anomalous_o_number_density_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               anomalous_o_number_density_dependent,
                           "")
                    .value("control_surface_deflection_dependent",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               control_surface_deflection_dependent,
                           "")
                    .value("undefined_independent_variable",
                           ta::AerodynamicCoefficientsIndependentVariables::
                               undefined_independent_variable,
                           "")
                    .export_values();


                py::enum_<ta::AerodynamicCoefficientFrames>(
                    m, "AerodynamicCoefficientFrames", "")
                    .value("positive_body_fixed_frame_coefficients",
                           ta::AerodynamicCoefficientFrames::
                               body_fixed_frame_coefficients,
                           "")
                    .value("negative_body_fixed_frame_coefficients",
                           ta::AerodynamicCoefficientFrames::
                               negative_body_fixed_frame_coefficients,
                           "")
                    .value("positive_aerodynamic_frame_coefficients",
                           ta::AerodynamicCoefficientFrames::
                               positive_aerodynamic_frame_coefficients,
                           "")
                    .value("negative_aerodynamic_frame_coefficients",
                           ta::AerodynamicCoefficientFrames::
                               negative_aerodynamic_frame_coefficients,
                           "")
                    .export_values();


                py::class_<
                    ta::AerodynamicCoefficientInterface,
                    std::shared_ptr<ta::AerodynamicCoefficientInterface>>(
                    m, "AerodynamicCoefficientInterface", "")
                    .def_property_readonly(
                        "reference_area",
                        &ta::AerodynamicCoefficientInterface::getReferenceArea,
                        "")
                    .def_property_readonly(
                        "current_force_coefficients",
                        &ta::AerodynamicCoefficientInterface::
                            getCurrentForceCoefficients,
                        "")
                    .def_property_readonly(
                        "current_moment_coefficients",
                        &ta::AerodynamicCoefficientInterface::
                            getCurrentMomentCoefficients,
                        "")
                    .def_property_readonly(
                        "current_coefficients",
                        &ta::AerodynamicCoefficientInterface::
                            getCurrentAerodynamicCoefficients,
                        "")
                    .def_property_readonly(
                        "force_coefficient_frame",
                        &ta::AerodynamicCoefficientInterface::
                            getForceCoefficientsFrame,
                        "")
                    .def_property_readonly(
                        "moment_coefficient_frame",
                        &ta::AerodynamicCoefficientInterface::
                            getMomentCoefficientsFrame,
                        "")
                    .def_property_readonly(
                        "independent_variable_names",
                        &ta::AerodynamicCoefficientInterface::
                            getIndependentVariableNames,
                        "")
                    .def_property_readonly(
                        "current_control_surface_free_force_coefficients",
                        &ta::AerodynamicCoefficientInterface::
                            getCurrentControlSurfaceFreeForceCoefficients,
                        "")
                    .def_property_readonly(
                        "current_control_surface_free_moment_coefficients",
                        &ta::AerodynamicCoefficientInterface::
                            getCurrentControlSurfaceFreeMomentCoefficients,
                        "")
                    .def_property_readonly(
                        "control_surface_independent_variable_names",
                        &ta::AerodynamicCoefficientInterface::
                            getControlSurfaceIndependentVariables,
                        "")
                    .def("current_control_surface_force_coefficient_increment",
                         &ta::AerodynamicCoefficientInterface::
                             getCurrentForceCoefficientIncrement,
                         py::arg("control_surface_name"), "")
                    .def("current_control_surface_moment_coefficient_increment",
                         &ta::AerodynamicCoefficientInterface::
                             getCurrentMomentCoefficientIncrement,
                         py::arg("control_surface_name"), "")
                    .def("set_control_surface_increments",
                         &ta::AerodynamicCoefficientInterface::
                             setControlSurfaceIncrements,
                         py::arg("control_surface_list"), "")
                    .def("update_coefficients",
                         &ta::AerodynamicCoefficientInterface::
                             updateCurrentCoefficients,
                         py::arg("independent_variables"), py::arg("time"), "")
                    .def("update_full_coefficients",
                         &ta::AerodynamicCoefficientInterface::
                             updateFullCurrentCoefficients,
                         py::arg("independent_variables"),
                         py::arg("control_surface_independent_variables"),
                         py::arg("time"),
                         py::arg("check_force_contribution") = true, "");

                py::class_<
                    ta::AerodynamicCoefficientGenerator<3, 6>,
                    std::shared_ptr<ta::AerodynamicCoefficientGenerator<3, 6>>,
                    ta::AerodynamicCoefficientInterface>(
                    m, "AerodynamicCoefficientGenerator36",
                    "<no_doc, only_dec>");

                py::class_<
                    ta::HypersonicLocalInclinationAnalysis,
                    std::shared_ptr<ta::HypersonicLocalInclinationAnalysis>,
                    ta::AerodynamicCoefficientGenerator<3, 6>>(
                    m, "HypersonicLocalInclinationAnalysis")
                    .def(
                        py::init<const std::vector<std::vector<double>> &,
                                 const std::shared_ptr<tudat::SurfaceGeometry>,
                                 const std::vector<int> &,
                                 const std::vector<int> &,
                                 const std::vector<bool> &,
                                 const std::vector<std::vector<int>> &,
                                 const double, const double,
                                 const Eigen::Vector3d &, const bool>(),
                        py::arg("independent_variable_points"),
                        py::arg("body_shape"), py::arg("number_of_lines"),
                        py::arg("number_of_points"), py::arg("invert_orders"),
                        py::arg("selected_methods"), py::arg("reference_area"),
                        py::arg("reference_length"),
                        py::arg("moment_reference_point"),
                        py::arg("save_pressure_coefficients") = false,
                        R"doc(Class constructor, taking the shape of the vehicle, and various analysis options as input.

	:param independent_variable_points:
		List containing three lists, with each sublist containing the data points of each of the
		independent variables for the coefficient generation. The physical meaning of each of the
		three independent variables is: 0 = mach number, 1 = angle of attack, 2 = angle of sideslip.
		Each of the subvectors must be sorted in ascending order.

	:param body_shape:
		Class that defines the shape of the vehicle as a continuous surface. The local inclination analysis
		discretizes the surface of the vehicle into quadrilateral panels, defined by the other inputs to
		this constructor. In case the :class:`tudat.geometry.SurfaceGeometry` object is made up of multiple
		sub-shapes, different settings may be used for each

	:param number_of_lines:
		Number of discretization points in the first independent surface variable of each of the subparts of body_shape.
		The size of this list should match the number of parts of which the body_shape is composed. The first independent
		variable of a subpart typically runs along the longitudinal vehicle direction

	:param number_of_points:
		Number of discretization points in the second independent surface variable of each of the subparts of body_shape.
		The size of this list should match the number of parts of which the body_shape is composed. The first independent
		variable of a subpart typically runs along the lateral vehicle direction

	:param invert_orders:
		Booleans to denote whether the surface normals of the panels of each discretized body_shape subpart are to be inverted
		(i.e. inward-facing->outward facing or vice versa). The size of this list should match the number of parts of which the body_shape is composed.

	:param selected_methods:
		Double list of selected local inclination methods, the first index (outer list) represents compression or expansion (0 and 1),
		the second index (inner list) denotes the vehicle part index. The size of this inner list should match the number of parts of which the body_shape is composed.
		The int defining the method type is interpreted as follows.
		For the compression methods, the following are available:
		*  0: Newtonian Method.
		*  1: Modified Newtonian.
		*  2 and 3: not available at this moment.
		*  4: Tangent-wedge method.
		*  5: Tangent-cone method.
		*  6: Modified Dahlem-Buck method.
		*  7: VanDyke unified pressure method.
		*  8: Smyth Delta Wing method.
		*  9: Hankey flat surface method
		The expansion method has the following options:
		*  0: Vacuum Pressure coefficient method.
		*  1: Zero Pressure function.
		*  4: High Mach base pressure method.
		*  3 or 5: Prandtl-Meyer method.
		*  6: ACM empirical pressure coefficient.

	:param reference_area:
		Reference area used to non-dimensionalize aerodynamic forces and moments.

	:param moment_reference_point:
		Reference point wrt which aerodynamic moments are calculated.

	:param save_pressure_coefficients:
		Boolean denoting whether to save the pressure coefficients that are computed to files

)doc")
                    .def("clear_data",
                         &ta::HypersonicLocalInclinationAnalysis::clearData);

                py::class_<
                    ta::ControlSurfaceIncrementAerodynamicInterface,
                    std::shared_ptr<
                        ta::ControlSurfaceIncrementAerodynamicInterface>>(
                    m, "ControlSurfaceIncrementAerodynamicInterface",
                    "<no_doc, only_dec>");

                py::class_<
                    ta::CustomControlSurfaceIncrementAerodynamicInterface,
                    std::shared_ptr<
                        ta::CustomControlSurfaceIncrementAerodynamicInterface>,
                    ta::ControlSurfaceIncrementAerodynamicInterface>(
                    m, "CustomControlSurfaceIncrementAerodynamicInterface",
                    "<no_doc, only_dec>")
                    .def(
                        py::init<
                            const std::function<Eigen::Vector6d(
                                const std::vector<double> &)>,
                            const std::vector<
                                ta::AerodynamicCoefficientsIndependentVariables>>(),
                        py::arg("coefficient_function"),
                        py::arg("independent_variable_names"));


                m.def("get_default_local_inclination_mach_points",
                      &ta::getDefaultHypersonicLocalInclinationMachPoints,
                      py::arg("mach_regime") = "Full");

                m.def(
                    "get_default_local_inclination_angle_of_attack_points",
                    &ta::
                        getDefaultHypersonicLocalInclinationAngleOfAttackPoints);

                m.def(
                    "get_default_local_inclination_sideslip_angle_points",
                    &ta::
                        getDefaultHypersonicLocalInclinationAngleOfSideslipPoints);


                m.def(
                    "save_vehicle_mesh_to_file", &ta::saveVehicleMeshToFile,
                    py::arg("local_inclination_analysis_object"),
                    py::arg("output_directory"),
                    py::arg("output_file_prefix") = "",
                    R"doc(Function to save the mesh used for a hypersonic local inclination analysis to a file.

	Function to save the mesh used for a hypersonic local inclination analysis to a file. This function saves
	two files to the specified directory, with filenames: "ShapeFile.dat" and "SurfaceNormalFile.dat", where these
	files names may be prefixed by an optional string (see below). The first of these files contains four columns defining
	the surface points that define mesh, with Column 0: point index; Column 1: x-position of point; Column 1: y-position of point;
	Column 2: z-position of point. The second file contains four columns with Column 0: point index; Column 1: x-component of surface normal;
	Column 1: y-position of surface normal; Column 2: z-position of surface normal.


	:param local_inclination_analysis_object:
		Object used to calculate the aerodynamics of the vehicle

	:param output_directory:
		Directory to which the files are to be saved

	:param output_file_prefix:
		Optional prefix of output file names

)doc");

                m.def("get_local_inclination_total_vehicle_area",
                      &ta::getTotalSurfaceArea,
                      py::arg("local_inclination_analysis_object"));

                m.def("get_local_inclination_mesh", &ta::getVehicleMesh,
                      py::arg("local_inclination_analysis_object"));


                py::class_<tsm::VehicleSystems,
                           std::shared_ptr<tsm::VehicleSystems>>(
                    m, "VehicleSystems", "")
                    .def(py::init<>())
                    .def("set_control_surface_deflection",
                         &tsm::VehicleSystems::
                             setCurrentControlSurfaceDeflection,
                         py::arg("control_surface_id"),
                         py::arg("deflection_angle"), "")
                    .def(
                        "set_transponder_turnaround_ratio",
                        py::overload_cast<std::map<
                            std::pair<tom::FrequencyBands, tom::FrequencyBands>,
                            double> &>(&tsm::VehicleSystems::
                                           setTransponderTurnaroundRatio),
                        py::arg("transponder_ratio_per_uplink_and_downlink_"
                                "frequency_band"),
                        "")
                    .def("get_control_surface_deflection",
                         &tsm::VehicleSystems::
                             getCurrentControlSurfaceDeflection,
                         py::arg("control_surface_id"), "")
                    .def("get_engine_model",
                         &tsm::VehicleSystems::getEngineModel,
                         py::arg("engine_name"), "");

                py::class_<tsm::EngineModel, std::shared_ptr<tsm::EngineModel>>(
                    m, "EngineModel")
                    .def_property_readonly(
                        "thrust_magnitude_calculator",
                        &tsm::EngineModel::getThrustMagnitudeWrapper);


                /*!
                 **************   FLIGHT CONDITIONS AND ASSOCIATED FUNCTIONALITY
                 *******************
                 */

                py::enum_<trf::AerodynamicsReferenceFrameAngles>(
                    m, "AerodynamicsReferenceFrameAngles")
                    .value(
                        "latitude_angle",
                        trf::AerodynamicsReferenceFrameAngles::latitude_angle)
                    .value(
                        "longitude_angle",
                        trf::AerodynamicsReferenceFrameAngles::longitude_angle)
                    .value("heading_angle",
                           trf::AerodynamicsReferenceFrameAngles::heading_angle)
                    .value("flight_path_angle",
                           trf::AerodynamicsReferenceFrameAngles::
                               flight_path_angle)
                    .value(
                        "angle_of_attack",
                        trf::AerodynamicsReferenceFrameAngles::angle_of_attack)
                    .value("angle_of_sideslip",
                           trf::AerodynamicsReferenceFrameAngles::
                               angle_of_sideslip)
                    .value("bank_angle",
                           trf::AerodynamicsReferenceFrameAngles::bank_angle)
                    .export_values();

                py::enum_<trf::AerodynamicsReferenceFrames>(
                    m, "AerodynamicsReferenceFrames",
                    R"doc(Enumeration of reference frame identifiers typical for aerodynamic calculations.

	Enumeration of reference frame identifiers typical for aerodynamic calculations. Note that the frames are also defined
	in the absence of any aerodynamic forces and/or atmosphere. They define frames of a body w.r.t. a central body, with
	the details given by Mooij (1994). The chain of frames starts from the inertial frame, to the frame fixed to the
	central body (corotating), to the vertical frame (defined by the body's relative position), the trajectory and aerodynamic frames
	(defined by the body's relative velocity) and finally the body's own body-fixed frame.


	:member inertial_frame:
	:member corotating_frame:
	:member vertical_frame:
	:member trajectory_frame:
	:member aerodynamic_frame:
	:member body_frame:
)doc")
                    .value("inertial_frame",
                           trf::AerodynamicsReferenceFrames::inertial_frame, "")
                    .value("corotating_frame",
                           trf::AerodynamicsReferenceFrames::corotating_frame,
                           "")
                    .value("vertical_frame",
                           trf::AerodynamicsReferenceFrames::vertical_frame, "")
                    .value("trajectory_frame",
                           trf::AerodynamicsReferenceFrames::trajectory_frame,
                           "")
                    .value("aerodynamic_frame",
                           trf::AerodynamicsReferenceFrames::aerodynamic_frame,
                           "")
                    .value("body_frame",
                           trf::AerodynamicsReferenceFrames::body_frame, "")
                    .export_values();

                py::class_<trf::AerodynamicAngleCalculator,
                           std::shared_ptr<trf::AerodynamicAngleCalculator>>(
                    m, "AerodynamicAngleCalculator",
                    R"doc(Object to calculate (aerodynamic) orientation angles, and frame transformations,
from current vehicle state.


	Object to calculate (aerodynamic) orientation angles (list given by the :class:`~AerodynamicsReferenceFrameAngles` enum)
	and transformations between frames (list given by the :class:`~AerodynamicsReferenceFrames` enum) from current vehicle state.

)doc")
                    .def(
                        "get_rotation_matrix_between_frames",
                        &trf::AerodynamicAngleCalculator::
                            getRotationMatrixBetweenFrames,
                        py::arg("original_frame"), py::arg("target_frame"),
                        R"doc(Function to get the rotation matrix between two frames.


	Function to get the rotation matrix between two frames. This function
	is meant to be used only *during* a numerical propagation, in particular
	for the definition of a custom (e.g. guidance) model.


	:param original_frame:
		The frame :math:`A` from which the rotation matrix is to be calculated

	:param target_frame:
		The frame :math:`B` to which the rotation matrix is to be calculated

	:return:
		Rotation matrix :math:`\mathbf{R}^{B/A}` from frame :math:`A` to frame `B`
)doc")
                    .def("get_angle",
                         &trf::AerodynamicAngleCalculator::getAerodynamicAngle,
                         py::arg("angle_type"),
                         R"doc(Function to get a single orientation angle


	Function to get a single orientation angle. This function
	is meant to be used only *during* a numerical propagation, in particular
	for the definition of a custom (e.g. guidance) model.


	:param original_frame:
		The identifier for the angle that is to be returnd

	:return:
		Value of requested angle
)doc")
                    // Function removed; error is shown
                    .def("set_body_orientation_angles",
                         &trf::AerodynamicAngleCalculator::
                             setOrientationAngleFunctionsRemoved2,
                         py::arg("angle_of_attack") = TUDAT_NAN,
                         py::arg("angle_of_sideslip") = TUDAT_NAN,
                         py::arg("bank_angle") = TUDAT_NAN,
                         py::arg("silence_warnings") = false)
                    // Function removed; error is shown
                    .def("set_body_orientation_angle_functions",
                         &trf::AerodynamicAngleCalculator::
                             setOrientationAngleFunctionsRemoved1,
                         py::arg("angle_of_attack_function") = std::function<
                             double()>(),  // <pybind11/functional.h>
                         py::arg("angle_of_sideslip_function") = std::function<
                             double()>(),  // <pybind11/functional.h>
                         py::arg("bank_angle_function") = std::function<
                             double()>(),  // <pybind11/functional.h>
                         py::arg("angle_update_function") =
                             std::function<void(const double)>(),
                         py::arg("silence_warnings") = false);


                py::class_<ta::FlightConditions,
                           std::shared_ptr<ta::FlightConditions>>(
                    m, "FlightConditions",
                    R"doc(Object that calculates various state-derived quantities typically
relevant for flight dynamics.


	Object that calculates various state-derived quantities typically
	relevant for flight dynamics, such as latitude, longitude,
	altitude, etc. It also contains an
	:py:class:`~AerodynamicAngleCalculator` that computes derived
	angles (flight path, heading angle, etc.). This object is limited
	to non-atmospheric flight. For flight through Body objects
	endowed with an atmosphere model, the derived class
	:py:class:`~AtmosphericFlightConditions` is used. This object is
	stored inside a Body object, and represents the flight conditions
	of a single body w.r.t. a single central body.

)doc")
                    //            .def(py::init<
                    //                 const
                    //                 std::shared_ptr<tudat::basic_astrodynamics::BodyShapeModel>,
                    //                 const
                    //                 std::shared_ptr<tudat::reference_frames::AerodynamicAngleCalculator>>(),
                    //                 py::arg("shape_model"),
                    //                 py::arg("aerodynamic_angle_calculator") =
                    //                 std::shared_ptr<
                    //                 tr::AerodynamicAngleCalculator>())
                    .def("update_conditions",
                         &ta::FlightConditions::updateConditions,
                         py::arg("current_time"))
                    .def_property_readonly(
                        "aerodynamic_angle_calculator",
                        &ta::FlightConditions::getAerodynamicAngleCalculator,
                        R"doc(The object that is responsible for computing various relevant
flight dynamics angles and frame rotations.

	)doc")
                    .def_property_readonly(
                        "longitude", &ta::FlightConditions::getCurrentLongitude,
                        R"doc(The body-fixed longitude of the body w.r.t. its central body.

	)doc")
                    .def_property_readonly(
                        "latitude", &ta::FlightConditions::getCurrentLatitude,
                        R"doc(The body-fixed geographic latitude of the body w.r.t. its
central body.

	)doc")
                    .def_property_readonly(
                        "geodetic_latitude",
                        &ta::FlightConditions::getCurrentGeodeticLatitude,
                        R"doc(The body-fixed geographic latitude of the body w.r.t. its
central body.

	)doc")
                    .def_property_readonly(
                        "time", &ta::FlightConditions::getCurrentTime,
                        R"doc(The current time, at which this object was last updated

	)doc")
                    .def_property_readonly(
                        "body_centered_body_fixed_state",
                        &ta::FlightConditions::
                            getCurrentBodyCenteredBodyFixedState,
                        R"doc(Cartesian translational state, expressed in a frame centered
at, and fixed to, the central body. Note that, due to the
rotation of the central body, the norm of the body-fixed,
body-centered, velocity differs from the norm of the inertial
body-centered velocity.

	)doc")
                    .def_property_readonly(
                        "altitude", &ta::FlightConditions::getCurrentAltitude,
                        R"doc(The current time, at which this object was last updated

	)doc");

                py::class_<ta::AtmosphericFlightConditions,
                           std::shared_ptr<ta::AtmosphericFlightConditions>,
                           ta::FlightConditions>(
                    m, "AtmosphericFlightConditions",
                    R"doc(Object that calculates various state-derived quantities typically
relevant for flight dynamics, for flight in an atmosphere.


	Object that calculates various state-derived quantities typically
	relevant for flight dynamics, for flight in an atmosphere, such
	as latitude,  longitude, altitude, density, Mach number etc. It
	also contains an ``AerodynamicAngleCalculator`` that computes
	derived angles (flight path, heading angle, etc.). This object is
	derived from ``FlightConditions``, which performs computations for
	non-atmospheric flight only. This object is stored inside a Body
	object, and represents the flight conditions of a single body
	w.r.t. a single central body.

)doc")
                    .def_property_readonly(
                        "density",
                        &ta::AtmosphericFlightConditions::getCurrentDensity,
                        R"doc(The freestream atmospheric density at the body's current
location.

	)doc")
                    .def_property_readonly(
                        "temperature",
                        &ta::AtmosphericFlightConditions::
                            getCurrentFreestreamTemperature,
                        R"doc(The freestream atmospheric temperature at the body's current
location.

	)doc")
                    .def_property_readonly(
                        "dynamic_pressure",
                        &ta::AtmosphericFlightConditions::
                            getCurrentDynamicPressure,
                        R"doc(The freestream atmospheric dynamic pressure at the body's
current location.

	)doc")
                    .def_property_readonly(
                        "pressure",
                        &ta::AtmosphericFlightConditions::getCurrentPressure,
                        R"doc(The freestream atmospheric static pressure at the body's
current location.

	)doc")
                    .def_property_readonly(
                        "airspeed",
                        &ta::AtmosphericFlightConditions::getCurrentAirspeed,
                        R"doc(The airspeed of the body w.r.t. the atmosphere.

	)doc")
                    .def_property_readonly(
                        "mach_number",
                        &ta::AtmosphericFlightConditions::getCurrentMachNumber,
                        R"doc(The freestream Mach number of the body.

	)doc")
                    .def_property_readonly(
                        "airspeed_velocity",
                        &ta::AtmosphericFlightConditions::
                            getCurrentAirspeedBasedVelocity,
                        R"doc(The velocity vector of the body w.r.t. the freestream
atmosphere (e.g. vectorial counterpart of airspeed).

	)doc")
                    .def_property_readonly(
                        "speed_of_sound",
                        &ta::AtmosphericFlightConditions::
                            getCurrentSpeedOfSound,
                        R"doc(The freestream atmospheric speed of sound at the body's current
location.

	)doc")
                    .def_property_readonly(
                        "aero_coefficient_independent_variables",
                        &ta::AtmosphericFlightConditions::
                            getAerodynamicCoefficientIndependentVariables,
                        R"doc(List of current values of independent variables of aerodynamic
coefficients. This list is only defined if the body has an
:py:class:`~AerodynamicCoefficientInterface` that has
dependencies on environmental variables (e.g. Mach number,
angle of attack, etc.).

	)doc")
                    .def_property_readonly(
                        "control_surface_aero_coefficient_independent_"
                        "variables",
                        &ta::AtmosphericFlightConditions::
                            getControlSurfaceAerodynamicCoefficientIndependentVariables,
                        R"doc(List of lists current values of independent variables of
aerodynamic coefficients for control surfaces. The outer list
defines the control surface, the inner list the values of the
independent variables. This list is only defined if the body
has an :py:class:`~AerodynamicCoefficientInterface` with
control surfaces that have dependencies on environmental
variables (e.g. Mach number, angle of attack, etc.).

	)doc")
                    .def_property_readonly(
                        "aerodynamic_coefficient_interface",
                        &ta::AtmosphericFlightConditions::
                            getAerodynamicCoefficientInterface,
                        R"doc(Object extracted from the same Body object as this
:py:class:`~AtmosphericFlightConditions` object, which defines
the aerodynamic coefficients.

	)doc");


                /*!
                 **************   EPHEMERIDES  ******************
                 */

                py::class_<te::Ephemeris, std::shared_ptr<te::Ephemeris>>(
                    m, "Ephemeris", "")
                    .def("cartesian_state", &te::Ephemeris::getCartesianState,
                         py::arg("current_time"), "")
                    .def("cartesian_position",
                         &te::Ephemeris::getCartesianPosition,
                         py::arg("current_time"), "")
                    .def("cartesian_velocity",
                         &te::Ephemeris::getCartesianVelocity,
                         py::arg("current_time"), "")
                    .def_property_readonly(
                        "frame_origin", &te::Ephemeris::getReferenceFrameOrigin,
                        "")
                    .def_property_readonly(
                        "frame_orientation",
                        &te::Ephemeris::getReferenceFrameOrientation, "");


                py::class_<te::ConstantEphemeris,
                           std::shared_ptr<te::ConstantEphemeris>,
                           te::Ephemeris>(m, "ConstantEphemeris",
                                          "")
                    .def(
                        py::init<
                            const std::function<
                                Eigen::
                                    Vector6d()>,  //<pybind11/functional.h>,<pybind11/eigen.h>
                            const std::string &, const std::string &>(),
                        py::arg("constant_state_function"),
                        py::arg("reference_frame_origin") = "SSB",
                        py::arg("reference_frame_orientation") = "ECLIPJ2000")
                    .def(py::init<const Eigen::Vector6d,  //<pybind11/eigen.h>
                                  const std::string &, const std::string &>(),
                         py::arg("constant_state"),
                         py::arg("reference_frame_origin") = "SSB",
                         py::arg("reference_frame_orientation") = "ECLIPJ2000")
                    .def("update_constant_state",
                         &te::ConstantEphemeris::updateConstantState,
                         py::arg("new_state"), "");


                py::class_<te::KeplerEphemeris,
                           std::shared_ptr<te::KeplerEphemeris>, te::Ephemeris>(
                    m, "KeplerEphemeris");


                py::class_<te::MultiArcEphemeris,
                           std::shared_ptr<te::MultiArcEphemeris>,
                           te::Ephemeris>(m, "MultiArcEphemeris")
                    .def(py::init<const std::map<
                                      double, std::shared_ptr<te::Ephemeris>> &,
                                  const std::string &, const std::string &>(),
                         py::arg("single_arc_ephemerides"),
                         py::arg("reference_frame_origin") = "SSB",
                         py::arg("reference_frame_orientation") = "ECLIPJ2000");


                py::class_<te::TabulatedCartesianEphemeris<double, double>,
                           std::shared_ptr<
                               te::TabulatedCartesianEphemeris<double, double>>,
                           te::Ephemeris>(m, "TabulatedEphemeris")
                    .def_property(
                        "interpolator",
                        &te::TabulatedCartesianEphemeris<
                            double, double>::getDynamicVectorSizeInterpolator,
                        py::overload_cast<const std::shared_ptr<
                            ti::OneDimensionalInterpolator<double,
                                                           Eigen::VectorXd>>>(
                            &te::TabulatedCartesianEphemeris<
                                double, double>::resetInterpolator));


                py::class_<te::Tle, std::shared_ptr<te::Tle>>(m, "Tle")
                    .def(py::init<  // ctor 1
                             const std::string &>(),
                         py::arg("lines"))
                    .def(py::init<  // ctor 2
                             const std::string &, const std::string &>(),
                         py::arg("line_1"), py::arg("line_2"))
                    .def("get_epoch", &te::Tle::getEpoch)
                    .def("get_b_star", &te::Tle::getBStar)
                    .def("get_epoch", &te::Tle::getEpoch)
                    .def("get_inclination", &te::Tle::getInclination)
                    .def("get_right_ascension", &te::Tle::getRightAscension)
                    .def("get_eccentricity", &te::Tle::getEccentricity)
                    .def("get_arg_of_perigee", &te::Tle::getArgOfPerigee)
                    .def("get_mean_anomaly", &te::Tle::getMeanAnomaly)
                    .def("get_mean_motion", &te::Tle::getMeanMotion);

                py::class_<te::TleEphemeris, std::shared_ptr<te::TleEphemeris>,
                           te::Ephemeris>(m, "TleEphemeris")
                    .def(py::init<const std::string &, const std::string &,
                                  const std::shared_ptr<te::Tle>, const bool>(),
                         py::arg("frame_origin") = "Earth",
                         py::arg("frame_orientation") = "J2000",
                         py::arg("tle") = nullptr, py::arg("use_sdp") = false);

                m.def(
                    "create_ground_station_ephemeris",
                    py::overload_cast<const std::shared_ptr<tss::Body>,
                                      const std::string &>(
                        &tss::createReferencePointEphemeris<TIME_TYPE, double>),
                    "body_with_ground_station", "station_name");

                /*!
                 **************   ROTATION MODELS  ******************
                 */


                py::class_<te::RotationalEphemeris,
                           std::shared_ptr<te::RotationalEphemeris>>(
                    m, "RotationalEphemeris", "")
                    .def("body_fixed_to_inertial_rotation",
                         &te::RotationalEphemeris::getRotationMatrixToBaseFrame,
                         py::arg("time"), "")
                    .def("time_derivative_body_fixed_to_inertial_rotation",
                         &te::RotationalEphemeris::
                             getDerivativeOfRotationToBaseFrame,
                         py::arg("time"), "")
                    .def("inertial_to_body_fixed_rotation",
                         &te::RotationalEphemeris::
                             getRotationMatrixToTargetFrame,
                         py::arg("time"), "")
                    .def("time_derivative_inertial_to_body_fixed_rotation",
                         &te::RotationalEphemeris::
                             getDerivativeOfRotationToTargetFrame,
                         py::arg("time"), "")
                    .def("angular_velocity_in_body_fixed_frame",
                         &te::RotationalEphemeris::
                             getRotationalVelocityVectorInTargetFrame,
                         py::arg("time"), "")
                    .def("angular_velocity_in_inertial_frame",
                         &te::RotationalEphemeris::
                             getRotationalVelocityVectorInBaseFrame,
                         py::arg("time"), "")
                    .def_property_readonly(
                        "body_fixed_frame_name",
                        &te::RotationalEphemeris::getTargetFrameOrientation, "")
                    .def_property_readonly(
                        "inertial_frame_name",
                        &te::RotationalEphemeris::getBaseFrameOrientation, "");

                m.def("transform_to_inertial_orientation",
                      &te::transformStateToInertialOrientation<double, double>,
                      py::arg("state_in_body_fixed_frame"),
                      py::arg("current_time"), py::arg("rotational_ephemeris"));


                py::class_<te::LongitudeLibrationCalculator,
                           std::shared_ptr<te::LongitudeLibrationCalculator>>(
                    m, "LongitudeLibrationCalculator");

                py::class_<
                    te::DirectLongitudeLibrationCalculator,
                    std::shared_ptr<te::DirectLongitudeLibrationCalculator>,
                    te::LongitudeLibrationCalculator>(
                    m, "DirectLongitudeLibrationCalculator")
                    .def(py::init<const double>(),
                         py::arg("scaled_libration_amplitude"));


                py::class_<te::SynchronousRotationalEphemeris,
                           std::shared_ptr<te::SynchronousRotationalEphemeris>,
                           te::RotationalEphemeris>(
                    m, "SynchronousRotationalEphemeris")
                    .def_property("libration_calculator",
                                  &te::SynchronousRotationalEphemeris::
                                      getLongitudeLibrationCalculator,
                                  &te::SynchronousRotationalEphemeris::
                                      setLibrationCalculation);

                py::class_<
                    te::AerodynamicAngleRotationalEphemeris,
                    std::shared_ptr<te::AerodynamicAngleRotationalEphemeris>,
                    te::RotationalEphemeris>(
                    m, "AerodynamicAngleRotationalEphemeris")
                    .def("reset_aerodynamic_angle_function",
                         &te::AerodynamicAngleRotationalEphemeris::
                             setAerodynamicAngleFunction);

                py::class_<te::GcrsToItrsRotationModel,
                           std::shared_ptr<te::GcrsToItrsRotationModel>,
                           te::RotationalEphemeris>(m,
                                                    "GcrsToItrsRotationModel");


                py::class_<
                    te::DirectionBasedRotationalEphemeris,
                    std::shared_ptr<te::DirectionBasedRotationalEphemeris>,
                    te::RotationalEphemeris>(
                    m, "CustomInertialDirectionBasedRotationalEphemeris")
                    .def_property_readonly(
                        "inertial_body_axis_calculator",
                        &te::DirectionBasedRotationalEphemeris::
                            getInertialBodyAxisDirectionCalculator);

                py::class_<
                    te::InertialBodyFixedDirectionCalculator,
                    std::shared_ptr<te::InertialBodyFixedDirectionCalculator>>(
                    m, "InertialBodyFixedDirectionCalculator");

                py::class_<
                    te::CustomBodyFixedDirectionCalculator,
                    std::shared_ptr<te::CustomBodyFixedDirectionCalculator>,
                    te::InertialBodyFixedDirectionCalculator>(
                    m, "CustomBodyFixedDirectionCalculator")
                    .def_property("inertial_body_axis_direction_function",
                                  &te::CustomBodyFixedDirectionCalculator::
                                      getInertialBodyAxisDirectionFunction,
                                  &te::CustomBodyFixedDirectionCalculator::
                                      resetInertialBodyAxisDirectionFunction);


                /*!
                 **************   GRAVITY FIELD  ******************
                 */

                py::class_<tg::GravityFieldModel,
                           std::shared_ptr<tg::GravityFieldModel>>(
                    m, "GravityFieldModel")
                    .def(py::init<const double, const std::function<void()>>(),
                         py::arg("gravitational_parameter"),
                         py::arg("update_inertia_tensor") =
                             std::function<void()>()  // <pybind11/functional.h>
                         )
                    .def("get_gravitational_parameter",
                         &tg::GravityFieldModel::getGravitationalParameter)
                    .def_property(
                        "gravitational_parameter",
                        &tg::GravityFieldModel::getGravitationalParameter,
                        &tg::GravityFieldModel::resetGravitationalParameter);

                py::class_<tg::SphericalHarmonicsGravityField,
                           std::shared_ptr<tg::SphericalHarmonicsGravityField>,
                           tg::GravityFieldModel>(
                    m, "SphericalHarmonicsGravityField")
                    .def_property_readonly(
                        "reference_radius",
                        &tg::SphericalHarmonicsGravityField::getReferenceRadius)
                    .def_property_readonly("maximum_degree",
                                           &tg::SphericalHarmonicsGravityField::
                                               getDegreeOfExpansion)
                    .def_property_readonly("maximum_order",
                                           &tg::SphericalHarmonicsGravityField::
                                               getOrderOfExpansion)
                    .def_property("cosine_coefficients",
                                  &tg::SphericalHarmonicsGravityField::
                                      getCosineCoefficients,
                                  &tg::SphericalHarmonicsGravityField::
                                      setCosineCoefficients)
                    .def_property("sine_coefficients",
                                  &tg::SphericalHarmonicsGravityField::
                                      getSineCoefficients,
                                  &tg::SphericalHarmonicsGravityField::
                                      setSineCoefficients);

                py::class_<tg::PolyhedronGravityField,
                           std::shared_ptr<tg::PolyhedronGravityField>,
                           tg::GravityFieldModel>(m, "PolyhedronGravityField")
                    .def_property_readonly(
                        "volume", &tg::PolyhedronGravityField::getVolume)
                    .def_property_readonly(
                        "vertices_coordinates",
                        &tg::PolyhedronGravityField::getVerticesCoordinates)
                    .def_property_readonly("vertices_defining_each_facet",
                                           &tg::PolyhedronGravityField::
                                               getVerticesDefiningEachFacet);

                /*!
                 **************   SHAPE MODELS  ******************
                 */

                py::class_<tba::BodyShapeModel,
                           std::shared_ptr<tba::BodyShapeModel>>(m,
                                                                 "ShapeModel")
                    .def("get_average_radius",
                         &tba::BodyShapeModel::getAverageRadius)
                    .def_property_readonly(
                        "average_radius",
                        &tba::BodyShapeModel::getAverageRadius);


                /*!
                 **************   GROUND STATION FUNCTIONALITY
                 *******************
                 */


                py::class_<tgs::GroundStationState,
                           std::shared_ptr<tgs::GroundStationState>>(
                    m, "GroundStationState")
                    .def("get_cartesian_state",
                         &tgs::GroundStationState::getCartesianStateInTime,
                         py::arg("seconds_since_epoch"),
                         py::arg("target_frame_origin"))
                    .def("get_cartesian_position",
                         &tgs::GroundStationState::getCartesianPositionInTime,
                         py::arg("seconds_since_epoch"),
                         py::arg("target_frame_origin"))
                    .def_property_readonly(
                        "cartesian_positon_at_reference_epoch",
                        &tgs::GroundStationState::getNominalCartesianPosition)
                    .def_property_readonly(
                        "spherical_positon_at_reference_epoch",
                        &tgs::GroundStationState::getNominalSphericalPosition)
                    .def_property_readonly(
                        "geodetic_positon_at_reference_epoch",
                        &tgs::GroundStationState::getNominalGeodeticPosition)
                    .def_property_readonly(
                        "rotation_matrix_body_fixed_to_topocentric",
                        &tgs::GroundStationState::
                            getRotationMatrixFromBodyFixedToTopocentricFrame);

                py::class_<tgs::GroundStation,
                           std::shared_ptr<tgs::GroundStation>>(m,
                                                                "GroundStation")
                    .def(
                        "set_transmitting_frequency_calculator",
                        &tgs::GroundStation::setTransmittingFrequencyCalculator,
                        py::arg("transmitting_frequency_calculator"))
                    .def("set_water_vapor_partial_pressure_function",
                         &tgs::GroundStation::
                             setWaterVaporPartialPressureFunction,
                         py::arg("water_vapor_partial_pressure_function"))
                    .def("set_temperature_function",
                         &tgs::GroundStation::setTemperatureFunction,
                         py::arg("temperature_function"))
                    .def("set_pressure_function",
                         &tgs::GroundStation::setPressureFunction,
                         py::arg("pressure_function"))
                    .def("set_relative_humidity_function",
                         &tgs::GroundStation::setRelativeHumidityFunction,
                         py::arg("relative_humidity_function"))
                    .def_property_readonly(
                        "temperature_function",
                        &tgs::GroundStation::getTemperatureFunction)
                    .def_property_readonly(
                        "pressure_function",
                        &tgs::GroundStation::getPressureFunction)
                    .def_property_readonly(
                        "relative_humidity_function",
                        &tgs::GroundStation::getRelativeHumidityFunction)
                    .def_property_readonly(
                        "pointing_angles_calculator",
                        &tgs::GroundStation::getPointingAnglesCalculator)
                    .def_property_readonly(
                        "station_state",
                        &tgs::GroundStation::getNominalStationState);


                py::class_<tgs::StationFrequencyInterpolator,
                           std::shared_ptr<tgs::StationFrequencyInterpolator>>(
                    m, "StationFrequencyInterpolator", "");

                py::class_<tgs::ConstantFrequencyInterpolator,
                           std::shared_ptr<tgs::ConstantFrequencyInterpolator>,
                           tgs::StationFrequencyInterpolator>(
                    m, "ConstantFrequencyInterpolator")
                    .def(py::init<double>(), py::arg("frequency"));


                py::class_<tgs::PointingAnglesCalculator,
                           std::shared_ptr<tgs::PointingAnglesCalculator>>(
                    m, "PointingAnglesCalculator")
                    .def("calculate_elevation_angle",
                         py::overload_cast<const Eigen::Vector3d &,
                                           const double>(
                             &tgs::PointingAnglesCalculator::
                                 calculateElevationAngleFromInertialVector),
                         py::arg("inertial_vector_to_target"), py::arg("time"))
                    .def("calculate_azimuth_angle",
                         py::overload_cast<const Eigen::Vector3d &,
                                           const double>(
                             &tgs::PointingAnglesCalculator::
                                 calculateAzimuthAngleFromInertialVector),
                         py::arg("inertial_vector_to_target"), py::arg("time"))
                    .def("convert_inertial_vector_to_topocentric",
                         &tgs::PointingAnglesCalculator::
                             convertVectorFromInertialToTopocentricFrame,
                         py::arg("inertial_vector"), py::arg("time"));


                /*!
                 **************   BODY OBJECTS AND ASSOCIATED FUNCTIONALITY
                 *******************
                 */

                py::class_<tss::Body, std::shared_ptr<tss::Body>>(
                    m, "Body",
                    R"doc(Object that stores the environment properties and current state of
a single body.


	Object that stores the environment properties and current state
	of a single celestial body (natural or artificial). Each separate
	environment model (gravity field, ephemeris, etc.) is stored as a
	member object in this class. During each time step, the Body gets
	updated to the current time/propagated state, and the current
	properties, in as much as they are time-dependent, can be
	extracted from this object

)doc")
                    .def_property("ephemeris_frame_to_base_frame",
                                  &tss::Body::getEphemerisFrameToBaseFrame,
                                  &tss::Body::setEphemerisFrameToBaseFrame)
                    .def_property_readonly(
                        "state", &tss::Body::getState,
                        R"doc(The translational state of the Body, as set during the current
step of the numerical propagation. The translational state
stored here is always in Cartesian elements, w.r.t. the global
frame origin, with axes along the global frame orientation. If
the body's translational state is numerically propagated, this
property gets extracted from the propagated state vector. If it
is not propagated, the state is extracted from this body's
ephemeris. In both cases, any required state transformations
are automatically applied. Note that this function  is *only*
valid during the numerical propagation if any aspects of the
dynamics or dependent variables require the body's state.

	)doc")
                    .def_property_readonly(
                        "position", &tss::Body::getPosition,
                        R"doc(The translational position of the Body, as set during the
current step of the numerical propagation
(see :py:attr:`~state`).

	)doc")
                    .def_property_readonly(
                        "velocity", &tss::Body::getVelocity,
                        R"doc(The translational velocity of the Body, as set during the
current step of the numerical propagation
(see :py:attr:`~state`).

	)doc")
                    .def_property_readonly(
                        "inertial_to_body_fixed_frame",
                        &tss::Body::getCurrentRotationMatrixToLocalFrame,
                        R"doc(The rotation from inertial frame (with global frame
orientation) to this Body's body-fixed frame. The rotation is
always returned here as a rotation matrix.  If the body's
rotational state is numerically propagated, this property gets
extracted from the propagated state vector. If it is not
propagated, the state is extracted from this body's rotational
ephemeris.

.. note:: This function is **only** valid during the
          numerical propagation if any aspects of the dynamics
          or dependent variables require the body's rotational
          state.

	)doc")
                    .def_property_readonly(
                        "body_fixed_to_inertial_frame",
                        &tss::Body::getCurrentRotationMatrixToGlobalFrame,
                        R"doc(The rotation from this Body's body-fixed frame to inertial
frame (see :py:attr:`~inertial_to_body_fixed_frame`).

	)doc")
                    .def_property_readonly(
                        "inertial_to_body_fixed_frame_derivative",
                        &tss::Body::
                            getCurrentRotationMatrixDerivativeToLocalFrame,
                        R"doc(Time derivative of rotation matrix from inertial frame to this
Body's body-fixed frame
(see :py:attr:`~inertial_to_body_fixed_frame`).

	)doc")
                    .def_property_readonly(
                        "body_fixed_to_inertial_frame_derivative",
                        &tss::Body::
                            getCurrentRotationMatrixDerivativeToGlobalFrame,
                        R"doc(Time derivative of rotation matrix from this Body's body-fixed
frame to inertial frame
(see :py:attr:`~inertial_to_body_fixed_frame`).

	)doc")
                    .def_property_readonly(
                        "inertial_angular_velocity",
                        &tss::Body::
                            getCurrentAngularVelocityVectorInGlobalFrame,
                        R"doc(Angular velocity vector of the body, expressed in inertial
frame (see :py:attr:`~inertial_to_body_fixed_frame`).

	)doc")
                    .def_property_readonly(
                        "body_fixed_angular_velocity",
                        &tss::Body::getCurrentAngularVelocityVectorInLocalFrame,
                        R"doc(Angular velocity vector of the body, expressed in body-fixed
frame (see :py:attr:`~inertial_to_body_fixed_frame`).

	)doc")
                    .def_property(
                        "mass", &tss::Body::getBodyMass,
                        &tss::Body::setConstantBodyMass,
                        R"doc(Denotes the current mass of the vehicle, as used in the calculation of
non-conservative acceleration. Note that defining a mass for a vehicle
does *not* define a gravity field (this is done through a gravity field model).
However, defining a gravity field model in a body automatically assigns it
a mass here.
Unlike the attributes containing the state, orientation, angular velocity
of the Body, this attribute may be used to retrieve the state during the
propagation *and* to define the mass of a vehicle

	)doc")
                    .def("set_constant_mass", &tss::Body::setConstantBodyMass,
                         py::arg("mass"))
                    .def_property("inertia_tensor",
                                  &tss::Body::getBodyInertiaTensor,
                                  py::overload_cast<const Eigen::Matrix3d &>(
                                      &tss::Body::setBodyInertiaTensor),
                                  "")
                    .def("state_in_base_frame_from_ephemeris",
                         &tss::Body::getStateInBaseFrameFromEphemeris<double,
                                                                      double>,
                         py::arg("time"))
                    .def_property_readonly(
                        "ephemeris", &tss::Body::getEphemeris,
                        R"doc(Ephemeris model of this body, used to calculate its current state as a function of time.
Depending on the selected type of model, the type of this attribute
is of type Ephemeris, or a derived class thereof.

	)doc")
                    .def_property(
                        "atmosphere_model", &tss::Body::getAtmosphereModel,
                        &tss::Body::setAtmosphereModel,
                        R"doc(Atmosphere model of this body, used to calculate density, temperature, etc. at a given
state/time. Depending on the selected type of model, the type of this attribute
is of type AtmosphereModel, or a derived class thereof.

	)doc")
                    .def_property(
                        "shape_model", &tss::Body::getShapeModel,
                        &tss::Body::setShapeModel,
                        R"doc(Shape model of this body, used to define the exterior shape of the body, for instance for
the calculation of vehicle's altitude. Depending on the selected type of model, the type of this attribute
is of type BodyShapeModel, or a derived class thereof.

	)doc")
                    .def_property(
                        "gravity_field_model", &tss::Body::getGravityFieldModel,
                        &tss::Body::setGravityFieldModel,
                        R"doc(Gravity field model of this body, used to define the exterior gravitational potential, and
its gradient(s). Depending on the selected type of model, the type of this attribute
is of type GravityFieldModel, or a derived class thereof.

	)doc")
                    .def_property(
                        "aerodynamic_coefficient_interface",
                        &tss::Body::getAerodynamicCoefficientInterface,
                        &tss::Body::setAerodynamicCoefficientInterface,
                        R"doc(Object defining the aerodynamic coefficients of a vehicle (force-only, or force and moment)
as a function of any number of independent variables. Depending on the selected type of model, the type of this attribute
is of type AerodynamicCoefficientInterface, or a derived class thereof.

	)doc")
                    .def_property(
                        "flight_conditions", &tss::Body::getFlightConditions,
                        &tss::Body::setFlightConditions,
                        R"doc(Current flight conditions of the body, which can be accessed during the propagation to get the current altitude, aerodynamic angle calculator, longitude, etc.

	)doc")
                    .def_property(
                        "rotation_model", &tss::Body::getRotationalEphemeris,
                        &tss::Body::setRotationalEphemeris,
                        R"doc(Object defining the orientation of a body, used to calculate the rotation to/from a body-fixed
frame (and its derivate). Depending on the selected type of model, the type of this attribute
is of type RotationalEphemeris, or a derived class thereof.

	)doc")
                    .def_property("system_models",
                                  &tss::Body::getVehicleSystems,
                                  &tss::Body::setVehicleSystems, "")
                    .def_property("rigid_body_properties",
                                  &tss::Body::getMassProperties,
                                  &tss::Body::setMassProperties, "")
                    .def_property_readonly(
                        "gravitational_parameter",
                        &tss::Body::getGravitationalParameter,
                        R"doc(Attribute of convenience, equivalent to ``.gravity_field_model.gravitational_parameter``

	)doc")
                    .def("get_ground_station", &tss::Body::getGroundStation,
                         py::arg("station_name"), "")
                    .def_property_readonly("ground_station_list",
                                           &tss::Body::getGroundStationMap, "");


                py::class_<tss::SystemOfBodies,
                           std::shared_ptr<tss::SystemOfBodies>>(
                    m, "SystemOfBodies",
                    R"doc(Object that contains a set of Body objects and associated frame
information.


	Object that contains a set of Body objects and associated frame
	information. This object stored the entire environment for a
	typical Tudat numerical simulation, and is fundamental for the
	overall Tudat architecture.

)doc")
                    .def(
                        "get", &tss::SystemOfBodies::getBody,
                        py::arg("body_name"),
                        R"doc(This function extracts a single Body object from the SystemOfBodies.

	:param body_name:
		Name of the Body that is to be retrieved.

	:return:
		Body object of the requested name
)doc")
                    .def("get_body", &tss::SystemOfBodies::getBody,
                         py::arg("body_name"),
                         R"doc(Deprecated version of :py:func:`~get`
)doc")
                    .def("create_empty_body",
                         &tss::SystemOfBodies::createEmptyBody<double,
                                                               TIME_TYPE>,
                         py::arg("body_name"), py::arg("process_body") = 1,
                         R"doc(This function creates a new empty body.

	This function creates a new empty body, and adds it to the
	:py:class:`~SystemOfBodies`. Since the body is empty, it will
	not have any environment models defined. These must all be
	added manually by a user.


	:param body_name:
		Name of the Body that is to be added

	:param process_body:
		Variable that defines whether this new Body will have its
		global frame origin/orientation set to conform to rest of
		the environment.

		.. warning:: Only in very rare cases should
		             this variable be anything other than ``True``.
		             Users are recommended to keep this default value
		             intact.

)doc")
                    .def("does_body_exist", &tss::SystemOfBodies::doesBodyExist,
                         py::arg("body_name"), "")
                    .def("list_of_bodies",
                         &tss::SystemOfBodies::getListOfBodies, "")
                    //            .def("get_body_dict",
                    //            &tss::SystemOfBodies::getMap,
                    .def(
                        "add_body",
                        &tss::SystemOfBodies::addBody<double, TIME_TYPE>,
                        py::arg("body_to_add"), py::arg("body_name"),
                        py::arg("process_body") = 1,
                        R"doc(This function adds an existing body, which the user has
separately created, to the :py:class:`~SystemOfBodies`.


	:param body_to_add:
		Body object that is to be added.

	:param body_name:
		Name of the Body that is to be added.

	:param process_body:
		Variable that defines whether this new Body will have its
		global frame origin/orientation set to conform to rest of
		the environment.

		.. warning:: Only in very rare cases should this variable be
		             anything other than ``True``. Users are
		             recommended to keep this default value intact.

)doc")
                    .def("remove_body", &tss::SystemOfBodies::deleteBody,
                         py::arg("body_name"),
                         R"doc(This function removes an existing body from the
:py:class:`~SystemOfBodies`.



	.. warning:: This function does *not* necessarily delete the
	             Body object, it only removes it from this object.
	             If any existing models in the simulation refer to
	             this Body, it will persist in memory.


	:param body_name:
		Name of the Body that is to be removed.

)doc")
                    .def("global_frame_orientation",
                         &tss::SystemOfBodies::getFrameOrientation, "")
                    .def("global_frame_origin",
                         &tss::SystemOfBodies::getFrameOrigin, "");

                //            .def_property_readonly("number_of_bodies",
                //            &tss::SystemOfBodies::getNumberOfBodies,
                //                                   );

                /*!
                 **************   SUPPORTING FUNCTIONS USED ENVIRONMENT MODELS
                 *******************
                 */
            }
        }  // namespace environment
    }  // namespace numerical_simulation
}  // namespace tudatpy
