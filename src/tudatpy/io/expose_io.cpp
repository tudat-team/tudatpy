/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define PYBIND11_DETAILED_ERROR_MESSAGES
#include <pybind11/stl.h>
#include <tudat/io/basicInputOutput.h>
#include <tudat/io/missileDatcomData.h>
#include <tudat/io/readHistoryFromFile.h>
#include <tudat/io/readOdfFile.h>
#include <tudat/io/readTabulatedWeatherData.h>
#include <tudat/io/readTrackingTxtFile.h>
#include <tudat/io/solarActivityData.h>

#include "tudatpy/docstrings.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy {

    namespace io {

        PYBIND11_MODULE(expose_io, m) {

            py::module_::import("tudatpy.math.interpolators");
            m.def("get_resource_path", &tudat::paths::get_resource_path,
R"doc(Get the path at which tudat resources are located.
	:return:
		Local path at which tudat resources are located.
)doc");
            m.def("get_ephemeris_path",
                  &tudat::paths::getEphemerisDataFilesPath,
R"doc(Get the path at which the ephemeris used by tudat are located.
	:return:
		Local path at which the tudat ephemeris resources are located.
)doc");
            m.def("get_earth_orientation_path",
                  &tudat::paths::getEarthOrientationDataFilesPath,
R"doc(Get the path at which the Earth orientation resources used by tudat are located.
	:return:
		Local path at which tudat Earth orientation resources are located.
)doc");
            m.def("get_quadrature_path", &tudat::paths::getQuadratureDataPath,
R"doc(Get the path at which the Gaussian quadrature resources are located.
	:return:
		Local path at which tudat Gaussian quadrature resources are located.
)doc");
            m.def("get_spice_kernel_path", &tudat::paths::getSpiceKernelPath,
R"doc(Get the path at which the SPICE kernel used by tudat is located.
	:return:
		Local path at which the SPICE kernel is located.
)doc");
            m.def("get_atmosphere_tables_path",
                  &tudat::paths::getAtmosphereTablesPath,
R"doc(Get the path at which tudat atmosphere tables are located.
	:return:
		Local path at which tudat atmosphere tables are located.
)doc");
            m.def("get_gravity_models_path",
                  &tudat::paths::getGravityModelsPath,
R"doc(Get the path at which tudat gravity models are located.
	:return:
		Local path at which tudat gravity models are located.
)doc");
            m.def("get_space_weather_path",
                  &tudat::paths::getSpaceWeatherDataPath,
R"doc(Get the path at which tudat space weather is located.
	:return:
		Local path at which tudat space weather is located.
)doc");

            m.def(
                "read_vector_history_from_file",
                &tudat::input_output::readVectorHistoryFromFile<double, double>,
                py::arg("vector_size"), py::arg("file_name"),
R"doc(Read a vector history from a file.

	:param vector_size:
		Size of the vector at each epoch.
	:param file_name:
		Name of the file containing the vector history.
	:return:
		Dictionary mapping epochs to the vector at the given epoch.
)doc");

            m.def(
                "read_matrix_history_from_file",
                &tudat::input_output::readMatrixHistoryFromFile<double, double>,
                py::arg("matrix_rows"), py::arg("matrix_columns"),
                py::arg("file_name"),
R"doc(Read a matrix history from a file.

	:param matrix_rows:
		Number of rows in the matrix at each epoch.
	:param matrix_columns:
		Number of columns in the matrix at each epoch.
	:param file_name:
		Name of the file containing the matrix history.
	:return:
		Dictionary mapping epochs to the matrix at the given epoch.
)doc");

            py::class_<tudat::input_output::MissileDatcomData,
                       std::shared_ptr<tudat::input_output::MissileDatcomData>>(
                m, "missile_DATCOM_data",
R"doc(Class containing data and methods interfacing the Missile DATCOM software.

	This class is the main method that can be used to interface tudat with the Missile DATCOM software.
	It can be initialised with the output file from Missile DATCOM, and provides methods to convert these results
	into tudat-compatible data.

	.. note:: The Missile DATCOM software from which outputs can be interfaced to TUDAT is an entirely separate software from Tudat(Py).
	          Please refer to Missile DATCOM user manuals for information on how to use it. These can be accessed on the US Defence Technical
	          Information Center at accession numbers `ADA267447 <https://apps.dtic.mil/sti/citations/ADA267447>`_ and
	          `ADA503576 <https://apps.dtic.mil/sti/citations/ADA503576>`_.

	.. note:: The interfacing of Missile DATCOM to tudat assumes that aerodynamic coefficients are computed as a function of both
	          Mach number and angle of attack.

)doc")
                .def(py::init<const std::string&>(),
                     py::arg("file_name_and_path"),
R"doc(Class constructor.

	Function used to construct and initialise the class. In essence, it can be used to read and extract the aerodynamic coefficients
	computed by Missile DATCOM, and save them in different formats.


	:param file_name_and_path:
		Full path and file name of the `for004.dat` Missile DATCOM results output file.
)doc")
                .def("get_static_coefficient",
                     &tudat::input_output::MissileDatcomData::
                         getStaticCoefficient,
                     py::arg("mach_index"), py::arg("angle_of_attack_index"),
                     py::arg("coefficient_index"),
R"doc(Get a specific static coefficient from the result database.

	:param mach_index:
		Index of the Mach number for which to get the static coefficient.
	:param angle_of_attack_index:
		Index of the angle of attack for which to get the static coefficient.
	:param coefficient_index:
		Type of the static aerodynamic coefficient.
	:return:
		Static aerodynamic coefficient.
)doc")
                .def(
                    "get_dynamic_coefficient",
                    &tudat::input_output::MissileDatcomData::
                        getDynamicCoefficient,
                    py::arg("mach_index"), py::arg("angle_of_attack_index"),
                    py::arg("coefficient_index"),
R"doc(Get a specific dynamic coefficient from the result database.

	:param mach_index:
		Index of the Mach number for which to get the static coefficient.
	:param angle_of_attack_index:
		Index of the angle of attack for which to get the static coefficient.
	:param coefficient_index:
		Type of the dynamic aerodynamic coefficient.
	:return:
		Dynamic aerodynamic coefficient.
)doc")
                .def("get_angle_of_attacks",
                     &tudat::input_output::MissileDatcomData::getAngleOfAttacks,
R"doc(Get the list of angle of attacks at which Missile DATCOM has been run.
	:return:
		List of angle of attacks.
)doc")
                .def("get_mach_numbers",
                     &tudat::input_output::MissileDatcomData::getMachNumbers,
R"doc(Get the list of Mach numbers at which Missile DATCOM has been run.
	:return:
		List of Mach numbers.
)doc")
                .def(
                    "get_Reynolds_numbers",
                    &tudat::input_output::MissileDatcomData::getReynoldsNumbers,
R"doc(Get the list of Reynolds numbers at which Missile DATCOM has been run.
	:return:
		List of Reynolds numbers.
)doc")
                .def("write_all_coefficients_to_files",
                     &tudat::input_output::MissileDatcomData::
                         writeAllCoefficientsToFiles,
                     py::arg("file_name_base"), py::arg("base_precision") = 15,
                     py::arg("exponent_width") = 2,
R"doc(Write all the aerodynamic coefficients to CSV files.

	:param file_name_base:
		Full base path and name of the file that will be saved. The name of each aerodynamic coefficient will be included at the end of the file name.
	:param base_precision:
		Number of digits to represent the base of the floating-point number.
	:param exponent_width:
		Number of digits to represent the exponent of the floating-point number.
)doc")
                .def("write_force_and_moment_coefficients_to_files",
                     &tudat::input_output::MissileDatcomData::
                         writeForceAndMomentCoefficientsToFiles,
                     py::arg("file_name_base"), py::arg("base_precision") = 15,
                     py::arg("exponent_width") = 2,
R"doc(Write the force and moment coefficients to a file in the format taken by the :func:`~tudatpy.numerical_simulation.environment_setup.aerodynamic_coefficients.tabulated_from_files` function.

	:param file_name_base:
		Full base path and name of the file that will be saved. The name of each aerodynamic coefficient will be included at the end of the file name.
	:param base_precision:
		Number of digits to represent the base of the floating-point number.
	:param exponent_width:
		Number of digits to represent the exponent of the floating-point number.
)doc");

            py::enum_<tudat::input_output::MissileDatcomData::
                          DynamicCoefficientNames>(
                m, "DynamicCoefficientNames",
R"doc(Enumeration of Missile DATCOM dynamic aerodynamic coefficient types.


	:member cnq:
	:member cmq:
	:member caq:
	:member cyq:
	:member clnq:
	:member cllq:
	:member cnr:
	:member cmr:
	:member car:
	:member cyr:
	:member clnr:
	:member cllr:
	:member cnp:
	:member cmp:
	:member cap:
	:member cyp:
	:member clnp:
	:member cllp:
	:member cnad:
	:member cmad:
)doc")
                .value("cnq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnq,
get_docstring("DynamicCoefficientNames.cnq").c_str())
                .value("cmq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmq,
get_docstring("DynamicCoefficientNames.cmq").c_str())
                .value("caq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::caq,
get_docstring("DynamicCoefficientNames.caq").c_str())
                .value("cyq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cyq,
get_docstring("DynamicCoefficientNames.cyq").c_str())
                .value("clnq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::clnq,
get_docstring("DynamicCoefficientNames.clnq").c_str())
                .value("cllq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cllq,
get_docstring("DynamicCoefficientNames.cllq").c_str())
                .value("cnr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnr,
get_docstring("DynamicCoefficientNames.cnr").c_str())
                .value("cmr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmr,
get_docstring("DynamicCoefficientNames.cmr").c_str())
                .value("car",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::car,
get_docstring("DynamicCoefficientNames.car").c_str())
                .value("cyr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cyr,
get_docstring("DynamicCoefficientNames.cyr").c_str())
                .value("clnr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::clnr,
get_docstring("DynamicCoefficientNames.clnr").c_str())
                .value("cllr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cllr,
get_docstring("DynamicCoefficientNames.cllr").c_str())
                .value("cnp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnp,
get_docstring("DynamicCoefficientNames.cnp").c_str())
                .value("cmp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmp,
get_docstring("DynamicCoefficientNames.cmp").c_str())
                .value("cap",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cap,
get_docstring("DynamicCoefficientNames.cap").c_str())
                .value("cyp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cyp,
get_docstring("DynamicCoefficientNames.cyp").c_str())
                .value("clnp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::clnp,
get_docstring("DynamicCoefficientNames.clnp").c_str())
                .value("cllp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cllp,
get_docstring("DynamicCoefficientNames.cllp").c_str())
                .value("cnad",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnad,
get_docstring("DynamicCoefficientNames.cnad").c_str())
                .value("cmad",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmad,
get_docstring("DynamicCoefficientNames.cmad").c_str())
                .export_values();

            py::enum_<
                tudat::input_output::MissileDatcomData::StaticCoefficientNames>(
                m, "StaticCoefficientNames",
R"doc(Enumeration of Missile DATCOM static aerodynamic coefficient types.


	:member cn:
	:member cm:
	:member ca:
	:member cy:
	:member cln:
	:member cll:
	:member cna:
	:member cma:
	:member cyb:
	:member cnb:
	:member clb:
)doc")
                .value("cn",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cn,
get_docstring("StaticCoefficientNames.cn").c_str())
                .value("cm",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cm,
get_docstring("StaticCoefficientNames.cm").c_str())
                .value("ca",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::ca,
get_docstring("StaticCoefficientNames.ca").c_str())
                .value("cy",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cy,
get_docstring("StaticCoefficientNames.cy").c_str())
                .value("cln",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cln,
get_docstring("StaticCoefficientNames.cln").c_str())
                .value("cll",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cll,
get_docstring("StaticCoefficientNames.cll").c_str())
                .value("cna",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cna,
get_docstring("StaticCoefficientNames.cna").c_str())
                .value("cma",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cma,
get_docstring("StaticCoefficientNames.cma").c_str())
                .value("cyb",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cyb,
get_docstring("StaticCoefficientNames.cyb").c_str())
                .value("cnb",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cnb,
get_docstring("StaticCoefficientNames.cnb").c_str())
                .value("clb",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::clb,
get_docstring("StaticCoefficientNames.clb").c_str())
                .export_values();

            py::class_<tio::solar_activity::SolarActivityData,
                       std::shared_ptr<tio::solar_activity::SolarActivityData>>(
                m, "SolarActivityData",
get_docstring("SolarActivityData").c_str());

            m.def("read_solar_activity_data",
                  &tio::solar_activity::readSolarActivityData,
                  py::arg("file_path"),
                  R"doc(
Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.

:param file_path: Path to the space weather data file.
)doc");

            py::class_<tio::OdfRawFileContents,
                       std::shared_ptr<tio::OdfRawFileContents>>(
                m, "OdfRawFileContents",
get_docstring("OdfRawFileContents").c_str())
                .def("write_to_text_file",
                     &tio::OdfRawFileContents::writeOdfToTextFile,
                     py::arg("output_file"),
get_docstring("OdfRawFileContents.write_to_text_file").c_str());

            m.def("read_odf_file", &tio::readOdfFile, py::arg("file_name"),
get_docstring("read_odf_file").c_str());

            m.def(
                "set_dsn_weather_data_in_ground_stations",
                py::overload_cast<
                    tudat::simulation_setup::SystemOfBodies&,
                    const std::vector<std::string>&,
                    std::shared_ptr<tudat::interpolators::InterpolatorSettings>,
                    const std::map<int, std::vector<std::string>>&,
                    const std::string&>(
                    &tio::setDsnWeatherDataInGroundStations),
                py::arg("bodies"), py::arg("weather_file_names"),
                py::arg("interpolator_settings") =
                    tudat::interpolators::linearInterpolation(),
                py::arg("ground_stations_per_complex") = tudat::
                    simulation_setup::getDefaultDsnStationNamesPerComplex(),
                py::arg("body_with_ground_stations_name") = "Earth",
get_docstring("set_dsn_weather_data_in_ground_stations").c_str());

            m.def("read_tracking_txt_file", &tio::createTrackingTxtFileContents,
                  py::arg("file_name"), py::arg("column_types"),
                  py::arg("comment_symbol") = '#',
                  py::arg("value_separators") = ",:\t ",
get_docstring("read_tracking_txt_file").c_str());

            py::class_<
                tudat::input_output::TrackingTxtFileContents,
                std::shared_ptr<tudat::input_output::TrackingTxtFileContents>>(
                m, "TrackingTxtFileContents",
get_docstring("TrackingTxtFileContents").c_str())
                .def(py::init<const std::string, const std::vector<std::string>,
                              const char, const std::string>(),
                     py::arg("file_name"), py::arg("column_types"),
                     py::arg("comment_symbol") = '#',
                     py::arg("value_separators") = ",:\t ",
get_docstring("TrackingTxtFileContents").c_str())
                .def_property_readonly(
                    "column_field_types",
                    &tio::TrackingTxtFileContents::getRawColumnTypes,
get_docstring("TrackingTxtFileContents.column_field_types").c_str())
                .def_property_readonly(
                    "double_datamap",
                    &tio::TrackingTxtFileContents::getDoubleDataMap,
get_docstring("TrackingTxtFileContents.double_datamap").c_str())
                .def_property_readonly(
                    "raw_datamap", &tio::TrackingTxtFileContents::getRawDataMap,
get_docstring("TrackingTxtFileContents.raw_datamap").c_str())
                .def_property_readonly(
                    "num_rows", &tio::TrackingTxtFileContents::getNumRows,
get_docstring("TrackingTxtFileContents.num_rows").c_str());

            py::enum_<tudat::input_output::TrackingDataType>(
                m, "TrackingDataType",
get_docstring("TrackingDataType").c_str())
                .value("year", tudat::input_output::TrackingDataType::year,
get_docstring("TrackingDataType.year").c_str())
                .value("month", tudat::input_output::TrackingDataType::month,
get_docstring("TrackingDataType.month").c_str())
                .value("day", tudat::input_output::TrackingDataType::day,
get_docstring("TrackingDataType.day").c_str())
                .value("hour", tudat::input_output::TrackingDataType::hour,
get_docstring("TrackingDataType.hour").c_str())
                .value("minute", tudat::input_output::TrackingDataType::minute,
get_docstring("TrackingDataType.minute").c_str())
                .value("second", tudat::input_output::TrackingDataType::second,
get_docstring("TrackingDataType.second").c_str())
                .value("time_tag_delay",
                       tudat::input_output::TrackingDataType::time_tag_delay,
get_docstring("TrackingDataType.time_tag_delay").c_str())
                .value("observation_time_scale",
                       tudat::input_output::TrackingDataType::
                           observation_time_scale,
get_docstring("TrackingDataType.observation_time_scale").c_str())
                .value("file_name",
                       tudat::input_output::TrackingDataType::file_name,
get_docstring("TrackingDataType.file_name").c_str())
                .value(
                    "n_way_light_time",
                    tudat::input_output::TrackingDataType::n_way_light_time,
get_docstring("TrackingDataType.n_way_light_time").c_str())
                .value("light_time_measurement_delay",
                       tudat::input_output::TrackingDataType::
                           light_time_measurement_delay,
get_docstring("TrackingDataType.light_time_measurement_delay").c_str())
                .value("light_time_measurement_accuracy",
                       tudat::input_output::TrackingDataType::
                           light_time_measurement_accuracy,
get_docstring("TrackingDataType.light_time_measurement_accuracy").c_str())
                .value("dsn_transmitting_station_nr",
                       tudat::input_output::TrackingDataType::
                           dsn_transmitting_station_nr,
get_docstring("TrackingDataType.dsn_transmitting_station_nr").c_str())
                .value(
                    "dsn_receiving_station_nr",
                    tudat::input_output::TrackingDataType::
                        dsn_receiving_station_nr,
get_docstring("TrackingDataType.dsn_receiving_station_nr").c_str())
                .value(
                    "observation_body",
                    tudat::input_output::TrackingDataType::observation_body,
get_docstring("TrackingDataType.observation_body").c_str())
                .value("observed_body",
                       tudat::input_output::TrackingDataType::observed_body,
get_docstring("TrackingDataType.observed_body").c_str())
                .value("spacecraft_id",
                       tudat::input_output::TrackingDataType::spacecraft_id,
get_docstring("TrackingDataType.spacecraft_id").c_str())
                .value("planet_nr",
                       tudat::input_output::TrackingDataType::planet_nr,
get_docstring("TrackingDataType.planet_nr").c_str())
                .value("tdb_time_j2000",
                       tudat::input_output::TrackingDataType::tdb_time_j2000,
get_docstring("TrackingDataType.tdb_time_j2000").c_str())
                .value(
                    "tdb_spacecraft_j2000",
                    tudat::input_output::TrackingDataType::tdb_spacecraft_j2000,
get_docstring("TrackingDataType.tdb_spacecraft_j2000").c_str())
                .value("x_planet_frame",
                       tudat::input_output::TrackingDataType::x_planet_frame,
get_docstring("TrackingDataType.x_planet_frame").c_str())
                .value("y_planet_frame",
                       tudat::input_output::TrackingDataType::y_planet_frame,
get_docstring("TrackingDataType.y_planet_frame").c_str())
                .value("z_planet_frame",
                       tudat::input_output::TrackingDataType::z_planet_frame,
get_docstring("TrackingDataType.z_planet_frame").c_str())
                .value(
                    "vx_planet_frame",
                    tudat::input_output::TrackingDataType::vx_planet_frame,
get_docstring("TrackingDataType.vx_planet_frame").c_str())
                .value(
                    "vy_planet_frame",
                    tudat::input_output::TrackingDataType::vy_planet_frame,
get_docstring("TrackingDataType.vy_planet_frame").c_str())
                .value(
                    "vz_planet_frame",
                    tudat::input_output::TrackingDataType::vz_planet_frame,
get_docstring("TrackingDataType.vz_planet_frame").c_str())
                .value("residual_de405",
                       tudat::input_output::TrackingDataType::residual_de405,
get_docstring("TrackingDataType.residual_de405").c_str())
                .value("spacecraft_transponder_delay",
                       tudat::input_output::TrackingDataType::
                           spacecraft_transponder_delay,
get_docstring("TrackingDataType.spacecraft_transponder_delay").c_str())
                .value(
                    "uplink_frequency",
                    tudat::input_output::TrackingDataType::uplink_frequency,
get_docstring("TrackingDataType.uplink_frequency").c_str())
                .value(
                    "downlink_frequency",
                    tudat::input_output::TrackingDataType::downlink_frequency,
get_docstring("TrackingDataType.downlink_frequency").c_str())
                .value(
                    "signal_to_noise",
                    tudat::input_output::TrackingDataType::signal_to_noise,
get_docstring("TrackingDataType.signal_to_noise").c_str())
                .value("spectral_max",
                       tudat::input_output::TrackingDataType::spectral_max,
get_docstring("TrackingDataType.spectral_max").c_str())
                .value(
                    "doppler_measured_frequency",
                    tudat::input_output::TrackingDataType::
                        doppler_measured_frequency,
get_docstring("TrackingDataType.doppler_measured_frequency").c_str())
                .value("doppler_base_frequency",
                       tudat::input_output::TrackingDataType::
                           doppler_base_frequency,
get_docstring("TrackingDataType.doppler_base_frequency").c_str())
                .value("doppler_noise",
                       tudat::input_output::TrackingDataType::doppler_noise,
get_docstring("TrackingDataType.doppler_noise").c_str())
                .value(
                    "doppler_bandwidth",
                    tudat::input_output::TrackingDataType::doppler_bandwidth,
get_docstring("TrackingDataType.doppler_bandwidth").c_str())
                .value(
                    "vlbi_station_name",
                    tudat::input_output::TrackingDataType::vlbi_station_name,
get_docstring("TrackingDataType.vlbi_station_name").c_str())
                .export_values();
        };

    }  // namespace io
}  // namespace tudatpy
