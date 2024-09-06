/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <pybind11/stl.h>
#include <tudat/io/basicInputOutput.h>
#include <tudat/io/missileDatcomData.h>
#include <tudat/io/readHistoryFromFile.h>
#include <tudat/io/readOdfFile.h>
#include <tudat/io/readTabulatedWeatherData.h>
#include <tudat/io/readTrackingTxtFile.h>
#include <tudat/io/solarActivityData.h>

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

            m.def(
                "get_ephemeris_path", &tudat::paths::getEphemerisDataFilesPath,
                R"doc(Get the path at which the ephemeris used by tudat are located.
	:return:
		Local path at which the tudat ephemeris resources are located.
)doc");
            m.def(
                "get_earth_orientation_path",
                &tudat::paths::getEarthOrientationDataFilesPath,
                R"doc(Get the path at which the Earth orientation resources used by tudat are located.
	:return:
		Local path at which tudat Earth orientation resources are located.
)doc");
            m.def(
                "get_quadrature_path", &tudat::paths::getQuadratureDataPath,
                R"doc(Get the path at which the Gaussian quadrature resources are located.
	:return:
		Local path at which tudat Gaussian quadrature resources are located.
)doc");
            m.def(
                "get_spice_kernel_path", &tudat::paths::getSpiceKernelPath,
                R"doc(Get the path at which the SPICE kernel used by tudat is located.
	:return:
		Local path at which the SPICE kernel is located.
)doc");
            m.def(
                "get_atmosphere_tables_path",
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
                .def(
                    "get_static_coefficient",
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
                .def(
                    "get_angle_of_attacks",
                    &tudat::input_output::MissileDatcomData::getAngleOfAttacks,
                    R"doc(Get the list of angle of attacks at which Missile DATCOM has been run.
	:return:
		List of angle of attacks.
)doc")
                .def(
                    "get_mach_numbers",
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
                .def(
                    "write_force_and_moment_coefficients_to_files",
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
"")
                .value("cmq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmq,
"")
                .value("caq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::caq,
"")
                .value("cyq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cyq,
"")
                .value("clnq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::clnq,
"")
                .value("cllq",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cllq,
"")
                .value("cnr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnr,
"")
                .value("cmr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmr,
"")
                .value("car",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::car,
"")
                .value("cyr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cyr,
"")
                .value("clnr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::clnr,
"")
                .value("cllr",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cllr,
"")
                .value("cnp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnp,
"")
                .value("cmp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmp,
"")
                .value("cap",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cap,
"")
                .value("cyp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cyp,
"")
                .value("clnp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::clnp,
"")
                .value("cllp",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cllp,
"")
                .value("cnad",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cnad,
"")
                .value("cmad",
                       tudat::input_output::MissileDatcomData::
                           DynamicCoefficientNames::cmad,
"")
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
"")
                .value("cm",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cm,
"")
                .value("ca",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::ca,
"")
                .value("cy",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cy,
"")
                .value("cln",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cln,
"")
                .value("cll",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cll,
"")
                .value("cna",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cna,
"")
                .value("cma",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cma,
"")
                .value("cyb",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cyb,
"")
                .value("cnb",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::cnb,
"")
                .value("clb",
                       tudat::input_output::MissileDatcomData::
                           StaticCoefficientNames::clb,
"")
                .export_values();

            py::class_<tio::solar_activity::SolarActivityData,
                       std::shared_ptr<tio::solar_activity::SolarActivityData>>(
                m, "SolarActivityData",
"");

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
"")
                .def("write_to_text_file",
                     &tio::OdfRawFileContents::writeOdfToTextFile,
                     py::arg("output_file"),
"");

            m.def("read_odf_file", &tio::readOdfFile, py::arg("file_name"),
"");

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
"");

            m.def("read_tracking_txt_file", &tio::createTrackingTxtFileContents,
                  py::arg("file_name"), py::arg("column_types"),
                  py::arg("comment_symbol") = '#',
                  py::arg("value_separators") = ",:\t ",
"");

            py::class_<
                tudat::input_output::TrackingTxtFileContents,
                std::shared_ptr<tudat::input_output::TrackingTxtFileContents>>(
                m, "TrackingTxtFileContents",
"")
                .def(py::init<const std::string, const std::vector<std::string>,
                              const char, const std::string>(),
                     py::arg("file_name"), py::arg("column_types"),
                     py::arg("comment_symbol") = '#',
                     py::arg("value_separators") = ",:\t ",
"")
                .def_property_readonly(
                    "column_field_types",
                    &tio::TrackingTxtFileContents::getRawColumnTypes,
"")
                .def_property_readonly(
                    "double_datamap",
                    &tio::TrackingTxtFileContents::getDoubleDataMap,
"")
                .def_property_readonly(
                    "raw_datamap", &tio::TrackingTxtFileContents::getRawDataMap,
"")
                .def_property_readonly(
                    "num_rows", &tio::TrackingTxtFileContents::getNumRows,
"");

            py::enum_<tudat::input_output::TrackingDataType>(
                m, "TrackingDataType",
"")
                .value("year", tudat::input_output::TrackingDataType::year,
"")
                .value("month", tudat::input_output::TrackingDataType::month,
"")
                .value("day", tudat::input_output::TrackingDataType::day,
"")
                .value("hour", tudat::input_output::TrackingDataType::hour,
"")
                .value("minute", tudat::input_output::TrackingDataType::minute,
"")
                .value("second", tudat::input_output::TrackingDataType::second,
"")
                .value("time_tag_delay",
                       tudat::input_output::TrackingDataType::time_tag_delay,
"")
                .value("observation_time_scale",
                       tudat::input_output::TrackingDataType::
                           observation_time_scale,
"")
                .value("file_name",
                       tudat::input_output::TrackingDataType::file_name,
"")
                .value(
                    "n_way_light_time",
                    tudat::input_output::TrackingDataType::n_way_light_time,
"")
                .value("light_time_measurement_delay",
                       tudat::input_output::TrackingDataType::
                           light_time_measurement_delay,
"")
                .value("light_time_measurement_accuracy",
                       tudat::input_output::TrackingDataType::
                           light_time_measurement_accuracy,
"")
                .value("dsn_transmitting_station_nr",
                       tudat::input_output::TrackingDataType::
                           dsn_transmitting_station_nr,
"")
                .value(
                    "dsn_receiving_station_nr",
                    tudat::input_output::TrackingDataType::
                        dsn_receiving_station_nr,
"")
                .value(
                    "observation_body",
                    tudat::input_output::TrackingDataType::observation_body,
"")
                .value("observed_body",
                       tudat::input_output::TrackingDataType::observed_body,
"")
                .value("spacecraft_id",
                       tudat::input_output::TrackingDataType::spacecraft_id,
"")
                .value("planet_nr",
                       tudat::input_output::TrackingDataType::planet_nr,
"")
                .value("tdb_time_j2000",
                       tudat::input_output::TrackingDataType::tdb_time_j2000,
"")
                .value(
                    "tdb_spacecraft_j2000",
                    tudat::input_output::TrackingDataType::tdb_spacecraft_j2000,
"")
                .value("x_planet_frame",
                       tudat::input_output::TrackingDataType::x_planet_frame,
"")
                .value("y_planet_frame",
                       tudat::input_output::TrackingDataType::y_planet_frame,
"")
                .value("z_planet_frame",
                       tudat::input_output::TrackingDataType::z_planet_frame,
"")
                .value(
                    "vx_planet_frame",
                    tudat::input_output::TrackingDataType::vx_planet_frame,
"")
                .value(
                    "vy_planet_frame",
                    tudat::input_output::TrackingDataType::vy_planet_frame,
"")
                .value(
                    "vz_planet_frame",
                    tudat::input_output::TrackingDataType::vz_planet_frame,
"")
                .value("residual_de405",
                       tudat::input_output::TrackingDataType::residual_de405,
"")
                .value("spacecraft_transponder_delay",
                       tudat::input_output::TrackingDataType::
                           spacecraft_transponder_delay,
"")
                .value(
                    "uplink_frequency",
                    tudat::input_output::TrackingDataType::uplink_frequency,
"")
                .value(
                    "downlink_frequency",
                    tudat::input_output::TrackingDataType::downlink_frequency,
"")
                .value(
                    "signal_to_noise",
                    tudat::input_output::TrackingDataType::signal_to_noise,
"")
                .value("spectral_max",
                       tudat::input_output::TrackingDataType::spectral_max,
"")
                .value(
                    "doppler_measured_frequency",
                    tudat::input_output::TrackingDataType::
                        doppler_measured_frequency,
"")
                .value("doppler_base_frequency",
                       tudat::input_output::TrackingDataType::
                           doppler_base_frequency,
"")
                .value("doppler_noise",
                       tudat::input_output::TrackingDataType::doppler_noise,
"")
                .value(
                    "doppler_bandwidth",
                    tudat::input_output::TrackingDataType::doppler_bandwidth,
"")
                .value(
                    "vlbi_station_name",
                    tudat::input_output::TrackingDataType::vlbi_station_name,
"")
                .export_values();
        };

    }  // namespace io
}  // namespace tudatpy
