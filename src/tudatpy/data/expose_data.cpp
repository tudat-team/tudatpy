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
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <tudat/io/basicInputOutput.h>

#include "tudat/io/missileDatcomData.h"
#include "tudat/io/readHistoryFromFile.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/io/readVariousPdsFiles.h"
#include "tudat/io/solarActivityData.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data
{

PYBIND11_MODULE( expose_data, m )
{
    py::module_::import( "tudatpy.math.interpolators" );

    m.def( "get_resource_path",
           &tudat::paths::get_resource_path,
           R"doc(

 Get the path at which tudat resources are located.

 Returns
 -------
 str
     Local path at which tudat resources are located.






     )doc" );
    m.def( "get_ephemeris_path",
           &tudat::paths::getEphemerisDataFilesPath,
           R"doc(

 Get the path at which the ephemeris used by tudat are located.

 Returns
 -------
 str
     Local path at which the tudat ephemeris resources are located.






     )doc" );
    m.def( "get_earth_orientation_path",
           &tudat::paths::getEarthOrientationDataFilesPath,
           R"doc(

 Get the path at which the Earth orientation resources used by tudat are located.

 Returns
 -------
 str
     Local path at which tudat Earth orientation resources are located.






     )doc" );
    m.def( "get_quadrature_path",
           &tudat::paths::getQuadratureDataPath,
           R"doc(

 Get the path at which the Gaussian quadrature resources are located.

 Returns
 -------
 str
     Local path at which tudat Gaussian quadrature resources are located.






     )doc" );
    m.def( "get_spice_kernel_path",
           &tudat::paths::getSpiceKernelPath,
           R"doc(

 Get the path at which the SPICE kernel used by tudat is located.

 Returns
 -------
 str
     Local path at which the SPICE kernel is located.






     )doc" );
    m.def( "get_atmosphere_tables_path",
           &tudat::paths::getAtmosphereTablesPath,
           R"doc(

 Get the path at which tudat atmosphere tables are located.

 Returns
 -------
 str
     Local path at which tudat atmosphere tables are located.






     )doc" );
    m.def( "get_gravity_models_path",
           &tudat::paths::getGravityModelsPath,
           R"doc(

 Get the path at which tudat gravity models are located.

 Returns
 -------
 str
     Local path at which tudat gravity models are located.






     )doc" );
    m.def( "get_space_weather_path",
           &tudat::paths::getSpaceWeatherDataPath,
           R"doc(

 Get the path at which tudat space weather is located.

 Returns
 -------
 str
     Local path at which tudat space weather is located.






     )doc" );

    m.def( "read_vector_history_from_file",
           &tudat::input_output::readVectorHistoryFromFile< double, double >,
           py::arg( "vector_size" ),
           py::arg( "file_name" ),
           R"doc(

 Read a vector history from a file.


 Parameters
 ----------
 vector_size : int
     Size of the vector at each epoch.
 file_name : str
     Name of the file containing the vector history.
 Returns
 -------
 Dict[float, numpy.ndarray]
     Dictionary mapping epochs to the vector at the given epoch.






     )doc" );

    m.def( "read_matrix_history_from_file",
           &tudat::input_output::readMatrixHistoryFromFile< double, double >,
           py::arg( "matrix_rows" ),
           py::arg( "matrix_columns" ),
           py::arg( "file_name" ),
           R"doc(

 Read a matrix history from a file.


 Parameters
 ----------
 matrix_rows : int
     Number of rows in the matrix at each epoch.
 matrix_columns : int
     Number of columns in the matrix at each epoch.
 file_name : str
     Name of the file containing the matrix history.
 Returns
 -------
 Dict[float, numpy.ndarray]
     Dictionary mapping epochs to the matrix at the given epoch.






     )doc" );

    py::class_< tudat::input_output::MissileDatcomData,
                std::shared_ptr< tudat::input_output::MissileDatcomData > >(
            m, "missile_DATCOM_data", R"doc(No documentation available.)doc" )
            .def( py::init< const std::string & >( ), py::arg( "file_name_and_path" ) )
            .def( "get_static_coefficient",
                  &tudat::input_output::MissileDatcomData::getStaticCoefficient,
                  py::arg( "mach_index" ),
                  py::arg( "angle_of_attack_index" ),
                  py::arg( "coefficient_index" ),
                  R"doc(No documentation available.)doc" )
            .def( "get_dynamic_coefficient",
                  &tudat::input_output::MissileDatcomData::getDynamicCoefficient,
                  py::arg( "mach_index" ),
                  py::arg( "angle_of_attack_index" ),
                  py::arg( "coefficient_index" ),
                  R"doc(No documentation available.)doc" )
            .def( "get_angle_of_attacks",
                  &tudat::input_output::MissileDatcomData::getAngleOfAttacks,
                  R"doc(No documentation available.)doc" )
            .def( "get_mach_numbers", &tudat::input_output::MissileDatcomData::getMachNumbers )
            .def( "get_Reynolds_numbers",
                  &tudat::input_output::MissileDatcomData::getReynoldsNumbers,
                  R"doc(No documentation available.)doc" )
            .def( "write_all_coefficients_to_files",
                  &tudat::input_output::MissileDatcomData::writeAllCoefficientsToFiles,
                  py::arg( "file_name_base" ),
                  py::arg( "base_precision" ) = 15,
                  py::arg( "exponent_width" ) = 2,
                  R"doc(No documentation available.)doc" )
            .def( "write_force_and_moment_coefficients_to_files",
                  &tudat::input_output::MissileDatcomData::writeForceAndMomentCoefficientsToFiles,
                  py::arg( "file_name_base" ),
                  py::arg( "base_precision" ) = 15,
                  py::arg( "exponent_width" ) = 2,
                  R"doc(No documentation available.)doc" );

    py::enum_< tudat::input_output::MissileDatcomData::DynamicCoefficientNames >(
            m, "DynamicCoefficientNames" )
            .value( "cnq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnq )
            .value( "cmq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmq )
            .value( "caq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::caq )
            .value( "cyq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cyq )
            .value( "clnq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::clnq )
            .value( "cllq", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cllq )
            .value( "cnr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnr )
            .value( "cmr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmr )
            .value( "car", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::car )
            .value( "cyr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cyr )
            .value( "clnr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::clnr )
            .value( "cllr", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cllr )
            .value( "cnp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnp )
            .value( "cmp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmp )
            .value( "cap", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cap )
            .value( "cyp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cyp )
            .value( "clnp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::clnp )
            .value( "cllp", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cllp )
            .value( "cnad", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cnad )
            .value( "cmad", tudat::input_output::MissileDatcomData::DynamicCoefficientNames::cmad )
            .export_values( );

    py::enum_< tudat::input_output::MissileDatcomData::StaticCoefficientNames >(
            m, "StaticCoefficientNames" )
            .value( "cn", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cn )
            .value( "cm", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cm )
            .value( "ca", tudat::input_output::MissileDatcomData::StaticCoefficientNames::ca )
            .value( "cy", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cy )
            .value( "cln", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cln )
            .value( "cll", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cll )
            .value( "cna", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cna )
            .value( "cma", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cma )
            .value( "cyb", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cyb )
            .value( "cnb", tudat::input_output::MissileDatcomData::StaticCoefficientNames::cnb )
            .value( "clb", tudat::input_output::MissileDatcomData::StaticCoefficientNames::clb )
            .export_values( );

    py::class_< tio::solar_activity::SolarActivityData,
                std::shared_ptr< tio::solar_activity::SolarActivityData > >(
            m, "SolarActivityData", R"doc(No documentation available.)doc" );

    // py::class_<std::map<
    //     double,
    //     std::shared_ptr<tio::solar_activity::SolarActivityData>>>(
    //     m, "SolarActivityDataMap");

    m.def( "read_solar_activity_data",
           &tio::solar_activity::readSolarActivityData,
           py::arg( "file_path" ),
           R"doc(
 Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.

 :param file_path: Path to the space weather data file.
 )doc" );

    py::class_< tio::OdfRawFileContents, std::shared_ptr< tio::OdfRawFileContents > >(
            m, "OdfRawFileContents", R"doc(No documentation available.)doc" )
            .def( "write_to_text_file",
                  &tio::OdfRawFileContents::writeOdfToTextFile,
                  py::arg( "output_file" ),
                  R"doc(No documentation available.)doc" );

    m.def( "read_odf_file",
           &tio::readOdfFile,
           py::arg( "file_name" ),
           R"doc(No documentation available.)doc" );

    m.def( "set_dsn_weather_data_in_ground_stations",
           py::overload_cast< tudat::simulation_setup::SystemOfBodies &,
                              const std::vector< std::string > &,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::map< int, std::vector< std::string > > &,
                              const std::string & >( &tio::setDsnWeatherDataInGroundStations ),
           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg( "interpolator_settings" ) = tudat::interpolators::linearInterpolation( ),
           py::arg( "ground_stations_per_complex" ) =
                   tudat::simulation_setup::getDefaultDsnStationNamesPerComplex( ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           R"doc(No documentation available.)doc" );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ) = '#',
           py::arg( "value_separators" ) = ",:\t ",
           py::arg( "ignore_omitted_columns" ) = false );

    py::class_< tudat::input_output::TrackingTxtFileContents,
                std::shared_ptr< tudat::input_output::TrackingTxtFileContents > >(
            m, "TrackingTxtFileContents", R"doc(No documentation available.)doc" )
            .def( py::init< const std::string,
                            const std::vector< std::string >,
                            const char,
                            const std::string >( ),
                  py::arg( "file_name" ),
                  py::arg( "column_types" ),
                  py::arg( "comment_symbol" ) = '#',
                  py::arg( "value_separators" ) = ",:\t ",
                  R"doc(No documentation available.)doc" )
            .def( "add_metadata_val",
                  py::overload_cast< tio::TrackingDataType, double >(
                          &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "value" ),
                  R"doc(No documentation available.)doc" )
            .def( "get_available_datatypes",
                  &tio::TrackingTxtFileContents::getAllAvailableDataTypes,
                  R"doc(No documentation available.)doc" )
            .def( "add_metadata_str",
                  py::overload_cast< tio::TrackingDataType, const std::string & >(
                          &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "str_value" ),
                  R"doc(No documentation available.)doc" )
            .def_property_readonly( "column_field_types",
                                    &tio::TrackingTxtFileContents::getRawColumnTypes,
                                    R"doc(No documentation available.)doc" )
            .def_property_readonly( "double_datamap",
                                    &tio::TrackingTxtFileContents::getDoubleDataMap,
                                    R"doc(No documentation available.)doc" )
            .def_property_readonly( "raw_datamap",
                                    &tio::TrackingTxtFileContents::getRawDataMap,
                                    R"doc(No documentation available.)doc" )
            .def_property_readonly( "num_rows",
                                    &tio::TrackingTxtFileContents::getNumRows,
                                    R"doc(No documentation available.)doc" );

    m.def( "grail_antenna_file_reader",
           &tio::grailAntennaFileReader,
           py::arg( "file_name" ),
           R"doc(No documentation available.)doc" );
    m.def( "grail_mass_level_0_file_reader",
           &tio::grailMassLevel0FileReader,
           py::arg( "file_name" ),
           R"doc(No documentation available.)doc" );
    m.def( "grail_mass_level_1_file_reader",
           &tio::grailMassLevel1FileReader,
           py::arg( "file_name" ),
           py::arg( "data_level" ) = "1b",
           R"doc(No documentation available.)doc" );

    py::enum_< tudat::input_output::TrackingDataType >(
            m, "TrackingDataType", R"doc(No documentation available.)doc" )
            .value( "year",
                    tudat::input_output::TrackingDataType::year,
                    R"doc(No documentation available.)doc" )
            .value( "month",
                    tudat::input_output::TrackingDataType::month,
                    R"doc(No documentation available.)doc" )
            .value( "day",
                    tudat::input_output::TrackingDataType::day,
                    R"doc(No documentation available.)doc" )
            .value( "hour",
                    tudat::input_output::TrackingDataType::hour,
                    R"doc(No documentation available.)doc" )
            .value( "minute",
                    tudat::input_output::TrackingDataType::minute,
                    R"doc(No documentation available.)doc" )
            .value( "second",
                    tudat::input_output::TrackingDataType::second,
                    R"doc(No documentation available.)doc" )
            .value( "time_tag_delay",
                    tudat::input_output::TrackingDataType::time_tag_delay,
                    R"doc(No documentation available.)doc" )
            .value( "observation_time_scale",
                    tudat::input_output::TrackingDataType::observation_time_scale,
                    R"doc(No documentation available.)doc" )
            .value( "file_name",
                    tudat::input_output::TrackingDataType::file_name,
                    R"doc(No documentation available.)doc" )
            .value( "n_way_light_time",
                    tudat::input_output::TrackingDataType::n_way_light_time,
                    R"doc(No documentation available.)doc" )
            .value( "light_time_measurement_delay",
                    tudat::input_output::TrackingDataType::light_time_measurement_delay,
                    R"doc(No documentation available.)doc" )
            .value( "light_time_measurement_accuracy",
                    tudat::input_output::TrackingDataType::light_time_measurement_accuracy,
                    R"doc(No documentation available.)doc" )
            .value( "dsn_transmitting_station_nr",
                    tudat::input_output::TrackingDataType::dsn_transmitting_station_nr,
                    R"doc(No documentation available.)doc" )
            .value( "dsn_receiving_station_nr",
                    tudat::input_output::TrackingDataType::dsn_receiving_station_nr,
                    R"doc(No documentation available.)doc" )
            .value( "observation_body",
                    tudat::input_output::TrackingDataType::observation_body,
                    R"doc(No documentation available.)doc" )
            .value( "observed_body",
                    tudat::input_output::TrackingDataType::observed_body,
                    R"doc(No documentation available.)doc" )
            .value( "spacecraft_id",
                    tudat::input_output::TrackingDataType::spacecraft_id,
                    R"doc(No documentation available.)doc" )
            .value( "planet_nr",
                    tudat::input_output::TrackingDataType::planet_nr,
                    R"doc(No documentation available.)doc" )
            .value( "tdb_spacecraft_j2000",
                    tudat::input_output::TrackingDataType::tdb_spacecraft_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "x_planet_frame",
                    tudat::input_output::TrackingDataType::x_planet_frame,
                    R"doc(No documentation available.)doc" )
            .value( "y_planet_frame",
                    tudat::input_output::TrackingDataType::y_planet_frame,
                    R"doc(No documentation available.)doc" )
            .value( "z_planet_frame",
                    tudat::input_output::TrackingDataType::z_planet_frame,
                    R"doc(No documentation available.)doc" )
            .value( "vx_planet_frame",
                    tudat::input_output::TrackingDataType::vx_planet_frame,
                    R"doc(No documentation available.)doc" )
            .value( "vy_planet_frame",
                    tudat::input_output::TrackingDataType::vy_planet_frame,
                    R"doc(No documentation available.)doc" )
            .value( "vz_planet_frame",
                    tudat::input_output::TrackingDataType::vz_planet_frame,
                    R"doc(No documentation available.)doc" )
            .value( "residual_de405",
                    tudat::input_output::TrackingDataType::residual_de405,
                    R"doc(No documentation available.)doc" )
            .value( "spacecraft_transponder_delay",
                    tudat::input_output::TrackingDataType::spacecraft_transponder_delay,
                    R"doc(No documentation available.)doc" )
            .value( "uplink_frequency",
                    tudat::input_output::TrackingDataType::uplink_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "downlink_frequency",
                    tudat::input_output::TrackingDataType::downlink_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "signal_to_noise",
                    tudat::input_output::TrackingDataType::signal_to_noise,
                    R"doc(No documentation available.)doc" )
            .value( "spectral_max",
                    tudat::input_output::TrackingDataType::spectral_max,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_measured_frequency",
                    tudat::input_output::TrackingDataType::doppler_measured_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_base_frequency",
                    tudat::input_output::TrackingDataType::doppler_base_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_noise",
                    tudat::input_output::TrackingDataType::doppler_noise,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_bandwidth",
                    tudat::input_output::TrackingDataType::doppler_bandwidth,
                    R"doc(No documentation available.)doc" )
            .value( "receiving_station_name",
                    tudat::input_output::TrackingDataType::receiving_station_name,
                    R"doc(No documentation available.)doc" )
            .value( "transmitting_station_name",
                    tudat::input_output::TrackingDataType::transmitting_station_name,
                    R"doc(No documentation available.)doc" )
            .export_values( );
};

}  // namespace data
}  // namespace tudatpy
