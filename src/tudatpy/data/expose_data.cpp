/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif
#include "expose_data.h"

#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
// #include <pybind11/native_enum.h>
#include <tudat/io/basicInputOutput.h>
#include <string>
#include <vector>

#include "tudat/io/missileDatcomData.h"
#include "tudat/io/readCrdFile.h"
#include "tudat/io/readHistoryFromFile.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/io/readSinexFile.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/io/readTrackingTxtFile.h"
#include "tudat/io/readVariousPdsFiles.h"
#include "tudat/io/preProcessFdetsFile.h"
#include "tudat/io/preProcessIfmsFile.h"
#include "tudat/io/solarActivityData.h"
#include "coma_model/expose_coma_model.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/trackingSupplementaryData.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data
{

void expose_data( py::module& m )
{
    auto comaModel = m.def_submodule( "coma_model" );
    tudatpy::data::coma_model::expose_coma_model( comaModel );

    py::module_::import( "tudatpy.kernel.math.interpolators" ).attr( "InterpolatorSettings" );
    // py::module_::import( "tudatpy.math.interpolators" ).attr( "cubic_spline_interpolation" );
    // py::object cubic_spline_interpolation =
    //         (py::object)py::module_::import( "tudatpy.math.interpolators" ).attr( "cubic_spline_interpolation" );

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

    py::class_< tio::CrdPassConfigurationData >( m, "CrdPassConfigurationData", R"doc(
Container for CRD pass-level metadata.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_name", &tio::CrdPassConfigurationData::stationName_ )
            .def_readwrite( "cdp_pad_id", &tio::CrdPassConfigurationData::cdpPadId_ )
            .def_readwrite( "target_name", &tio::CrdPassConfigurationData::targetName_ )
            .def_readwrite( "start_year", &tio::CrdPassConfigurationData::startYear_ )
            .def_readwrite( "start_month", &tio::CrdPassConfigurationData::startMonth_ )
            .def_readwrite( "start_day", &tio::CrdPassConfigurationData::startDay_ )
            .def_readwrite( "end_year", &tio::CrdPassConfigurationData::endYear_ )
            .def_readwrite( "end_month", &tio::CrdPassConfigurationData::endMonth_ )
            .def_readwrite( "end_day", &tio::CrdPassConfigurationData::endDay_ )
            .def_readwrite( "transmit_wavelength_nm", &tio::CrdPassConfigurationData::transmitWavelengthNm_ );

    py::class_< tio::CrdNormalPointRecord >( m, "CrdNormalPointRecord", R"doc(
Container for a CRD normal-point (record ``11``).
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdNormalPointRecord::secondOfDay_ )
            .def_readwrite( "two_way_time_of_flight", &tio::CrdNormalPointRecord::twoWayTimeOfFlight_ )
            .def_readwrite( "one_way_range", &tio::CrdNormalPointRecord::oneWayRange_ )
            .def_readwrite( "system_configuration_id", &tio::CrdNormalPointRecord::systemConfigurationId_ )
            .def_readwrite( "epoch_event", &tio::CrdNormalPointRecord::epochEvent_ )
            .def_readwrite( "normal_point_window_length", &tio::CrdNormalPointRecord::normalPointWindowLength_ )
            .def_readwrite( "number_of_returns", &tio::CrdNormalPointRecord::numberOfReturns_ )
            .def_readwrite( "bin_rms", &tio::CrdNormalPointRecord::binRms_ );

    py::class_< tio::CrdFullRateRecord >( m, "CrdFullRateRecord", R"doc(
Container for a CRD full-rate observation (record ``10``).
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdFullRateRecord::secondOfDay_ )
            .def_readwrite( "two_way_time_of_flight", &tio::CrdFullRateRecord::twoWayTimeOfFlight_ )
            .def_readwrite( "one_way_range", &tio::CrdFullRateRecord::oneWayRange_ )
            .def_readwrite( "system_configuration_id", &tio::CrdFullRateRecord::systemConfigurationId_ )
            .def_readwrite( "epoch_event", &tio::CrdFullRateRecord::epochEvent_ )
            .def_readwrite( "filter_flag", &tio::CrdFullRateRecord::filterFlag_ )
            .def_readwrite( "detector_channel", &tio::CrdFullRateRecord::detectorChannel_ )
            .def_readwrite( "stop_number", &tio::CrdFullRateRecord::stopNumber_ );

    py::class_< tio::CrdMeteoRecord >( m, "CrdMeteoRecord", R"doc(
Container for a CRD meteorological record (record ``20``).
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdMeteoRecord::secondOfDay_ )
            .def_readwrite( "pressure", &tio::CrdMeteoRecord::pressure_ )
            .def_readwrite( "temperature", &tio::CrdMeteoRecord::temperature_ )
            .def_readwrite( "humidity", &tio::CrdMeteoRecord::humidity_ );

    py::class_< tio::CrdPassData >( m, "CrdPassData", R"doc(
Container for the measurement and meteorological data of a CRD pass.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "full_rate_data", &tio::CrdPassData::fullRateData_ )
            .def_readwrite( "normal_point_data", &tio::CrdPassData::normalPointData_ )
            .def_readwrite( "meteorological_data", &tio::CrdPassData::meteorologicalData_ );

    py::class_< tio::CrdPass >( m, "CrdPass", R"doc(
Container for a CRD pass, including configuration and data records.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "configuration", &tio::CrdPass::configuration_ )
            .def_readwrite( "data", &tio::CrdPass::data_ );

    m.def( "convert_crd_two_way_time_of_flight_to_slr_range",
           &tio::convertCrdTwoWayTimeOfFlightToSlrRange,
           py::arg( "two_way_time_of_flight" ),
           R"doc(
Convert a CRD two-way light time (seconds) into one-way SLR range (meters).
           )doc" );

    m.def( "read_crd_file",
           &tio::readCrdFile,
           py::arg( "file_name" ),
           R"doc(
Read a CRD file and return the parsed pass list.
           )doc" );

    m.def( "read_crd_files",
           &tio::readCrdFiles,
           py::arg( "file_names" ),
           R"doc(
Read multiple CRD files and concatenate all parsed passes.
           )doc" );

    m.def( "group_crd_data_per_target",
           &tio::groupCrdDataPerTarget,
           py::arg( "crd_passes" ),
           R"doc(
Group CRD passes by target name.
           )doc" );

    m.def( "group_crd_data_per_station",
           &tio::groupCrdDataPerStation,
           py::arg( "crd_passes" ),
           py::arg( "monument_id_to_ground_station_name_map" ),
           R"doc(
Group CRD passes by station name using a monument-ID to station-name map.
           )doc" );

    m.def( "extract_normal_point_measurements",
           py::overload_cast< const tio::CrdPass& >( &tio::extractNormalPointMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract normal-point measurements from one CRD pass as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "extract_normal_point_measurements_from_passes",
           py::overload_cast< const std::vector< tio::CrdPass >& >( &tio::extractNormalPointMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract normal-point measurements from multiple CRD passes as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "extract_full_rate_measurements",
           py::overload_cast< const tio::CrdPass& >( &tio::extractFullRateMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract full-rate measurements from one CRD pass as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "extract_full_rate_measurements_from_passes",
           py::overload_cast< const std::vector< tio::CrdPass >& >( &tio::extractFullRateMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
Extract full-rate measurements from multiple CRD passes as ``{epoch_utc_seconds_since_j2000: one_way_range_m}``.
           )doc" );

    m.def( "get_station_wavelengths",
           &tio::getStationWavelengths,
           py::arg( "grouped_data" ),
           R"doc(
Extract station transmit wavelengths from grouped CRD pass data.
           )doc" );

    py::class_< tio::SinexStationState >( m, "SinexStationState", R"doc(
Container for station state data parsed from a SINEX file.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "site_code", &tio::SinexStationState::siteCode_ )
            .def_readwrite( "domes_id", &tio::SinexStationState::domesId_ )
            .def_readwrite( "position", &tio::SinexStationState::position_ )
            .def_readwrite( "velocity", &tio::SinexStationState::velocity_ )
            .def_readwrite( "reference_epoch", &tio::SinexStationState::referenceEpoch_ );

    py::class_< tio::SinexStationEccentricity >( m, "SinexStationEccentricity", R"doc(
Container for SINEX station eccentricity data.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "domes_id", &tio::SinexStationEccentricity::domesId_ )
            .def_readwrite( "station_code", &tio::SinexStationEccentricity::stationCode_ )
            .def_readwrite( "eccentricity", &tio::SinexStationEccentricity::eccentricity_ )
            .def_readwrite( "start_epoch", &tio::SinexStationEccentricity::startEpoch_ )
            .def_readwrite( "end_epoch", &tio::SinexStationEccentricity::endEpoch_ )
            .def_readwrite( "has_open_end", &tio::SinexStationEccentricity::hasOpenEnd_ );

    py::class_< tio::IlrsStationRegistryEntry >( m, "IlrsStationRegistryEntry", R"doc(
Container for one ILRS station registry entry parsed from SINEX ``SITE/ID``.
    )doc" )
            .def( py::init<>( ) )
            .def_readwrite( "station_code", &tio::IlrsStationRegistryEntry::stationCode_ )
            .def_readwrite( "station_name", &tio::IlrsStationRegistryEntry::stationName_ )
            .def_readwrite( "domes_id", &tio::IlrsStationRegistryEntry::domesId_ )
            .def_readwrite( "approximate_longitude", &tio::IlrsStationRegistryEntry::approximateLongitude_ )
            .def_readwrite( "approximate_latitude", &tio::IlrsStationRegistryEntry::approximateLatitude_ )
            .def_readwrite( "approximate_height", &tio::IlrsStationRegistryEntry::approximateHeight_ );

    m.def( "convert_sinex_datetime_to_seconds_since_epoch",
           &tio::convertSinexDateTimeToSecondsSinceEpoch,
           py::arg( "date_time" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Convert a SINEX epoch string (``YY:DDD:SSSSS``) to seconds since a reference Julian day.
           )doc" );

    m.def( "read_sinex_station_data",
           &tio::readSinexStationData,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read station positions, velocities and reference epochs from a SINEX file.
           )doc" );

    m.def( "read_sinex_station_eccentricities",
           &tio::readSinexStationEccentricities,
           py::arg( "file_name" ),
           py::arg( "reference_julian_day" ) = tudat::basic_astrodynamics::JULIAN_DAY_ON_J2000,
           R"doc(
Read station eccentricity history from a SINEX file.
           )doc" );

    m.def( "read_ilrs_station_registry_from_sinex_site_id",
           &tio::readIlrsStationRegistryFromSinexSiteId,
           py::arg( "file_name" ),
           R"doc(
Read the ILRS station registry from a SINEX ``SITE/ID`` block as ``{station_code: entry}``.
           )doc" );

    m.def( "read_domes_id_numbers",
           &tio::readDomesIdNumbers,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from station name to DOMES id.
           )doc" );

    m.def( "read_monument_numbers",
           &tio::readMonumentNumbers,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from monument/station code to station name.
           )doc" );

    m.def( "read_ground_station_names",
           &tio::readGroundStationNames,
           py::arg( "file_name" ),
           R"doc(
Read a mapping from DOMES id to station name.
           )doc" );

    py::enum_< tudat::input_output::TrackingDataType >( m, "TrackingDataType", R"doc(No documentation available.)doc" )
            .value( "year", tudat::input_output::TrackingDataType::year, R"doc(No documentation available.)doc" )
            .value( "month", tudat::input_output::TrackingDataType::month, R"doc(No documentation available.)doc" )
            .value( "day", tudat::input_output::TrackingDataType::day, R"doc(No documentation available.)doc" )
            .value( "hour", tudat::input_output::TrackingDataType::hour, R"doc(No documentation available.)doc" )
            .value( "minute", tudat::input_output::TrackingDataType::minute, R"doc(No documentation available.)doc" )
            .value( "second", tudat::input_output::TrackingDataType::second, R"doc(No documentation available.)doc" )
            .value( "observation_time_scale",
                    tudat::input_output::TrackingDataType::observation_time_scale,
                    R"doc(No documentation available.)doc" )
            .value( "file_name", tudat::input_output::TrackingDataType::file_name, R"doc(No documentation available.)doc" )
            .value( "n_way_light_time", tudat::input_output::TrackingDataType::n_way_light_time, R"doc(No documentation available.)doc" )
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
            .value( "observation_body", tudat::input_output::TrackingDataType::observation_body, R"doc(No documentation available.)doc" )
            .value( "observed_body", tudat::input_output::TrackingDataType::observed_body, R"doc(No documentation available.)doc" )
            .value( "spacecraft_id", tudat::input_output::TrackingDataType::spacecraft_id, R"doc(No documentation available.)doc" )
            .value( "spacecraft_name", tudat::input_output::TrackingDataType::spacecraft_name, R"doc(No documentation available.)doc" )
            .value( "planet_nr", tudat::input_output::TrackingDataType::planet_nr, R"doc(No documentation available.)doc" )
            .value( "tdb_reception_time_j2000",
                    tudat::input_output::TrackingDataType::tdb_reception_time_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "utc_reception_time_j2000",
                    tudat::input_output::TrackingDataType::utc_reception_time_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "utc_ramp_referencee_j2000",
                    tudat::input_output::TrackingDataType::utc_ramp_referencee_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "tdb_spacecraft_j2000",
                    tudat::input_output::TrackingDataType::tdb_spacecraft_j2000,
                    R"doc(No documentation available.)doc" )
            .value( "x_planet_frame", tudat::input_output::TrackingDataType::x_planet_frame, R"doc(No documentation available.)doc" )
            .value( "y_planet_frame", tudat::input_output::TrackingDataType::y_planet_frame, R"doc(No documentation available.)doc" )
            .value( "z_planet_frame", tudat::input_output::TrackingDataType::z_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vx_planet_frame", tudat::input_output::TrackingDataType::vx_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vy_planet_frame", tudat::input_output::TrackingDataType::vy_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vz_planet_frame", tudat::input_output::TrackingDataType::vz_planet_frame, R"doc(No documentation available.)doc" )
            .value( "residual_de405", tudat::input_output::TrackingDataType::residual_de405, R"doc(No documentation available.)doc" )
            .value( "spacecraft_transponder_delay",
                    tudat::input_output::TrackingDataType::spacecraft_transponder_delay,
                    R"doc(No documentation available.)doc" )
            .value( "uplink_frequency", tudat::input_output::TrackingDataType::uplink_frequency, R"doc(No documentation available.)doc" )
            .value( "downlink_frequency",
                    tudat::input_output::TrackingDataType::downlink_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "signal_to_noise", tudat::input_output::TrackingDataType::signal_to_noise, R"doc(No documentation available.)doc" )
            .value( "spectral_max", tudat::input_output::TrackingDataType::spectral_max, R"doc(No documentation available.)doc" )
            .value( "doppler_measured_frequency",
                    tudat::input_output::TrackingDataType::doppler_measured_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_averaged_frequency",
                    tudat::input_output::TrackingDataType::doppler_averaged_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_integration_time",
                    tudat::input_output::TrackingDataType::doppler_integration_time,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_base_frequency",
                    tudat::input_output::TrackingDataType::doppler_base_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_noise", tudat::input_output::TrackingDataType::doppler_noise, R"doc(No documentation available.)doc" )
            .value( "doppler_bandwidth", tudat::input_output::TrackingDataType::doppler_bandwidth, R"doc(No documentation available.)doc" )
            .value( "receiving_station_name",
                    tudat::input_output::TrackingDataType::receiving_station_name,
                    R"doc(No documentation available.)doc" )
            .value( "transmitting_station_name",
                    tudat::input_output::TrackingDataType::transmitting_station_name,
                    R"doc(No documentation available.)doc" )
            .value( "time_tag_delay", tudat::input_output::TrackingDataType::time_tag_delay )
            .value( "sample_number", tudat::input_output::TrackingDataType::sample_number )
            .value( "utc_day_of_year", tudat::input_output::TrackingDataType::utc_day_of_year )
            .value( "reference_body_distance", tudat::input_output::TrackingDataType::reference_body_distance )
            .value( "transmission_frequency_constant_term", tudat::input_output::TrackingDataType::transmission_frequency_constant_term )
            .value( "transmission_frequency_linear_term", tudat::input_output::TrackingDataType::transmission_frequency_linear_term )
            .value( "doppler_predicted_frequency_hz", tudat::input_output::TrackingDataType::doppler_predicted_frequency_hz )
            .value( "doppler_troposphere_correction", tudat::input_output::TrackingDataType::doppler_troposphere_correction )
            .value( "scan_nr", tudat::input_output::TrackingDataType::scan_nr )
            .export_values( );

    py::enum_< tudat::input_output::TrackingTxtFileReadFilterType >(
            m, "TrackingTxtFileReadFilterType", R"doc(No documentation available.)doc" )
            .value( "no_tracking_txt_file_filter",
                    tudat::input_output::TrackingTxtFileReadFilterType::no_tracking_txt_file_filter,
                    R"doc(No documentation available.)doc" )
            .value( "ifms_tracking_txt_file_filter",
                    tudat::input_output::TrackingTxtFileReadFilterType::ifms_tracking_txt_file_filter,
                    R"doc(No documentation available.)doc" )
            .export_values( );

    py::enum_< tudat::input_output::FdetDateFormat >( m, "FdetDateFormat", R"doc(Date format used in an Fdets file.)doc" )
            .value( "datetime_string",
                    tudat::input_output::FdetDateFormat::datetime_string,
                    R"doc(Date is provided in a single UTC datetime string column.)doc" )
            .value( "pair_of_numbers",
                    tudat::input_output::FdetDateFormat::pair_of_numbers,
                    R"doc(Date is provided as a pair of numeric columns. This format is not currently supported.)doc" )
            .export_values( );

    py::class_< tio::solar_activity::SolarActivityData, std::shared_ptr< tio::solar_activity::SolarActivityData > >(
            m, "SolarActivityData", R"doc(No documentation available.)doc" )
            .def_readonly( "solar_radio_flux_107_observed", &tio::solar_activity::SolarActivityData::solarRadioFlux107Observed );

    // py::class_<std::map<
    //     double,
    //     std::shared_ptr<tio::solar_activity::SolarActivityData>>>(
    //     m, "SolarActivityDataMap");

    m.def( "read_solar_activity_data",
           &tio::solar_activity::readSolarActivityData,
           py::arg( "file_path" ),
           R"doc(
 Reads a space weather data file and produces a dictionary with solar activity data for a range of epochs. Data files can be obtained from http://celestrak.com/SpaceData and should follow the legacy format.

 Parameters
 ----------
 file_path : str
     Path to the space weather data file.
 )doc" );

    py::class_< tio::solar_activity::SolarActivityContainer, std::shared_ptr< tio::solar_activity::SolarActivityContainer > >(
            m, "SolarActivityContainer" )

            .def( py::init< const std::map< double, std::shared_ptr< tio::solar_activity::SolarActivityData > >& >( ),
                  py::arg( "solar_activity_data_map" ) )

            .def( "get_solar_activity_data",
                  &tio::solar_activity::SolarActivityContainer::getSolarActivityData,
                  py::arg( "time" ),
                  R"doc(
        Returns the nearest SolarActivityData (in UTC Julian days) for the given time in seconds since J2000.
        )doc" )

            .def( "get_solar_activity_data_map",
                  &tio::solar_activity::SolarActivityContainer::getSolarActivityDataMap,
                  R"doc(Returns the full map of SolarActivityData.)doc" );

    py::enum_< tio::OdfDataType >( m, "OdfDataType", R"doc(Possible data types in orbit section of ODF file)doc" )
            .value( "narrowband_spacecraft_vlbi_doppler_mode", tio::OdfDataType::narrowband_spacecraft_vlbi_doppler_mode )
            .value( "narrowband_spacecraft_vlbi_phase_mode", tio::OdfDataType::narrowband_spacecraft_vlbi_phase_mode )
            .value( "narrowband_quasar_vlbi_doppler_mode", tio::OdfDataType::narrowband_quasar_vlbi_doppler_mode )
            .value( "narrowband_quasar_vlbi_phase_mode", tio::OdfDataType::narrowband_quasar_vlbi_phase_mode )
            .value( "wideband_spacecraft_vlbi", tio::OdfDataType::wideband_spacecraft_vlbi )
            .value( "wideband_quasar_vlbi", tio::OdfDataType::wideband_quasar_vlbi )
            .value( "one_way_doppler", tio::OdfDataType::one_way_doppler )
            .value( "two_way_doppler", tio::OdfDataType::two_way_doppler )
            .value( "three_way_doppler", tio::OdfDataType::three_way_doppler )
            .value( "one_way_total_count_phase", tio::OdfDataType::one_way_total_count_phase )
            .value( "two_way_total_count_phase", tio::OdfDataType::two_way_total_count_phase )
            .value( "three_way_total_count_phase", tio::OdfDataType::three_way_total_count_phase )
            .value( "pra_planetary_operational_discrete_spectrum_range",
                    tio::OdfDataType::pra_planetary_operational_discrete_spectrum_range )
            .value( "sra_planetary_operational_discrete_spectrum_range",
                    tio::OdfDataType::sra_planetary_operational_discrete_spectrum_range )
            .value( "re_range", tio::OdfDataType::re_range )
            .value( "azimuth_angle", tio::OdfDataType::azimuth_angle )
            .value( "elevation_angle", tio::OdfDataType::elevation_angle )
            .value( "hour_angle", tio::OdfDataType::hour_angle )
            .value( "declination_angle", tio::OdfDataType::declination_angle )
            .value( "x_angle_east", tio::OdfDataType::x_angle_east )
            .value( "y_angle_east", tio::OdfDataType::y_angle_east )
            .value( "x_angle_south", tio::OdfDataType::x_angle_south )
            .value( "y_angle_south", tio::OdfDataType::y_angle_south );

    py::class_< tio::OdfCommonDataBlock, std::shared_ptr< tio::OdfCommonDataBlock > >(
            m, "OdfCommonDataBlock", R"doc(Base class observable-independent ODF data containers

        The data section of an ODF is split into blocks (lines), each associated with
        an observation epoch. The first elements of each block (e.g. observation epoch,
        value of the observable) are common for all the observable types, but the
        values in the remaining columns will have different meanings for different
        types of observations. This class serves as interface to the observable-independent
        part of an ODF data block. The different classes inheriting from OdfDataSpecificBlock
        provide interfaces to the observable-specific part of the blocks.
        )doc" )
            .def_property_readonly( "observable_time", &tio::OdfCommonDataBlock::getObservableTime )
            .def_property_readonly( "observable_value", &tio::OdfCommonDataBlock::getObservableValue )
            .def_property_readonly( "receiving_station_downlink_delay", &tio::OdfCommonDataBlock::getReceivingStationDownlinkDelay )
            .def_readonly( "format_id", &tio::OdfCommonDataBlock::formatId_ )
            .def_readonly( "receiving_station_id", &tio::OdfCommonDataBlock::receivingStationId_ )
            .def_readonly( "transmitting_station_id", &tio::OdfCommonDataBlock::transmittingStationId_ )
            .def_readonly( "transmitting_station_network_id", &tio::OdfCommonDataBlock::transmittingStationNetworkId_ )
            .def_readonly( "data_type", &tio::OdfCommonDataBlock::dataType_ )
            .def_readonly( "downlink_band_id", &tio::OdfCommonDataBlock::downlinkBandId_ )
            .def_readonly( "uplink_band_id", &tio::OdfCommonDataBlock::uplinkBandId_ )
            .def_readonly( "reference_band_id", &tio::OdfCommonDataBlock::referenceBandId_ )
            .def_readonly( "is_invalid", &tio::OdfCommonDataBlock::validity_ )
            .def( "print_data_block",
                  &tio::OdfCommonDataBlock::printDataBlock,
                  py::arg( "output_file" ),
                  R"doc(Write the contents of the data block to a text file

                  The file is created if it does not exist, and it can have, for example, txt extension

                  Parameters
                  ----------
                  output_file : str
                      Contents will be written to the file defined by this path
                  )doc" );

    py::class_< tio::OdfDataSpecificBlock, std::shared_ptr< tio::OdfDataSpecificBlock > >(
            m, "OdfDataSpecificBlock", R"doc(Base class observable-dependent ODF data containers

        The data section of an ODF is split into blocks (lines), each associated with
        an observation epoch. The first elements of each block (e.g. observation epoch,
        value of the observable) are common for all the observable types, but the
        values in the remaining columns will have different meanings for different
        types of observations. This base class serves as parent for interfaces to the
        observable-specific part of an ODF data block. The interface to the common
        part of the blocks is provided by the OdfCommonDataBlock class.
        )doc" );

    py::class_< tio::OdfDopplerDataBlock, std::shared_ptr< tio::OdfDopplerDataBlock >, tio::OdfDataSpecificBlock >(
            m, "OdfDopplerDataBlock", R"doc(Container for ODF Doppler-specific data)doc" )
            .def_property_readonly( "receiver_channel", &tio::OdfDopplerDataBlock::getReceiverChannel )
            .def_property_readonly( "spacecraft_id", &tio::OdfDopplerDataBlock::getSpacecraftId )
            .def_property_readonly( "receiver_exciter_flag", &tio::OdfDopplerDataBlock::getReceiverExciterFlag )
            .def_property_readonly( "reference_frequency", &tio::OdfDopplerDataBlock::getReferenceFrequency )
            .def_property_readonly( "compression_time", &tio::OdfDopplerDataBlock::getCompressionTime )
            .def_property_readonly( "transmitting_station_uplink_delay", &tio::OdfDopplerDataBlock::getTransmittingStationUplinkDelay );

    py::class_< tio::OdfDataBlock, std::shared_ptr< tio::OdfDataBlock > >(
            m, "OdfDataBlock", R"doc(Contents of a line of the data section of an ODF)doc" )
            .def_property_readonly( "observable_specific_data_block", &tio::OdfDataBlock::getObservableSpecificDataBlock )
            .def_property_readonly( "common_data_block", &tio::OdfDataBlock::getCommonDataBlock )
            .def( "print_data_block",
                  &tio::OdfDataBlock::printDataBlock,
                  py::arg( "output_file" ),
                  R"doc(Write the contents of the data block to a text file

                  The file is created if it does not exist, and it can have, for example, txt extension

                  Parameters
                  ----------
                  output_file : str
                      Contents will be written to the file defined by this path
                  )doc" );

    py::class_< tio::OdfRampBlock, std::shared_ptr< tio::OdfRampBlock > >(
            m, "OdfRampBlock", R"doc(Contents of a line of the ramp section of an ODF)doc" )
            .def_property_readonly( "ramp_start_frequency", &tio::OdfRampBlock::getRampStartFrequency )
            .def_property_readonly( "ramp_rate", &tio::OdfRampBlock::getRampRate )
            .def_property_readonly( "ramp_start_epoch", &tio::OdfRampBlock::getRampStartTime )
            .def_property_readonly( "ramp_end_epoch", &tio::OdfRampBlock::getRampEndTime )
            .def_property_readonly( "transmitting_station_id", &tio::OdfRampBlock::getTransmittingStationId );

    py::class_< tio::OdfRawFileContents, std::shared_ptr< tio::OdfRawFileContents > >(
            m, "OdfRawFileContents", R"doc(No documentation available.)doc" )
            .def_property_readonly( "data_blocks", &tio::OdfRawFileContents::getDataBlocks )
            .def_property_readonly( "ramp_blocks", &tio::OdfRawFileContents::getRampBlocks )
            .def_property_readonly( "clock_offset_blocks", &tio::OdfRawFileContents::getClockOffsetBlocks )
            .def_readonly( "file_reference_date", &tio::OdfRawFileContents::fileReferenceDate_ )
            .def_readonly( "file_reference_time", &tio::OdfRawFileContents::fileReferenceTime_ )
            .def( "write_to_text_file",
                  &tio::OdfRawFileContents::writeOdfToTextFile,
                  py::arg( "output_file" ),
                  R"doc(No documentation available.)doc" );

    m.def( "read_odf_file", &tio::readOdfFile, py::arg( "file_name" ), R"doc(No documentation available.)doc" );

    m.def( "set_dsn_weather_data_in_ground_stations",
           py::overload_cast< tudat::simulation_setup::SystemOfBodies&,
                              const std::vector< std::string >&,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::map< int, std::vector< std::string > >&,
                              const std::string& >( &tio::setDsnWeatherDataInGroundStations ),
           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg( "interpolator_settings" ) = tudat::interpolators::linearInterpolation( ),
           py::arg( "ground_stations_per_complex" ) = tudat::simulation_setup::getDefaultDsnStationNamesPerComplex( ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           R"doc(No documentation available.)doc" );

    py::class_< tudat::input_output::TrackingTxtFileContents, std::shared_ptr< tudat::input_output::TrackingTxtFileContents > >(
            m, "TrackingTxtFileContents", R"doc(No documentation available.)doc" )
            .def( py::init< const std::string, const std::vector< std::string >, const char, const std::string >( ),
                  py::arg( "file_name" ),
                  py::arg( "column_types" ),
                  py::arg( "comment_symbol" ) = '#',
                  py::arg( "value_separators" ) = ",:\t ",
                  R"doc(No documentation available.)doc" )
            .def( "add_metadata_val",
                  py::overload_cast< tio::TrackingDataType, double >( &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "value" ),
                  R"doc(No documentation available.)doc" )
            .def( "get_available_datatypes",
                  &tio::TrackingTxtFileContents::getAllAvailableDataTypes,
                  R"doc(No documentation available.)doc" )
            .def( "add_metadata_str",
                  py::overload_cast< tio::TrackingDataType, const std::string& >( &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "str_value" ),
                  R"doc(No documentation available.)doc" )
            .def_property_readonly(
                    "column_field_types", &tio::TrackingTxtFileContents::getRawColumnTypes, R"doc(No documentation available.)doc" )
            .def_property_readonly(
                    "double_datamap", &tio::TrackingTxtFileContents::getDoubleDataMap, R"doc(No documentation available.)doc" )
            .def_property_readonly( "raw_datamap", &tio::TrackingTxtFileContents::getRawDataMap, R"doc(No documentation available.)doc" )
            .def_property_readonly( "num_rows", &tio::TrackingTxtFileContents::getNumRows, R"doc(No documentation available.)doc" );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ) = '#',
           py::arg( "value_separators" ) = ",:\t ",
           py::arg( "ignore_omitted_columns" ) = false,
           py::arg( "data_filter_method" ) = tio::no_tracking_txt_file_filter );

    m.def( "grail_antenna_file_reader", &tio::grailAntennaFileReader, py::arg( "file_name" ), R"doc(No documentation available.)doc" );
    m.def( "grail_mass_level_0_file_reader",
           &tio::grailMassLevel0FileReader,
           py::arg( "file_name" ),
           R"doc(No documentation available.)doc" );
    m.def( "grail_mass_level_1_file_reader",
           &tio::grailMassLevel1FileReader,
           py::arg( "file_name" ),
           py::arg( "data_level" ) = "1b",
           R"doc(No documentation available.)doc" );

    m.def( "read_ifms_file",
           &tio::readIfmsFile,
           py::arg( "file_name" ),
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           R"doc(Load contents of IFMS file into object

           The keys of the dictionary represent the different columns of the IFMS file, and their values are lists with all the values in the associated column as strings.

           Two of the columns of an IFMS file contain, respectively, the Doppler averaged frequency and a tropospheric correction for the station. When the `apply_tropospheric_correction` option is set to true, the content of the first column is modified by subtracting the values in the second.

           Parameters
           ----------
           file_name : str
               String representing the path to the file to be loaded
           apply_tropospheric_correction : bool
               Whether to modify the averaged Doppler frequency as described above (Default: True)
           remove_invalid_lines : bool
               Boolean defining whether a line is skipped if the transmit frequency, observed frequency, or troposphere correction is undefined (Default: True)

           Returns
           -------
           ifms_contents : TrackingTxtFileContents
               Dictionary with contents of the IFMS file as lists of strings
           )doc" );

    m.def( "read_ifms_files",
           static_cast< std::pair< std::vector< std::shared_ptr< tudat::data::TrackingData< double, double > > >,
                                   std::vector< std::shared_ptr< tudat::data::TrackingSupplementaryData > > > ( * )(
                   const std::vector< std::string >&,
                   const std::string&,
                   const std::vector< std::string >&,
                   const std::string&,
                   const bool,
                   const bool ) >( &tio::readIfmsFiles< double, double > ),
           py::arg( "ifms_file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           py::arg( "apply_tropospheric_correction" ) = true,
           py::arg( "remove_invalid_lines" ) = true,
           R"doc(Load IFMS files into tracking data and supplementary data objects.)doc" );

    m.def( "set_estrack_weather_data_in_ground_stations",
           py::overload_cast< tudat::simulation_setup::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::string,
                              std::shared_ptr< tudat::interpolators::InterpolatorSettings >,
                              const std::string& >( &tio::setEstrackWeatherDataInGroundStation ),

           py::arg( "bodies" ),
           py::arg( "weather_file_names" ),
           py::arg( "ground_station_name" ),
           py::arg( "interpolator_settings" ) = tudat::interpolators::cubicSplineInterpolation( ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           R"doc(No documentation available.)doc" );

    m.def( "read_fdets_file",
           py::overload_cast< const std::string&, tio::FdetDateFormat >( &tio::readFdetsFile ),
           py::arg( "file_name" ),
           py::arg( "date_format" ) = tio::FdetDateFormat::datetime_string,
           R"doc(Load contents of Fdets file into object

           This function loads the contents of an Fdets file. For the currently supported `datetime_string` date format, the following structure is assumed:

           - UTC datetime string
           - Signal to noise ratio [-]
           - Normalized spectral max (Not sure what this is)
           - Doppler measured frequency [Hz]
           - Doppler noise [Hz]

           If the file contains an additional leading scan-number column, this column is detected automatically.

           Parameters
           ----------
           file_name : str
               String representing the path to the file to be loaded
           date_format : FdetDateFormat
               Date format used in the Fdets file

           Returns
           -------
           fdets_contents : TrackingTxtFileContents
               Dictionary with contents of the Fdets file as lists of strings
           )doc" );

    m.def( "read_fdets_file",
           py::overload_cast< const std::string&, const std::vector< std::string >& >( &tio::readFdetsFile ),
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           R"doc(Load contents of Fdets file into object using explicit column identifiers.

           .. deprecated::
              Passing explicit column identifiers is deprecated. Use the `date_format` argument instead.

           Parameters
           ----------
           file_name : str
               String representing the path to the file to be loaded
           column_types : List[str]
               Identifiers of the columns present in the Fdets file

           Returns
           -------
           fdets_contents : TrackingTxtFileContents
               Dictionary with contents of the Fdets file as lists of strings
           )doc" );

    m.def( "read_fdets_files",
           static_cast< std::pair< std::vector< std::shared_ptr< tudat::data::TrackingData< double, double > > >,
                                   std::vector< std::shared_ptr< tudat::data::TrackingSupplementaryData > > > ( * )(
                   const std::vector< std::string >&,
                   const std::vector< double >&,
                   const tio::FdetDateFormat,
                   const std::string&,
                   const std::vector< std::string >&,
                   const std::vector< std::string >&,
                   const std::string& ) >( &tio::readFdetsFiles< double, double > ),
           py::arg( "fdets_file_names" ),
           py::arg( "base_frequencies" ),
           py::arg( "date_format" ),
           py::arg( "spacecraft_name" ),
           py::arg( "transmitting_station_names" ),
           py::arg( "receiving_station_names" ),
           py::arg( "earth_name" ) = "Earth",
           R"doc(Load FDETS files into tracking data and supplementary data objects.)doc" );

    py::class_< tudat::data::TrackingData<>, std::shared_ptr< tudat::data::TrackingData<> > >( m, "TrackingData", R"doc(
 TrackingData Class container.
    )doc" )
            .def( py::init< const std::string,
                            const tudat::data::PlainLinkDefinition&,
                            const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >&,
                            const std::vector< double >,
                            const std::string,
                            const std::string >( ),
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "observations" ),
                  py::arg( "epochs" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "time_scale" ) = "TDB",
                  R"doc(Creates a TrackingData object.)doc" )
            .def_property_readonly( "observable_type", &tudat::data::TrackingData<>::getObservableType )
            .def_property_readonly( "link_ends", &tudat::data::TrackingData<>::getLinkEnds )
            .def_property_readonly( "observations", &tudat::data::TrackingData<>::getObservations )
            .def_property_readonly( "epochs", &tudat::data::TrackingData<>::getObservationEpochs )
            .def_property_readonly( "reference_link_end", &tudat::data::TrackingData<>::getReferenceLinkEnd )
            .def_property_readonly( "time_scale", &tudat::data::TrackingData<>::getTimeScale )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::vector< std::string > >(
                          &tudat::data::TrackingData<>::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (optional).)doc" )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const double >( &tudat::data::TrackingData<>::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (double type), (optional).)doc" )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::vector< double > >( &tudat::data::TrackingData<>::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (double vector type), (optional).)doc" )
            .def( "get_ancillary_settings_string_vector",
                  &tudat::data::TrackingData<>::getAncillarySettingsStringVector,
                  R"doc(Returns the ancillary settings in the TrackingData object as a vector of strings.)doc" )
            .def( "get_ancillary_settings_double",
                  &tudat::data::TrackingData<>::getAncillarySettingsDouble,
                  R"doc(Returns the ancillary settings in the TrackingData object as a double.)doc" )
            .def( "get_ancillary_settings_double_vector",
                  &tudat::data::TrackingData<>::getAncillarySettingsDoubleVector,
                  R"doc(Returns the ancillary settings in the TrackingData object as a double vector.)doc" )
            .def( "get_ancillary_settings_double_vector",
                  &tudat::data::TrackingData<>::getAncillarySettingsDoubleVector,
                  R"doc(Returns the ancillary settings in the TrackingData object as a double vector.)doc" )
            .def( "set_observation_weights",
                  ( &tudat::data::TrackingData<>::setObservationWeights ),
                  py::arg( "observationWeights" ),
                  R"doc(Adds (sets) observation weights to the TrackingData object (optional).
						It also perfroms some necessary consistency checks  (size, pre-existence) on the weights.)doc" )
            .def( "reset_single_observation_weight",
                  ( &tudat::data::TrackingData<>::setSingleObservationWeight ),
                  py::arg( "index" ),
                  py::arg( "observationWeight" ),
                  R"doc(Allows resetting of a single observation weight (specified by index i)  into the TrackingData object.
						It also perfroms some necessary consistency checks.)doc" )
            .def( "get_observation_weights",
                  &tudat::data::TrackingData<>::getObservationWeights,
                  R"doc(Returns the list (of arrays) of observation weights stored in the TrackingData object.)doc" )
            .def( "get_concatenated_observation_weights",
                  &tudat::data::TrackingData<>::getObservationWeightsVector,
                  R"doc(Returns the concatenated list of observation weights stored in the TrackingData object.)doc" )
            .def( "set_observation_corrections",
                  &tudat::data::TrackingData<>::setObservationCorrections,
                  py::arg( "observationCorrections" ),
                  R"doc(Adds corrections to the TrackingData object (optional).
					It also perfroms some necessary consistency checks.)doc" )
            .def( "reset_single_observation_correction",
                  ( &tudat::data::TrackingData<>::setSingleObservationCorrection ),
                  py::arg( "index" ),
                  py::arg( "observationCorrection" ),
                  R"doc(Allows resetting of a single observation correction (specified by index i)  into the TrackingData object.
						It also perfroms some necessary consistency checks.)doc" )
            .def( "remove_single_observation_entry",
                  ( &tudat::data::TrackingData<>::removeSingleObservationEntry ),
                  py::arg( "index" ),
                  R"doc(Removes a single observation entry from the TrackingData object.)doc" );

    py::class_< tudat::data::TrackingSupplementaryData, std::shared_ptr< tudat::data::TrackingSupplementaryData > >(
            m, "TrackingSupplementaryData", R"doc(Tracking supplementary data container.)doc" )
            .def_property_readonly( "body_name", &tudat::data::TrackingSupplementaryData::getBodyName )
            .def_property_readonly( "reference_point_name", &tudat::data::TrackingSupplementaryData::getReferencePointName );
};

}  // namespace data
}  // namespace tudatpy
