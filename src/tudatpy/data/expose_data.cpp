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
#include "tudat/io/preProcessIfmsFile.h"
#include "tudat/io/solarActivityData.h"
#include "tracking/expose_tracking.h"
#include "environment/expose_environment.h"
#include "paths/expose_paths.h"
#include "data_retrieval/expose_data_retrieval.h"
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
    auto tracking = m.def_submodule( "tracking" );
    tudatpy::data::tracking::expose_tracking( tracking );

    auto environment = m.def_submodule( "environment" );
    tudatpy::data::environment::expose_environment( environment );

    auto paths = m.def_submodule( "paths" );
    tudatpy::data::paths::expose_paths( paths );

    auto data_retrieval = m.def_submodule( "data_retrieval" );
    tudatpy::data::data_retrieval::expose_data_retrieval( data_retrieval );

    py::module_::import( "tudatpy.kernel.math.interpolators" ).attr( "InterpolatorSettings" );
    // py::module_::import( "tudatpy.math.interpolators" ).attr( "cubic_spline_interpolation" );
    // py::object cubic_spline_interpolation =
    //         (py::object)py::module_::import( "tudatpy.math.interpolators" ).attr( "cubic_spline_interpolation" );

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
};

}  // namespace data
}  // namespace tudatpy
