#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_slr.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/readCrdFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace slr
{

void expose_slr( py::module& m )
{
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
}

}  // namespace slr

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
