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

namespace data_input
{

namespace tracking_data
{

namespace slr
{

//! Register raw CRD records and ordered measurement extraction without applying physical corrections.
void expose_slr( py::module& m )
{
    py::class_< tio::CrdSystemConfigurationRecord >(
            m,
            "CrdSystemConfigurationRecord",
            "C0 system identifier, transmit wavelength in nanometres and component configuration identifiers." )
            .def( py::init<>( ) )
            .def_readwrite( "system_configuration_id", &tio::CrdSystemConfigurationRecord::systemConfigurationId_ )
            .def_readwrite( "transmit_wavelength_nm", &tio::CrdSystemConfigurationRecord::transmitWavelengthNm_ )
            .def_readwrite( "component_configuration_ids", &tio::CrdSystemConfigurationRecord::componentConfigurationIds_ );

    py::class_< tio::CrdPassConfigurationData >( m,
                                                 "CrdPassConfigurationData",
                                                 "Station, target and session flags. Applied flags use 0/1, with -1 for unavailable "
                                                 "fields; station time codes 3/4/7 denote UTC references." )
            .def( py::init<>( ) )
            .def_readwrite( "source", &tio::CrdPassConfigurationData::source_ )
            .def_readwrite( "format_version", &tio::CrdPassConfigurationData::formatVersion_ )
            .def_readwrite( "station_name", &tio::CrdPassConfigurationData::stationName_ )
            .def_readwrite( "cdp_pad_id", &tio::CrdPassConfigurationData::cdpPadId_ )
            .def_readwrite( "system_number", &tio::CrdPassConfigurationData::systemNumber_ )
            .def_readwrite( "occupancy", &tio::CrdPassConfigurationData::occupancy_ )
            .def_readwrite( "station_epoch_time_scale", &tio::CrdPassConfigurationData::stationEpochTimeScale_ )
            .def_readwrite( "station_network", &tio::CrdPassConfigurationData::stationNetwork_ )
            .def_readwrite( "target_name", &tio::CrdPassConfigurationData::targetName_ )
            .def_readwrite( "ilrs_satellite_id", &tio::CrdPassConfigurationData::ilrsSatelliteId_ )
            .def_readwrite( "data_type", &tio::CrdPassConfigurationData::dataType_ )
            .def_readwrite( "start_year", &tio::CrdPassConfigurationData::startYear_ )
            .def_readwrite( "start_month", &tio::CrdPassConfigurationData::startMonth_ )
            .def_readwrite( "start_day", &tio::CrdPassConfigurationData::startDay_ )
            .def_readwrite( "start_hour", &tio::CrdPassConfigurationData::startHour_ )
            .def_readwrite( "start_minute", &tio::CrdPassConfigurationData::startMinute_ )
            .def_readwrite( "start_second", &tio::CrdPassConfigurationData::startSecond_ )
            .def_readwrite( "end_year", &tio::CrdPassConfigurationData::endYear_ )
            .def_readwrite( "end_month", &tio::CrdPassConfigurationData::endMonth_ )
            .def_readwrite( "end_day", &tio::CrdPassConfigurationData::endDay_ )
            .def_readwrite( "end_hour", &tio::CrdPassConfigurationData::endHour_ )
            .def_readwrite( "end_minute", &tio::CrdPassConfigurationData::endMinute_ )
            .def_readwrite( "end_second", &tio::CrdPassConfigurationData::endSecond_ )
            .def_readwrite( "data_release", &tio::CrdPassConfigurationData::dataRelease_ )
            .def_readwrite( "tropospheric_correction_applied", &tio::CrdPassConfigurationData::troposphericCorrectionApplied_ )
            .def_readwrite( "center_of_mass_correction_applied", &tio::CrdPassConfigurationData::centerOfMassCorrectionApplied_ )
            .def_readwrite( "receive_amplitude_correction_applied", &tio::CrdPassConfigurationData::receiveAmplitudeCorrectionApplied_ )
            .def_readwrite( "station_system_delay_applied", &tio::CrdPassConfigurationData::stationSystemDelayApplied_ )
            .def_readwrite( "spacecraft_system_delay_applied", &tio::CrdPassConfigurationData::spacecraftSystemDelayApplied_ )
            .def_readwrite( "range_type", &tio::CrdPassConfigurationData::rangeType_ )
            .def_readwrite( "data_quality_alert", &tio::CrdPassConfigurationData::dataQualityAlert_ )
            .def_readwrite( "system_configurations", &tio::CrdPassConfigurationData::systemConfigurations_ );

    py::class_< tio::CrdNormalPointRecord >(
            m,
            "CrdNormalPointRecord",
            "A normal point retaining original event epoch, time of flight in seconds and distribution statistics in picoseconds." )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdNormalPointRecord::secondOfDay_ )
            .def_readwrite( "two_way_time_of_flight", &tio::CrdNormalPointRecord::twoWayTimeOfFlight_ )
            .def_readwrite( "one_way_range", &tio::CrdNormalPointRecord::oneWayRange_ )
            .def_readwrite( "system_configuration_id", &tio::CrdNormalPointRecord::systemConfigurationId_ )
            .def_readwrite( "epoch_event", &tio::CrdNormalPointRecord::epochEvent_ )
            .def_readwrite( "normal_point_window_length", &tio::CrdNormalPointRecord::normalPointWindowLength_ )
            .def_readwrite( "number_of_returns", &tio::CrdNormalPointRecord::numberOfReturns_ )
            .def_readwrite( "bin_rms", &tio::CrdNormalPointRecord::binRms_ )
            .def_readwrite( "bin_skew", &tio::CrdNormalPointRecord::binSkew_ )
            .def_readwrite( "bin_kurtosis", &tio::CrdNormalPointRecord::binKurtosis_ )
            .def_readwrite( "bin_peak_minus_mean", &tio::CrdNormalPointRecord::binPeakMinusMean_ )
            .def_readwrite( "return_rate", &tio::CrdNormalPointRecord::returnRate_ )
            .def_readwrite( "detector_channel", &tio::CrdNormalPointRecord::detectorChannel_ )
            .def_readwrite( "signal_to_noise_ratio", &tio::CrdNormalPointRecord::signalToNoiseRatio_ );

    py::class_< tio::CrdFullRateRecord >(
            m, "CrdFullRateRecord", "A full-rate measurement with its original event epoch and detector identifiers." )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdFullRateRecord::secondOfDay_ )
            .def_readwrite( "two_way_time_of_flight", &tio::CrdFullRateRecord::twoWayTimeOfFlight_ )
            .def_readwrite( "one_way_range", &tio::CrdFullRateRecord::oneWayRange_ )
            .def_readwrite( "system_configuration_id", &tio::CrdFullRateRecord::systemConfigurationId_ )
            .def_readwrite( "epoch_event", &tio::CrdFullRateRecord::epochEvent_ )
            .def_readwrite( "filter_flag", &tio::CrdFullRateRecord::filterFlag_ )
            .def_readwrite( "detector_channel", &tio::CrdFullRateRecord::detectorChannel_ )
            .def_readwrite( "stop_number", &tio::CrdFullRateRecord::stopNumber_ );

    py::class_< tio::CrdMeteoRecord >(
            m, "CrdMeteoRecord", "Weather: hPa pressure, kelvin temperature, percent relative humidity and measured/interpolated origin." )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdMeteoRecord::secondOfDay_ )
            .def_readwrite( "pressure", &tio::CrdMeteoRecord::pressure_ )
            .def_readwrite( "temperature", &tio::CrdMeteoRecord::temperature_ )
            .def_readwrite( "humidity", &tio::CrdMeteoRecord::humidity_ )
            .def_readwrite( "origin_of_values", &tio::CrdMeteoRecord::originOfValues_ );

    py::class_< tio::CrdRangeSupplementRecord >( m,
                                                 "CrdRangeSupplementRecord",
                                                 "Optional record12 corrections. Atmospheric delay is one-way ps; centre-of-mass length is "
                                                 "one-way metres; time bias is seconds. H4 flags indicate which corrections were applied." )
            .def( py::init<>( ) )
            .def_readwrite( "second_of_day", &tio::CrdRangeSupplementRecord::secondOfDay_ )
            .def_readwrite( "system_configuration_id", &tio::CrdRangeSupplementRecord::systemConfigurationId_ )
            .def_readwrite( "tropospheric_correction_ps", &tio::CrdRangeSupplementRecord::troposphericCorrectionPs_ )
            .def_readwrite( "center_of_mass_correction", &tio::CrdRangeSupplementRecord::centerOfMassCorrection_ )
            .def_readwrite( "neutral_density_filter", &tio::CrdRangeSupplementRecord::neutralDensityFilter_ )
            .def_readwrite( "time_bias_applied", &tio::CrdRangeSupplementRecord::timeBiasApplied_ )
            .def_readwrite( "range_rate", &tio::CrdRangeSupplementRecord::rangeRate_ );

    py::class_< tio::CrdCalibrationRecord >( m,
                                             "CrdCalibrationRecord",
                                             "Calibration summary40 or detail41, with delays and RMS in ps and target distance in metres; "
                                             "these values are metadata, not reapplied corrections." )
            .def( py::init<>( ) )
            .def_readwrite( "record_type", &tio::CrdCalibrationRecord::recordType_ )
            .def_readwrite( "second_of_day", &tio::CrdCalibrationRecord::secondOfDay_ )
            .def_readwrite( "data_type", &tio::CrdCalibrationRecord::dataType_ )
            .def_readwrite( "system_configuration_id", &tio::CrdCalibrationRecord::systemConfigurationId_ )
            .def_readwrite( "number_recorded", &tio::CrdCalibrationRecord::numberRecorded_ )
            .def_readwrite( "number_used", &tio::CrdCalibrationRecord::numberUsed_ )
            .def_readwrite( "target_distance", &tio::CrdCalibrationRecord::targetDistance_ )
            .def_readwrite( "system_delay_ps", &tio::CrdCalibrationRecord::systemDelayPs_ )
            .def_readwrite( "delay_shift_ps", &tio::CrdCalibrationRecord::delayShiftPs_ )
            .def_readwrite( "rms_ps", &tio::CrdCalibrationRecord::rmsPs_ )
            .def_readwrite( "skew", &tio::CrdCalibrationRecord::skew_ )
            .def_readwrite( "kurtosis", &tio::CrdCalibrationRecord::kurtosis_ )
            .def_readwrite( "peak_minus_mean_ps", &tio::CrdCalibrationRecord::peakMinusMeanPs_ )
            .def_readwrite( "calibration_type", &tio::CrdCalibrationRecord::calibrationType_ )
            .def_readwrite( "shift_type", &tio::CrdCalibrationRecord::shiftType_ )
            .def_readwrite( "detector_channel", &tio::CrdCalibrationRecord::detectorChannel_ )
            .def_readwrite( "span", &tio::CrdCalibrationRecord::span_ )
            .def_readwrite( "return_rate", &tio::CrdCalibrationRecord::returnRate_ );

    py::class_< tio::CrdPassData >(
            m, "CrdPassData", "Ordered measurement, weather, supplement and calibration records; duplicate epochs are retained." )
            .def( py::init<>( ) )
            .def_readwrite( "full_rate_data", &tio::CrdPassData::fullRateData_ )
            .def_readwrite( "normal_point_data", &tio::CrdPassData::normalPointData_ )
            .def_readwrite( "meteorological_data", &tio::CrdPassData::meteorologicalData_ )
            .def_readwrite( "range_supplement_data", &tio::CrdPassData::rangeSupplementData_ )
            .def_readwrite( "calibration_data", &tio::CrdPassData::calibrationData_ );

    py::class_< tio::CrdPass >( m, "CrdPass", "One session with its station, target, system and correction metadata." )
            .def( py::init<>( ) )
            .def_readwrite( "configuration", &tio::CrdPass::configuration_ )
            .def_readwrite( "data", &tio::CrdPass::data_ );

    m.def( "convert_crd_two_way_time_of_flight_to_slr_range",
           &tio::convertCrdTwoWayTimeOfFlightToSlrRange,
           py::arg( "two_way_time_of_flight" ),
           R"doc(
         Convert a CRD two-way light time (seconds) into one-way SLR range (meters).
      )doc" );

    m.def( "read_slr_data",
           &tio::readCrdFiles,
           py::arg( "file_names" ),
           R"doc(
         Read SLR CRD files and concatenate all parsed passes.

         The Consolidated Laser Ranging Data Format (CRD) is the ILRS format for
         full-rate, sampled engineering, and normal-point laser-ranging data:
         https://ilrs.gsfc.nasa.gov/data_and_products/formats/crd.html.
         This reader returns parsed ``CrdPass`` objects. The supporting
         functions in this module can then group the passes by target or station
         and extract normal-point or full-rate range measurements.

         Parameters
         ----------
         file_names : list[str]
             Paths to CRD files.

         Returns
         -------
         list[CrdPass]
             Parsed CRD passes from all input files.
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
         Extract normal-point measurements from one CRD pass as ordered ``(Time, one_way_range_m)`` pairs. Epochs retain the original event and station UTC reference; duplicate epochs are preserved.
      )doc" );

    m.def( "extract_normal_point_measurements_from_passes",
           py::overload_cast< const std::vector< tio::CrdPass >& >( &tio::extractNormalPointMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
         Extract normal-point measurements from multiple CRD passes as ordered ``(Time, one_way_range_m)`` pairs. Epochs retain the original event and station UTC reference; duplicate epochs are preserved.
      )doc" );

    m.def( "extract_full_rate_measurements",
           py::overload_cast< const tio::CrdPass& >( &tio::extractFullRateMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
         Extract full-rate measurements from one CRD pass as ordered ``(Time, one_way_range_m)`` pairs. Epochs retain the original event and station UTC reference; duplicate epochs are preserved.
      )doc" );

    m.def( "extract_full_rate_measurements_from_passes",
           py::overload_cast< const std::vector< tio::CrdPass >& >( &tio::extractFullRateMeasurements ),
           py::arg( "pass_data" ),
           R"doc(
         Extract full-rate measurements from multiple CRD passes as ordered ``(Time, one_way_range_m)`` pairs. Epochs retain the original event and station UTC reference; duplicate epochs are preserved.
      )doc" );

    m.def( "get_station_wavelengths",
           &tio::getStationWavelengths,
           py::arg( "grouped_data" ),
           R"doc(
         Extract station transmit wavelengths in metres; reject differing wavelengths across any C0 systems.
      )doc" );
}

}  // namespace slr

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
