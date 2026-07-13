#if TUDATPY_ENABLE_DETAILED_PYBIND11_ERRORS
#define PYBIND11_DETAILED_ERROR_MESSAGES
#endif

#include "expose_generic_text_file.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "tudat/io/readTrackingTxtFile.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{

namespace data_access
{

namespace tracking
{

namespace generic_text_file
{

void expose_generic_text_file( py::module& m )
{
    py::enum_< tio::TrackingDataType >( m, "TrackingDataType", R"doc(No documentation available.)doc" )
            .value( "year", tio::TrackingDataType::year, R"doc(No documentation available.)doc" )
            .value( "month", tio::TrackingDataType::month, R"doc(No documentation available.)doc" )
            .value( "day", tio::TrackingDataType::day, R"doc(No documentation available.)doc" )
            .value( "hour", tio::TrackingDataType::hour, R"doc(No documentation available.)doc" )
            .value( "minute", tio::TrackingDataType::minute, R"doc(No documentation available.)doc" )
            .value( "second", tio::TrackingDataType::second, R"doc(No documentation available.)doc" )
            .value( "observation_time_scale", tio::TrackingDataType::observation_time_scale, R"doc(No documentation available.)doc" )
            .value( "file_name", tio::TrackingDataType::file_name, R"doc(No documentation available.)doc" )
            .value( "n_way_light_time", tio::TrackingDataType::n_way_light_time, R"doc(No documentation available.)doc" )
            .value( "light_time_measurement_delay",
                    tio::TrackingDataType::light_time_measurement_delay,
                    R"doc(No documentation available.)doc" )
            .value( "light_time_measurement_accuracy",
                    tio::TrackingDataType::light_time_measurement_accuracy,
                    R"doc(No documentation available.)doc" )
            .value( "dsn_transmitting_station_nr",
                    tio::TrackingDataType::dsn_transmitting_station_nr,
                    R"doc(No documentation available.)doc" )
            .value( "dsn_receiving_station_nr", tio::TrackingDataType::dsn_receiving_station_nr, R"doc(No documentation available.)doc" )
            .value( "observation_body", tio::TrackingDataType::observation_body, R"doc(No documentation available.)doc" )
            .value( "observed_body", tio::TrackingDataType::observed_body, R"doc(No documentation available.)doc" )
            .value( "spacecraft_id", tio::TrackingDataType::spacecraft_id, R"doc(No documentation available.)doc" )
            .value( "spacecraft_name", tio::TrackingDataType::spacecraft_name, R"doc(No documentation available.)doc" )
            .value( "planet_nr", tio::TrackingDataType::planet_nr, R"doc(No documentation available.)doc" )
            .value( "tdb_reception_time_j2000", tio::TrackingDataType::tdb_reception_time_j2000, R"doc(No documentation available.)doc" )
            .value( "utc_reception_time_j2000", tio::TrackingDataType::utc_reception_time_j2000, R"doc(No documentation available.)doc" )
            .value( "utc_ramp_referencee_j2000", tio::TrackingDataType::utc_ramp_referencee_j2000, R"doc(No documentation available.)doc" )
            .value( "tdb_spacecraft_j2000", tio::TrackingDataType::tdb_spacecraft_j2000, R"doc(No documentation available.)doc" )
            .value( "x_planet_frame", tio::TrackingDataType::x_planet_frame, R"doc(No documentation available.)doc" )
            .value( "y_planet_frame", tio::TrackingDataType::y_planet_frame, R"doc(No documentation available.)doc" )
            .value( "z_planet_frame", tio::TrackingDataType::z_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vx_planet_frame", tio::TrackingDataType::vx_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vy_planet_frame", tio::TrackingDataType::vy_planet_frame, R"doc(No documentation available.)doc" )
            .value( "vz_planet_frame", tio::TrackingDataType::vz_planet_frame, R"doc(No documentation available.)doc" )
            .value( "residual_de405", tio::TrackingDataType::residual_de405, R"doc(No documentation available.)doc" )
            .value( "spacecraft_transponder_delay",
                    tio::TrackingDataType::spacecraft_transponder_delay,
                    R"doc(No documentation available.)doc" )
            .value( "uplink_frequency", tio::TrackingDataType::uplink_frequency, R"doc(No documentation available.)doc" )
            .value( "downlink_frequency", tio::TrackingDataType::downlink_frequency, R"doc(No documentation available.)doc" )
            .value( "signal_to_noise", tio::TrackingDataType::signal_to_noise, R"doc(No documentation available.)doc" )
            .value( "spectral_max", tio::TrackingDataType::spectral_max, R"doc(No documentation available.)doc" )
            .value( "doppler_measured_frequency",
                    tio::TrackingDataType::doppler_measured_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_averaged_frequency",
                    tio::TrackingDataType::doppler_averaged_frequency,
                    R"doc(No documentation available.)doc" )
            .value( "doppler_integration_time", tio::TrackingDataType::doppler_integration_time, R"doc(No documentation available.)doc" )
            .value( "doppler_base_frequency", tio::TrackingDataType::doppler_base_frequency, R"doc(No documentation available.)doc" )
            .value( "doppler_noise", tio::TrackingDataType::doppler_noise, R"doc(No documentation available.)doc" )
            .value( "doppler_bandwidth", tio::TrackingDataType::doppler_bandwidth, R"doc(No documentation available.)doc" )
            .value( "receiving_station_name", tio::TrackingDataType::receiving_station_name, R"doc(No documentation available.)doc" )
            .value( "transmitting_station_name", tio::TrackingDataType::transmitting_station_name, R"doc(No documentation available.)doc" )
            .value( "time_tag_delay", tio::TrackingDataType::time_tag_delay )
            .value( "sample_number", tio::TrackingDataType::sample_number )
            .value( "utc_day_of_year", tio::TrackingDataType::utc_day_of_year )
            .value( "reference_body_distance", tio::TrackingDataType::reference_body_distance )
            .value( "transmission_frequency_constant_term", tio::TrackingDataType::transmission_frequency_constant_term )
            .value( "transmission_frequency_linear_term", tio::TrackingDataType::transmission_frequency_linear_term )
            .value( "doppler_predicted_frequency_hz", tio::TrackingDataType::doppler_predicted_frequency_hz )
            .value( "doppler_troposphere_correction", tio::TrackingDataType::doppler_troposphere_correction )
            .value( "scan_nr", tio::TrackingDataType::scan_nr )
            .export_values( );

    py::enum_< tio::TrackingTxtFileReadFilterType >( m, "TrackingTxtFileReadFilterType" )
            .value( "no_tracking_txt_file_filter", tio::TrackingTxtFileReadFilterType::no_tracking_txt_file_filter )
            .value( "ifms_tracking_txt_file_filter", tio::TrackingTxtFileReadFilterType::ifms_tracking_txt_file_filter )
            .export_values( );

    py::class_< tio::TrackingTxtFileContents, std::shared_ptr< tio::TrackingTxtFileContents > >( m, "TrackingTxtFileContents" )
            .def( py::init< const std::string, const std::vector< std::string >, const char, const std::string >( ),
                  py::arg( "file_name" ),
                  py::arg( "column_types" ),
                  py::arg( "comment_symbol" ) = '#',
                  py::arg( "value_separators" ) = ",:\t " )
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
            .def_property_readonly( "column_field_types", &tio::TrackingTxtFileContents::getRawColumnTypes )
            .def_property_readonly( "double_datamap", &tio::TrackingTxtFileContents::getDoubleDataMap )
            .def_property_readonly( "raw_datamap", &tio::TrackingTxtFileContents::getRawDataMap )
            .def_property_readonly( "num_rows", &tio::TrackingTxtFileContents::getNumRows );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ) = '#',
           py::arg( "value_separators" ) = ",:\t ",
           py::arg( "ignore_omitted_columns" ) = false,
           py::arg( "data_filter_method" ) = tio::no_tracking_txt_file_filter );
}

}  // namespace generic_text_file

}  // namespace tracking

}  // namespace data_access

}  // namespace tudatpy
