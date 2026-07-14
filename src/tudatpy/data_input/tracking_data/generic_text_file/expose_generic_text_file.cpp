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

namespace data_input
{

namespace tracking_data
{

namespace generic_text_file
{

void expose_generic_text_file( py::module& m )
{
    py::enum_< tio::TrackingDataType >( m, "TrackingDataType", R"doc(Recognized data fields in generic tracking text files.)doc" )
            .value( "year", tio::TrackingDataType::year, R"doc(Generic tracking text-file data field.)doc" )
            .value( "month", tio::TrackingDataType::month, R"doc(Generic tracking text-file data field.)doc" )
            .value( "day", tio::TrackingDataType::day, R"doc(Generic tracking text-file data field.)doc" )
            .value( "hour", tio::TrackingDataType::hour, R"doc(Generic tracking text-file data field.)doc" )
            .value( "minute", tio::TrackingDataType::minute, R"doc(Generic tracking text-file data field.)doc" )
            .value( "second", tio::TrackingDataType::second, R"doc(Generic tracking text-file data field.)doc" )
            .value( "observation_time_scale",
                    tio::TrackingDataType::observation_time_scale,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "file_name", tio::TrackingDataType::file_name, R"doc(Generic tracking text-file data field.)doc" )
            .value( "n_way_light_time", tio::TrackingDataType::n_way_light_time, R"doc(Generic tracking text-file data field.)doc" )
            .value( "light_time_measurement_delay",
                    tio::TrackingDataType::light_time_measurement_delay,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "light_time_measurement_accuracy",
                    tio::TrackingDataType::light_time_measurement_accuracy,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "dsn_transmitting_station_nr",
                    tio::TrackingDataType::dsn_transmitting_station_nr,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "dsn_receiving_station_nr",
                    tio::TrackingDataType::dsn_receiving_station_nr,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "observation_body", tio::TrackingDataType::observation_body, R"doc(Generic tracking text-file data field.)doc" )
            .value( "observed_body", tio::TrackingDataType::observed_body, R"doc(Generic tracking text-file data field.)doc" )
            .value( "spacecraft_id", tio::TrackingDataType::spacecraft_id, R"doc(Generic tracking text-file data field.)doc" )
            .value( "spacecraft_name", tio::TrackingDataType::spacecraft_name, R"doc(Generic tracking text-file data field.)doc" )
            .value( "planet_nr", tio::TrackingDataType::planet_nr, R"doc(Generic tracking text-file data field.)doc" )
            .value( "tdb_reception_time_j2000",
                    tio::TrackingDataType::tdb_reception_time_j2000,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "utc_reception_time_j2000",
                    tio::TrackingDataType::utc_reception_time_j2000,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "utc_ramp_referencee_j2000",
                    tio::TrackingDataType::utc_ramp_referencee_j2000,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "tdb_spacecraft_j2000", tio::TrackingDataType::tdb_spacecraft_j2000, R"doc(Generic tracking text-file data field.)doc" )
            .value( "x_planet_frame", tio::TrackingDataType::x_planet_frame, R"doc(Generic tracking text-file data field.)doc" )
            .value( "y_planet_frame", tio::TrackingDataType::y_planet_frame, R"doc(Generic tracking text-file data field.)doc" )
            .value( "z_planet_frame", tio::TrackingDataType::z_planet_frame, R"doc(Generic tracking text-file data field.)doc" )
            .value( "vx_planet_frame", tio::TrackingDataType::vx_planet_frame, R"doc(Generic tracking text-file data field.)doc" )
            .value( "vy_planet_frame", tio::TrackingDataType::vy_planet_frame, R"doc(Generic tracking text-file data field.)doc" )
            .value( "vz_planet_frame", tio::TrackingDataType::vz_planet_frame, R"doc(Generic tracking text-file data field.)doc" )
            .value( "residual_de405", tio::TrackingDataType::residual_de405, R"doc(Generic tracking text-file data field.)doc" )
            .value( "spacecraft_transponder_delay",
                    tio::TrackingDataType::spacecraft_transponder_delay,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "uplink_frequency", tio::TrackingDataType::uplink_frequency, R"doc(Generic tracking text-file data field.)doc" )
            .value( "downlink_frequency", tio::TrackingDataType::downlink_frequency, R"doc(Generic tracking text-file data field.)doc" )
            .value( "signal_to_noise", tio::TrackingDataType::signal_to_noise, R"doc(Generic tracking text-file data field.)doc" )
            .value( "spectral_max", tio::TrackingDataType::spectral_max, R"doc(Generic tracking text-file data field.)doc" )
            .value( "doppler_measured_frequency",
                    tio::TrackingDataType::doppler_measured_frequency,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "doppler_averaged_frequency",
                    tio::TrackingDataType::doppler_averaged_frequency,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "doppler_integration_time",
                    tio::TrackingDataType::doppler_integration_time,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "doppler_base_frequency",
                    tio::TrackingDataType::doppler_base_frequency,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "doppler_noise", tio::TrackingDataType::doppler_noise, R"doc(Generic tracking text-file data field.)doc" )
            .value( "doppler_bandwidth", tio::TrackingDataType::doppler_bandwidth, R"doc(Generic tracking text-file data field.)doc" )
            .value( "receiving_station_name",
                    tio::TrackingDataType::receiving_station_name,
                    R"doc(Generic tracking text-file data field.)doc" )
            .value( "transmitting_station_name",
                    tio::TrackingDataType::transmitting_station_name,
                    R"doc(Generic tracking text-file data field.)doc" )
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

    py::enum_< tio::TrackingTxtFileReadFilterType >(
            m, "TrackingTxtFileReadFilterType", R"doc(Filter mode used when reading generic tracking text files.)doc" )
            .value( "no_tracking_txt_file_filter", tio::TrackingTxtFileReadFilterType::no_tracking_txt_file_filter )
            .value( "ifms_tracking_txt_file_filter", tio::TrackingTxtFileReadFilterType::ifms_tracking_txt_file_filter )
            .export_values( );

    py::class_< tio::TrackingTxtFileContents, std::shared_ptr< tio::TrackingTxtFileContents > >( m,
                                                                                                 "TrackingTxtFileContents",
                                                                                                 R"doc(
Container for parsed generic tracking text-file contents.

Parameters
----------
file_name : str
    Path to the tracking text file.
column_types : list[str]
    Names of the data fields represented by each column.
comment_symbol : str, default="#"
    Character used to mark comment lines.
value_separators : str, default=",:\t "
    Characters interpreted as value separators.
)doc" )
            .def( py::init< const std::string, const std::vector< std::string >, const char, const std::string >( ),
                  py::arg( "file_name" ),
                  py::arg( "column_types" ),
                  py::arg( "comment_symbol" ) = '#',
                  py::arg( "value_separators" ) = ",:\t ",
                  R"doc(
Read a generic tracking text file.

Parameters
----------
file_name : str
    Path to the tracking text file.
column_types : list[str]
    Names of the data fields represented by each column.
comment_symbol : str, default="#"
    Character used to mark comment lines.
value_separators : str, default=",:\t "
    Characters interpreted as value separators.
)doc" )
            .def( "add_metadata_val",
                  py::overload_cast< tio::TrackingDataType, double >( &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "value" ),
                  R"doc(
Add a scalar metadata value to the parsed file contents.

Parameters
----------
tracking_data_type : TrackingDataType
    Metadata field to set.
value : float
    Metadata value.

Returns
-------
None
)doc" )
            .def( "get_available_datatypes",
                  &tio::TrackingTxtFileContents::getAllAvailableDataTypes,
                  R"doc(
Return all data fields available in the parsed contents.

Returns
-------
list[TrackingDataType]
    Data fields for which data or metadata is available.
)doc" )
            .def( "add_metadata_str",
                  py::overload_cast< tio::TrackingDataType, const std::string& >( &tio::TrackingTxtFileContents::addMetaData ),
                  py::arg( "tracking_data_type" ),
                  py::arg( "str_value" ),
                  R"doc(
Add a string metadata value to the parsed file contents.

Parameters
----------
tracking_data_type : TrackingDataType
    Metadata field to set.
str_value : str
    Metadata value.

Returns
-------
None
)doc" )
            .def_property_readonly( "column_field_types",
                                    &tio::TrackingTxtFileContents::getRawColumnTypes,
                                    R"doc(list[TrackingDataType]: Data field assigned to each input column.)doc" )
            .def_property_readonly( "double_datamap",
                                    &tio::TrackingTxtFileContents::getDoubleDataMap,
                                    R"doc(dict[TrackingDataType, list[float]]: Numeric parsed data by field.)doc" )
            .def_property_readonly( "raw_datamap",
                                    &tio::TrackingTxtFileContents::getRawDataMap,
                                    R"doc(dict[TrackingDataType, list[str]]: Raw string parsed data by field.)doc" )
            .def_property_readonly( "num_rows", &tio::TrackingTxtFileContents::getNumRows, R"doc(int: Number of parsed data rows.)doc" );

    m.def( "read_tracking_txt_file",
           &tio::createTrackingTxtFileContents,
           py::arg( "file_name" ),
           py::arg( "column_types" ),
           py::arg( "comment_symbol" ) = '#',
           py::arg( "value_separators" ) = ",:\t ",
           py::arg( "ignore_omitted_columns" ) = false,
           py::arg( "data_filter_method" ) = tio::no_tracking_txt_file_filter,
           R"doc(
Read a generic tracking text file.

Parameters
----------
file_name : str
    Path to the tracking text file.
column_types : list[str]
    Names of the data fields represented by each column.
comment_symbol : str, default="#"
    Character used to mark comment lines.
value_separators : str, default=",:\t "
    Characters interpreted as value separators.
ignore_omitted_columns : bool, default=False
    Whether columns omitted from ``column_types`` are ignored.
data_filter_method : TrackingTxtFileReadFilterType, default=no_tracking_txt_file_filter
    Optional reader-specific filter mode.

Returns
-------
TrackingTxtFileContents
    Parsed tracking text-file contents.
)doc" );
}

}  // namespace generic_text_file

}  // namespace tracking_data

}  // namespace data_input

}  // namespace tudatpy
