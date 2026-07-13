/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_tracking.h"

#include "scalarTypes.h"
#include "tudat/io/trackingData.h"
#include "tudat/io/readTrackingTxtFile.h"

#include "odf/expose_odf.h"
#include "fdets/expose_fdets.h"
#include "ifms/expose_ifms.h"

namespace py = pybind11;
namespace tio = tudat::input_output;

namespace tudatpy
{
namespace data
{
namespace tracking
{

void expose_tracking( py::module& m )
{
    auto odf_submodule = m.def_submodule( "odf" );
    odf::expose_odf( odf_submodule );

    auto fdets_submodule = m.def_submodule( "fdets" );
    fdets::expose_fdets( fdets_submodule );

    auto ifms_submodule = m.def_submodule( "ifms" );
    ifms::expose_ifms( ifms_submodule );

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

    py::class_< tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE > > >( m, "TrackingData", R"doc(
 TrackingData Class container.
    )doc" )
            .def( py::init< const std::string,
                            const tudat::data::PlainLinkDefinition&,
                            const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                            const std::vector< TIME_TYPE >,
                            const std::string,
                            const std::string >( ),
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "observations" ),
                  py::arg( "epochs" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "time_scale" ) = "TDB",
                  R"doc(Creates a TrackingData object.)doc" )
            .def_property_readonly( "observable_type", &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType )
            .def_property_readonly( "link_ends", &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds )
            .def_property_readonly( "observations", &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations )
            .def_property_readonly( "epochs", &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationEpochs )
            .def_property_readonly( "reference_link_end", &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd )
            .def_property_readonly( "time_scale", &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeScale )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::string >(
                          &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (string type), (optional).)doc" )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::vector< std::string > >(
                          &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (optional).)doc" )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const double >(
                          &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (double type), (optional).)doc" )
            .def( "add_ancillary_settings",
                  py::overload_cast< const std::string, const std::vector< double > >(
                          &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::addAncillarySettings ),
                  py::arg( "ancillarySettingsType" ),
                  py::arg( "ancillarySettingsValue" ),
                  R"doc(Adds ancillary settings to the TrackingData object (double vector type), (optional).)doc" )
            .def( "get_ancillary_settings_string_vector",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsStringVector,
                  R"doc(Returns the ancillary settings in the TrackingData object as a vector of strings.)doc" )
            .def( "get_ancillary_settings_double",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDouble,
                  R"doc(Returns the ancillary settings in the TrackingData object as a double.)doc" )
            .def( "get_ancillary_settings_double_vector",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDoubleVector,
                  R"doc(Returns the ancillary settings in the TrackingData object as a double vector.)doc" )
            .def( "get_ancillary_settings_double_vector",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getAncillarySettingsDoubleVector,
                  R"doc(Returns the ancillary settings in the TrackingData object as a double vector.)doc" )
            .def( "set_observation_weights",
                  ( &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationWeights ),
                  py::arg( "observationWeights" ),
                  R"doc(Adds (sets) observation weights to the TrackingData object (optional).
						It also perfroms some necessary consistency checks  (size, pre-existence) on the weights.)doc" )
            .def( "reset_single_observation_weight",
                  ( &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationWeight ),
                  py::arg( "index" ),
                  py::arg( "observationWeight" ),
                  R"doc(Allows resetting of a single observation weight (specified by index i)  into the TrackingData object.
						It also perfroms some necessary consistency checks.)doc" )
            .def( "get_observation_weights",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeights,
                  R"doc(Returns the list (of arrays) of observation weights stored in the TrackingData object.)doc" )
            .def( "get_concatenated_observation_weights",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationWeightsVector,
                  R"doc(Returns the concatenated list of observation weights stored in the TrackingData object.)doc" )
            .def( "set_observation_corrections",
                  &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationCorrections,
                  py::arg( "observationCorrections" ),
                  R"doc(Adds corrections to the TrackingData object (optional).
					It also perfroms some necessary consistency checks.)doc" )
            .def( "reset_single_observation_correction",
                  ( &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::setSingleObservationCorrection ),
                  py::arg( "index" ),
                  py::arg( "observationCorrection" ),
                  R"doc(Allows resetting of a single observation correction (specified by index i)  into the TrackingData object.
						It also perfroms some necessary consistency checks.)doc" )
            .def( "remove_single_observation_entry",
                  ( &tudat::data::TrackingData< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationEntry ),
                  py::arg( "index" ),
                  R"doc(Removes a single observation entry from the TrackingData object.)doc" );
}

}  // namespace tracking
}  // namespace data
}  // namespace tudatpy
