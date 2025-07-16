/*    Copyright (c) 2010-2021, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_observation.h"

#include <pybind11/functional.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace tuc = tudat::unit_conversions;
namespace ti = tudat::interpolators;



namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation_setup
{
namespace observation
{

void expose_observation_setup( py::module& m )
{

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // FREQUENCIES
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::enum_< tom::FrequencyBands >( m, "FrequencyBands", R"doc(No documentation found.)doc" )
            .value( "s_band", tom::FrequencyBands::s_band )
            .value( "x_band", tom::FrequencyBands::x_band )
            .value( "ka_band", tom::FrequencyBands::ka_band )
            .value( "ku_band", tom::FrequencyBands::ku_band );

    // Create wrapper function
    py::cpp_function getDsnDefaultTurnaroundRatios_wrapper = []( tudat::observation_models::FrequencyBands band1,
                                                                 tudat::observation_models::FrequencyBands band2 ) {
        return tom::getDsnDefaultTurnaroundRatios( band1, band2 );
    };

    m.def( "dsn_default_turnaround_ratios",
           &tom::getDsnDefaultTurnaroundRatios,
           py::arg( "uplink_band" ),
           py::arg( "downlink_band" ),
           R"doc(No documentation found.)doc" );

    m.def( "cassini_turnaround_ratios",
           &tom::getCassiniTurnaroundRatio,
           py::arg( "uplink_band" ),
           py::arg( "downlink_band" ),
           R"doc(No documentation found.)doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // ODF OBSERVATIONS
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::class_< tom::ProcessedOdfFileContents< TIME_TYPE >, std::shared_ptr< tom::ProcessedOdfFileContents< TIME_TYPE > > >(
            m, "ProcessedOdfFileContents", R"doc(No documentation found.)doc" )
            .def_property_readonly( "ground_station_names",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getGroundStationsNames,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "processed_observable_types",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getProcessedObservableTypes,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "start_and_end_time",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getStartAndEndTime,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "ignored_odf_observable_types",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getIgnoredRawOdfObservableTypes,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "ignored_ground_stations",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getIgnoredGroundStations,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly(
                    "raw_odf_data", &tom::ProcessedOdfFileContents< TIME_TYPE >::getRawOdfData, R"doc(No documentation found.)doc" )
            .def( "define_antenna_id",
                  py::overload_cast< const std::string&, const std::string& >(
                          &tom::ProcessedOdfFileContents< TIME_TYPE >::defineSpacecraftAntennaId ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "antenna_name" ),
                  R"doc(No documentation found.)doc" );

    m.def( "process_odf_data_multiple_files",
           py::overload_cast< const std::vector< std::string >&,
                              const std::string&,
                              const bool,
                              const std::map< std::string, Eigen::Vector3d >& >( &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg( "earth_fixed_ground_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "process_odf_data_single_file",
           py::overload_cast< const std::string&, const std::string&, const bool, const std::map< std::string, Eigen::Vector3d >& >(
                   &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_name" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg( "earth_fixed_ground_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
           R"doc(No documentation found.)doc" );


    m.def( "set_odf_information_in_bodies",
           &tom::setOdfInformationInBodies< TIME_TYPE >,
           py::arg( "processed_odf_file" ),
           py::arg( "bodies" ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           py::arg( "turnaround_ratio_function" ) = getDsnDefaultTurnaroundRatios_wrapper,
           R"doc(No documentation found.)doc" );

    m.def( "create_odf_observed_observation_collection",
           &tom::createOdfObservedObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "processed_odf_file" ),
           py::arg( "observable_types_to_process" ),
           py::arg( "start_and_end_times_to_process" ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_odf_files",
           &tom::createOdfObservedObservationCollectionFromFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "odf_file_names" ),
           py::arg( "target_name" ),
           py::arg( "verbose_output" ) = true,
           py::arg( "earth_fixed_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_ifms_files",
           &tom::createIfmsObservedObservationCollectionFromFiles< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_names" ),
           py::arg( "bodies" ),
           py::arg( "target_name" ),
           py::arg( "ground_station_name" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "apply_troposphere_correction" ) = true,
           py::arg( "earth_fixed_station_positions" ) = tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_multi_station_ifms_files",
           &tom::createMultiStationIfmsObservedObservationCollectionFromFiles< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_names" ),
           py::arg( "bodies" ),
           py::arg( "target_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "apply_troposphere_correction" ) = true,
           py::arg( "earth_fixed_station_positions" ) = tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "observations_from_fdets_files",
           &tom::createFdetsObservedObservationCollectionFromFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_name" ),
           py::arg( "base_frequency" ),
           py::arg( "column_types" ),
           py::arg( "target_name" ),
           py::arg( "transmitting_station_name" ),
           py::arg( "receiving_station_name" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "earth_fixed_station_positions" ) = tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(No documentation found.)doc" );

    m.def( "create_compressed_doppler_collection",
           &tom::createCompressedDopplerCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "compression_ratio" ),
           py::arg( "minimum_number_of_observations" ) = 10,
           R"doc(No documentation found.)doc" );

    //    m.def("create_odf_observation_simulation_settings_list",
    //          &tom::createOdfObservationSimulationSettingsList<
    //          STATE_SCALAR_TYPE, TIME_TYPE >,
    //          py::arg("observed_observation_collection"),
    //          get_docstring("create_odf_observation_simulation_settings_list").c_str()
    //          );

    

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // Tracking Txt OBSERVATIONS
    /////////////////////////////////////////////////////////////////////////////////////////////////

    m.def( "create_tracking_txtfile_observation_collection",
           py::overload_cast< const std::shared_ptr< tudat::input_output::TrackingTxtFileContents >,
                              const std::string,
                              const std::vector< tom::ObservableType >,
                              const std::map< std::string, Eigen::Vector3d >,
                              const tom::ObservationAncilliarySimulationSettings& >(
                   &tom::createTrackingTxtFileObservationCollection< double, TIME_TYPE > ),
           py::arg( "raw_tracking_txtfile_contents" ),
           py::arg( "spacecraft_name" ),
           py::arg( "observable_types_to_process" ) = std::vector< tom::ObservableType >( ),
           py::arg( "earth_fixed_ground_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
           py::arg( "ancillary_settings" ) = tom::ObservationAncilliarySimulationSettings( ),
           R"doc(No documentation found.)doc" );


}


}  // namespace observation
}  // namespace estimation_setup
}  // namespace numerical_simulation
}  // namespace tudatpy