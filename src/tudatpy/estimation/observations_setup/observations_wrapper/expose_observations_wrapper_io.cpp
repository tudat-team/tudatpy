/*    Copyright (c) 2010-2021, Delft University of Technology
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
#include "expose_observations_wrapper_bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"

#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/simulation/environment_setup/defaultGroundStationSettings.h"

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_wrapper
{

void expose_observations_wrapper_io_bindings( py::module &m )
{
    py::cpp_function getDsnDefaultTurnaroundRatios_wrapper = []( tudat::observation_models::FrequencyBands band1,
                                                                 tudat::observation_models::FrequencyBands band2 ) {
        return tom::getDsnDefaultTurnaroundRatios( band1, band2 );
    };

    py::class_< tom::ProcessedOdfFileContents< TIME_TYPE >, std::shared_ptr< tom::ProcessedOdfFileContents< TIME_TYPE > > >(
            m, "ProcessedOdfFileContents", R"doc(
        Class containing processed ODF data.
        )doc" )
            .def_property_readonly( "ground_station_names",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getGroundStationsNames,
                                    R"doc(
        Get the names of the ground stations included in the ODF files.

        Returns
        -------
        list[str]
            List of ground station names.
        )doc" )
            .def_property_readonly( "processed_observable_types",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getProcessedObservableTypes,
                                    R"doc(
        Get the observable types in the ODF files.

        Returns
        -------
        list[tudatpy.estimation.observable_models_setup.model_settings.ObservableType]
            List of observable types.
        )doc" )
            .def_property_readonly( "start_and_end_time",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getStartAndEndTime,
                                    R"doc(
        Get pair of < start time, end time > of the data contained in the ODF files.

        Returns
        -------
        tuple[float, float]
            Start and end time of the data.
        )doc" )
            .def_property_readonly( "ignored_odf_observable_types",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getIgnoredRawOdfObservableTypes,
                                    R"doc(
        Return ODF observable types IDs (as per TRK-2-18) that were not included in the processed data.

        Returns
        -------
        list[tudatpy.io.odf.OdfDataType]
            List of ignored ODF observable type IDs.
        )doc" )
            .def_property_readonly( "ignored_ground_stations",
                                    &tom::ProcessedOdfFileContents< TIME_TYPE >::getIgnoredGroundStations,
                                    R"doc(
        Return ground stations for which observations were not included in the processed data (due to absence of ramp tables).

        Returns
        -------
        list[str]
            List of ignored ground station names.
        )doc" )
            .def_property_readonly( "raw_odf_data", &tom::ProcessedOdfFileContents< TIME_TYPE >::getRawOdfData, R"doc(
        Return the raw ODF data.

        Returns
        -------
        list[tudatpy.io.OdfRawFileContents]
            List of raw ODF data objects.
        )doc" )
            .def( "define_antenna_id",
                  py::overload_cast< const std::string &, const std::string & >(
                          &tom::ProcessedOdfFileContents< TIME_TYPE >::defineSpacecraftAntennaId ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "antenna_name" ),
                  R"doc(
        Define the antenna ID for a given spacecraft.
        )doc" );

    m.def( "process_odf_data_multiple_files",
           py::overload_cast< const std::vector< std::string > &,
                              const std::string &,
                              const bool,
                              const std::map< std::string, Eigen::Vector3d > & >( &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg_v( "earth_fixed_ground_station_positions",
                      tss::getApproximateDsnGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions()" ),
           R"doc(
        Process multiple ODF files.

        Parameters
        ----------
        file_names : list[str]
            List of ODF file names.
        spacecraft_name : str
            Name of the spacecraft.
        verbose : bool, optional
            Whether to print warnings, by default True.
        earth_fixed_ground_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions of DSN ground stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions`) will be used.

        Returns
        -------
        tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
            Processed ODF file contents.
        )doc" );

    m.def( "process_odf_data_single_file",
           py::overload_cast< const std::string &, const std::string &, const bool, const std::map< std::string, Eigen::Vector3d > & >(
                   &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_name" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg_v( "earth_fixed_ground_station_positions",
                      tss::getApproximateDsnGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions()" ),
           R"doc(
        Process a single ODF file.

        Parameters
        ----------
        file_name : str
            ODF file name.
        spacecraft_name : str
            Name of the spacecraft.
        verbose : bool, optional
            Whether to print warnings, by default True.
        earth_fixed_ground_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions of DSN ground stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions`) will be used.

        Returns
        -------
        tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
            Processed ODF file contents.
        )doc" );

    m.def( "set_odf_information_in_bodies",
           &tom::setOdfInformationInBodies< TIME_TYPE >,
           py::arg( "processed_odf_file" ),
           py::arg( "bodies" ),
           py::arg( "body_with_ground_stations_name" ) = "Earth",
           py::arg( "turnaround_ratio_function" ) = getDsnDefaultTurnaroundRatios_wrapper,
           R"doc(
        Sets the ODF information required for simulating observations into the system of bodies.

        This includes:
        - Setting the transmitting frequencies objects in the ground stations
        - Setting the turnaround ratio in the spacecraft

        Parameters
        ----------
        processed_odf_file : tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
            Processed ODF file contents.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies.
        body_with_ground_stations_name : str, optional
            Name of the body in which the ground stations are located, by default "Earth".
        turnaround_ratio_function : function, optional
            Function returning the turnaround ratio as a function of the uplink and downlink bands.
        )doc" );

    m.def( "create_odf_observed_observation_collection",
           &tom::createOdfObservedObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "processed_odf_file" ),
           py::arg( "observable_types_to_process" ),
           py::arg( "start_and_end_times_to_process" ),
           py::arg( "allow_duplicate_observations_within_single_set" ) = true,
           R"doc(
        Creates an observation collection containing the provided ODF data.

        Only the specified observable types are loaded from the processed ODF data into the observation collection.

        Parameters
        ----------
        processed_odf_file : tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
            Processed ODF data.
        observable_types_to_process : list[tudatpy.estimation.observable_models_setup.model_settings.ObservableType]
            Observable types to process.
        start_and_end_times_to_process : tuple[float, float]
            Start and end times of the data to process.
        allow_duplicate_observations_within_single_set: bool
            Determines if duplicate observations should be erased on SingleObservationSet level before ObservationCollection creation, default is True.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );

    m.def( "observations_from_odf_files",
           &tom::createOdfObservedObservationCollectionFromFile< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "bodies" ),
           py::arg( "odf_file_names" ),
           py::arg( "target_name" ),
           py::arg( "verbose_output" ) = true,
           py::arg_v( "earth_fixed_station_positions",
                      tss::getApproximateDsnGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions()" ),
           py::arg( "allow_duplicate_observations_within_single_set" ) = true,
           R"doc(
        Create an observation collection from ODF files.

        This function processes ODF files, sets the required information in the bodies, and creates an observation collection.

        Parameters
        ----------
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies.
        odf_file_names : list[str]
            List of ODF file names.
        target_name : str
            Name of the target spacecraft.
        verbose_output : bool, optional
            Whether to print verbose output during processing, by default True.
        earth_fixed_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions of DSN ground stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions`) will be used.
        allow_duplicate_observations_within_single_set: bool
            Determines if duplicate observations should be erased on SingleObservationSet level before ObservationCollection creation, default is True.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );

    m.def( "observations_from_ifms_files",
           &tom::createIfmsObservedObservationCollectionFromFiles< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_names" ),
           py::arg( "bodies" ),
           py::arg( "target_name" ),
           py::arg( "ground_station_name" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "apply_troposphere_correction" ) = true,
           py::arg_v( "earth_fixed_station_positions",
                      tss::getCombinedApproximateGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions()" ),
           py::arg( "remove_invalid_lines" ) = true,
           R"doc(
        Create an observation collection from IFMS files for a single station.

        This function processes IFMS files, sets the required information in the bodies, and creates an observation collection.

        Parameters
        ----------
        ifms_file_names : list[str]
            List of IFMS file names.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies.
        target_name : str
            Name of the target spacecraft.
        ground_station_name : str
            Name of the ground station.
        reception_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Reception frequency band.
        transmission_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Transmission frequency band.
        apply_troposphere_correction : bool, optional
            Whether to apply troposphere correction, by default True.
        earth_fixed_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions` will be used.
        remove_invalid_lines : bool, optional
            Boolean (default true) defining whether a line is skipped if the transmit frequency, observed frequency, or troposphere correction is undefined

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );

    m.def( "observations_from_multi_station_ifms_files",
           &tom::createMultiStationIfmsObservedObservationCollectionFromFiles< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "ifms_file_names" ),
           py::arg( "bodies" ),
           py::arg( "target_name" ),
           py::arg( "ground_station_names" ),
           py::arg( "reception_band" ),
           py::arg( "transmission_band" ),
           py::arg( "apply_troposphere_correction" ) = true,
           py::arg_v( "earth_fixed_station_positions",
                      tss::getCombinedApproximateGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions()" ),
           py::arg( "remove_invalid_lines" ) = true,
           R"doc(
        Create an observation collection from IFMS files for multiple stations.

        This function processes IFMS files, sets the required information in the bodies, and creates an observation collection.

        Parameters
        ----------
        ifms_file_names : list[str]
            List of IFMS file names.
        bodies : tudatpy.dynamics.environment.SystemOfBodies
            System of bodies.
        target_name : str
            Name of the target spacecraft.
        ground_station_names : list[str]
            List of ground station names, must be same size as ifms_file_names.
        reception_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Reception frequency band.
        transmission_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Transmission frequency band.
        apply_troposphere_correction : bool, optional
            Whether to apply troposphere correction, by default True.
        earth_fixed_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions` will be used.
        remove_invalid_lines : bool, optional
            Boolean (default true) defining whether a line is skipped if the transmit frequency, observed frequency, or troposphere correction is undefined

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );

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
           py::arg_v( "earth_fixed_station_positions",
                      tss::getCombinedApproximateGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions()" ),
           R"doc(
        Create an observation collection from an FDETS file.

        This function processes an FDETS file and creates an observation collection.

        Parameters
        ----------
        ifms_file_name : str
            FDETS file name.
        base_frequency : float
            Base frequency for Doppler observables.
        column_types : list[str]
            List of column types in the FDETS file.
        target_name : str
            Name of the target spacecraft.
        transmitting_station_name : str
            Name of the transmitting station.
        receiving_station_name : str
            Name of the receiving station.
        reception_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Reception frequency band.
        transmission_band : tudatpy.estimation.observations_setup.ancillary_settings.FrequencyBands
            Transmission frequency band.
        earth_fixed_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_radio_telescope_positions` will be used.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );

    m.def( "create_compressed_doppler_collection",
           &tom::createCompressedDopplerCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "compression_ratio" ),
           py::arg( "minimum_number_of_observations" ) = 10,
           py::arg( "max_arc_gap" ) = 300.0,
           R"doc(
        Create a compressed Doppler observation collection.

        This function takes a collection of Doppler observations, splits them into arcs, and compresses each arc by averaging observations.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The original observation collection containing Doppler data.
        compression_ratio : int
            The number of observations to average into a single compressed observation.
        minimum_number_of_observations : int, optional
            The minimum number of observations required in an arc to be considered for compression, by default 10.
        max_arc_gap : float, optional
            Maximum time gap (in seconds) between consecutive observations before splitting into a new arc, by default 300.0.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with compressed Doppler data.
        )doc" );

    m.def( "create_tracking_txtfile_observation_collection",
           py::overload_cast< const std::shared_ptr< tudat::input_output::TrackingTxtFileContents >,
                              const std::string,
                              const std::vector< tom::ObservableType >,
                              const std::map< std::string, Eigen::Vector3d >,
                              const tom::ObservationAncillarySimulationSettings & >(
                   &tom::createTrackingTxtFileObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "raw_tracking_txtfile_contents" ),
           py::arg( "spacecraft_name" ),
           py::arg( "observable_types_to_process" ) = std::vector< tom::ObservableType >( ),
           py::arg_v( "earth_fixed_ground_station_positions",
                      tss::getApproximateDsnGroundStationPositions( ),
                      "tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions()" ),
           py::arg( "ancillary_settings" ) = tom::ObservationAncillarySimulationSettings( ),
           R"doc(
        Create an observation collection from raw tracking file data.

        Parameters
        ----------
        raw_tracking_txtfile_contents : tudatpy.io.TrackingTxtFileContents
            The raw tracking file contents.
        spacecraft_name : str
            Name of the spacecraft.
        observable_types_to_process : list[tudatpy.estimation.observable_models_setup.model_settings.ObservableType], optional
            List of observable types to process. If empty, all available types are processed.
        earth_fixed_ground_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame. If none is provided, the approximate positions of DSN ground stations (as given by :func:`~tudatpy.dynamics.environment_setup.ground_station.get_approximate_dsn_ground_station_positions`) will be used.
        ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings, optional
            Ancillary settings for the observations.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );
}

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
