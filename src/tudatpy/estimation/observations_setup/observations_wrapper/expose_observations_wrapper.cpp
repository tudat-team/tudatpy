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
#include "expose_observations_wrapper.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"

#include "tudat/simulation/estimation_setup/processOdfFile.h"
#include "tudat/simulation/estimation_setup/processTrackingTxtFile.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"

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

void expose_observations_wrapper( py::module& m )
{
    py::module_::import( "tudatpy.estimation.observations" ).attr( "ObservationCollection" );

    // Create wrapper function
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
        list[tudatpy.astro.ObservableType]
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
                  py::overload_cast< const std::string&, const std::string& >(
                          &tom::ProcessedOdfFileContents< TIME_TYPE >::defineSpacecraftAntennaId ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "antenna_name" ),
                  R"doc(
        Define the antenna ID for a given spacecraft.
        )doc" );

    m.def( "process_odf_data_multiple_files",
           py::overload_cast< const std::vector< std::string >&,
                              const std::string&,
                              const bool,
                              const std::map< std::string, Eigen::Vector3d >& >( &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_names" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg( "earth_fixed_ground_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
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
            Map with approximate positions of ground stations in Earth-fixed frame.

        Returns
        -------
        tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
            Processed ODF file contents.
        )doc" );

    m.def( "process_odf_data_single_file",
           py::overload_cast< const std::string&, const std::string&, const bool, const std::map< std::string, Eigen::Vector3d >& >(
                   &tom::processOdfData< TIME_TYPE > ),
           py::arg( "file_name" ),
           py::arg( "spacecraft_name" ),
           py::arg( "verbose" ) = true,
           py::arg( "earth_fixed_ground_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
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
            Map with approximate positions of ground stations in Earth-fixed frame.

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
        bodies : tudatpy.dynamics.SystemOfBodies
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
           R"doc(
        Creates an observation collection containing the provided ODF data.

        Only the specified observable types are loaded from the processed ODF data into the observation collection.

        Parameters
        ----------
        processed_odf_file : tudatpy.estimation.observations_setup.observations_wrapper.ProcessedOdfFileContents
            Processed ODF data.
        observable_types_to_process : list[tudatpy.astro.ObservableType]
            Observable types to process.
        start_and_end_times_to_process : tuple[float, float]
            Start and end times of the data to process.

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
           py::arg( "earth_fixed_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
           R"doc(
        Create an observation collection from ODF files.

        This function processes ODF files, sets the required information in the bodies, and creates an observation collection.

        Parameters
        ----------
        bodies : tudatpy.dynamics.SystemOfBodies
            System of bodies.
        odf_file_names : list[str]
            List of ODF file names.
        target_name : str
            Name of the target spacecraft.
        verbose_output : bool, optional
            Whether to print verbose output during processing, by default True.
        earth_fixed_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame.

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
           py::arg( "earth_fixed_station_positions" ) = tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(
        Create an observation collection from IFMS files for a single station.

        This function processes IFMS files, sets the required information in the bodies, and creates an observation collection.

        Parameters
        ----------
        ifms_file_names : list[str]
            List of IFMS file names.
        bodies : tudatpy.dynamics.SystemOfBodies
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
            Map with approximate positions of ground stations in Earth-fixed frame.

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
           py::arg( "earth_fixed_station_positions" ) = tss::getCombinedApproximateGroundStationPositions( ),
           R"doc(
        Create an observation collection from IFMS files for multiple stations.

        This function processes IFMS files, sets the required information in the bodies, and creates an observation collection.

        Parameters
        ----------
        ifms_file_names : list[str]
            List of IFMS file names.
        bodies : tudatpy.dynamics.SystemOfBodies
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
            Map with approximate positions of ground stations in Earth-fixed frame.

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
           py::arg( "earth_fixed_station_positions" ) = tss::getCombinedApproximateGroundStationPositions( ),
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
            Map with approximate positions of ground stations in Earth-fixed frame.

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
                              const tom::ObservationAncilliarySimulationSettings& >(
                   &tom::createTrackingTxtFileObservationCollection< double, TIME_TYPE > ),
           py::arg( "raw_tracking_txtfile_contents" ),
           py::arg( "spacecraft_name" ),
           py::arg( "observable_types_to_process" ) = std::vector< tom::ObservableType >( ),
           py::arg( "earth_fixed_ground_station_positions" ) = tss::getApproximateDsnGroundStationPositions( ),
           py::arg( "ancillary_settings" ) = tom::ObservationAncilliarySimulationSettings( ),
           R"doc(
        Create an observation collection from raw tracking file data.

        Parameters
        ----------
        raw_tracking_txtfile_contents : tudatpy.io.TrackingTxtFileContents
            The raw tracking file contents.
        spacecraft_name : str
            Name of the spacecraft.
        observable_types_to_process : list[tudatpy.astro.ObservableType], optional
            List of observable types to process. If empty, all available types are processed.
        earth_fixed_ground_station_positions : dict[str, numpy.ndarray[3]], optional
            Map with approximate positions of ground stations in Earth-fixed frame.
        ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncilliarySimulationSettings, optional
            Ancillary settings for the observations.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            Observation collection.
        )doc" );

    m.def( "create_pseudo_observations_and_models",
    py::overload_cast<  const tss::SystemOfBodies &,
                        const std::vector< std::string > &,
                        const std::vector< std::string > &,
                        const TIME_TYPE,
                        const TIME_TYPE,
                        const TIME_TYPE > ( &tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "time_step" ),
           R"doc(No documentation found.)doc" );


    m.def( "create_pseudo_observations_and_models_from_observation_times",
            py::overload_cast<  const tss::SystemOfBodies &,
                    const std::vector< std::string > &,
                    const std::vector< std::string > &,
                    const std::vector< TIME_TYPE > >
                    ( &tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE > ),
       py::arg( "bodies" ),
       py::arg( "observed_bodies" ),
       py::arg( "central_bodies" ),
       py::arg( "observation_times" ),
       R"doc(No documentation found.)doc" );


    m.def( "set_existing_observations",
           &tss::setExistingObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observations" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancilliary_settings_per_observatble" ) =
                   std::map< tom::ObservableType, std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >( ) );

    m.def( "simulate_observations",
           &tss::simulateObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "simulation_settings" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(

 Function to simulate observations.

 Function to simulate observations from set observation simulators and observation simulator settings.
 Automatically iterates over all provided observation simulators, generating the full set of simulated observations.


 Parameters
 ----------
 observation_to_simulate : List[ :class:`ObservationSimulationSettings` ]
     List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.

 observation_simulators : List[ :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` ]
     List of :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.

 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

 Returns
 -------
 :class:`~tudatpy.estimation.observations.ObservationCollection`
     Object collecting all products of the observation simulation.






     )doc" );

    m.def( "single_type_observation_collection",
           py::overload_cast< const tom::ObservableType,
                              const tom::LinkDefinition&,
                              const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                              const std::vector< TIME_TYPE >,
                              const tom::LinkEndType,
                              const std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
                   &tom::createManualObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations_list" ),
           py::arg( "times_list" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancilliary_settings" ) = nullptr,
           R"doc(No documentation found.)doc" );
}

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
