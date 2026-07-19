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

#include "tudat/simulation/estimation_setup/simulateObservationsLegacy.h"
#include "tudat/simulation/estimation_setup/simulatePseudoObservations.h"

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace
{

const char* legacyObservationWrapperDeprecationGuide =
        "https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-dataset-deprecation.html";

void warnLegacyObservationWrapperInterface( const std::string& interfaceName, const std::string& replacementApi )
{
    const std::string message = interfaceName + " is deprecated and kept only for backwards compatibility. Use " + replacementApi +
            " instead. API reference: https://py.api.tudat.space/en/latest/estimation/observations_setup/observations_wrapper.html#"
            "tudatpy.estimation.observations_setup.observations_wrapper." +
            replacementApi + ". Migration guide: " + legacyObservationWrapperDeprecationGuide;
    if( PyErr_WarnEx( PyExc_DeprecationWarning, message.c_str( ), 1 ) < 0 )
    {
        throw pybind11::error_already_set( );
    }
}

}  // namespace

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_wrapper
{

void expose_observations_wrapper_simulation_bindings( py::module& m )
{
    m.def(
            "create_pseudo_observations_and_models",
            []( const tss::SystemOfBodies& bodies,
                const std::vector< std::string >& observedBodies,
                const std::vector< std::string >& centralBodies,
                const TIME_TYPE initialTime,
                const TIME_TYPE finalTime,
                const TIME_TYPE timeStep ) {
                warnLegacyObservationWrapperInterface( "create_pseudo_observations_and_models",
                                                       "create_pseudo_observation_dataset_and_models" );
                return tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE >(
                        bodies, observedBodies, centralBodies, initialTime, finalTime, timeStep );
            },
            py::arg( "bodies" ),
            py::arg( "observed_bodies" ),
            py::arg( "central_bodies" ),
            py::arg( "initial_time" ),
            py::arg( "final_time" ),
            py::arg( "time_step" ),
            R"doc(No documentation found.)doc" );

    m.def(
            "create_pseudo_observations_and_models_from_observation_times",
            []( const tss::SystemOfBodies& bodies,
                const std::vector< std::string >& observedBodies,
                const std::vector< std::string >& centralBodies,
                const std::vector< TIME_TYPE >& observationTimes ) {
                warnLegacyObservationWrapperInterface( "create_pseudo_observations_and_models_from_observation_times",
                                                       "create_pseudo_observation_dataset_and_models_from_observation_times" );
                return tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE >(
                        bodies, observedBodies, centralBodies, observationTimes );
            },
            py::arg( "bodies" ),
            py::arg( "observed_bodies" ),
            py::arg( "central_bodies" ),
            py::arg( "observation_times" ),
            R"doc(No documentation found.)doc" );

    m.def( "create_pseudo_observation_dataset_and_models",
           py::overload_cast< const tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const TIME_TYPE,
                              const TIME_TYPE,
                              const TIME_TYPE >( &tss::simulatePseudoObservationDataset< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "time_step" ),
           R"doc(
Create pseudo-observation simulators and an :class:`~tudatpy.estimation.observations.ObservationDataset`.

This is the dataset-backed counterpart of
``create_pseudo_observations_and_models``. It creates the observation models for
the requested observed bodies and returns the generated pseudo-observations in
the new dataset representation.

Returns
-------
tuple[list[ObservationSimulator], tudatpy.estimation.observations.ObservationDataset]
    Observation simulators and the generated dataset.
)doc" );

    m.def( "create_pseudo_observation_dataset_and_models_from_observation_times",
           py::overload_cast< const tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const std::vector< TIME_TYPE >& >( &tss::simulatePseudoObservationDataset< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "observation_times" ),
           R"doc(
Create pseudo-observation simulators and a dataset at explicit observation times.

This is the dataset-backed counterpart of
``create_pseudo_observations_and_models_from_observation_times``.

Returns
-------
tuple[list[ObservationSimulator], tudatpy.estimation.observations.ObservationDataset]
    Observation simulators and the generated dataset.
)doc" );

    m.def(
            "set_existing_observations",
            []( const std::map< tom::ObservableType,
                                std::pair< tom::LinkEnds,
                                           std::pair< std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >,
                                                      std::vector< TIME_TYPE > > > >& observations,
                const tom::LinkEndType referenceLinkEnd,
                const std::map< tom::ObservableType, std::shared_ptr< tom::ObservationAncillarySimulationSettings > >&
                        ancillarySettingsPerObservable ) {
                warnLegacyObservationWrapperInterface( "set_existing_observations", "set_existing_observation_dataset" );
                return tss::setExistingObservations< STATE_SCALAR_TYPE, TIME_TYPE >(
                        observations, referenceLinkEnd, ancillarySettingsPerObservable );
            },
            py::arg( "observations" ),
            py::arg( "reference_link_end" ),
            py::arg( "ancillary_settings_per_observatble" ) =
                    std::map< tom::ObservableType, std::shared_ptr< tom::ObservationAncillarySimulationSettings > >( ) );

    m.def( "set_existing_observation_dataset",
           &tss::setExistingObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observations" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings_per_observatble" ) =
                   std::map< tom::ObservableType, std::shared_ptr< tom::ObservationAncillarySimulationSettings > >( ),
           R"doc(
Create an observation dataset from existing observation values.

Parameters
----------
observations : dict
    Existing observations grouped by observable type and link definition.
reference_link_end : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Link end used as the reference for the observation times.
ancillary_settings_per_observatble : dict, optional
    Ancillary settings to attach to each observable type.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    Dataset containing the provided observations.
)doc" );

    m.def(
            "simulate_observations",
            []( const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >& simulationSettings,
                const std::vector< std::shared_ptr< tom::ObservationSimulatorBase< STATE_SCALAR_TYPE, TIME_TYPE > > >&
                        observationSimulators,
                const tss::SystemOfBodies& bodies ) {
                warnLegacyObservationWrapperInterface( "simulate_observations", "simulate_observation_dataset" );
                return tss::simulateObservations< STATE_SCALAR_TYPE, TIME_TYPE >( simulationSettings, observationSimulators, bodies );
            },
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

    m.def( "simulate_observation_dataset",
           &tss::simulateObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "simulation_settings" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(
Simulate observations and return the dataset representation.

This function uses the supplied simulation settings, observation simulators and
bodies to create an :class:`~tudatpy.estimation.observations.ObservationDataset`.

Parameters
----------
simulation_settings : list[ObservationSimulationSettings]
    Settings defining what observations should be simulated.
observation_simulators : list[ObservationSimulator]
    Observation simulators that evaluate the requested observables.
bodies : tudatpy.dynamics.environment.SystemOfBodies
    Physical environment used by the observation models.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    Dataset containing all simulated observations, residual placeholders,
    weights, dependent variables and metadata.
)doc" );

    m.def(
            "single_type_observation_collection",
            []( const tom::ObservableType observableType,
                const tom::LinkDefinition& linkEnds,
                const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& observationsList,
                const std::vector< TIME_TYPE > timesList,
                const tom::LinkEndType referenceLinkEnd,
                const std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings ) {
                warnLegacyObservationWrapperInterface( "single_type_observation_collection", "single_type_observation_dataset" );
                return tom::createManualObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >(
                        observableType, linkEnds, observationsList, timesList, referenceLinkEnd, ancillarySettings );
            },
            py::arg( "observable_type" ),
            py::arg( "link_ends" ),
            py::arg( "observations_list" ),
            py::arg( "times_list" ),
            py::arg( "reference_link_end" ),
            py::arg( "ancillary_settings" ) = nullptr,
            R"doc(No documentation found.)doc" );

    m.def(
            "single_type_observation_dataset",
            []( const tom::ObservableType observableType,
                const tom::LinkDefinition& linkEnds,
                const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& observationsList,
                const std::vector< TIME_TYPE > timesList,
                const tom::LinkEndType referenceLinkEnd,
                const std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings ) {
                std::shared_ptr< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > > dataset =
                        std::make_shared< tom::ObservationDataset< STATE_SCALAR_TYPE, TIME_TYPE > >( );
                dataset->addObservationSet( observableType,
                                            linkEnds,
                                            observationsList,
                                            timesList,
                                            referenceLinkEnd,
                                            std::vector< Eigen::VectorXd >( ),
                                            nullptr,
                                            ancillarySettings );
                return dataset;
            },
            py::arg( "observable_type" ),
            py::arg( "link_ends" ),
            py::arg( "observations_list" ),
            py::arg( "times_list" ),
            py::arg( "reference_link_end" ),
            py::arg( "ancillary_settings" ) = nullptr,
            R"doc(
Create a single-set observation dataset from existing values.

Parameters
----------
observable_type : tudatpy.estimation.observable_models_setup.model_settings.ObservableType
    Observable type of the provided observations.
link_ends : tudatpy.estimation.observable_models_setup.links.LinkDefinition
    Link definition for the observation set.
observations_list : list[numpy.ndarray]
    Vector-valued observations to store.
times_list : list[float]
    Observation time for each observation.
reference_link_end : tudatpy.estimation.observable_models_setup.links.LinkEndType
    Link end used as the reference for the observation times.
ancillary_settings : tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings, optional
    Ancillary settings to attach to the observation set.

Returns
-------
tudatpy.estimation.observations.ObservationDataset
    Dataset containing one observation set.
)doc" );
}

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
