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

#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/simulatePseudoObservations.h"

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

void expose_observations_wrapper_simulation_bindings( py::module& m )
{
    m.def( "create_pseudo_observations_and_models",
           py::overload_cast< const tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const TIME_TYPE,
                              const TIME_TYPE,
                              const TIME_TYPE >( &tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "initial_time" ),
           py::arg( "final_time" ),
           py::arg( "time_step" ),
           R"doc(

         Create relative-position model settings and pseudo-observations from the ephemerides in ``bodies``.

         ``observed_bodies`` are paired with ``central_bodies``. Epochs start one hour after ``initial_time`` and stop
         before one hour before ``final_time``, with spacing ``time_step``; times are in seconds since J2000.

         Parameters
         ----------
         bodies : SystemOfBodies
             System of bodies defining the physical environment.
         observed_bodies : list[str]
             Names of the bodies whose ephemerides provide the observations.
         central_bodies : list[str]
             Names of the reference bodies, paired with the propagated or observed bodies.
         initial_time : Time
             Initial epoch, in seconds since J2000.
         final_time : Time
             Final epoch, in seconds since J2000.
         time_step : Time
             Sampling interval, in seconds.

         Returns
         -------
         tuple[list[ObservationModelSettings], ObservationCollection]
             Pair containing the relative-position model settings and the simulated observation collection.

      )doc" );

    m.def( "create_pseudo_observations_and_models_from_observation_times",
           py::overload_cast< const tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const std::vector< TIME_TYPE > >( &tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "observation_times" ),
           R"doc(

         Create relative-position model settings and pseudo-observations from ephemerides at ``observation_times``.

         ``observed_bodies`` are paired with ``central_bodies`` in the supplied system of ``bodies``. Epochs are in
         seconds since J2000.

         Parameters
         ----------
         bodies : SystemOfBodies
             System of bodies defining the physical environment.
         observed_bodies : list[str]
             Names of the bodies whose ephemerides provide the observations.
         central_bodies : list[str]
             Names of the reference bodies, paired with the propagated or observed bodies.
         observation_times : list[Time]
             Observation epochs, in seconds since J2000.

         Returns
         -------
         tuple[list[ObservationModelSettings], ObservationCollection]
             Pair containing the relative-position model settings and the simulated observation collection.

      )doc" );

    m.def( "set_existing_observations",
           &tss::setExistingObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observations" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings_per_observatble" ) =
                   std::map< tom::ObservableType, std::shared_ptr< tom::ObservationAncillarySimulationSettings > >( ),
           R"doc(

         Create an ObservationCollection from existing measurements grouped by observable type.

         ``observations`` maps each observable type to its link ends and a pair of measurement-vector and epoch lists.
         ``reference_link_end`` identifies the time reference. ``ancillary_settings_per_observatble`` optionally
         supplies ancillary data for each observable type.

         Parameters
         ----------
         observations : dict[ObservableType, tuple[dict[LinkEndType, LinkEndId], tuple[list[numpy.ndarray[numpy.float64[m, 1]]], list[Time]]]]
             Measurements grouped by observable type; each entry contains the link ends and paired lists of measurement vectors and epochs.
         reference_link_end : LinkEndType
             Link end at which the observation epochs are defined.
         ancillary_settings_per_observatble : dict[ObservableType, ObservationAncillarySimulationSettings], optional
             Ancillary data for each observable type.

         Returns
         -------
         ObservationCollection
             Collection containing the supplied measurements and epochs.

      )doc" );

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
                              const std::shared_ptr< tom::ObservationAncillarySimulationSettings > >(
                   &tom::createManualObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations_list" ),
           py::arg( "times_list" ),
           py::arg( "reference_link_end" ),
           py::arg_v( "ancillary_settings", std::shared_ptr< tom::ObservationAncillarySimulationSettings >( ), "None" ),
           R"doc(

         Create an ObservationCollection for one ``observable_type`` and ``link_ends`` definition.

         ``observations_list`` contains measurement vectors and ``times_list`` the corresponding epochs in seconds
         since J2000, referenced to ``reference_link_end``. ``ancillary_settings`` optionally supplies observable-
         specific supporting data.

         Parameters
         ----------
         observable_type : ObservableType
             Type of observable to which the settings or measurements apply.
         link_ends : LinkDefinition
             Definition of the bodies and reference points participating in the observation link.
         observations_list : list[numpy.ndarray[numpy.float64[m, 1]]]
             Measurement vectors, with one vector for each epoch.
         times_list : list[Time]
             Observation epochs in seconds since J2000, paired with the measurement vectors.
         reference_link_end : LinkEndType
             Link end at which the observation epochs are defined.
         ancillary_settings : ObservationAncillarySimulationSettings, optional
             Optional observable-specific ancillary data.

         Returns
         -------
         ObservationCollection
             Collection containing the supplied measurements for the single observable and link definition.

      )doc" );
}

}  // namespace observations_wrapper
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
