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
#include "expose_observations_bindings.h"

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
namespace observations
{

void expose_observations_simulation_bindings( py::module& m )
{
    m.def( "simulate_pseudo_observations",
           py::overload_cast< const tss::SystemOfBodies&,
                              const std::vector< std::string >&,
                              const std::vector< std::string >&,
                              const std::vector< TIME_TYPE > >( &tss::simulatePseudoObservations< TIME_TYPE, STATE_SCALAR_TYPE > ),
           py::arg( "bodies" ),
           py::arg( "observed_bodies" ),
           py::arg( "central_bodies" ),
           py::arg( "observation_times" ),
           R"doc(No documentation found.)doc" );

    m.def( "create_observation_collection_from_arrays",
           &tss::setExistingObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observations" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings_per_observable" ) =
                   std::map< tom::ObservableType, std::shared_ptr< tom::ObservationAncillarySimulationSettings > >( ) );

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

    m.def( "create_single_type_observation_collection_from_arrays",
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
           py::arg( "ancillary_settings" ) = nullptr,
           R"doc(No documentation found.)doc" );
}

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
