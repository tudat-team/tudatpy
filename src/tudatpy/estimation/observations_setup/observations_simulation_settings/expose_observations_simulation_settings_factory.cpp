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
#include "expose_observations_simulation_settings_bindings.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{
namespace observations_simulation_settings
{

void expose_observation_simulation_settings_factory_bindings( py::module& m )
{
    m.def( "observation_settings_from_collection",
           &tss::getObservationSimulationSettingsFromObservations< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "bodies" ),
           R"doc(No documentation found.)doc" );

    m.def( "change_simulation_settings_observable_types",
           &tom::changeObservableTypesOfObservationSimulationSettings< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_simulation_settings" ),
           py::arg( "replacement_observable_types" ) =
                   std::map< tom::ObservableType, tom::ObservableType >{
                           { tom::dsn_n_way_averaged_doppler, tom::n_way_differenced_range },
                           { tom::dsn_one_way_averaged_doppler, tom::one_way_differenced_range } },
           R"doc(No documentation found.)doc" );

    //    m.def("create_odf_observation_simulation_settings_list",
    //          &tom::createOdfObservationSimulationSettingsList<
    //          STATE_SCALAR_TYPE, TIME_TYPE >,
    //          py::arg("observed_observation_collection"),
    //          get_docstring("create_odf_observation_simulation_settings_list").c_str()
    //          );

    // #   Observation Model Settings --> Observation Simulator #
    m.def( "create_observation_simulators",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&, const tss::SystemOfBodies& >(
                   &tom::createObservationSimulators< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function for creating observation simulator objects.

 Function for creating observation simulator objects from observation settings.
 Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.


 Parameters
 ----------
 observation_settings : List[ ObservationModelSettings ]
     List of settings objects, each object defining the observation model settings for one combination of observable and link geometry that is to be simulated.

 bodies : :class:`~tudatpy.dynamics.environment.SystemOfBodies`
     Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

 Returns
 -------
 List[ :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` ]
     List of :class:`~tudatpy.estimation.observable_models.observables_simulation.ObservationSimulator` objects, each object hosting the functionality for simulating one combination of observable type and link geometry.

 Examples
 --------
 .. code-block:: python

     from tudatpy.estimation.observations_setup import observations_simulation_settings

     # Create bodies
     bodies = ...
     # Define parameters settings
     observation_settings = ...
     # Create observation simulators
     observation_simulators = observations_simulation_settings.create_observation_simulators(observation_settings, bodies)

 This code snippet closely follows what is done in: The following snippet closely follows what is done in: `Galilean Moons State Estimation Example <https://github.com/tudat-team/tudatpy-examples/blob/master/estimation/galilean_moons_state_estimation.ipynb>`_.



     )doc" );
}

}  // namespace observations_simulation_settings
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
