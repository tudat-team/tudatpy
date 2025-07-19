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
#include "expose_observations_simulation_settings.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"
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

void expose_observations_simulation_settings( py::module& m )
{
    py::class_< tss::ObservationSimulationSettings< TIME_TYPE >, std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >(
            m,
            "ObservationSimulationSettings",
            R"doc(
         Base class for defining settings for simulated observations.
      )doc" )
            .def_property( "viability_settings_list",
                           &tss::ObservationSimulationSettings< TIME_TYPE >::getViabilitySettingsList,
                           &tss::ObservationSimulationSettings< TIME_TYPE >::setViabilitySettingsList,
                           R"doc(
         viability_settings_list : List[ ObservationViabilitySettings ], default = [ ]) -
         Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
      )doc" )
            .def_property( "noise_function",
                           &tss::ObservationSimulationSettings< TIME_TYPE >::getObservationNoiseFunction,
                           py::overload_cast< const std::function< double( const double ) >& >(
                                   &tss::ObservationSimulationSettings< TIME_TYPE >::setObservationNoiseFunction ),
                           R"doc(
         noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None -
         Function providing the observation noise as a function of observation time (can be constant or time-dependent), default is None.
      )doc" );
    //            .def_property("observable_type",
    //                         &tss::ObservationSimulationSettings<double>::getObservableType,
    //                         &tss::ObservationSimulationSettings<double>::setObservableType,
    //                         get_docstring("ObservationSimulationSettings.observable_type").c_str()
    //                         )
    //            .def_property_readonly("link_ends",
    //                         &tss::ObservationSimulationSettings<double>::getLinkEnds,
    //                         get_docstring("ObservationSimulationSettings.link_ends").c_str()
    //                         );

    py::class_< tss::TabulatedObservationSimulationSettings< TIME_TYPE >,
                std::shared_ptr< tss::TabulatedObservationSimulationSettings< TIME_TYPE > >,
                tss::ObservationSimulationSettings< TIME_TYPE > >( m,
                                                                   "TabulatedObservationSimulationSettings",
                                                                   R"doc(

         Class for defining settings for simulating observations at a predefined set of times.

         Class for defining settings for simulating observations at a predefined set of times.
         This type defines predefined time epochs at which applicable observations are to be simulated, stored in a rigid, "tabulated" form.
         Some observation times may be discarded due to the use of viability settings.
         Instances of this class are typically created via the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings`
         and :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings_list` functions.

         Associated base class: :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservationSettings`

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of a TabulatedObservationSimulationSettings object
             import numpy as np
             from tudatpy.astro.time_conversion import DateTime
             from tudatpy.estimation.observable_models_setup import links, model_settings
             from tudatpy.estimation.observations_setup import observations_simulation_settings

             # Set simulation start and end epochs
             simulation_start_epoch = DateTime(2000, 1, 1).epoch()
             simulation_end_epoch   = DateTime(2000, 1, 4).epoch()

             # Define the uplink link ends for one-way observable
             link_ends = dict()
             link_ends[links.transmitter] = links.body_origin_link_end_id("Earth")
             link_ends[links.receiver] = links.body_origin_link_end_id("Delfi-C3")

             # Create LinkDefinition Object and set observation settings for each link/observable
             link_definition = links.LinkDefinition(link_ends)
             observation_settings_list = [model_settings.one_way_doppler_instantaneous(link_definition)]

             # Define observation simulation times (separated by steps of 1 minute)
             observation_times = np.arange(simulation_start_epoch, simulation_end_epoch, 60.0)

             # Create TabulatedObservationSimulationSettings object
             tabulated_observation_simulation_settings = observations_simulation_settings.tabulated_simulation_settings(
                 model_settings.one_way_instantaneous_doppler_type,
                 link_definition,
                 observation_times
             )

             # Show that this is indeed a TabulatedObservationSimulationSettings object
             print(tabulated_observation_simulation_settings)



      )doc" );
    

    m.def( "tabulated_simulation_settings",
           &tss::tabulatedObservationSimulationSettings< TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "simulation_times" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           py::arg( "noise_function" ) = nullptr,
           py::arg( "ancilliary_settings" ) = nullptr,
           R"doc(

 Function for creating settings object for observation simulation, using a predefined list of observation times.

 Function for creating single simulation settings object, using a predefined list of observation times.
 The list of resulting observations may be reduced compared to the ``simulation_times`` should be provided here, as
 only observations that meet the viability settings are retained during observation simulation (these may be
 provide directly here through the ``viability_settings`` input, or added later to the resulting settings object).


 Parameters
 ----------
 observable_type : :class:`ObservableType`
     Observable type of which observations are to be simulated.
 link_ends : LinkDefinition
     Link ends for which observations are to be simulated.
 simulation_times : List[float]
     List of times at which to perform the observation simulation.
 reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType.receiver`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference time for the observation.
 viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
     Function providing the observation noise factors as a function of observation time.
 Returns
 -------
 :class:`TabulatedObservationSimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` class.






     )doc" );

    m.def( "tabulated_simulation_settings_list",
           &tss::createTabulatedObservationSimulationSettingsList< TIME_TYPE >,
           py::arg( "link_ends_per_observable" ),
           py::arg( "simulation_times" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           R"doc(

 Function for creating a list of settings object for observation simulation, using a predefined list of observation times.

 Function for creating multiple tabulated observation simulation settings objects in a list. This function is
 equivalent to calling the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings` repeatedly, with the different
 observables and link definition provided here through `link_ends_per_observable`.
 During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


 Parameters
 ----------
 link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
     Link geometry per observable type of which observations are to be simulated.
 simulation_times : List[ float ]
     List of times at which to perform the observation simulation.
 reference_link_end_type : :class:`LinkEndType`, default = :class:`LinkEndType.receiver`
     Defines the link end (via the :class:`LinkEndType`) which is used as a reference for observation times.
     The single link end specified here will be considered as the reference link end for all simulation settings object created in the function call.

 viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
     The single settings list given here will be considered as potential viability settings for all simulation settings object created in the function call.

 Returns
 -------
 List[ TabulatedObservationSimulationSettings ]
     List of :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` objects.






     )doc" );

    m.def( "continuous_arc_simulation_settings",
           &tss::perArcObservationSimulationSettings< TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "interval_between_observations" ),
           py::arg( "arc_limiting_constraints" ),
           py::arg( "minimum_arc_duration" ),
           py::arg( "maximum_arc_duration" ),
           py::arg( "minimum_time_between_arcs" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "additional_viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           py::arg( "noise_function" ) = nullptr,
           R"doc(

 Function for creating settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

 Function for creating settings object for observation simulation. Unlike the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.tabulated_simulation_settings`
 function, the resulting settings do not define the observation times explicitly. Instead, this settings object determines the observation times adaptively during the
 simulation of the observation, with the requirement that observations should be simulated over a set of contiguous arcs (if possible). The exact algorithm meets the following conditions:

 * Observations are only simulated within the time span of ``start_time`` and ``end_time``
 * A contiguous tracking arc has simulated observations separated by ``interval_between_observations``
 * Starting from ``start_time``, an observation is simulated each ``interval_between_observations``. Once an observation is unviable, as defined by
   the ``arc_limiting_constraints`` input, it is checked whether the arc up until that point
   is longer in duration than ``minimum_arc_duration``. If it is, the arc is added to the simulated observations. If not, the arc is discarded. In either case, a new arc is started once a
   viable is observation is encountered
 * If the current arc reaches a duration greater than ``maximum_arc_duration``, the arc is added to the existing observations, and a new arc is started
 * If defined (e.g. if not NaN), the current observation time is incremented by ``minimum_time_between_arcs`` when an arc has been added to the observations.

 Nominally, this algorithm ensures that any arc of observations has a minimum and maximum duration. In addition, it ensures that (if desired) there is a minimum time interval
 between two tracking arcs. This behaviour can be modified by adding ``additional_viability_settings``, which are *not* used when computing the tracking arcs, but which are instead only used
 to reduce the set of simulated observations afterwards.


 Parameters
 ----------
 observable_type : :class:`ObservableType`
     Observable type of which observations are to be simulated.
 link_ends : LinkDefinition
     Link ends for which observations are to be simulated.
 start_time : float
     First time at which an observation is to be simulated (and checked for viability).
 end_time : float
     Maximum time at which an observation is to be simulated (and checked for viability).
 interval_between_observations : float
     Cadence (in seconds) of subsequent observations in an arc
 arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
     whether an arc should be terminated.

 minimum_arc_duration : float
     Minimum permissible time for a tracking arc
 maximum_arc_duration : float
     Maximum permissible time for a tracking arc
 minimum_time_between_arc : float, default = NaN
     Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
 additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
     These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ], default = None
     Function providing the observation noise factors as a function of observation time.
 Returns
 -------
 :class:`TabulatedObservationSimulationSettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` class.






     )doc" );

    m.def( "continuous_arc_simulation_settings_list",
           &tss::perArcObservationSimulationSettingsList< TIME_TYPE >,
           py::arg( "link_ends_per_observable" ),
           py::arg( "start_time" ),
           py::arg( "end_time" ),
           py::arg( "interval_between_observations" ),
           py::arg( "arc_limiting_constraints" ),
           py::arg( "minimum_arc_duration" ),
           py::arg( "maximum_arc_duration" ),
           py::arg( "minimum_time_between_arcs" ),
           py::arg( "reference_link_end_type" ) = tom::receiver,
           py::arg( "additional_viability_settings" ) = std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >( ),
           R"doc(

 Function for creating a list of settings object for observation simulation, using observation times according to a requirement for a continuous tracking arc.

 Function for creating multiple settings objects for observation simulation in a list. This function is
 equivalent to calling the :func:`~tudatpy.estimation.observations_setup.observations_simulation_settings.continuous_arc_simulation_settings` repeatedly, with the different
 observables and link definition provided here through `link_ends_per_observable`.
 During a single call to this function, one simulation settings object is created for each combination of observable type and link geometry given by the `link_ends_per_observable` parameter.


 Parameters
 ----------
 link_ends_per_observable : Dict[:class:`ObservableType`, List[LinkDefinition]]]
     Link geometry per observable type of which observations are to be simulated.
 start_time : float
     First time at which an observation is to be simulated (and checked for viability).
 end_time : float
     Maximum time at which an observation is to be simulated (and checked for viability).
 interval_between_observations : float
     Cadence (in seconds) of subsequent observations in an arc
 arc_limiting_constraints : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     List of settings for the creation of the viability criteria calculators, which are used to check if an observation is viable, and define
     whether an arc should be terminated.

 minimum_arc_duration : float
     Minimum permissible time for a tracking arc
 maximum_arc_duration : float
     Maximum permissible time for a tracking arc
 minimum_time_between_arc : float, default = NaN
     Minimum time between two tracking arcs. If NaN, this is effectively set to the ``interval_between_observations``
 additional_viability_settings : List[ :class:`ObservationViabilitySettings` ], default = [ ]
     Settings for the creation of the viability criteria calculators, which conduct viability checks on the simulated observations.
     These settings are *not* used to determine whether an arc is to be terminated, but are instead applied after the arcs have been computed.

 Returns
 -------
 List[ :class:`TabulatedObservationSimulationSettings` ]
     List of :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` derived :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.TabulatedObservationSimulationSettings` objects.






     )doc" );

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
           py::overload_cast<
               const std::vector< std::shared_ptr< tom::ObservationModelSettings > >&,
               const tss::SystemOfBodies& >(
               &tom::createObservationSimulators< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observation_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function for creating observation simulator objects.

 Function for creating observation simulator objects from observation settings.
 Note that each observation (i.e. combination of observable and link geometry) requires its own observation simulator object.


 Parameters
 ----------
 observation_settings : List[ ObservationSettings ]
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

}
}
}
}