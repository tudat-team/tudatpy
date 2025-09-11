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
#include "expose_viability.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace simulation_setup
{

    void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList )
{
    tss::addViabilityToObservationSimulationSettings< TIME_TYPE >( observationSimulationSettings, viabilitySettingsList );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList,
        const tom::ObservableType observableType )
{
    tss::addViabilityToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType >(
            observationSimulationSettings, viabilitySettingsList, observableType );
}

void addViabilityToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >& viabilitySettingsList,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addViabilityToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType, const tom::LinkDefinition& >(
            observationSimulationSettings, viabilitySettingsList, observableType, linkEnds );
}

} // namespace simulation_setup

} // namespace tudat


namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

namespace viability
{

void expose_viability( py::module& m )
{
    py::class_< tom::ObservationViabilitySettings, std::shared_ptr< tom::ObservationViabilitySettings > >( m,
                                                                                                           "ObservationViabilitySettings",
                                                                                                           R"doc(

         Class for defining observation viability calculator settings.

         Class for defining the settings for observation viability calculator creation.
         Instances of this class are typically be created through various dedicated functions,such as :func:`~tudatpy.estimation.observations_setup.viability.elevation_angle_viability`, :func:`~tudatpy.estimation.observations_setup.viability.body_avoidance_viability` and :func:`~tudatpy.estimation.observations_setup.viability.body_occultation_viability`

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationViabilitySettings object
             import numpy as np
             from tudatpy.estimation.observations_setup import viability

             # Create ObservationViabilitySettings object
             # In this case, we exclude observations for which the local elevation angle at link end is less 15 degrees.
             min_elevation = np.deg2rad(15)
             # We apply these settings to every ground station on Earth using the following link_end_id: [“Earth”, “”]
             viability_settings = viability.elevation_angle_viability(["Earth", ""], min_elevation)

             # Show that this is indeed an ObservationViabilitySettings object
             print(viability_settings)




      )doc" );

    py::class_< tom::CustomObservationViabilitySettings,
                std::shared_ptr< tom::CustomObservationViabilitySettings >,
                tom::ObservationViabilitySettings >( m, "CustomObservationViabilitySettings", R"doc()doc" );

    py::enum_< tom::ObservationViabilityType >( m, "ObservationViabilityType", R"doc(

Enumeration of observation viability criterion types.

Examples
--------
.. code-block:: python

    # Code snippet to print all available Observation Viability Types
    from tudatpy.estimation.observations_setup import viability

    num_observation_viability_types = len(viability.ObservationViabilityType.__members__)
    print(f'The length of all available Tudatpy Observation Viability Types is: {num_observation_viability_types}')

    # Print all available Observation Viability Types using the "name" property
    for i in range(num_observation_viability_types):
        print(i, viability.ObservationViabilityType(i).name)




      )doc" )
            .value( "minimum_elevation_angle", tom::ObservationViabilityType::minimum_elevation_angle )
            .value( "body_avoidance_angle", tom::ObservationViabilityType::body_avoidance_angle )
            .value( "body_occultation", tom::ObservationViabilityType::body_occultation )
            .value( "custom_viability", tom::ObservationViabilityType::custom_viability )
            .export_values( );


    m.def( "elevation_angle_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const double >( &tom::elevationAngleViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "elevation_angle" ),
           R"doc(

 Function for defining single elevation angle viability setting.

 Function for defining elevation angle viability settings for single link end.
 When simulating observations, this setting ensures that any applicable observations, for which the local elevation angle at link end is less than some limit value, will be omitted.


 Parameters
 ----------
 link_end_id : Tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

 elevation_angle : float
     Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` class, defining the settings for observation viability






     )doc" );

    m.def( "body_avoidance_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const std::string, const double >(
                   &tom::bodyAvoidanceAngleViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "body_to_avoid" ),
           py::arg( "avoidance_angle" ),
           R"doc(

 Function for defining body avoidance observation viability settings.

 Function for defining body avoidance observation viability settings for single link ends.
 When simulating observations, this settings ensures that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
 The definition of 'too close' is computed as the angle between:

 * The line-of-sight vector from a link end to a given third body
 * The line-of-sight between two link ends

 This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
 a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


 Parameters
 ----------
 link_end_id : Tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` ), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end passes too close to the specified body.

 body_to_avoid : str
     Name of the body which the signal path should not pass 'too close' to.

 avoidance_angle : float
     Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`, defining the settings for observation viability.






     )doc" );

    m.def( "body_occultation_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const std::string >( &tom::bodyOccultationViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "occulting_body" ),
           R"doc(

 Function for defining body occultation viability settings.

 Function for defining body occultation viability settings for single link ends.
 When simulating observations, this setting ensures that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
 The occultation is computed using the shape model of the specified body.


 Parameters
 ----------
 link_end_id : Tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.

 body_to_avoid : str
     Name of the body which the signal path should not be occulted by.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`, defining the settings for observation viability.






     )doc" );

    m.def( "custom_observation_viability",
           py::overload_cast< const std::pair< std::string, std::string >,
                              const std::function< bool( const std::vector< Eigen::Vector6d >, const std::vector< double > ) > >(
                   &tom::customObservationViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "custom_viability_function" ),
           R"doc(No documentation found.)doc" );

    m.def( "elevation_angle_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >, const double >(
                   &tom::elevationAngleViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "elevation_angle" ),
           R"doc(

 Function for defining list of elevation angle viability settings.

 Function for defining elevation angle viability settings for multiple link ends.
 Each entry in the returned list contains the observation viability settings for one link end.
 When simulating observations, these settings ensure that any applicable observations, for which the local elevation angle at a link end is less than some limit value, will be omitted.


 Parameters
 ----------
 link_end_ids : List[ Tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end violates the minimum elevation angle constraint.

 elevation_angle : float
     Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "body_avoidance_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >, const std::string, const double >(
                   &tom::bodyAvoidanceAngleViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "body_to_avoid" ),
           py::arg( "avoidance_angle" ),
           R"doc(

 Function for defining list of body avoidance viability settings.

 Function for defining body avoidance viability settings for multiple link ends.
 Each entry in the returned list contains the observation viability settings for one link end.
 When simulating observations, these settings ensure that any applicable observations, for which the signal path passes 'too close' to a body, will be omitted.
 The definition of 'too close' is computed as the angle between:

 * The line-of-sight vector from a link end to a given third body
 * The line-of-sight between two link ends

 This constraint is typically used to prevent the Sun from being too close to the field-of-view of the telescope(s), as defined by
 a so-called 'SPE' (Sun-Probe-Earth) angle constraint. The present viability setting generalizes this constraint.


 Parameters
 ----------
 link_end_ids : List[ Tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

 body_to_avoid : str
     Name of the body which the signal path should not pass 'too close' to.

 avoidance_angle : float
     Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "body_occultation_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >, const std::string >(
                   &tom::bodyOccultationViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "occulting_body" ),
           R"doc(

 Function for defining body occultation viability settings.

 Function for defining body occultation viability settings for multiple link ends.
 Each entry in the returned list contains the observation viability settings for one link end.
 When simulating observations, these settings ensure that any applicable observations, for which the signal path is occulted by a given body, will be omitted.
 The occultation is computed using the shape model of the specified body.


 Parameters
 ----------
 link_end_ids : List[ Tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end is occulted by the specified body.

 body_to_avoid : str
     Name of the body which the signal path should not be occulted by.

 Returns
 -------
 :class:`ObservationViabilitySettings`
     List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "custom_observation_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >,
                              const std::function< bool( const std::vector< Eigen::Vector6d >, const std::vector< double > ) > >(
                   &tom::customObservationViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "custom_viability_function" ),
           R"doc(No documentation found.)doc" );

    m.def( "add_viability_check_to_all",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >& >(
                   &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
 The noise settings are added to all :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
 list.
 Note: the :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects.
 viability_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` ]
     List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.

 Returns
 -------
 None
     The :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_viability_check_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
                              const tom::ObservableType >( &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 As :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, except that the function only adds viabilitt settings to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects.
 viability_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` ]
     List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.

 observable_type : :class:`ObservableType`
     Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

 Returns
 -------
 None
     The :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_viability_check_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >&,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >( &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 As :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.viability.ObservationSimulationSettings` objects.
 viability_settings : List[ :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` ]
     List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.

 observable_type : :class:`ObservableType`
     Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

 link_definition : :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.

 Returns
 -------
 None







     )doc" );
}

}
}
}
}