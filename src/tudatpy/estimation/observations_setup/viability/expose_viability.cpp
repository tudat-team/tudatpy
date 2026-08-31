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
#include "expose_viability.h"

#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/simulation/estimation_setup/createObservationModelSettings.h"
#include "tudat/simulation/estimation_setup/createObservationViability.h"
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

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

namespace viability
{

void expose_observation_viability_settings_type( py::module& m )
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
}

void expose_viability( py::module& m )
{
    py::class_< tom::ObservationBoundariesViabilitySettings,
                std::shared_ptr< tom::ObservationBoundariesViabilitySettings >,
                tom::ObservationViabilitySettings >( m,
                                                     "ObservationBoundariesViabilitySettings",
                                                     R"doc(

         Class for defining observation boundaries viability settings.

         Class for defining the settings for observation boundaries viability calculator creation.
         Instances of this class are typically created through the :func:`~tudatpy.estimation.observations_setup.viability.observation_boundaries_viability` function.
        )doc" );

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
            .value( "observation_boundaries", tom::ObservationViabilityType::observation_boundaries )
            .value( "ground_station_darkness", tom::ObservationViabilityType::ground_station_darkness )
            .value( "body_in_sunlight", tom::ObservationViabilityType::body_in_sunlight )
            .export_values( );

    m.def( "observation_boundaries_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const std::vector< std::pair< double, double > > >(
                   &tom::observationBoundariesViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "boundaries" ),
           R"doc(

    Function for defining observation boundaries viability settings.

    Function for defining observation boundaries viability settings for single link ends.
    When simulating observations, this setting ensures that any applicable observations, for which the observed value is outside of given boundaries, will be omitted.

    Examples
    --------
    .. code-block:: python

        # Code snippet to show the creation of an ObservationBoundariesViabilitySettings object
        import numpy as np
        from tudatpy.estimation.observations_setup import viability

        # Create ObservationBoundariesViabilitySettings object
        # In this case, we exclude observations for which the observed value is outside of the following boundaries: [0, 100] for the first entry of the observation vector, and [-50, 50] for the second entry of the observation vector.
        boundaries = [(0, 100), (-50, 50)]
        viability_settings = viability.observation_boundaries_viability(["Earth", ""], boundaries)

        # Show that this is indeed an ObservationBoundariesViabilitySettings object
        print(viability_settings)
    
    Parameters
    ----------
    link_end_id : tuple[str,str]
    Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId` ), for which the viability settings are to be created.
    
    boundaries : list[tuple[float, float]]
    List of pairs of minimum and maximum allowed values for the observation. Each entry on the list corresponds to minimum and maximum allowed for each entry in the observation vector.

    Returns
    -------
    ObservationBoundariesViabilitySettings
        Observation-boundary viability settings for the link end.

     )doc" );

    m.def( "elevation_angle_viability",
           py::overload_cast< const std::pair< std::string, std::string >, const double >( &tom::elevationAngleViabilitySettings ),
           py::arg( "link_end_id" ),
           py::arg( "elevation_angle" ),
           R"doc(

 Function for defining single elevation angle viability setting.

 Function for defining elevation angle viability settings for single link end.
 When simulating observations, this setting ensures that any applicable observations for which the local elevation angle at link end
 ``link_end_id`` is less than some limit value, will be omitted. Note that for (for instance) a two-way observable where the given ``link_end_id``
 acts as both receiver and transmitter, the check will be performed both at reception and transmission time, and the observation will only
 be accepted if the elevation angle is sufficient at both epochs.

 The elevation angle used by this functionality can also be computed manually using this :meth:`~tudatpy.dynamics.environment.PointingAnglesCalculator.calculate_elevation_angle`
 function of the :class:~tudatpy.dynamics.environment.PointingAnglesCalculator` class, which can be extracted from a
 :class:`~tudatpy.dynamics.environment.GroundStation` object using :func:`~tudatpy.dynamics.environment.GroundStation.pointing_angles_calculator`


 Parameters
 ----------
 link_end_id : tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

 elevation_angle : float
     Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 ObservationViabilitySettings
     Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` class, defining the settings for observation viability

     )doc" );

    m.def( "ground_station_darkness_viability",
           &tom::groundStationDarknessViabilitySettings,
           py::arg( "link_end_id" ),
           py::arg( "maximum_sun_elevation" ) = -12.0 * tudat::mathematical_constants::PI / 180.0,
           R"doc(

 Function for defining a ground-station darkness viability setting.

 Observations are retained only when the Sun elevation at every epoch at which the specified ground station participates
 in the observation is at or below ``maximum_sun_elevation``. The Sun state is obtained from the environment in the global
 frame, independently of the observation link-end states.

 Parameters
 ----------
 link_end_id : tuple[str,str]
     Ground-station link end for which darkness is required.
 maximum_sun_elevation : float, default=-0.20943951023931953
     Maximum allowed Sun elevation in radians. The default is -12 degrees.

 Returns
 -------
 ObservationViabilitySettings
     Settings defining the ground-station darkness viability condition.

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

 The epoch at which the avoided body is evaluated is computed in an identical way as the epoch at which the occulted body is evaluated
 (see :func:`~body_occultation_viability`)

 Parameters
 ----------
 link_end_id : tuple[str,str]
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
 ObservationViabilitySettings
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
 The occultation is computed using the shape model of the specified body, using a spherical body approximation
 (using :attr:`~tudatpy.dynamics.environment.BodyShapeModel.average_radius` attribute of the shape model).

 The epoch :math:`t_{\text{occ}}` a which the ephemeris :math:`\mathbf{r}_{\text{occ}}` of the ``occulting_body`` is evaluated to compute the occultation is defined as follows:

 .. math::
   {r}_{\text{occ},1}&=||r_{\text{occ}}(t_{1})-\mathbf{r}_{1}(t_{1})||\\
   {r}_{\text{occ},2}&=||r_{\text{occ}}(t_{2})-\mathbf{r}_{2}(t_{2})||\\
   t_{\text{occ}}&=t_{1}\frac{{r}_{\text{occ},2}}{{r}_{\text{occ},1}+{r}_{\text{occ},2}}+t_{2}\frac{{r}_{\text{occ},1}}{{r}_{\text{occ},1}+{r}_{\text{occ},2}}

 for a link between link ends 1 and 2 (and transmission/reception epochs :math:`t_{1}` and :math:`t_{2}`)

 Parameters
 ----------
 link_end_id : tuple[str,str]
     Link end (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.

 body_to_avoid : str
     Name of the body which the signal path should not be occulted by.

 Returns
 -------
 ObservationViabilitySettings
     Instance of the :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`, defining the settings for observation viability.






     )doc" );

    m.def( "body_in_sunlight_viability",
           &tom::bodyInSunlightViabilitySettings,
           py::arg( "link_end_id" ),
           py::arg( "occulting_bodies" ),
           R"doc(

 Function for defining a body-in-sunlight viability setting.

 Observations are retained only when the Sun is not occulted by any of ``occulting_bodies`` at every epoch at which the
 specified link end participates in the observation. Sun and occulting-body states are obtained from the environment in
 the global frame, independently of the other observation link-end states.

 Parameters
 ----------
 link_end_id : tuple[str,str]
     Link end representing the body that must be illuminated by the Sun.
 occulting_bodies : list[str]
     Bodies that may occult the Sun, using their spherical average radii.

 Returns
 -------
 ObservationViabilitySettings
     Settings defining the body-in-sunlight viability condition.

     )doc" );

    m.def( "observation_boundaries_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >,
                              const std::vector< std::pair< double, double > > >( &tom::observationBoundariesViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "boundaries" ),
           R"doc(

 Function for defining list of observation boundaries viability settings, equivalent to a series of calls to :func:`~observation_boundaries_viability`.

 Parameters
 ----------
 link_end_ids : List[ tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the observation boundaries viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

    boundaries : list[tuple[float, float]]
    List of pairs of minimum and maximum allowed values for the observation. Each entry on the list corresponds to minimum and maximum allowed for each entry in the observation vector.

 Returns
 -------
 list[ObservationBoundariesViabilitySettings]
     List of observation-boundary viability settings, one for each link end.
    )doc" );

    m.def( "elevation_angle_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >, const double >(
                   &tom::elevationAngleViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "elevation_angle" ),
           R"doc(

 Function for defining list of elevation angle viability settings, equivalent to a series of calls to :func:`~elevation_angle_viability`.

 Parameters
 ----------
 link_end_ids : List[ tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end violates the minimum elevation angle constraint.

 elevation_angle : float
     Limit elevation angle, below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 list[ObservationViabilitySettings]
     List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "body_avoidance_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >, const std::string, const double >(
                   &tom::bodyAvoidanceAngleViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "body_to_avoid" ),
           py::arg( "avoidance_angle" ),
           R"doc(

 Function for defining list of body avoidance viability settings, equivalent to a series of calls to :func:`~body_avoidance_viability`.


 Parameters
 ----------
 link_end_ids : List[ tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the elevation angle viability setting is to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""].

 body_to_avoid : str
     Name of the body which the signal path should not pass 'too close' to.

 avoidance_angle : float
     Limit angle (generalization of SPE angle), below which no observations are produced when using the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. Note: this
     value must be in radians.

 Returns
 -------
 list[ObservationViabilitySettings]
     List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "body_occultation_viability_list",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >, const std::string >(
                   &tom::bodyOccultationViabilitySettings ),
           py::arg( "link_end_ids" ),
           py::arg( "occulting_body" ),
           R"doc(

 Function for defining body occultation viability settingsFunction for defining list of body avoidance viability settings, equivalent to a series of calls to :func:`~body_occultation_viability`.


 Parameters
 ----------
 link_end_ids : List[ tuple[str,str] ]
     List of individual link ends (as defined by body/reference point pair, see :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`), for which the viability settings are to be created.
     To apply these settings to *all* ground station on a given body (such as "Earth"), use ["Earth", ""] is entry in this list.
     For each link end included in this list, it will be checked if a signal received by and/or transmitted (or reflected) by this
     link end is occulted by the specified body.

 body_to_avoid : str
     Name of the body which the signal path should not be occulted by.

 Returns
 -------
 list[ObservationViabilitySettings]
     List of :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, each defining the settings for observation viability of one link end.






     )doc" );

    m.def( "add_viability_check_to_all",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tom::ObservationViabilitySettings > >& >(
                   &tss::addViabilityToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "viability_settings" ),
           R"doc(

 Function for including viability checks into existing observation simulation settings.

 Function for adding viability checks to the observation simulation settings, such that only observations meeting certain conditions are retained.
 The viability settings are added to all :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) in the ``observation_simulation_settings_list``
 list.
 Note: the :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings_list : list[:class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings`]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
 viability_settings : list[:class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`]
     List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.

 Returns
 -------
 None
     The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place.







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

 As :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, except that the function only adds viability settings to entries of the
 ``observation_simulation_settings_list`` list that matches the specified `observable_type`.

 Parameters
 ----------
 observation_simulation_settings_list : list[:class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings`]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
 viability_settings : list[:class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`]
     List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.
 observable_type : :class:`tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
     Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.

 Returns
 -------
 None
     The :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` object(s) are changed in-place.







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

 As :func:`~tudatpy.estimation.observations_setup.viability.add_viability_check_to_all`, except that the function only adds viability settings to entries of the
 ``observation_simulation_settings_list`` list that matches the specified `observable_type` and `link_ends`.


 Parameters
 ----------
 observation_simulation_settings_list : list[:class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings`]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.estimation.observations_setup.observations_simulation_settings.ObservationSimulationSettings` objects.
 viability_settings : list[:class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings`]
     List of one or more :class:`~tudatpy.estimation.observations_setup.viability.ObservationViabilitySettings` objects, defining the viability checks to be included.
 observable_type : :class:`tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
     Identifies the observable type in the observation simulation settings for which the viability checks are to be considered.
 link_ends : :class:`~tudatpy.estimation.observable_models_setup.links.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the viability checks are to be considered.

 Returns
 -------
 None







     )doc" );
}

}  // namespace viability
}  // namespace observations_setup
}  // namespace estimation
}  // namespace tudatpy
