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
#include "expose_observations_dependent_variables.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

namespace tudat
{

namespace simulation_setup
{

void addDependentVariablesToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableList,
        const SystemOfBodies& bodies )
{
    tss::addDependentVariablesToObservationSimulationSettings< TIME_TYPE >( observationSimulationSettings, dependentVariableList, bodies );
}

void addDependentVariablesToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableList,
        const SystemOfBodies& bodies,
        const tom::ObservableType observableType )
{
    tss::addDependentVariablesToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType >(
            observationSimulationSettings, dependentVariableList, bodies, observableType );
}

void addDependentVariablesToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableList,
        const SystemOfBodies& bodies,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addDependentVariablesToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType, const tom::LinkDefinition& >(
            observationSimulationSettings, dependentVariableList, bodies, observableType, linkEnds );
}

} // namespace simulation_setup

} // namespace tudat


namespace tudatpy
{
namespace estimation
{
namespace observations_setup
{

namespace observations_dependent_variables
{

void expose_observations_dependent_variables( py::module& m )
{
    py::class_< tss::ObservationDependentVariableSettings, std::shared_ptr< tss::ObservationDependentVariableSettings > >(
            m,
            "ObservationDependentVariableSettings",
            R"doc(

         Base class for setting observation dependent variables as part of the observation output.

         Base class for setting observation dependent variables as part of the observation output.
         The user can create instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.observation.elevation_angle_dependent_variable` function.
         Note: The associated functionality is not yet mature enough for the end user. Class is exposed for development purposes only.

         Examples
         --------
         .. code-block:: python

             # Code snippet to show the creation of an ObservationDependentVariableSettings object
             from tudatpy.numerical_simulation.estimation_setup import observation

             # Create ObservationDependentVariableSettings object
             elevation_angle_settings = observation.elevation_angle_dependent_variable(observation.receiver)

             # Show that this is indeed an ObservationDependentVariableSettings object
             print(elevation_angle_settings)



      )doc" );

    m.def( "add_dependent_variables_to_all",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
                              const tss::SystemOfBodies& >( &tss::addDependentVariablesToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings" ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           R"doc(

 Function for including dependent variables into all existing observation simulation settings.

 Function for including the computation and reporting of dependent variables into the observation simulation settings of all observables.
 Note: The associated functionality is not yet mature enough for the end user. Function is exposed for development purposes only.

 Modifications are applied to all given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s),
 matching each :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object with the corresponding :class:`ObservationDependentVariableSettings` entry in the `dependent_variable_settings` parameter.
 Note that the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

 bodies : :class:`~tudatpy.dynamics.environment_setup.SystemOfBodies`
     Object consolidating all bodies and environment models that constitute the physical environment.






     )doc" );

    m.def( "add_dependent_variables_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
                              const tss::SystemOfBodies&,
                              const tom::ObservableType >( &tss::addDependentVariablesToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings" ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for including dependent variables into selected existing observation simulation settings.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_dependent_variables_to_all`, except that the function only adds includes the
 computation and reporting of dependent variables to entries of the `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 dependent_variable_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` ]
     List of one or more :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationDependentVariableSettings` objects, defining the dependent variables to be considered.

 bodies : :class:`~tudatpy.dynamics.environment_setup.SystemOfBodies`
     Object consolidating all bodies and environment models that constitute the physical environment.

 observable_type : :class:`ObservableType`
     Identifies the observable type in the observation simulation settings for which the dependent variables are to be included.






     )doc" );

    m.def( "add_dependent_variables_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::vector< std::shared_ptr< tss::ObservationDependentVariableSettings > >&,
                              const tss::SystemOfBodies&,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >( &tss::addDependentVariablesToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings" ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "bodies" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(No documentation found.)doc" );

    /////////////////////////////////////////////////////////////////////////////////////////////////
    // DEPENDENT VARIABLES
    /////////////////////////////////////////////////////////////////////////////////////////////////

    py::enum_< tss::IntegratedObservationPropertyHandling >(
            m, "IntegratedObservationPropertyHandling", R"doc(No documentation found.)doc" )
            .value( "interval_start", tss::IntegratedObservationPropertyHandling::interval_start )
            .value( "interval_end", tss::IntegratedObservationPropertyHandling::interval_end )
            .value( "interval_undefined", tss::IntegratedObservationPropertyHandling::interval_undefined )
            .export_values( );

    m.def( "elevation_angle_dependent_variable",
           &tss::elevationAngleDependentVariable,
           py::arg( "link_end_type" ) = tom::unidentified_link_end,
           py::arg( "link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "originating_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "originating_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "azimuth_angle_dependent_variable",
           &tss::azimuthAngleDependentVariable,
           py::arg( "link_end_type" ) = tom::unidentified_link_end,
           py::arg( "link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "originating_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "originating_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "target_range_between_link_ends_dependent_variable",
           &tss::targetRangeBetweenLinkEndsDependentVariable,
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "avoidance_angle_dependent_variable",
           &tss::bodyAvoidanceAngleDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "body_center_distance_dependent_variable",
           &tss::linkBodyCenterDistanceDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "body_limb_distance_dependent_variable",
           &tss::linkLimbDistanceDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "angle_wrt_orbital_plane_dependent_variable",
           &tss::linkAngleWrtOrbitalPlaneDependentVariable,
           py::arg( "body_name" ),
           py::arg( "start_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "end_link_end_type" ) = tom::unidentified_link_end,
           py::arg( "start_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "end_link_end_id" ) = tom::LinkEndId( "", "" ),
           py::arg( "integrated_observation_handling" ) = tss::interval_start,
           R"doc(No documentation found.)doc" );

    m.def( "integration_time_dependent_variable",
           &tss::integrationTimeDependentVariable,
           py::arg( "observable_type" ) = tom::undefined_observation_model,
           R"doc(No documentation found.)doc" );

    m.def( "retransmission_delays_dependent_variable",
           &tss::retransmissionDelaysDependentVariable,
           py::arg( "observable_type" ) = tom::undefined_observation_model,
           R"doc(No documentation found.)doc" );
}

}
}
}
}