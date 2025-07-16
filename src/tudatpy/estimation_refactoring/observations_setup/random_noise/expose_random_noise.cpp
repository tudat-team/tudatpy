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
#include "expose_random_noise.h"
#include <pybind11/functional.h>
#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/observationSimulationSettings.h"

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace tudat
{

namespace simulation_setup
{

void addNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction )
{
    tss::addNoiseFunctionToObservationSimulationSettings< TIME_TYPE, Eigen::VectorXd >( observationSimulationSettings,
                                                                                        observationNoiseFunction );
}

void addNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction,
        const tom::ObservableType observableType )
{
    tss::addNoiseFunctionToObservationSimulationSettings< TIME_TYPE, Eigen::VectorXd, const tom::ObservableType >(
            observationSimulationSettings, observationNoiseFunction, observableType );
}

void addNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const std::function< Eigen::VectorXd( const double ) > observationNoiseFunction,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addNoiseFunctionToObservationSimulationSettings< TIME_TYPE,
                                                          Eigen::VectorXd,
                                                          const tom::ObservableType,
                                                          const tom::LinkDefinition& >(
            observationSimulationSettings, observationNoiseFunction, observableType, linkEnds );
}

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const double observationNoiseAmplitude )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< TIME_TYPE >( observationSimulationSettings, observationNoiseAmplitude );
}

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const double observationNoiseAmplitude,
        const tom::ObservableType observableType )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType >(
            observationSimulationSettings, observationNoiseAmplitude, observableType );
}

void addGaussianNoiseFunctionToObservationSimulationSettingsPy(
        const std::vector< std::shared_ptr< ObservationSimulationSettings< TIME_TYPE > > >& observationSimulationSettings,
        const double observationNoiseAmplitude,
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds )
{
    tss::addGaussianNoiseFunctionToObservationSimulationSettings< TIME_TYPE, const tom::ObservableType, const tom::LinkDefinition& >(
            observationSimulationSettings, observationNoiseAmplitude, observableType, linkEnds );
}

} // namespace simulation_setup

} // namespace tudat


namespace tudatpy
{
namespace estimation_refactoring
{
namespace observations_setup
{

namespace random_noise
{

void expose_random_noise( py::module& m )
{
    m.def( "add_noise_function_to_all",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::function< Eigen::VectorXd( const double ) > >(
                   &tss::addNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           R"doc(

 Function for adding a custom noise function to all existing observation simulation settings.

 Function for including a custom noise function to the simulation settings of all observables.
 The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
 list.

 Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings_list : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
     Function providing the observation noise factors as a function of observation time.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_noise_function_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::function< Eigen::VectorXd( const double ) >,
                              const tom::ObservableType >( &tss::addNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for adding a custom noise function to selected existing observation simulation settings of a given observable type.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings_list : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
     Function providing the observation noise factors as a function of observation time.

 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_noise_function_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const std::function< Eigen::VectorXd( const double ) >,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >( &tss::addNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           R"doc(

 Function for adding a custom noise function to existing observation simulation settings of a given observable type and link definition.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_noise_function_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.

 noise_function : Callable[ [float], numpy.ndarray[numpy.float64[m, 1]] ]
     Function providing the observation noise factors as a function of observation time.

 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_gaussian_noise_to_all",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&, const double >(
                   &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           R"doc(

 Function for adding gaussian noise function to all existing observation simulation settings.

 Function for including simple time-independent and time-uncorrelated Gaussian noise function to the simulation settings of one or more observable(s).
 The noise settings are added to all :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) in the `observation_simulation_settings`
 list.

 Note: the :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects are modified in-place by this function,
 and thus the function does not return anything.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 noise_amplitude : float
     Standard deviation defining the un-biased Gaussian distribution for the noise.
 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_gaussian_noise_to_observable",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const double,
                              const tom::ObservableType >( &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           R"doc(

 Function for adding gaussian noise function to existing observation simulation settings of a given observable type.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 noise_amplitude : float
     Standard deviation defining the un-biased Gaussian distribution for the noise.
 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );

    m.def( "add_gaussian_noise_to_observable_for_link_ends",
           py::overload_cast< const std::vector< std::shared_ptr< tss::ObservationSimulationSettings< TIME_TYPE > > >&,
                              const double,
                              const tom::ObservableType,
                              const tom::LinkDefinition& >( &tss::addGaussianNoiseFunctionToObservationSimulationSettingsPy ),
           py::arg( "observation_simulation_settings_list" ),
           py::arg( "noise_amplitude" ),
           py::arg( "observable_type" ),
           py::arg( "link_definition" ),
           R"doc(

 Function for adding gaussian noise function to existing observation simulation settings of a given observable type and link definition.

 As :func:`~tudatpy.numerical_simulation.estimation_setup.observation.add_gaussian_noise_to_all`, except that the function only adds noise to entries of the
 `observation_simulation_settings` list that matches the specified `observable_type` and `link_definition`.


 Parameters
 ----------
 observation_simulation_settings : List[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` ]
     Observation simulation settings, given by a list of one or more existing :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects.
 noise_amplitude : float
     Standard deviation defining the un-biased Gaussian distribution for the noise.
 observable_type : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`
     Identifies the observable type in the observation simulation settings to which the noise is to be added.

 link_definition : :class:`~tudatpy.numerical_simulation.estimation_setup.observation.LinkDefinition`
     Identifies the link definition in the observation simulation settings for which the noise is to be added.

 Returns
 -------
 None
     The :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` object(s) are changed in-place.







     )doc" );
}

}
}
}
}