/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#define PYBIND11_DETAILED_ERROR_MESSAGES
#include "expose_observations.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"

namespace py = pybind11;
namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace te = tudat::ephemerides;

namespace tudat
{

namespace simulation_setup
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSetWithoutDependentVariables(
        const tom::ObservableType observableType,
        const tom::LinkDefinition& linkEnds,
        const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
        const std::vector< TimeType > observationTimes,
        const tom::LinkEndType referenceLinkEnd,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr )
{
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >( observableType,
                                                                                             linkEnds,
                                                                                             observations,
                                                                                             observationTimes,
                                                                                             referenceLinkEnd,
                                                                                             std::vector< Eigen::VectorXd >( ),
                                                                                             nullptr,
                                                                                             ancilliarySettings );
}

}  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy
{
namespace estimation
{
namespace observations
{

void expose_observations( py::module& m )
{
    auto observations_processing = m.def_submodule( "observations_processing" );
    observations_processing::expose_observations_processing( observations_processing );

    auto observations_geometry = m.def_submodule( "observations_geometry" );
    observations_geometry::expose_observations_geometry( observations_geometry );

    // SINGLE OBSERVATION SET

    py::class_< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                "SingleObservationSet",
                                                                                                R"doc(
        Class for storing and handling a single set of observations.

        This class stores a set of observations of a single observable type and a single link definition.
        The observations are stored as a vector of Eigen matrices, where each entry represents a single
        observation at a given time. The class also stores the observation times, reference link end,
        and other metadata.

        The pybind ``SingleObservationSet`` object is a list-like object, which can be iterated over, and from which specific
        observation data can be retrieved by index.

        When iterating, a tuple is returned with the following entries:
        - ``int``: index of the observation.
        - ``float``: time of the observation.
        - ``numpy.ndarray``: value of the observation.

        When using the ``[]`` operator, a tuple is returned with the following entries:
        - ``float``: time of the observation.
        - ``numpy.ndarray``: value of the observation.

        Parameters
        ----------
        observable_type : tudatpy.astro.observation_models.ObservableType
            Type of observable.
        link_ends : tudatpy.astro.observation_models.LinkDefinition
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_epochs : list[float]
            List of observation times.
        reference_link_end : tudatpy.astro.observation_models.LinkEndType
            Reference link end for the observation.
        observation_dependent_variables : list[numpy.ndarray], optional
            List of dependent variables for each observation.
        dependent_variable_calculator : tudatpy.estimation.observations.ObservationDependentVariableCalculator, optional
            Calculator for dependent variables.
        ancilliary_settings : tudatpy.astro.observation_models.ObservationAncilliarySimulationSettings, optional
            Ancillary settings for the observation.
      )doc" )
            .def( py::init< const tom::ObservableType,
                            const tom::LinkDefinition,
                            const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >,
                            const std::vector< TIME_TYPE >,
                            const tom::LinkEndType,
                            const std::vector< Eigen::VectorXd >,
                            const std::shared_ptr< tss::ObservationDependentVariableBookkeeping >,
                            const std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >( ),
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "observations" ),
                  py::arg( "observation_epochs" ),
                  py::arg( "reference_link_end" ),
                  py::arg( "observation_dependent_variables" ) = std::vector< Eigen::VectorXd >( ),
                  py::arg( "dependent_variable_bookkeeping" ) = nullptr,
                  py::arg( "ancilliary_settings" ) = nullptr )
            .def( "set_observations",
                  py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(
Sets the observation values from a list of vectors.

Parameters
----------
observations : list[numpy.ndarray]
    The new list of observations.
)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(
Sets the observation values from a single concatenated vector.

Parameters
----------
observations : numpy.ndarray
    A single vector containing all observations concatenated.
)doc" )
            .def( "set_residuals",
                  py::overload_cast< const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(
Sets the residuals for all observations from a list of vectors.

Parameters
----------
residuals : list[numpy.ndarray]
    The new list of residuals.
)doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(
Sets the residuals for all observations from a single concatenated vector.

Parameters
----------
residuals : numpy.ndarray
    A single vector containing all residuals concatenated.
)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const double >( &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  R"doc(
Sets a constant scalar weight for all observations.

Parameters
----------
weight : float
    The constant weight to apply.
)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  R"doc(
Sets a constant weight vector for all observations.

Parameters
----------
weight : numpy.ndarray
    The constant weight vector to apply to each observation.
)doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< const Eigen::VectorXd& >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "weights" ),
                  R"doc(
Sets weights for all observations from a single concatenated vector.

Parameters
----------
weights : numpy.ndarray
    A single vector containing all weights concatenated.
)doc" )
            .def( "filter_observations",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations,
                  py::arg( "filter" ),
                  py::arg( "save_filtered_obs" ) = true,
                  R"doc(
Filters observations based on a given filter criterion.

Parameters
----------
filter : tudatpy.numerical_simulation.estimation.ObservationFilterBase
    The filter to apply.
save_filtered_obs : bool, optional
    If true, the filtered observations are stored in a separate set. Defaults to true.
)doc" )
            .def_property_readonly( "observable_type",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableType,
                                    R"doc(

         **read-only**

         Type of observable for which the object stores observations

         :type: ObservableType
      )doc" )
            .def_property( "link_definition",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEnds,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setLinkEnds,
                           R"doc(

         Definition of the link ends for which the object stores observations

         :type: LinkDefinition
      )doc" )
            .def_property_readonly( "reference_link_end",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getReferenceLinkEnd,
                                    R"doc(

         **read-only**

         Reference link end for all stored observations

         :type: LinkEndType
      )doc" )
            .def_property_readonly( "number_of_observables",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfObservables,
                                    R"doc(
Returns the number of observations in the set.

Returns
-------
int
    The number of observations.
)doc" )
            .def_property_readonly( "single_observable_size",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservableSize,
                                    R"doc(
Returns the size of a single observation value (e.g., 1 for range, 2 for angular position).

Returns
-------
int
    The size of a single observation.
)doc" )
            .def_property_readonly( "total_observation_set_size",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalObservationSetSize,
                                    R"doc(
Returns the total size of all observation values in the set.

Returns
-------
int
    The total size of the observation set.
)doc" )
            .def_property_readonly( "time_bounds",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                    R"doc(
Returns the start and end time of the observation set.

Returns
-------
tuple[tudatpy.astro.time_representation.Time, tudatpy.astro.time_representation.Time]
    The start and end time of the observations.
)doc" )
            .def_property_readonly( "list_of_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                                    R"doc(

         **read-only**

         List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1

         :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
            .def_property_readonly( "observation_times",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                                    R"doc(

         **read-only**

         List of reference times for each of the observations in ``list_of_observations``

         :type: list[ tudatpy.astro.time_representation.Time ]
      )doc" )
            .def_property_readonly( "concatenated_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsVector,
                                    R"doc(

         **read-only**

         Concatenated vector of all stored observations

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "computed_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                                    R"doc(
Returns the vector of computed observation values (observed - residual).

Returns
-------
list[numpy.ndarray]
    The list of computed observations.
)doc" )
            .def_property_readonly( "concatenated_computed_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservationsVector,
                                    R"doc(
Returns all computed observations concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all computed observations.
)doc" )
            .def_property_readonly( "residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                                    R"doc(
Returns the vector of residuals for all observations.

Returns
-------
list[numpy.ndarray]
    The list of residuals.
)doc" )
            .def_property_readonly( "concatenated_residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getResidualsVector,
                                    R"doc(
Returns all residuals concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all residuals.
)doc" )
            .def_property_readonly( "rms_residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                                    R"doc(
Returns the Root Mean Square (RMS) of the residuals.

Returns
-------
numpy.ndarray
    A vector containing the RMS of residuals for each component of the observable.
)doc" )
            .def_property_readonly( "mean_residuals",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                                    R"doc(
Returns the mean of the residuals.

Returns
-------
numpy.ndarray
    A vector containing the mean of residuals for each component of the observable.
)doc" )
            .def_property_readonly( "weights", &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights, R"doc(
Returns the vector of weights for all observations.

Returns
-------
list[numpy.ndarray]
    The list of weights.
)doc" )
            .def_property_readonly( "concatenad_weights",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                                    R"doc(
Returns all weights concatenated into a single vector.

Returns
-------
numpy.ndarray
    A single vector containing all weights.
)doc" )
            .def_property( "dependent_variables",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariables,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setObservationsDependentVariables,
                           R"doc(List of dependent variables for all observations.)doc" )
            .def_property_readonly( "dependent_variables_history",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                                    R"doc(
Time history of all observation dependent variables.

Returns
-------
dict[tudatpy.astro.time_representation.Time, numpy.ndarray]
    A map from observation time to the vector of dependent variables.
)doc" )
            .def_property_readonly( "observations_history",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsHistory,
                                    R"doc(

         **read-only**

         Dictionary of observations sorted by time. Created by making a dictionary with ``observation_times`` as keys and ``list_of_observations`` as values.

         :type: dict[ tudatpy.astro.time_representation.Time, numpy.ndarray[numpy.float64[m, 1]] ]
      )doc" )
            .def_property_readonly( "ancilliary_settings",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAncilliarySettings,
                                    R"doc(

         **read-only**

         Ancilliary settings for all stored observations

         :type: ObservationAncilliarySimulationSettings
      )doc" )
            .def_property( "weights_vector",
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getWeightsVector,
                           &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights,
                           R"doc(Concatenated vector of weights for all observations.)doc" )
            .def_property_readonly( "filtered_observation_set",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getFilteredObservationSet,
                                    R"doc(
Returns the set of filtered observations.

Returns
-------
SingleObservationSet
    The observation set containing filtered observations.
)doc" )
            .def_property_readonly( "number_filtered_observations",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getNumberOfFilteredObservations,
                                    R"doc(
Returns the number of observations that have been filtered out.

Returns
-------
int
    The number of filtered observations.
)doc" )
            .def( "single_dependent_variable",
                  py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >, const bool >(
                          &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariable ),
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "return_first_compatible_settings" ) = false,
                  R"doc(
Returns the values of a single dependent variable (specified by dependent variable settings).

Parameters
----------
dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
    Settings for the dependent variable to retrieve.
return_first_compatible_settings : bool, optional
    If true, returns the first compatible variable found. Defaults to false.

Returns
-------
numpy.ndarray
    A matrix of the dependent variable values over all observation times.
)doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  R"doc(
Returns the list of all dependent variable settings compatible with the settings provided as inputs.

Parameters
----------
dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
    The settings to check for compatibility.

Returns
-------
list[tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings]
    A list of compatible dependent variable settings.
)doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  R"doc(
Returns a vector containing the values of all dependent variables compatible with the settings provided as input.

Parameters
----------
dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
    The settings to check for compatibility.

Returns
-------
list[numpy.ndarray]
    A list of matrices, each containing values of a compatible dependent variable.
)doc" )
            .def( "single_dependent_variable_history",
                  &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleDependentVariableHistory,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "return_first_compatible_settings" ) = false,
                  R"doc(
Returns the time history of a single dependent variable (specified by settings).

Parameters
----------
dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
    Settings for the dependent variable to retrieve.
return_first_compatible_settings : bool, optional
    If true, returns the first compatible variable found. Defaults to false.

Returns
-------
dict[float, numpy.ndarray]
    A map from observation time to the value of the specified dependent variable.
)doc" )
            .def_property_readonly( "dependent_variables_matrix",
                                    &tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsDependentVariablesMatrix,
                                    R"doc(
Returns the dependent variables for all observations as a matrix.

Returns
-------
numpy.ndarray
    A matrix where each row corresponds to an observation and columns to dependent variables.
)doc" );

    m.def( "single_observation_set",
           &tss::singleObservationSetWithoutDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observable_type" ),
           py::arg( "link_definition" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancilliary_settings" ) = nullptr,
           R"doc(
        Creates a `SingleObservationSet` object without dependent variables.

        This function is a factory function that simplifies the creation of a `SingleObservationSet`
        when no dependent variables are associated with the observations.

        Parameters
        ----------
        observable_type : tudatpy.astro.observation_models.ObservableType
            Type of observable.
        link_definition : tudatpy.astro.observation_models.LinkDefinition
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_times : list[float]
            List of observation times.
        reference_link_end : tudatpy.astro.observation_models.LinkEndType
            Reference link end for the observation.
        ancilliary_settings : tudatpy.astro.observation_models.ObservationAncilliarySimulationSettings, optional
            Ancillary settings for the observation.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A `SingleObservationSet` object.
        )doc" );

    m.def( "create_single_observation_set",
           py::overload_cast< const tom::ObservableType,
                              const tom::LinkEnds&,
                              const std::vector< Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >&,
                              const std::vector< TIME_TYPE >,
                              const tom::LinkEndType,
                              const std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
                   &tom::createSingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings" ),
           R"doc(
        Factory function to create a `SingleObservationSet` object.

        This function creates a `SingleObservationSet` object from a list of observations and their corresponding times.
        This is a simplified factory function that does not include dependent variables.

        Parameters
        ----------
        observable_type : tudatpy.astro.observation_models.ObservableType
            Type of observable.
        link_ends : tudatpy.astro.observation_models.LinkEnds
            Definition of the link ends for the observation.
        observations : list[numpy.ndarray]
            List of observations. Each entry is a vector representing a single observation.
        observation_times : list[float]
            List of observation times.
        reference_link_end : tudatpy.astro.observation_models.LinkEndType
            Reference link end for the observation.
        ancillary_settings : tudatpy.astro.observation_models.ObservationAncilliarySimulationSettings, optional
            Ancillary settings for the observation.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A `SingleObservationSet` object.
        )doc" );

    // OBSERVATION COLLECTION

    py::class_< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                 "ObservationCollection",
                                                                                                 R"doc(

         Class collecting all observations and associated data for use in an estimation.

         Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
         this class is instantiated via a call to the :func:`~tudatpy.estimation.observations_setup.observations_wrapper.simulate_observations` function. More information is provided
         on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_





      )doc" )
            .def( py::init< std::vector< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > > >( ),
                  py::arg( "observation_sets" ),
                  R"doc(
Constructor for the ObservationCollection class.

Parameters
----------
observation_sets : list[tudatpy.numerical_simulation.estimation.SingleObservationSet]
    List of single observation sets to be included in the collection.
)doc" )
            .def_property_readonly( "concatenated_times",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleTimeVector,
                                    R"doc(

         **read-only**

         Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "concatenated_times_objects",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedTimeVector,
                                    R"doc(

         **read-only**

         Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly( "concatenated_weights",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getUnparsedConcatenatedWeights,
                                    R"doc(
         **read-only**

         Vector containing concatenated observation weights.

         :type: numpy.ndarray[numpy.float64[m, 1]]
)doc" )
            .def_property_readonly( "concatenated_observations",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationVector,
                                    R"doc(

         **read-only**

         Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

         :type: numpy.ndarray[numpy.float64[m, 1]]
      )doc" )
            .def_property_readonly(
                    "concatenated_link_definition_ids",
                    py::overload_cast<>( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                    R"doc(

         **read-only**

         Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order of the present vector.

         :type: numpy.ndarray[ int ]
      )doc" )
            .def_property_readonly( "link_definition_ids",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getInverseLinkEndIdentifierMap,
                                    R"doc(

         **read-only**

         Dictionaty mapping a link end integer identifier to the specific link ends

         :type: dict[ int, dict[ LinkEndType, LinkEndId ] ]
      )doc" )
            .def_property_readonly( "observable_type_start_index_and_size",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTypeStartAndSize,
                                    R"doc(

         **read-only**

         Dictionary defining per obervable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type

         :type: dict[ ObservableType, [ int, int ] ]
      )doc" )
            .def_property_readonly(
                    "observation_set_start_index_and_size",
                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationSetStartAndSizePerLinkEndIndex,
                    R"doc(

         **read-only**

         The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the correspoding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per obervable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.

         :type: dict[ ObservableType, dict[ int, list[ int, int ] ] ]
      )doc" )
            .def_property_readonly( "observation_vector_size",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTotalObservableSize,
                                    R"doc(

         **read-only**

         Length of the total vector of observations

         :type: int
      )doc" )
            .def_property_readonly( "sorted_observation_sets",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSets,
                                    R"doc(

         **read-only**

         The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_

         :type: dict[ ObservableType, dict[ int, list[ SingleObservationSet ] ] ]
      )doc" )
            .def_property_readonly( "link_ends_per_observable_type",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkEndsPerObservableType,
                                    R"doc(
         **read-only**

         Dictionary mapping each observable type to a list of link ends for which observations are available.

         :type: dict[ ObservableType, list[LinkEnds] ]
)doc" )
            .def_property_readonly( "link_definitions_per_observable",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsPerObservable,
                                    R"doc(
         **read-only**

         Dictionary mapping each observable type to a list of link definitions.

         :type: dict[ ObservableType, list[LinkDefinition] ]
)doc" )
            .def_property_readonly( "time_bounds",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsDouble,
                                    R"doc(
         **read-only**

         Pair of minimum and maximum observation time in the collection.

         :type: tuple[float, float]
)doc" )
            .def_property_readonly( "time_bounds_time_object",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                    R"doc(
         **read-only**

         Pair of minimum and maximum observation time in the collection, as Time objects.

         :type: tuple[Time, Time]
)doc" )
            .def_property_readonly( "sorted_per_set_time_bounds",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSetsTimeBounds,
                                    R"doc(
         **read-only**

         Nested dictionary providing the time bounds for each single observation set, sorted by observable type and link end ID.

         :type: dict[ ObservableType, dict[ int, list[tuple[float, float]] ] ]
)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(
Function to reset the full list of observable values. The order of the observations must be the same as for :attr:`~ObservationCollection.concatenated_observations`

Parameters
----------
observations : np.ndarray
    New list of observable values
     )doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  py::arg( "observation_parser" ),
                  R"doc(
Function to reset a subset of all observable values, with this subset defined by the ``observation_parser`` input.
The order and size of the new observation vector must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
an ``ObservationCollection`` containing only the parsed observation.

Parameters
----------
observations : np.ndarray
    New list of observable values
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of observations that is to be reset
     )doc" )
            .def( "set_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations_per_parser" ),
                  R"doc(
Function to reset a subset of all observable values, with this subset defined by a list of observation parsers input.
Each observation parser is associated with a new set of observable values.
The order and size of the new observation vector for each parser must be the same as when calling :attr:`~ObservationCollection.concatenated_observations` on
an ``ObservationCollection`` containing only the parsed observation (from a single parser). NOTE: since the multiple parsers
are handled in order (iterating over the keys of ``observations_per_parser``) some observations may be reset several times,
in case.

Parameters
----------
observations : np.ndarray
    New list of observable values
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of observations that is to be reset
     )doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(
Function to reset the full list of observation residuals. The order of the residuals must be the same as for :attr:`~ObservationCollection.concatenated_observations`

Parameters
----------
residuals : np.ndarray
    New list of observation residuals
     )doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  py::arg( "observation_parser" ),
                  R"doc(
Function to reset a subset of all observation residuals, with this subset defined by the ``observation_parser`` input.

Parameters
----------
residuals : np.ndarray
    New list of observation residuals
observation_parser : ObservationCollectionParser
    Observation parser with which to select the subset of residuals that is to be reset
     )doc" )
            .def( "set_residuals",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals_per_parser" ),
                  R"doc(
Function to reset a subset of all observation residuals, with this subset defined by a list of observation parsers input.

Parameters
----------
residuals_per_parser : dict[ObservationCollectionParser, np.ndarray]
    Dictionary mapping observation parsers to new residual values.
     )doc" )
            .def( "get_link_definitions_for_observables",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsForSingleObservable,
                  py::arg( "observable_type" ),
                  R"doc(
         Function to get all link definitions for a given observable type.

         Parameters
         ----------
         observable_type : :class:`ObservableType`
             Observable type for which link definitions are to be retrieved.
         Returns
         -------
         list[ LinkDefinition ]
             List of link definitions for the given observable type.
     )doc" )
            .def( "get_single_link_and_type_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleLinkAndTypeObservationSets,
                  py::arg( "observable_type" ),
                  py::arg( "link_definition" ),
                  R"doc(
         Function to get all observation sets for a given observable type and link definition.


         Parameters
         ----------
         observable_type : :class:`ObservableType`
             Observable type of which observations are to be simulated.
         link_ends : LinkDefinition
             Link ends for which observations are to be simulated.
         Returns
         -------
         list[ SingleObservationSet ]
             List of observation sets for given observable type and link definition.
     )doc" )
            .def( "get_observable_types",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservableTypes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observable types for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[ObservableType]
             A list of observable types present in the selected subset.
     )doc" )
            .def( "get_bodies_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getBodiesInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the names of bodies present in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of body names.
     )doc" )
            .def( "get_reference_points_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getReferencePointsInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the names of reference points (e.g., ground stations) in the link ends of a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[str]
             A list of reference point names.
     )doc" )
            .def( "get_time_bounds_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsListDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs (start_time, end_time).
     )doc" )
            .def( "get_time_bounds_list_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsList,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs (start_time, end_time) as Time objects.
     )doc" )
            .def( "get_time_bounds_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSetDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for each set in a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[float, float]]
             A list of time bound pairs for each observation set.
     )doc" )
            .def( "get_time_bounds_per_set_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSet,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the time bounds for each set in a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tuple[Time, Time]]
             A list of time bound pairs for each observation set as Time objects.
     )doc" )
            .def( "get_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of observation vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of observation values.
     )doc" )
            .def( "get_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimesDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[float]]
             A list of lists of observation times, one for each selected observation set.
     )doc" )
            .def( "get_observation_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[list[Time]]
             A list of lists of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observation_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observation times for a subset of observation sets as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[Time]
             A concatenated list of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observation times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[float]
             A concatenated list of observation times.
     )doc" )
            .def( "get_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimesDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[float]]]
             A pair containing a list of observation vectors and a list of lists of observation times.
     )doc" )
            .def( "get_observations_and_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], list[list[Time]]]
             A pair containing a list of observation vectors and a list of lists of observation times as Time objects.
     )doc" )
            .def( "get_concatenated_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimesDouble,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observations and times for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[float]]
             A pair containing the concatenated observation vector and list of times.
     )doc" )
            .def( "get_concatenated_observations_and_times_objects",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated observations and times for a subset of observation sets, with times as Time objects.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         tuple[numpy.ndarray, list[Time]]
             A pair containing the concatenated observation vector and list of times as Time objects.
     )doc" )
            .def( "get_concatenated_link_definition_ids",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                  py::arg( "observation_parser" ),
                  R"doc(
         Get the concatenated link definition IDs for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser
             Object that is used to select a subset of the observation sets.

         Returns
         -------
         list[int]
             A list of concatenated link end IDs.
     )doc" )
            .def( "get_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of weight vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated weights for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of weights.
     )doc" )
            .def( "get_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of residuals.
     )doc" )
            .def( "get_rms_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the RMS of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of RMS residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_mean_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the mean of residuals for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of mean residual vectors, one for each selected observation set.
     )doc" )
            .def( "get_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[numpy.ndarray]
             A list of computed observation vectors, one for each selected observation set.
     )doc" )
            .def( "get_concatenated_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get the concatenated computed observations for a subset of observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         numpy.ndarray
             The concatenated vector of computed observation values.
     )doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const double, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a constant weight for a subset of observation sets.

         Parameters
         ----------
         weight : float
             The constant weight to set.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a constant weight vector for a subset of observation sets.

         Parameters
         ----------
         weight : numpy.ndarray
             The constant weight vector to set.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, double > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  R"doc(
         Set a constant weight for multiple subsets of observation sets.

         Parameters
         ----------
         weights_per_observation_parser : dict[tudatpy.numerical_simulation.estimation.ObservationCollectionParser, float]
             A dictionary mapping observation parsers to constant weights.
     )doc" )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  R"doc(
         Set a constant weight vector for multiple subsets of observation sets.

         Parameters
         ----------
         weights_per_observation_parser : dict[tudatpy.numerical_simulation.estimation.ObservationCollectionParser, numpy.ndarray]
             A dictionary mapping observation parsers to constant weight vectors.
     )doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set tabulated weights for a subset of observation sets.

         Parameters
         ----------
         tabulated_weights : numpy.ndarray
             The vector of tabulated weights to set.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  R"doc(
         Set tabulated weights for multiple subsets of observation sets.

         Parameters
         ----------
         tabulated_weights : dict[tudatpy.numerical_simulation.estimation.ObservationCollectionParser, numpy.ndarray]
             A dictionary mapping observation parsers to tabulated weight vectors.
     )doc" )
            .def( "append",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::appendObservationCollection,
                  py::arg( "observation_collection_to_append" ) )
            .def( "filter_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     std::shared_ptr< tom::ObservationFilterBase > >&,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(
         Filter observations using a set of filters.

         This function filters the observations in the collection based on a map of observation filters, each associated with an observation parser.

         Parameters
         ----------
         observation_filters : dict[tudatpy.numerical_simulation.estimation.ObservationCollectionParser, tudatpy.numerical_simulation.estimation.ObservationFilterBase]
             A dictionary mapping observation parsers to observation filters.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
            .def( "filter_observations",
                  py::overload_cast< std::shared_ptr< tom::ObservationFilterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser >,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(
         Filter observations using a single filter.

         This function filters a subset of observations (or all) using a single observation filter.

         Parameters
         ----------
         observation_filter : tudatpy.numerical_simulation.estimation.ObservationFilterBase
             The observation filter to apply.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
         save_filtered_observations : bool, optional
             If true, the filtered-out observations are saved within each observation set, by default True.
     )doc" )
            .def( "split_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationSetSplitterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::splitObservationSets ),
                  py::arg( "observation_set_splitter" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Split observation sets based on a splitter.

         This function splits a subset of observation sets (or all) into smaller sets based on the criteria defined by the splitter.

         Parameters
         ----------
         observation_set_splitter : tudatpy.numerical_simulation.estimation.ObservationSetSplitterBase
             The splitter to use for splitting the observation sets.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, applying to all observation sets.
     )doc" )
            .def( "get_single_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservationSets,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get a subset of single observation sets.

         This function returns a list of the single observation sets that are selected by the provided observation parser.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Object that is used to select a subset of the observation sets, by default an empty parser, retrieving all observation sets.

         Returns
         -------
         list[tudatpy.numerical_simulation.estimation.SingleObservationSet]
             A list of the selected single observation sets.
     )doc" )
            .def( "print_observation_sets_start_and_size",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::printObservationSetsStartAndSize,
                  R"doc(Prints the structure of the observation collection, showing the start index and size of each individual observation set.)doc" )
            .def( "remove_single_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationSets ),
                  py::arg( "observation_parser" ),
                  R"doc(
         Remove a subset of single observation sets.

         Parameters
         ----------
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser
             Object that is used to select the subset of observation sets to remove.
     )doc" )
            .def( "set_reference_point",
                  py::overload_cast< tss::SystemOfBodies&,
                                     const Eigen::Vector3d&,
                                     const std::string&,
                                     const std::string&,
                                     const tom::LinkEndType,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                  py::arg( "bodies" ),
                  py::arg( "antenna_position" ),
                  py::arg( "antenna_name" ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "link_end_type" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a fixed reference point for a subset of observations.

         Parameters
         ----------
         bodies : tudatpy.numerical_simulation.environment.SystemOfBodies
             System of bodies.
         antenna_position : numpy.ndarray
             Position of the antenna in the spacecraft body-fixed frame.
         antenna_name : str
             Name of the antenna/reference point.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference point should be applied.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select the observation sets to which the reference point should be applied.
     )doc" )
            .def( "set_reference_points",
                  py::overload_cast< tss::SystemOfBodies&,
                                     const std::map< double, Eigen::Vector3d >&,
                                     const std::string&,
                                     const tom::LinkEndType,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoints ),
                  py::arg( "bodies" ),
                  py::arg( "antenna_switch_history" ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "link_end_type" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set multiple reference points based on an antenna switch history.

         Parameters
         ----------
         bodies : tudatpy.numerical_simulation.environment.SystemOfBodies
             System of bodies.
         antenna_switch_history : dict[float, numpy.ndarray]
             Dictionary mapping time to antenna position.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference points should be applied.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
            .def( "set_reference_point",
                  py::overload_cast< tss::SystemOfBodies&,
                                     const std::shared_ptr< te::Ephemeris >,
                                     const std::string&,
                                     const std::string&,
                                     const tom::LinkEndType,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setReferencePoint ),
                  py::arg( "bodies" ),
                  py::arg( "antenna_body_fixed_ephemeris" ),
                  py::arg( "antenna_name" ),
                  py::arg( "spacecraft_name" ),
                  py::arg( "link_end_type" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set a time-varying reference point for a subset of observations.

         Parameters
         ----------
         bodies : tudatpy.numerical_simulation.environment.SystemOfBodies
             System of bodies.
         antenna_body_fixed_ephemeris : tudatpy.numerical_simulation.ephemerides.Ephemeris
             Ephemeris of the antenna in the spacecraft body-fixed frame.
         antenna_name : str
             Name of the antenna/reference point.
         spacecraft_name : str
             Name of the spacecraft body.
         link_end_type : LinkEndType
             Link end type to which the reference point should be applied.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
            .def( "set_transponder_delay",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTransponderDelay,
                  py::arg( "spacecraft_name" ),
                  py::arg( "transponder_delay" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Set the transponder delay for a subset of observations.

         Parameters
         ----------
         spacecraft_name : str
             Name of the spacecraft with the transponder.
         transponder_delay : float
             The transponder delay in seconds.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select the observation sets.
     )doc" )
            .def( "remove_empty_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeEmptySingleObservationSets,
                  R"doc(Remove all single observation sets that contain no observations.)doc" )
            .def( "add_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariable,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Add an observation dependent variable to a subset of the single observation sets.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to add.
         bodies : tudatpy.numerical_simulation.environment.SystemOfBodies
             System of bodies containing the environment.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select the observation sets to which the variable should be added.
         
         Returns
         -------
         tudatpy.numerical_simulation.estimation.ObservationCollectionParser
             A parser that can be used to retrieve the added dependent variable.
     )doc" )
            .def( "dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the values of a given dependent variable, sorted per single observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[numpy.ndarray], tudatpy.numerical_simulation.estimation.ObservationCollectionParser]
             A pair containing a list of matrices with the dependent variable values and the parser used.
     )doc" )
            .def( "concatenated_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the concatenated values of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[numpy.ndarray, tudatpy.numerical_simulation.estimation.ObservationCollectionParser]
             A pair containing a matrix with the concatenated dependent variable values and the parser used.
     )doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get a list of all dependent variable settings compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings]], tudatpy.numerical_simulation.estimation.ObservationCollectionParser]
             A pair containing a list of lists of compatible settings and the parser used.
     )doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Get all dependent variables compatible with the input settings.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         tuple[list[list[numpy.ndarray]], tudatpy.numerical_simulation.estimation.ObservationCollectionParser]
             A pair containing a list of lists of dependent variable values and the parser used.
     )doc" )
            .def( "dependent_variable_history_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSetDouble,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         list[dict[float, numpy.ndarray]]
             A list of maps from time to dependent variable value, one for each set.
     )doc" )
            .def( "dependent_variable_history_per_set_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSet,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the time history of a given dependent variable, sorted per observation set, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         list[dict[Time, numpy.ndarray]]
             A list of maps from time to dependent variable value, one for each set, with times as Time objects.
     )doc" )
            .def( "dependent_variable_history",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryDouble,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the concatenated time history of a given dependent variable.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[float, numpy.ndarray]
             A map from time to dependent variable value.
     )doc" )
            .def( "dependent_variable_history_time_object",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(
         Retrieve the concatenated time history of a given dependent variable, with times as Time objects.

         Parameters
         ----------
         dependent_variable_settings : tudatpy.numerical_simulation.estimation.ObservationDependentVariableSettings
             Settings for the dependent variable to retrieve.
         first_compatible_settings : bool, optional
             If true, returns the first compatible variable found, by default False.
         observation_parser : tudatpy.numerical_simulation.estimation.ObservationCollectionParser, optional
             Parser to select a subset of observation sets.

         Returns
         -------
         dict[Time, numpy.ndarray]
             A map from time to dependent variable value, with times as Time objects.
     )doc" );

    m.def( "compute_residuals_and_dependent_variables",
           &tss::computeResidualsAndDependentVariables< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection" ),
           py::arg( "observation_simulators" ),
           py::arg( "bodies" ),
           R"doc(
        Computes residuals and dependent variables for a given observation collection.

        This function simulates observations based on the settings of the input `observation_collection`
        (which typically contains real data). It then computes the residuals by subtracting the
        simulated observations from the original observations. The computed residuals and any
        observation-dependent variables are then stored in the input `observation_collection`.

        Parameters
        ----------
        observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection of observations for which to compute residuals and dependent variables.
            This object is modified in-place.
        observation_simulators : list[tudatpy.astro.observation_models.ObservationSimulatorBase]
            List of observation simulators to be used for computing the observations.
        bodies : tudatpy.numerical_simulation.environment.SystemOfBodies
            The system of bodies required for the observation simulation.
        )doc" );

    m.def( "filter_observations",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const bool >( &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_filter" ),
           py::arg( "save_filtered_observations" ) = false,
           R"doc(
        Filters a single observation set and returns a new set containing the filtered observations.

        This function creates a copy of the input observation set and applies the given filter.
        The original observation set is not modified.

        Parameters
        ----------
        original_observation_set : tudatpy.estimation.observations.SingleObservationSet
            The observation set to filter.
        observation_filter : tudatpy.estimation.observations.ObservationFilterBase
            The filter to apply.
        save_filtered_observations : bool, optional
            If true, the observations that are filtered out are stored within the new observation set. Defaults to false.

        Returns
        -------
        tudatpy.estimation.observations.SingleObservationSet
            A new observation set containing only the observations that passed the filter.
        )doc" );

    m.def( "split_observation_set",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationSetSplitterBase >,
                              const bool >( &tom::splitObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_splitter" ),
           py::arg( "print_warning" ) = true,
           R"doc(
        Splits a single observation set into multiple sets based on a splitter.

        This function takes an observation set and divides it into a list of smaller observation sets
        according to the criteria defined in the observation splitter.

        Parameters
        ----------
        original_observation_set : tudatpy.estimation.observations.SingleObservationSet
            The observation set to split.
        observation_splitter : tudatpy.estimation.observations.ObservationSetSplitterBase
            The splitter defining the splitting criteria.
        print_warning : bool, optional
            If true, a warning is printed if the original set contains filtered observations that will be lost. Defaults to true.

        Returns
        -------
        list[tudatpy.estimation.observations.SingleObservationSet]
            A list of new observation sets resulting from the split.
        )doc" );

    m.def( "merge_observation_collections",
           &tss::mergeObservationCollections< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection_list" ) );

    // The following functions create a new ObservationCollection object from an existing one

    m.def( "create_filtered_observation_collection",
           py::overload_cast<
                   const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                   const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >& >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filters_map" ),
           R"doc(

        Creates a new, filtered observation collection from an existing one using multiple filters.

        This function applies a set of filters, each associated with a parser, to an observation collection
        and returns a new collection with the filtered observations. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection to filter.
        observation_filters_map : dict[tudatpy.estimation.observations.ObservationCollectionParser, tudatpy.estimation.observations.ObservationFilterBase]
            A dictionary mapping parsers to filters.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

    m.def( "create_filtered_observation_collection",
           py::overload_cast< const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const std::shared_ptr< tom::ObservationCollectionParser > >(
                   &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_collection" ),
           py::arg( "observation_filter" ),
           py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(

        Creates a new, filtered observation collection from an existing one using a single filter.

        This function applies a single filter to a subset of an observation collection (selected by a parser)
        and returns a new collection with the filtered observations. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection to filter.
        observation_filter : tudatpy.estimation.observations.ObservationFilterBase
            The filter to apply.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select the subset of observations to filter. Defaults to an empty parser (all observations).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the filtered observations.

        )doc" );

    m.def( "split_observation_collection",
           &tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_set_splitter" ),
           py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(
        Creates a new observation collection by splitting sets from an existing collection.

        This function splits observation sets from the original collection based on a splitter and
        returns a new collection containing the resulting sets. The original collection is not modified.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The observation collection from which to split sets.
        observation_set_splitter : tudatpy.estimation.observations.ObservationSetSplitterBase
            The splitter defining how to split the sets.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select which observation sets to split. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection with the split observation sets.
        )doc" );

    m.def( "create_new_observation_collection",
           &tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "original_observation_collection" ),
           py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
           R"doc(
        Creates a new observation collection containing a subset of an existing collection.

        This function selects a subset of observation sets from an original collection using a parser
        and creates a new collection containing only those sets.

        Parameters
        ----------
        original_observation_collection : tudatpy.estimation.observations.ObservationCollection
            The collection from which to extract a subset.
        observation_parser : tudatpy.estimation.observations.ObservationCollectionParser, optional
            Parser to select the observation sets to include in the new collection. Defaults to an empty parser (all sets).

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollection
            A new observation collection containing the selected subset of observation sets.
        )doc" );
}

}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
