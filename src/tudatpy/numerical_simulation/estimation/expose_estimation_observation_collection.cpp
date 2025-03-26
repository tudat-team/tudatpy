/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "expose_estimation.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"
#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;

namespace tudatpy
{
namespace numerical_simulation
{
namespace estimation
{

void expose_estimation_observation_collection( py::module& m )
{
    py::class_< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
                std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > > >( m,
                                                                                                 "ObservationCollection",
                                                                                                 R"doc(

        Class collecting all observations and associated data for use in an estimation.

        Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
        this class is instantiated via a call to the :func:`~tudatpy.numerical_simulation.estimation.simulate_observations` function. More information is provided
        on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_





     )doc" )
            .def( py::init< std::vector< std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > > > >( ),
                  py::arg( "observation_sets" ) )
            .def_property_readonly( "concatenated_times",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedTimeVector,
                                    R"doc(

        **read-only**

        Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc" )
            .def_property_readonly( "concatenated_float_times",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleTimeVector,
                                    R"doc(

        **read-only**

        Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc" )
            .def_property_readonly( "concatenated_weights",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getUnparsedConcatenatedWeights,
                                    R"doc(No documentation found.)doc" )
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
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "link_definitions_per_observable",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsPerObservable,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "time_bounds",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBounds,
                                    R"doc(No documentation found.)doc" )
            .def_property_readonly( "sorted_per_set_time_bounds",
                                    &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSortedObservationSetsTimeBounds,
                                    R"doc(No documentation found.)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_observations",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations" ),
                  py::arg( "observation_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setObservations ),
                  py::arg( "observations_per_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_residuals",
                  py::overload_cast< const Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 >&,
                                     const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals" ),
                  py::arg( "observation_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_residuals",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     Eigen::Matrix< STATE_SCALAR_TYPE, Eigen::Dynamic, 1 > >& >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setResiduals ),
                  py::arg( "residuals_per_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "get_link_definitions_for_observables",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getLinkDefinitionsForSingleObservable,
                  py::arg( "observable_type" ),
                  R"doc(No documentation found.)doc" )
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
                  R"doc(No documentation found.)doc" )
            .def( "get_bodies_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getBodiesInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_reference_points_in_link_ends",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getReferencePointsInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_time_bounds_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsList,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_time_bounds_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getTimeBoundsPerSet,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_float_observation_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDoubleObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_observations_and_times",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_link_definition_ids",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedLinkEndIds ),
                  py::arg( "observation_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "get_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_weights",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_rms_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getRmsResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_mean_residuals",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getMeanResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_concatenated_computed_observations",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const double, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, double > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  R"doc(No documentation found.)doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< const Eigen::VectorXd, const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "set_tabulated_weights",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  R"doc(No documentation found.)doc" )
            .def( "append",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::appendObservationCollection,
                  py::arg( "observation_collection_to_append" ) )
            .def( "filter_observations",
                  py::overload_cast< const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                                                     std::shared_ptr< tom::ObservationFilterBase > >&,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(No documentation found.)doc" )
            .def( "filter_observations",
                  py::overload_cast< std::shared_ptr< tom::ObservationFilterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser >,
                                     const bool >( &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  py::arg( "save_filtered_observations" ) = true,
                  R"doc(No documentation found.)doc" )
            .def( "split_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationSetSplitterBase >,
                                     std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::splitObservationSets ),
                  py::arg( "observation_set_splitter" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "get_single_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getSingleObservationSets,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "print_observation_sets_start_and_size",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::printObservationSetsStartAndSize,
                  R"doc(No documentation found.)doc" )
            .def( "remove_single_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeSingleObservationSets ),
                  py::arg( "observation_parser" ),
                  R"doc(No documentation found.)doc" )
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
                  R"doc(No documentation found.)doc" )
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
                  R"doc(No documentation found.)doc" )
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
                  R"doc(No documentation found.)doc" )
            .def( "set_transponder_delay",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::setTransponderDelay,
                  py::arg( "spacecraft_name" ),
                  py::arg( "transponder_delay" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "remove_empty_observation_sets",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::removeEmptySingleObservationSets,
                  R"doc(No documentation found.)doc" )
            .def( "add_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::addDependentVariable,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "bodies" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "concatenated_dependent_variable",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getConcatenatedDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "compatible_dependent_variable_settings",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "compatible_dependent_variables_list",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getAllCompatibleDependentVariables,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "dependent_variable_history_per_set",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistoryPerObservationSet,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" )
            .def( "dependent_variable_history",
                  &tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >::getDependentVariableHistory,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "first_compatible_settings" ) = false,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  R"doc(No documentation found.)doc" );

    m.def( "merge_observation_collections",
           &tss::mergeObservationCollections< STATE_SCALAR_TYPE, TIME_TYPE >,
           py::arg( "observation_collection_list" ) );

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
           R"doc(No documentation found.)doc" );

    m.def( "filter_observations",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationFilterBase >,
                              const bool >( &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_filter" ),
           py::arg( "save_filtered_observations" ) = false,
           R"doc(No documentation found.)doc" );

    m.def( "split_observation_set",
           py::overload_cast< const std::shared_ptr< tom::SingleObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > >,
                              const std::shared_ptr< tom::ObservationSetSplitterBase >,
                              const bool >( &tom::splitObservationSet< STATE_SCALAR_TYPE, TIME_TYPE > ),
           py::arg( "original_observation_set" ),
           py::arg( "observation_splitter" ),
           py::arg( "print_warning" ) = true,
           R"doc(No documentation found.)doc" );
}

}  // namespace estimation
}  // namespace numerical_simulation
}  // namespace tudatpy
