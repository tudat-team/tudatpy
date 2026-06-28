/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "expose_observations_processing.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/astro/observation_models/observationAncillarySettings.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"

namespace py = pybind11;
namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

namespace
{

const char* legacyObservationProcessingDeprecationGuide =
        "https://docs.tudat.space/en/latest/user-guide/state-estimation/observation-dataset-deprecation.html";

void warnLegacyObservationProcessingInterface( const std::string& interfaceName )
{
    const std::string message = interfaceName +
            " is deprecated and kept only for backwards compatibility with legacy observation processing. Use ObservationDataset "
            "conditions/viewers instead. Migration guide: " +
            legacyObservationProcessingDeprecationGuide;
    if( PyErr_WarnEx( PyExc_DeprecationWarning, message.c_str( ), 1 ) < 0 )
    {
        throw py::error_already_set( );
    }
}

}  // namespace

namespace tudatpy
{
namespace estimation
{
namespace observations
{
namespace observations_processing
{

void expose_observations_processing( py::module& m )
{
    py::enum_< tom::ObservationFilterType >( m, "ObservationFilterType", R"doc(
        Enum for types of observation filters.

        This enum defines the available types of observation filters that can be used to reject observations from a collection.
        )doc" )
            .value( "residual_filtering", tom::ObservationFilterType::residual_filtering )
            .value( "absolute_value_filtering", tom::ObservationFilterType::absolute_value_filtering )
            .value( "epochs_filtering", tom::ObservationFilterType::epochs_filtering )
            .value( "time_bounds_filtering", tom::ObservationFilterType::time_bounds_filtering )
            .value( "dependent_variable_filtering", tom::ObservationFilterType::dependent_variable_filtering )
            .export_values( );

    py::enum_< tom::ObservationSetSplitterType >( m, "ObservationSetSplitterType", R"doc(
        Enum for types of observation set splitters.

        This enum defines the available types of observation set splitters that can be used to divide a collection of observations into multiple sets.
        )doc" )
            .value( "time_tags_splitter", tom::ObservationSetSplitterType::time_tags_splitter )
            .value( "time_interval_splitter", tom::ObservationSetSplitterType::time_interval_splitter )
            .value( "time_span_splitter", tom::ObservationSetSplitterType::time_span_splitter )
            .value( "nb_observations_splitter", tom::ObservationSetSplitterType::nb_observations_splitter )
            .export_values( );

    py::enum_< tom::ObservationParserType >( m, "ObservationParserType", R"doc(
        Enum for types of observation parsers.

        This enum defines the available types of observation parsers that can be used to select observations from a collection based on various criteria.
        )doc" )
            .value( "empty_parser", tom::ObservationParserType::empty_parser )
            .value( "observable_type_parser", tom::ObservationParserType::observable_type_parser )
            .value( "link_ends_parser", tom::ObservationParserType::link_ends_parser )
            .value( "link_end_str_parser", tom::ObservationParserType::link_end_string_parser )
            .value( "link_end_id_parser", tom::ObservationParserType::link_end_id_parser )
            .value( "link_end_type_parser", tom::ObservationParserType::link_end_type_parser )
            .value( "single_link_end_parser", tom::ObservationParserType::single_link_end_parser )
            .value( "time_bounds_parser", tom::ObservationParserType::time_bounds_parser )
            .value( "ancillary_settings_parser", tom::ObservationParserType::ancillary_settings_parser )
            .value( "multi_type_parser", tom::ObservationParserType::multi_type_parser )
            .export_values( );

    py::class_< tom::ObservationFilterBase, std::shared_ptr< tom::ObservationFilterBase > >( m, "ObservationFilterBase", R"doc(
        Base class for observation filters.

        This is the base class from which all observation filters are derived. It is not intended to be instantiated directly.
        )doc" );

    py::class_< tom::ObservationCollectionParser, std::shared_ptr< tom::ObservationCollectionParser > >(
            m, "ObservationCollectionParser", R"doc(
        Base class for observation collection parsers.

        This is the base class from which all observation collection parsers are derived. It is not intended to be instantiated directly.
        )doc" );

    m.def(
            "observation_filter",
            []( const tom::ObservationFilterType filterType,
                const double filterValue,
                const bool filterOut,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_filter" );
                return tom::observationFilter( filterType, filterValue, filterOut, useOppositeCondition );
            },
            py::arg( "filter_type" ),
            py::arg( "filter_value" ),
            py::arg( "filter_out" ) = true,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation filter with a single double value.

        Parameters
        ----------
        filter_type : tudatpy.estimation.observations.ObservationFilterType
            Type of observation filter.
        filter_value : float
            Value to be used by the filter.
        filter_out : bool, optional
            Whether to filter out observations that satisfy the condition (True) or keep them (False). Default is True.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            An observation filter object.
        )doc" );

    m.def(
            "observation_filter",
            []( const tom::ObservationFilterType filterType,
                const std::vector< double > filterValue,
                const bool filterOut,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_filter" );
                return tom::observationFilter( filterType, filterValue, filterOut, useOppositeCondition );
            },
            py::arg( "filter_type" ),
            py::arg( "filter_value" ),
            py::arg( "filter_out" ) = true,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation filter with a list of double values.

        Parameters
        ----------
        filter_type : tudatpy.estimation.observations.ObservationFilterType
            Type of observation filter.
        filter_value : list[float]
            List of values to be used by the filter.
        filter_out : bool, optional
            Whether to filter out observations that satisfy the condition (True) or keep them (False). Default is True.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            An observation filter object.
        )doc" );

    //    m.def("observation_filter",
    //          py::overload_cast< tom::ObservationFilterType,
    //          const std::pair< double, double >, const bool, const
    //          bool >( &tom::observationFilter ),
    //          py::arg("filter_type"),
    //          py::arg("filter_value"),
    //          py::arg("filter_out") = true,
    //          py::arg("use_opposite_condition") = false,
    //          get_docstring("observation_filter").c_str() );

    m.def(
            "observation_filter",
            []( const tom::ObservationFilterType filterType,
                const double firstFilterValue,
                const double secondFilterValue,
                const bool filterOut,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_filter" );
                return tom::observationFilter( filterType, firstFilterValue, secondFilterValue, filterOut, useOppositeCondition );
            },
            py::arg( "filter_type" ),
            py::arg( "first_filter_value" ),
            py::arg( "second_filter_value" ),
            py::arg( "filter_out" ) = true,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation filter with two double values (e.g., for a time range).

        Parameters
        ----------
        filter_type : tudatpy.estimation.observations.ObservationFilterType
            Type of observation filter.
        first_filter_value : float
            First value to be used by the filter (e.g., start time).
        second_filter_value : float
            Second value to be used by the filter (e.g., end time).
        filter_out : bool, optional
            Whether to filter out observations that satisfy the condition (True) or keep them (False). Default is True.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            An observation filter object.
        )doc" );

    m.def(
            "observation_filter",
            []( const tom::ObservationFilterType filterType,
                const Eigen::VectorXd& filterValue,
                const bool filterOut,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_filter" );
                return tom::observationFilter( filterType, filterValue, filterOut, useOppositeCondition );
            },
            py::arg( "filter_type" ),
            py::arg( "filter_value" ),
            py::arg( "filter_out" ) = true,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation filter with a numpy array.

        Parameters
        ----------
        filter_type : tudatpy.estimation.observations.ObservationFilterType
            Type of observation filter.
        filter_value : numpy.ndarray
            Numpy array to be used by the filter.
        filter_out : bool, optional
            Whether to filter out observations that satisfy the condition (True) or keep them (False). Default is True.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            An observation filter object.
        )doc" );

    m.def(
            "observation_filter",
            []( const std::shared_ptr< tss::ObservationDependentVariableSettings > dependentVariableSettings,
                const Eigen::VectorXd& filterValue,
                const bool filterOut,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_filter" );
                return tom::observationFilter( dependentVariableSettings, filterValue, filterOut, useOppositeCondition );
            },
            py::arg( "dependent_variable_settings" ),
            py::arg( "filter_value" ),
            py::arg( "filter_out" ) = true,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create a dependent variable observation filter.

        Parameters
        ----------
        dependent_variable_settings : :class:`~tudatpy.estimation.observations_setup.observations_dependent_variables.ObservationDependentVariableSettings`
            Settings for the dependent variable to be used for filtering.
        filter_value : numpy.ndarray
            Numpy array to be used by the filter.
        filter_out : bool, optional
            Whether to filter out observations that satisfy the condition (True) or keep them (False). Default is True.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationFilterBase
            An observation filter object.
        )doc" );

    py::class_< tom::ObservationSetSplitterBase, std::shared_ptr< tom::ObservationSetSplitterBase > >(
            m, "ObservationSetSplitterBase", R"doc(
        Base class for observation set splitters.

        This is the base class from which all observation set splitters are derived. It is not intended to be instantiated directly.
        )doc" );

    m.def(
            "observation_set_splitter",
            []( const tom::ObservationSetSplitterType splitterType,
                const std::vector< double > splitterValue,
                const int minNumberObservations ) {
                warnLegacyObservationProcessingInterface( "observation_set_splitter" );
                return tom::observationSetSplitter( splitterType, splitterValue, minNumberObservations );
            },
            py::arg( "splitter_type" ),
            py::arg( "splitter_value" ),
            py::arg( "min_number_observations" ) = 0,
            R"doc(
        Create an observation set splitter with a list of double values.

        Parameters
        ----------
        splitter_type : tudatpy.estimation.observations.ObservationSetSplitterType
            Type of observation set splitter.
        splitter_value : list[float]
            List of values to be used by the splitter.
        min_number_observations : int, optional
            Minimum number of observations per split set. Default is 0.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            An observation set splitter object.
        )doc" );

    m.def(
            "observation_set_splitter",
            []( const tom::ObservationSetSplitterType splitterType, const double splitterValue, const int minNumberObservations ) {
                warnLegacyObservationProcessingInterface( "observation_set_splitter" );
                return tom::observationSetSplitter( splitterType, splitterValue, minNumberObservations );
            },
            py::arg( "splitter_type" ),
            py::arg( "splitter_value" ),
            py::arg( "min_number_observations" ) = 0,
            R"doc(
        Create an observation set splitter with a single double value.

        Parameters
        ----------
        splitter_type : tudatpy.estimation.observations.ObservationSetSplitterType
            Type of observation set splitter.
        splitter_value : float
            Value to be used by the splitter.
        min_number_observations : int, optional
            Minimum number of observations per split set. Default is 0.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            An observation set splitter object.
        )doc" );

    m.def(
            "observation_set_splitter",
            []( const tom::ObservationSetSplitterType splitterType, const int splitterValue, const int minNumberObservations ) {
                warnLegacyObservationProcessingInterface( "observation_set_splitter" );
                return tom::observationSetSplitter( splitterType, splitterValue, minNumberObservations );
            },
            py::arg( "splitter_type" ),
            py::arg( "splitter_value" ),
            py::arg( "min_number_observations" ) = 0,
            R"doc(
        Create an observation set splitter with a single integer value.

        Parameters
        ----------
        splitter_type : tudatpy.estimation.observations.ObservationSetSplitterType
            Type of observation set splitter.
        splitter_value : int
            Value to be used by the splitter.
        min_number_observations : int, optional
            Minimum number of observations per split set. Default is 0.

        Returns
        -------
        tudatpy.estimation.observations.observations_processing.ObservationSetSplitterBase
            An observation set splitter object.
        )doc" );

    m.def(
            "observation_parser",
            []( ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( );
            },
            R"doc(
        Create an empty observation parser.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An empty observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const tom::ObservableType observableType, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( observableType, useOppositeCondition );
            },
            py::arg( "observable_type" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on observable type.

        Parameters
        ----------
        observable_type : :class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`
            Observable type to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< tom::ObservableType >& observableTypeVector, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( observableTypeVector, useOppositeCondition );
            },
            py::arg( "observable_type_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of observable types.

        Parameters
        ----------
        observable_type_vector : list[:class:`~tudatpy.estimation.observable_models_setup.model_settings.ObservableType`]
            List of observable types to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const tom::LinkEnds linkEnds, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEnds, useOppositeCondition );
            },
            py::arg( "link_ends" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on link ends.

        Parameters
        ----------
        link_ends : dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]
            Link ends to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< tom::LinkEnds >& linkEndsVector, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndsVector, useOppositeCondition );
            },
            py::arg( "link_ends_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of link ends.

        Parameters
        ----------
        link_ends_vector : list[dict[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]]
            List of link ends to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::string linkEndsString, const bool isReferencePoint, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndsString, isReferencePoint, useOppositeCondition );
            },
            py::arg( "link_ends_str" ),
            py::arg( "is_reference_point" ) = false,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a link end string (body name).

        Parameters
        ----------
        link_ends_str : str
            Name of the body involved in the link end.
        is_reference_point : bool, optional
            Whether the body is a reference point. Default is False.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< std::string >& linkEndsStringVector, const bool isReferencePoint, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndsStringVector, isReferencePoint, useOppositeCondition );
            },
            py::arg( "link_ends_str_vector" ),
            py::arg( "is_reference_point" ) = false,
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of link end strings (body names).

        Parameters
        ----------
        link_ends_str_vector : list[str]
            List of names of bodies involved in the link ends.
        is_reference_point : bool, optional
            Whether the bodies are reference points. Default is False.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::pair< std::string, std::string >& linkEndId, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndId, useOppositeCondition );
            },
            py::arg( "link_end_id" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a link end ID.

        Parameters
        ----------
        link_end_id : tuple[str, str]
            Link end ID, as a tuple of (body name, station name).
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< std::pair< std::string, std::string > >& linkEndIdsVector, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndIdsVector, useOppositeCondition );
            },
            py::arg( "link_end_ids_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of link end IDs.

        Parameters
        ----------
        link_end_ids_vector : list[tuple[str, str]]
            List of link end IDs, each as a tuple of (body name, station name).
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const tom::LinkEndType& linkEndType, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndType, useOppositeCondition );
            },
            py::arg( "link_end_type" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a link end type.

        Parameters
        ----------
        link_end_type : :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`
            Link end type to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< tom::LinkEndType >& linkEndTypesVector, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( linkEndTypesVector, useOppositeCondition );
            },
            py::arg( "link_end_types_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of link end types.

        Parameters
        ----------
        link_end_types_vector : list[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`]
            List of link end types to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::pair< tom::LinkEndType, tom::LinkEndId >& singleLinkEnd, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( singleLinkEnd, useOppositeCondition );
            },
            py::arg( "single_link_end" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a single link end (type and ID).

        Parameters
        ----------
        single_link_end : tuple[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]
            A single link end, specified by its type and ID.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< std::pair< tom::LinkEndType, tom::LinkEndId > >& singleLinkEndsVector,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( singleLinkEndsVector, useOppositeCondition );
            },
            py::arg( "single_link_ends_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of single link ends (type and ID).

        Parameters
        ----------
        single_link_ends_vector : list[tuple[:class:`~tudatpy.estimation.observable_models_setup.links.LinkEndType`, :class:`~tudatpy.estimation.observable_models_setup.links.LinkEndId`]]
            A list of single link ends, each specified by its type and ID.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::pair< double, double >& timeBounds, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( timeBounds, useOppositeCondition );
            },
            py::arg( "time_bounds" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on time bounds.

        Parameters
        ----------
        time_bounds : tuple[float, float]
            Time bounds (start and end time) for parsing.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< std::pair< double, double > >& timeBoundsVector, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( timeBoundsVector, useOppositeCondition );
            },
            py::arg( "time_bounds_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of time bounds.

        Parameters
        ----------
        time_bounds_vector : list[tuple[float, float]]
            List of time bounds (start and end time) for parsing.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::shared_ptr< tom::ObservationAncillarySimulationSettings > ancillarySettings, const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( ancillarySettings, useOppositeCondition );
            },
            py::arg( "ancillary_settings" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on ancillary settings.

        Parameters
        ----------
        ancillary_settings : :class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings`
            Ancillary settings for parsing.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< std::shared_ptr< tom::ObservationAncillarySimulationSettings > >& ancillarySettingsVector,
                const bool useOppositeCondition ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( ancillarySettingsVector, useOppositeCondition );
            },
            py::arg( "ancillary_settings_vector" ),
            py::arg( "use_opposite_condition" ) = false,
            R"doc(
        Create an observation parser based on a list of ancillary settings.

        Parameters
        ----------
        ancillary_settings_vector : list[:class:`~tudatpy.estimation.observations_setup.ancillary_settings.ObservationAncillarySimulationSettings`]
            List of ancillary settings for parsing.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def(
            "observation_parser",
            []( const std::vector< std::shared_ptr< tom::ObservationCollectionParser > >& observationParsers,
                const bool combineConditions ) {
                warnLegacyObservationProcessingInterface( "observation_parser" );
                return tom::observationParser( observationParsers, combineConditions );
            },
            py::arg( "observation_parsers" ),
            py::arg( "combine_conditions" ) = false,
            R"doc(
        Create a multi-type observation parser from a list of other parsers.

        Parameters
        ----------
        observation_parsers : list[tudatpy.estimation.observations.ObservationCollectionParser]
            List of observation parsers to combine.
        combine_conditions : bool, optional
            If True, conditions are combined with AND (intersection). If False, with OR (union). Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            A multi-type observation parser object.
        )doc" );
}

}  // namespace observations_processing
}  // namespace observations
}  // namespace estimation
}  // namespace tudatpy
