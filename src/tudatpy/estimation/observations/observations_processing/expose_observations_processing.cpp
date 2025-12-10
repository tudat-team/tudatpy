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
#include "expose_observations_processing.h"

#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "scalarTypes.h"
#include "tudat/simulation/estimation_setup/observationsProcessing.h"

namespace py = pybind11;
namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

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

    m.def( "observation_filter",
           py::overload_cast< tom::ObservationFilterType, const double, const bool, const bool >( &tom::observationFilter ),
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
        tudatpy.estimation.observations.ObservationFilterBase
            An observation filter object.
        )doc" );

    m.def( "observation_filter",
           py::overload_cast< tom::ObservationFilterType, const std::vector< double >, const bool, const bool >( &tom::observationFilter ),
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
        tudatpy.estimation.observations.ObservationFilterBase
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

    m.def( "observation_filter",
           py::overload_cast< tom::ObservationFilterType, const double, const double, const bool, const bool >( &tom::observationFilter ),
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
        tudatpy.estimation.observations.ObservationFilterBase
            An observation filter object.
        )doc" );

    m.def( "observation_filter",
           py::overload_cast< tom::ObservationFilterType, const Eigen::VectorXd&, const bool, const bool >( &tom::observationFilter ),
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
        tudatpy.estimation.observations.ObservationFilterBase
            An observation filter object.
        )doc" );

    m.def( "observation_filter",
           py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >,
                              const Eigen::VectorXd&,
                              const bool,
                              const bool >( &tom::observationFilter ),
           py::arg( "dependent_variable_settings" ),
           py::arg( "filter_value" ),
           py::arg( "filter_out" ) = true,
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create a dependent variable observation filter.

        Parameters
        ----------
        dependent_variable_settings : tudatpy.kernel.simulation.estimation_setup.ObservationDependentVariableSettings
            Settings for the dependent variable to be used for filtering.
        filter_value : numpy.ndarray
            Numpy array to be used by the filter.
        filter_out : bool, optional
            Whether to filter out observations that satisfy the condition (True) or keep them (False). Default is True.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationFilterBase
            An observation filter object.
        )doc" );

    py::class_< tom::ObservationSetSplitterBase, std::shared_ptr< tom::ObservationSetSplitterBase > >(
            m, "ObservationSetSplitterBase", R"doc(
        Base class for observation set splitters.

        This is the base class from which all observation set splitters are derived. It is not intended to be instantiated directly.
        )doc" );

    m.def( "observation_set_splitter",
           py::overload_cast< tom::ObservationSetSplitterType, const std::vector< double >, const int >( &tom::observationSetSplitter ),
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
        tudatpy.estimation.observations.ObservationSetSplitterBase
            An observation set splitter object.
        )doc" );

    m.def( "observation_set_splitter",
           py::overload_cast< tom::ObservationSetSplitterType, const double, const int >( &tom::observationSetSplitter ),
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
        tudatpy.estimation.observations.ObservationSetSplitterBase
            An observation set splitter object.
        )doc" );

    m.def( "observation_set_splitter",
           py::overload_cast< tom::ObservationSetSplitterType, const int, const int >( &tom::observationSetSplitter ),
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
        tudatpy.estimation.observations.ObservationSetSplitterBase
            An observation set splitter object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast<>( &tom::observationParser ),
           R"doc(
        Create an empty observation parser.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An empty observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const tom::ObservableType, const bool >( &tom::observationParser ),
           py::arg( "observable_type" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on observable type.

        Parameters
        ----------
        observable_type : tudatpy.kernel.astro.ObservableType
            Observable type to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::vector< tom::ObservableType >&, const bool >( &tom::observationParser ),
           py::arg( "observable_type_vector" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a list of observable types.

        Parameters
        ----------
        observable_type_vector : list[tudatpy.kernel.astro.ObservableType]
            List of observable types to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const tom::LinkEnds, const bool >( &tom::observationParser ),
           py::arg( "link_ends" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on link ends.

        Parameters
        ----------
        link_ends : tudatpy.kernel.astro.LinkEnds
            Link ends to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::vector< tom::LinkEnds >&, const bool >( &tom::observationParser ),
           py::arg( "link_ends_vector" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a list of link ends.

        Parameters
        ----------
        link_ends_vector : list[tudatpy.kernel.astro.LinkEnds]
            List of link ends to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::string, const bool, const bool >( &tom::observationParser ),
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

    m.def( "observation_parser",
           py::overload_cast< const std::vector< std::string >&, const bool, const bool >( &tom::observationParser ),
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

    m.def( "observation_parser",
           py::overload_cast< const std::pair< std::string, std::string >&, const bool >( &tom::observationParser ),
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

    m.def( "observation_parser",
           py::overload_cast< const std::vector< std::pair< std::string, std::string > >&, const bool >( &tom::observationParser ),
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

    m.def( "observation_parser",
           py::overload_cast< const tom::LinkEndType&, const bool >( &tom::observationParser ),
           py::arg( "link_end_type" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a link end type.

        Parameters
        ----------
        link_end_type : tudatpy.kernel.astro.LinkEndType
            Link end type to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::vector< tom::LinkEndType >&, const bool >( &tom::observationParser ),
           py::arg( "link_end_types_vector" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a list of link end types.

        Parameters
        ----------
        link_end_types_vector : list[tudatpy.kernel.astro.LinkEndType]
            List of link end types to parse.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::pair< tom::LinkEndType, tom::LinkEndId >&, const bool >( &tom::observationParser ),
           py::arg( "single_link_end" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a single link end (type and ID).

        Parameters
        ----------
        single_link_end : tuple[tudatpy.kernel.astro.LinkEndType, tudatpy.kernel.astro.LinkEndId]
            A single link end, specified by its type and ID.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::vector< std::pair< tom::LinkEndType, tom::LinkEndId > >&, const bool >( &tom::observationParser ),
           py::arg( "single_link_ends_vector" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a list of single link ends (type and ID).

        Parameters
        ----------
        single_link_ends_vector : list[tuple[tudatpy.kernel.astro.LinkEndType, tudatpy.kernel.astro.LinkEndId]]
            A list of single link ends, each specified by its type and ID.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::pair< double, double >&, const bool >( &tom::observationParser ),
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

    m.def( "observation_parser",
           py::overload_cast< const std::vector< std::pair< double, double > >&, const bool >( &tom::observationParser ),
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

    m.def( "observation_parser",
           py::overload_cast< const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >, const bool >(
                   &tom::observationParser ),
           py::arg( "ancillary_settings" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on ancillary settings.

        Parameters
        ----------
        ancillary_settings : tudatpy.kernel.astro.ObservationAncilliarySimulationSettings
            Ancillary settings for parsing.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >&, const bool >(
                   &tom::observationParser ),
           py::arg( "ancillary_settings_vector" ),
           py::arg( "use_opposite_condition" ) = false,
           R"doc(
        Create an observation parser based on a list of ancillary settings.

        Parameters
        ----------
        ancillary_settings_vector : list[tudatpy.kernel.astro.ObservationAncilliarySimulationSettings]
            List of ancillary settings for parsing.
        use_opposite_condition : bool, optional
            Whether to use the opposite of the default condition. Default is False.

        Returns
        -------
        tudatpy.estimation.observations.ObservationCollectionParser
            An observation parser object.
        )doc" );

    m.def( "observation_parser",
           py::overload_cast< const std::vector< std::shared_ptr< tom::ObservationCollectionParser > >&, const bool >(
                   &tom::observationParser ),
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
