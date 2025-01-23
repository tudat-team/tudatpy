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

#include "tudat/simulation/estimation_setup/fitOrbitToEphemeris.h"
#include "tudat/astro/propagators/propagateCovariance.h"
#include "tudat/basics/utilities.h"

#include "docstrings.h"
#include "scalarTypes.h"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>
#include <pybind11/chrono.h>

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;
namespace te = tudat::ephemerides;



namespace tudatpy {
namespace numerical_simulation {
namespace estimation {



void expose_estimation2(py::module &m) {

    /*!
     *************** PARAMETERS ***************
     */


    py::class_<tep::EstimatableParameterSet<STATE_SCALAR_TYPE>,
            std::shared_ptr<tep::EstimatableParameterSet<STATE_SCALAR_TYPE>>>(m, "EstimatableParameterSet",
                                                                   get_docstring("EstimatableParameterSet").c_str() )
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getEstimatedParameterSetSize,
                                    get_docstring("EstimatableParameterSet.parameter_set_size").c_str() )
            .def_property_readonly( "initial_states_size",
                                    &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getInitialDynamicalStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_states_size").c_str() )
            .def_property_readonly( "initial_single_arc_states_size",
                                    &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getInitialDynamicalSingleArcStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_single_arc_states_size").c_str() )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getInitialDynamicalMultiArcStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_multi_arc_states_size").c_str() )
            .def_property_readonly( "constraints_size",
                                    &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getConstraintSize,
                                    get_docstring("EstimatableParameterSet.constraints_size").c_str() )
            .def_property( "parameter_vector",
                           &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getFullParameterValues< double >,
                           &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::resetParameterValues< double >,
                           get_docstring("EstimatableParameterSet.parameter_vector").c_str() )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::getIndicesForParameterType,
                  py::arg("parameter_type"),
                  get_docstring("EstimatableParameterSet.indices_for_parameter_type").c_str() );

    /*!
     *************** OBSERVATIONS ***************
     */

    py::enum_< tom::ObservationFilterType >(m, "ObservationFilterType",
                                               get_docstring("ObservationFilterType").c_str() )
            .value("residual_filtering", tom::ObservationFilterType::residual_filtering )
            .value("absolute_value_filtering", tom::ObservationFilterType::absolute_value_filtering )
            .value("epochs_filtering", tom::ObservationFilterType::epochs_filtering )
            .value("time_bounds_filtering", tom::ObservationFilterType::time_bounds_filtering )
            .value("dependent_variable_filtering", tom::ObservationFilterType::dependent_variable_filtering )
            .export_values();

    py::enum_< tom::ObservationSetSplitterType >(m, "ObservationSetSplitterType",
                                               get_docstring("ObservationSetSplitterType").c_str() )
            .value("time_tags_splitter", tom::ObservationSetSplitterType::time_tags_splitter )
            .value("time_interval_splitter", tom::ObservationSetSplitterType::time_interval_splitter )
            .value("time_span_splitter", tom::ObservationSetSplitterType::time_span_splitter )
            .value("nb_observations_splitter", tom::ObservationSetSplitterType::nb_observations_splitter )
            .export_values();

    py::enum_< tom::ObservationParserType >(m, "ObservationParserType",
                                               get_docstring("ObservationParserType").c_str() )
            .value("empty_parser", tom::ObservationParserType::empty_parser )
            .value("observable_type_parser", tom::ObservationParserType::observable_type_parser )
            .value("link_ends_parser", tom::ObservationParserType::link_ends_parser )
            .value("link_end_str_parser", tom::ObservationParserType::link_end_string_parser )
            .value("link_end_id_parser", tom::ObservationParserType::link_end_id_parser )
            .value("link_end_type_parser", tom::ObservationParserType::link_end_type_parser )
            .value("single_link_end_parser", tom::ObservationParserType::single_link_end_parser )
            .value("time_bounds_parser", tom::ObservationParserType::time_bounds_parser )
            .value("ancillary_settings_parser", tom::ObservationParserType::ancillary_settings_parser )
            .value("multi_type_parser", tom::ObservationParserType::multi_type_parser )
            .export_values();

    py::class_<tom::ObservationFilterBase,
            std::shared_ptr<tom::ObservationFilterBase>>(m, "ObservationFilterBase",
                    get_docstring("ObservationFilterBase").c_str() );

    m.def("observation_filter",
          py::overload_cast< tom::ObservationFilterType,
          const double, const bool, const bool >( &tom::observationFilter ),
          py::arg("filter_type"),
          py::arg("filter_value"),
          py::arg("filter_out") = true,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_filter").c_str() );

    m.def("observation_filter",
          py::overload_cast< tom::ObservationFilterType,
          const std::vector< double >, const bool, const bool >( &tom::observationFilter ),
          py::arg("filter_type"),
          py::arg("filter_value"),
          py::arg("filter_out") = true,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_filter").c_str() );

//    m.def("observation_filter",
//          py::overload_cast< tom::ObservationFilterType,
//          const std::pair< double, double >, const bool, const bool >( &tom::observationFilter ),
//          py::arg("filter_type"),
//          py::arg("filter_value"),
//          py::arg("filter_out") = true,
//          py::arg("use_opposite_condition") = false,
//          get_docstring("observation_filter").c_str() );

    m.def("observation_filter",
          py::overload_cast< tom::ObservationFilterType, const double, const double, const bool, const bool >( &tom::observationFilter ),
          py::arg("filter_type"),
          py::arg("first_filter_value"),
          py::arg("second_filter_value"),
          py::arg("filter_out") = true,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_filter").c_str() );

    m.def("observation_filter",
          py::overload_cast< tom::ObservationFilterType, const Eigen::VectorXd&, const bool, const bool >( &tom::observationFilter ),
          py::arg("filter_type"),
          py::arg("filter_value"),
          py::arg("filter_out") = true,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_filter").c_str() );

    m.def("observation_filter",
          py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >, const Eigen::VectorXd&, const bool, const bool >(
                  &tom::observationFilter ),
          py::arg("dependent_variable_settings"),
          py::arg("filter_value"),
          py::arg("filter_out") = true,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_filter").c_str() );

    py::class_<tom::ObservationSetSplitterBase,
            std::shared_ptr<tom::ObservationSetSplitterBase>>(m, "ObservationSetSplitterBase",
                    get_docstring("ObservationSetSplitterBase").c_str() );

    m.def("observation_set_splitter",
      py::overload_cast< tom::ObservationSetSplitterType,
        const std::vector< double >, const int >( &tom::observationSetSplitter ),
      py::arg("splitter_type"),
      py::arg("splitter_value"),
      py::arg("min_number_observations") = 0,
      get_docstring("observation_set_splitter").c_str() );

    m.def("observation_set_splitter",
          py::overload_cast< tom::ObservationSetSplitterType,
                  const double, const int >( &tom::observationSetSplitter),
          py::arg("splitter_type"),
          py::arg("splitter_value"),
          py::arg("min_number_observations") = 0,
          get_docstring("observation_set_splitter").c_str() );

    m.def("observation_set_splitter",
          py::overload_cast< tom::ObservationSetSplitterType,
                  const int, const int >( &tom::observationSetSplitter ),
          py::arg("splitter_type"),
          py::arg("splitter_value"),
          py::arg("min_number_observations") = 0,
          get_docstring("observation_set_splitter").c_str() );

    py::class_<tom::ObservationCollectionParser,
                std::shared_ptr<tom::ObservationCollectionParser>>(m, "ObservationCollectionParser",
                        get_docstring("ObservationCollectionParser").c_str() );

    m.def("observation_parser",
          py::overload_cast<  >( &tom::observationParser ),
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const tom::ObservableType, const bool >( &tom::observationParser ),
          py::arg("observable_type"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< tom::ObservableType >&, const bool >( &tom::observationParser ),
          py::arg("observable_type_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const tom::LinkEnds, const bool >( &tom::observationParser ),
          py::arg("link_ends"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< tom::LinkEnds >&, const bool >( &tom::observationParser ),
          py::arg("link_ends_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::string, const bool, const bool >( &tom::observationParser ),
          py::arg("link_ends_str"),
          py::arg( "is_reference_point") = false,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< std::string >&, const bool, const bool >( &tom::observationParser ),
          py::arg("link_ends_str_vector"),
          py::arg( "is_reference_point") = false,
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::pair< std::string, std::string >&, const bool >( &tom::observationParser ),
          py::arg("link_end_id"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< std::pair< std::string, std::string > >&, const bool >( &tom::observationParser ),
          py::arg("link_end_ids_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const tom::LinkEndType&, const bool >( &tom::observationParser ),
          py::arg("link_end_type"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< tom::LinkEndType >&, const bool >( &tom::observationParser ),
          py::arg("link_end_types_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::pair< tom::LinkEndType, tom::LinkEndId >&, const bool >( &tom::observationParser ),
          py::arg("single_link_end"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< std::pair< tom::LinkEndType, tom::LinkEndId > >&, const bool >( &tom::observationParser ),
          py::arg("single_link_ends_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::pair< double, double >&, const bool >( &tom::observationParser ),
          py::arg("time_bounds"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast< const std::vector< std::pair< double, double > >&, const bool >( &tom::observationParser ),
          py::arg("time_bounds_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationAncilliarySimulationSettings >,
                  const bool >( &tom::observationParser ),
          py::arg("ancillary_settings"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast<
                  const std::vector< std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >&,
                  const bool >( &tom::observationParser ),
          py::arg("ancillary_settings_vector"),
          py::arg("use_opposite_condition") = false,
          get_docstring("observation_parser").c_str() );

    m.def("observation_parser",
          py::overload_cast<
                  const std::vector< std::shared_ptr< tom::ObservationCollectionParser > >&,
                  const bool >( &tom::observationParser ),
          py::arg("observation_parsers"),
          py::arg("combine_conditions") = false,
          get_docstring("observation_parser").c_str() );







}

}
}
}// namespace tudatpy
