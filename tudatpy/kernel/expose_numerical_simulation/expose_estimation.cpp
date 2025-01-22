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



void expose_estimation(py::module &m) {

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


    py::class_<tom::ObservationViabilityCalculator,
            std::shared_ptr<tom::ObservationViabilityCalculator>>(m, "ObservationViabilityCalculator",
                                                                  get_docstring("ObservationViabilityCalculator").c_str() )
            .def("is_observation_viable", &tom::ObservationViabilityCalculator::isObservationViable,
                 py::arg( "link_end_states" ),
                 py::arg( "link_end_times" ),
                 get_docstring("ObservationViabilityCalculator.is_observation_viable").c_str() );

    py::class_<tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>>(m, "ObservationSimulator",
                                                                           get_docstring("ObservationSimulator").c_str() );

    py::class_<tom::ObservationSimulator<1,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<1,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_1",
                                                          get_docstring("ObservationSimulator_1").c_str() );

    py::class_<tom::ObservationSimulator<2,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<2,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_2",
                                                          get_docstring("ObservationSimulator_2").c_str() );

    py::class_<tom::ObservationSimulator<3,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<3,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_3",
                                                          get_docstring("ObservationSimulator_3").c_str() );

    py::class_<tom::ObservationSimulator<6,STATE_SCALAR_TYPE, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<6,STATE_SCALAR_TYPE, TIME_TYPE>>,
            tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>>(m, "ObservationSimulator_6",
                                                          get_docstring("ObservationSimulator_6").c_str() );

    m.def("simulate_observations",
          &tss::simulateObservations<STATE_SCALAR_TYPE, TIME_TYPE>,
          py::arg("simulation_settings"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("simulate_observations").c_str() );

    m.def("compute_residuals_and_dependent_variables",
          &tss::computeResidualsAndDependentVariables<STATE_SCALAR_TYPE, TIME_TYPE>,
          py::arg("observation_collection"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("compute_and_set_residuals").c_str() );


    m.def("create_pseudo_observations_and_models",
          &tss::simulatePseudoObservations<TIME_TYPE, STATE_SCALAR_TYPE>,
          py::arg("bodies"),
          py::arg("observed_bodies" ),
          py::arg("central_bodies" ),
          py::arg("initial_time"),
          py::arg("final_time"),
          py::arg("time_step"),
          get_docstring("create_pseudo_observations_and_models").c_str() );

    m.def("create_best_fit_to_ephemeris",
          &tss::createBestFitToCurrentEphemeris<TIME_TYPE, STATE_SCALAR_TYPE>,
          py::arg("bodies"),
          py::arg("acceleration_models"),
          py::arg("observed_bodies" ),
          py::arg("central_bodies" ),
          py::arg("integrator_settings" ),
          py::arg("initial_time"),
          py::arg("final_time"),
          py::arg("data_point_interval"),
          py::arg("additional_parameter_names") = std::vector< std::shared_ptr< tep::EstimatableParameterSettings > >( ),
          py::arg("number_of_iterations") = 3,
          py::arg("reintegrate_variational_equations") = true,
          py::arg("results_print_frequency") = 0.0,
          get_docstring("create_best_fit_to_ephemeris").c_str() );

    m.def("set_existing_observations",
          &tss::setExistingObservations<STATE_SCALAR_TYPE, TIME_TYPE>,
          py::arg("observations"),
          py::arg("reference_link_end" ),
          py::arg("ancilliary_settings_per_observatble" ) = std::map< tom::ObservableType,
          std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >( ) );

    m.def("compute_target_angles_and_range",
          &tss::getTargetAnglesAndRange,
          py::arg("bodies"),
          py::arg("station_id" ),
          py::arg("target_body" ),
          py::arg("observation_times"),
          py::arg("is_station_transmitting"),
          get_docstring("compute_target_angles_and_range").c_str() );

    
    m.def("create_filtered_observation_collection",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                  const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >& >(
          &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
          py::arg("original_observation_collection"),
          py::arg("observation_filters_map"),
          get_docstring("create_filtered_observation_collection").c_str() );

    m.def("create_filtered_observation_collection",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationFilterBase >,
                  const std::shared_ptr< tom::ObservationCollectionParser > >(
          &tom::filterObservations< STATE_SCALAR_TYPE, TIME_TYPE > ),
          py::arg("original_observation_collection"),
          py::arg("observation_filter"),
          py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
          get_docstring("create_filtered_observation_collection").c_str() );

    m.def("split_observation_collection",
          &tom::splitObservationSets< STATE_SCALAR_TYPE, TIME_TYPE >,
          py::arg("original_observation_collection"),
          py::arg("observation_set_splitter"),
          py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
          get_docstring("split_observation_collection").c_str() );

    m.def("create_new_observation_collection",
         &tom::createNewObservationCollection< STATE_SCALAR_TYPE, TIME_TYPE >,
         py::arg("original_observation_collection"),
         py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
         get_docstring("create_new_observation_collection ").c_str() );





}

}
}
}// namespace tudatpy
