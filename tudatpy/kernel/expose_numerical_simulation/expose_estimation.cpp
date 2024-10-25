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


namespace tudat
{

namespace propagators
{

std::map< double, Eigen::MatrixXd > propagateCovarianceRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput< double, TIME_TYPE > > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<double, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance(
                propagatedCovariance, estimationOutput->getUnnormalizedCovarianceMatrix( ),
                orbitDeterminationManager->getStateTransitionAndSensitivityMatrixInterface( ), evaluationTimes );

    tss::SystemOfBodies bodies = orbitDeterminationManager->getBodies( );

    std::shared_ptr< tep::EstimatableParameterSet<double> > parameterSet = orbitDeterminationManager->getParametersToEstimate( );

    std::map< int, std::shared_ptr< tep::EstimatableParameter< Eigen::VectorXd > > > initialStates = parameterSet->getInitialStateParameters( );
    std::map< std::pair< std::string, std::string >, std::vector< int > > transformationList;
    for( auto it : initialStates )
    {
        if( std::dynamic_pointer_cast< tep::InitialTranslationalStateParameter< double > >( it.second ) )
        {
            std::shared_ptr< tep::InitialTranslationalStateParameter< double > > currentInitialState =
                    std::dynamic_pointer_cast< tep::InitialTranslationalStateParameter< double > >( it.second );
            transformationList[ std::make_pair( currentInitialState->getParameterName( ).second.first, currentInitialState->getCentralBody( ) ) ].push_back(
                        it.first );

        }
        else if( std::dynamic_pointer_cast< tep::ArcWiseInitialTranslationalStateParameter< double > >( it.second ) )
        {
            throw std::runtime_error( "Error, multi-arc not yet supported in automatic covariance conversion" );
        }
    }

    Eigen::Matrix3d currentInertialToRswPosition;
    Eigen::Matrix6d currentInertialToRswState;
    Eigen::MatrixXd currentFullInertialToRswState = Eigen::MatrixXd::Zero( 6, 6 );

    std::map< double, Eigen::MatrixXd > propagatedRswCovariance;
    for( auto it : propagatedCovariance )
    {
        double currentTime = it.first;
        Eigen::MatrixXd currentCovariance = it.second;
        currentFullInertialToRswState.setZero( );

        for( auto it_body : transformationList )
        {
            Eigen::Vector6d relativeState =
                    bodies.getBody( it_body.first.first )->getStateInBaseFrameFromEphemeris( currentTime ) -
                    bodies.getBody( it_body.first.second )->getStateInBaseFrameFromEphemeris( currentTime );
            currentInertialToRswPosition = trf::getInertialToRswSatelliteCenteredFrameRotationMatrix( relativeState );
            currentInertialToRswState.block( 0, 0, 3, 3 ) = currentInertialToRswPosition;
            currentInertialToRswState.block( 3, 3, 3, 3 ) = currentInertialToRswPosition;
            for( unsigned int j = 0; j < it_body.second.size( ); j++ )
            {
                int currentStartIndex = it_body.second.at( j );
                currentFullInertialToRswState.block( currentStartIndex, currentStartIndex, 6, 6 ) =  currentInertialToRswState;
            }
        }
        propagatedRswCovariance[ currentTime ] = currentFullInertialToRswState * currentCovariance * currentFullInertialToRswState.transpose( );

    }
    return propagatedRswCovariance;
}


std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput<double, TIME_TYPE> > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<double, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedRswCovariance = propagateCovarianceRsw(
            estimationOutput, orbitDeterminationManager, evaluationTimes );

    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedRswCovariance ),
                           utilities::createVectorFromMapValues(
                               propagatedRswCovariance ) );
}

std::map< double, Eigen::VectorXd > propagateFormalErrorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput<double, TIME_TYPE> > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<double, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;

    propagatedCovariance = propagateCovarianceRsw(
                estimationOutput, orbitDeterminationManager, evaluationTimes );
    tp::convertCovarianceHistoryToFormalErrorHistory( propagatedFormalErrors, propagatedCovariance );

    return propagatedFormalErrors;
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectorsRsw(
        const std::shared_ptr< tss::CovarianceAnalysisOutput<double, TIME_TYPE> > estimationOutput,
        const std::shared_ptr< tss::OrbitDeterminationManager<double, TIME_TYPE> > orbitDeterminationManager,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors =
            propagateFormalErrorsRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );
    tp::propagateFormalErrorsRsw( estimationOutput, orbitDeterminationManager, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedFormalErrors ),
                           utilities::createVectorFromMapValues(
                               propagatedFormalErrors ) );
}



std::pair< std::vector< double >, std::vector< Eigen::MatrixXd > > propagateCovarianceVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::MatrixXd > propagatedCovariance;
    tp::propagateCovariance(
                propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedCovariance ),
                           utilities::createVectorFromMapValues(
                               propagatedCovariance ) );
}

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > propagateFormalErrorVectors(
        const Eigen::MatrixXd initialCovariance,
        const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface > stateTransitionInterface,
        const std::vector< double > evaluationTimes )
{
    std::map< double, Eigen::VectorXd > propagatedFormalErrors;
    tp::propagateFormalErrors( propagatedFormalErrors, initialCovariance, stateTransitionInterface, evaluationTimes );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               propagatedFormalErrors ),
                           utilities::createVectorFromMapValues(
                               propagatedFormalErrors ) );
}

}

namespace simulation_setup
{

std::pair< std::vector< double >, std::vector< Eigen::VectorXd > > getTargetAnglesAndRangeVector(
        const simulation_setup::SystemOfBodies& bodies,
        const std::pair< std::string, std::string > groundStationId,
        const std::string& targetBody,
        const std::vector< double > times,
        const bool transmittingToTarget )
{
    std::map< double, Eigen::VectorXd > targetAnglesAndRange = getTargetAnglesAndRange(
                bodies, groundStationId, targetBody, times, transmittingToTarget );
    return std::make_pair( utilities::createVectorFromMapKeys(
                               targetAnglesAndRange ),
                           utilities::createVectorFromMapValues(
                               targetAnglesAndRange ) );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< tom::SingleObservationSet< ObservationScalarType, TimeType > > singleObservationSetWithoutDependentVariables(
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds,
            const std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >& observations,
            const std::vector< TimeType > observationTimes,
            const tom::LinkEndType referenceLinkEnd,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr )
{
    return std::make_shared< tom::SingleObservationSet< ObservationScalarType, TimeType > >(
            observableType, linkEnds, observations, observationTimes, referenceLinkEnd,
            std::vector< Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > >( ), nullptr, ancilliarySettings );
}


}

}

namespace tudatpy {
namespace numerical_simulation {
namespace estimation {



void expose_estimation(py::module &m) {

    /*!
     *************** PARAMETERS ***************
     */


    py::class_<tep::EstimatableParameterSet<double>,
            std::shared_ptr<tep::EstimatableParameterSet<double>>>(m, "EstimatableParameterSet",
                                                                   get_docstring("EstimatableParameterSet").c_str() )
            .def_property_readonly( "parameter_set_size",
                                    &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize,
                                    get_docstring("EstimatableParameterSet.parameter_set_size").c_str() )
            .def_property_readonly( "initial_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_states_size").c_str() )
            .def_property_readonly( "initial_single_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalSingleArcStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_single_arc_states_size").c_str() )
            .def_property_readonly( "initial_multi_arc_states_size",
                                    &tep::EstimatableParameterSet<double>::getInitialDynamicalMultiArcStateParameterSize,
                                    get_docstring("EstimatableParameterSet.initial_multi_arc_states_size").c_str() )
            .def_property_readonly( "constraints_size",
                                    &tep::EstimatableParameterSet<double>::getConstraintSize,
                                    get_docstring("EstimatableParameterSet.constraints_size").c_str() )
            .def_property( "parameter_vector",
                           &tep::EstimatableParameterSet<double>::getFullParameterValues< double >,
                           &tep::EstimatableParameterSet<double>::resetParameterValues< double >,
                           get_docstring("EstimatableParameterSet.parameter_vector").c_str() )
            .def( "indices_for_parameter_type",
                  &tep::EstimatableParameterSet<double>::getIndicesForParameterType,
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

    m.def("observation_filter",
          py::overload_cast< tom::ObservationFilterType,
          const std::pair< double, double >, const bool, const bool >( &tom::observationFilter ),
          py::arg("filter_type"),
          py::arg("filter_value"),
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

    py::class_<tom::ObservationSimulatorBase<double, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulatorBase<double, TIME_TYPE>>>(m, "ObservationSimulator",
                                                                           get_docstring("ObservationSimulator").c_str() );

    py::class_<tom::ObservationSimulator<1,double, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<1,double, TIME_TYPE>>,
            tom::ObservationSimulatorBase<double, TIME_TYPE>>(m, "ObservationSimulator_1",
                                                          get_docstring("ObservationSimulator_1").c_str() );

    py::class_<tom::ObservationSimulator<2,double, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<2,double, TIME_TYPE>>,
            tom::ObservationSimulatorBase<double, TIME_TYPE>>(m, "ObservationSimulator_2",
                                                          get_docstring("ObservationSimulator_2").c_str() );

    py::class_<tom::ObservationSimulator<3,double, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<3,double, TIME_TYPE>>,
            tom::ObservationSimulatorBase<double, TIME_TYPE>>(m, "ObservationSimulator_3",
                                                          get_docstring("ObservationSimulator_3").c_str() );

    py::class_<tom::ObservationSimulator<6,double, TIME_TYPE>,
            std::shared_ptr<tom::ObservationSimulator<6,double, TIME_TYPE>>,
            tom::ObservationSimulatorBase<double, TIME_TYPE>>(m, "ObservationSimulator_6",
                                                          get_docstring("ObservationSimulator_6").c_str() );

    m.def("simulate_observations",
          &tss::simulateObservations<double, TIME_TYPE>,
          py::arg("simulation_settings"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("simulate_observations").c_str() );

    m.def("compute_residuals_and_dependent_variables",
          &tss::computeResidualsAndDependentVariables<double, TIME_TYPE>,
          py::arg("observation_collection"),
          py::arg("observation_simulators" ),
          py::arg("bodies"),
          get_docstring("compute_and_set_residuals").c_str() );


    m.def("create_pseudo_observations_and_models",
          &tss::simulatePseudoObservations<TIME_TYPE, double>,
          py::arg("bodies"),
          py::arg("observed_bodies" ),
          py::arg("central_bodies" ),
          py::arg("initial_time"),
          py::arg("final_time"),
          py::arg("time_step"),
          get_docstring("create_pseudo_observations_and_models").c_str() );

    m.def("create_best_fit_to_ephemeris",
          &tss::createBestFitToCurrentEphemeris<TIME_TYPE, double>,
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
          get_docstring("create_best_fit_to_ephemeris").c_str() );

    m.def("set_existing_observations",
          &tss::setExistingObservations<double, TIME_TYPE>,
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

    m.def("compute_target_angles_and_range_vectors",
          &tss::getTargetAnglesAndRangeVector,
          py::arg("bodies"),
          py::arg("station_id" ),
          py::arg("target_body" ),
          py::arg("observation_times"),
          py::arg("is_station_transmitting"),
          get_docstring("compute_target_angles_and_range").c_str() );

    m.def("create_filtered_observation_collection",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationCollection< double, TIME_TYPE > >,
                  const std::map< std::shared_ptr< tom::ObservationCollectionParser >, std::shared_ptr< tom::ObservationFilterBase > >& >(
          &tom::filterObservations< double, TIME_TYPE > ),
          py::arg("original_observation_collection"),
          py::arg("observation_filters_map"),
          get_docstring("create_filtered_observation_collection").c_str() );

    m.def("create_filtered_observation_collection",
          py::overload_cast<
                  const std::shared_ptr< tom::ObservationCollection< double, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationFilterBase >,
                  const std::shared_ptr< tom::ObservationCollectionParser > >(
          &tom::filterObservations< double, TIME_TYPE > ),
          py::arg("original_observation_collection"),
          py::arg("observation_filter"),
          py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
          get_docstring("create_filtered_observation_collection").c_str() );

    m.def("split_observation_collection",
          &tom::splitObservationSets< double, TIME_TYPE >,
          py::arg("original_observation_collection"),
          py::arg("observation_set_splitter"),
          py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
          get_docstring("split_observation_collection").c_str() );

    m.def("create_new_observation_collection",
         &tom::createNewObservationCollection< double, TIME_TYPE >,
         py::arg("original_observation_collection"),
         py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
         get_docstring("create_new_observation_collection ").c_str() );


    py::class_< tom::ObservationCollection<double, TIME_TYPE>,
            std::shared_ptr<tom::ObservationCollection<double, TIME_TYPE>>>(m, "ObservationCollection",
                                                           get_docstring("ObservationCollection").c_str() )
            .def(py::init< std::vector< std::shared_ptr< tom::SingleObservationSet< double, TIME_TYPE > > > >(),
                 py::arg("observation_sets") )
            .def_property_readonly("concatenated_times", &tom::ObservationCollection<double, TIME_TYPE>::getConcatenatedTimeVector,
                                   get_docstring("ObservationCollection.concatenated_times").c_str() )
            .def_property_readonly("concatenated_float_times", &tom::ObservationCollection<double, TIME_TYPE>::getConcatenatedDoubleTimeVector,
                                   get_docstring("ObservationCollection.concatenated_times").c_str() )
            .def_property_readonly("concatenated_weights", &tom::ObservationCollection<double, TIME_TYPE>::getConcatenatedWeights,
                                   get_docstring("ObservationCollection.concatenated_weights").c_str() )
            .def_property_readonly("concatenated_observations", &tom::ObservationCollection<double, TIME_TYPE>::getObservationVector,
                                   get_docstring("ObservationCollection.concatenated_observations").c_str() )
            .def_property_readonly("concatenated_link_definition_ids", &tom::ObservationCollection<double, TIME_TYPE>::getConcatenatedLinkEndIds,
                                   get_docstring("ObservationCollection.concatenated_link_definition_ids").c_str() )
            .def_property_readonly("link_definition_ids", &tom::ObservationCollection<double, TIME_TYPE>::getInverseLinkEndIdentifierMap,
                                   get_docstring("ObservationCollection.link_definition_ids").c_str() )
            .def_property_readonly("observable_type_start_index_and_size", &tom::ObservationCollection<double, TIME_TYPE>::getObservationTypeStartAndSize,
                                   get_docstring("ObservationCollection.observable_type_start_index_and_size").c_str() )
            .def_property_readonly("observation_set_start_index_and_size", &tom::ObservationCollection<double, TIME_TYPE>::getObservationSetStartAndSizePerLinkEndIndex,
                                   get_docstring("ObservationCollection.observation_set_start_index_and_size").c_str() )
            .def_property_readonly("observation_vector_size", &tom::ObservationCollection<double, TIME_TYPE>::getTotalObservableSize,
                                   get_docstring("ObservationCollection.observation_vector_size").c_str() )
            .def_property_readonly("sorted_observation_sets", &tom::ObservationCollection<double, TIME_TYPE>::getSortedObservationSets,
                                   get_docstring("ObservationCollection.sorted_observation_sets").c_str() )
            .def_property_readonly("link_ends_per_observable_type", &tom::ObservationCollection<double, TIME_TYPE>::getLinkEndsPerObservableType,
                                   get_docstring("ObservationCollection.link_ends_per_observable_type").c_str() )
            .def_property_readonly("link_definitions_per_observable", &tom::ObservationCollection<double, TIME_TYPE>::getLinkDefinitionsPerObservable,
                                   get_docstring("ObservationCollection.link_definitions_per_observable").c_str() )
            .def_property_readonly("time_bounds", &tom::ObservationCollection<double, TIME_TYPE>::getTimeBounds,
                                   get_docstring("ObservationCollection.time_bounds").c_str() )
            .def_property_readonly("sorted_per_set_time_bounds", &tom::ObservationCollection<double, TIME_TYPE>::getSortedObservationSetsTimeBounds,
                                   get_docstring("ObservationCollection.time_bounds").c_str() )
            .def("set_observations", py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                 &tom::ObservationCollection<double, TIME_TYPE>::setObservations ),
                 py::arg("observations"), get_docstring("set_observations").c_str( ) )
            .def("set_observations", py::overload_cast<
                const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
                const std::shared_ptr< tom::ObservationCollectionParser > >(
                &tom::ObservationCollection<double, TIME_TYPE>::setObservations ),
                py::arg("observations"), py::arg("observation_parser"), get_docstring("set_observations").c_str( ) )
            .def("set_observations", py::overload_cast<
                    const std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& >(
                    &tom::ObservationCollection<double, TIME_TYPE>::setObservations ),
                 py::arg("observations_per_parser"), get_docstring("set_observations").c_str( ) )
            .def("set_residuals", py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                 &tom::ObservationCollection<double, TIME_TYPE>::setResiduals ),
                 py::arg("residuals"), get_docstring("set_residuals").c_str( ) )
            .def("set_residuals", py::overload_cast<
                 const Eigen::Matrix< double, Eigen::Dynamic, 1 >&,
                 const std::shared_ptr< tom::ObservationCollectionParser > >(
                 &tom::ObservationCollection<double, TIME_TYPE>::setResiduals ),
                 py::arg("residuals"), py::arg("observation_parser"), get_docstring("set_residuals").c_str( ) )
            .def("set_residuals", py::overload_cast<
                    const std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::Matrix< double, Eigen::Dynamic, 1 > >& >(
                    &tom::ObservationCollection<double, TIME_TYPE>::setResiduals ),
                 py::arg("residuals_per_parser"), get_docstring("set_residuals").c_str( ) )
            .def("get_link_definitions_for_observables", &tom::ObservationCollection<double, TIME_TYPE>::getLinkDefinitionsForSingleObservable,
                 py::arg( "observable_type" ),
                 get_docstring("ObservationCollection.get_link_definitions_for_observables").c_str() )
            .def("get_single_link_and_type_observations", &tom::ObservationCollection<double, TIME_TYPE>::getSingleLinkAndTypeObservationSets,
                                   py::arg( "observable_type" ),
                                   py::arg( "link_definition" ),
                                   get_docstring("ObservationCollection.get_single_link_and_type_observations").c_str() )
            .def( "get_observable_types", &tom::ObservationCollection< double, TIME_TYPE >::getObservableTypes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observable_types" ).c_str() )
            .def( "get_bodies_in_link_ends", &tom::ObservationCollection< double, TIME_TYPE >::getBodiesInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_bodies_in_link_ends" ).c_str() )
            .def( "get_reference_points_in_link_ends", &tom::ObservationCollection< double, TIME_TYPE >::getReferencePointsInLinkEnds,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_reference_points_in_link_ends" ).c_str() )
            .def( "get_time_bounds_list", &tom::ObservationCollection< double, TIME_TYPE >::getTimeBoundsList,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_time_bounds_list" ).c_str() )
            .def( "get_observations", &tom::ObservationCollection< double, TIME_TYPE >::getObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observations" ).c_str() )
            .def( "get_concatenated_observations", &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_observations" ).c_str() )
            .def( "get_observation_times", &tom::ObservationCollection< double, TIME_TYPE >::getObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observation_times" ).c_str() )
            .def( "get_concatenated_observation_times", &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedObservationTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_observation_times" ).c_str() )
            .def( "get_observations_and_times", &tom::ObservationCollection< double, TIME_TYPE >::getObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_observations_and_times" ).c_str() )
            .def( "get_concatenated_observations_and_times", &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedObservationsAndTimes,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_observations_and_times" ).c_str() )
            .def( "get_weights", &tom::ObservationCollection< double, TIME_TYPE >::getWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_weights" ).c_str() )
            .def( "get_concatenated_weights", &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedWeights,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_weights" ).c_str() )
            .def( "get_residuals", &tom::ObservationCollection< double, TIME_TYPE >::getResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_residuals" ).c_str() )
            .def( "get_concatenated_residuals", &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedResiduals,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_residuals" ).c_str() )
            .def( "get_computed_observations", &tom::ObservationCollection< double, TIME_TYPE >::getComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_computed_observations" ).c_str() )
            .def( "get_concatenated_computed_observations",
                  &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedComputedObservations,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_concatenated_computed_observations" ).c_str() )
            .def( "set_constant_weight", py::overload_cast<
                    const double,
                    const std::shared_ptr< tom::ObservationCollectionParser > >(
                    &tom::ObservationCollection< double, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_constant_weight" ).c_str() )
            .def( "set_constant_weight", py::overload_cast<
                    const Eigen::VectorXd,
                    const std::shared_ptr< tom::ObservationCollectionParser > >(
                    &tom::ObservationCollection< double, TIME_TYPE >::setConstantWeight ),
                  py::arg( "weight" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_constant_weight" ).c_str() )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, double > >(
                          &tom::ObservationCollection< double, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  get_docstring( "set_constant_weight_per_observation_parser" ).c_str() )
            .def( "set_constant_weight_per_observation_parser",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< double, TIME_TYPE >::setConstantWeightPerObservable ),
                  py::arg( "weights_per_observation_parser" ),
                  get_docstring( "set_constant_weight_per_observation_parser" ).c_str() )
            .def( "set_tabulated_weights",
                  py::overload_cast<
                  const Eigen::VectorXd,
                  const std::shared_ptr< tom::ObservationCollectionParser > >(
                          &tom::ObservationCollection< double, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_tabulated_weights" ).c_str() )
            .def( "set_tabulated_weights",
                  py::overload_cast< std::map< std::shared_ptr< tom::ObservationCollectionParser >, Eigen::VectorXd > >(
                          &tom::ObservationCollection< double, TIME_TYPE >::setTabulatedWeights ),
                  py::arg( "tabulated_weights" ),
                  get_docstring( "set_tabulated_weights" ).c_str() )
            .def( "filter_observations",
                  py::overload_cast<
                  const std::map< std::shared_ptr< tom::ObservationCollectionParser >,
                  std::shared_ptr< tom::ObservationFilterBase > >&,
                  const bool >( &tom::ObservationCollection< double, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg("save_filtered_observations") = true,
                  get_docstring( "filter_observations" ).c_str( ) )
            .def( "filter_observations",
                  py::overload_cast<
                  std::shared_ptr< tom::ObservationFilterBase >,
                  std::shared_ptr< tom::ObservationCollectionParser >,
                  const bool >( &tom::ObservationCollection< double, TIME_TYPE >::filterObservations ),
                  py::arg( "observation_filters" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  py::arg("save_filtered_observations") = true,
                  get_docstring( "filter_observations" ).c_str(  ) )
            .def( "split_observation_sets",
                  py::overload_cast<
                  std::shared_ptr< tom::ObservationSetSplitterBase >,
                  std::shared_ptr< tom::ObservationCollectionParser > >(
                  &tom::ObservationCollection< double, TIME_TYPE >::splitObservationSets ),
                  py::arg( "observation_set_splitter" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "split_observation_sets" ).c_str( ) )
            .def( "get_single_observation_sets",
                  &tom::ObservationCollection< double, TIME_TYPE >::getSingleObservationSets,
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "get_single_observation_sets" ).c_str( ) )
            .def( "print_observation_sets_start_and_size",
                  &tom::ObservationCollection< double, TIME_TYPE >::printObservationSetsStartAndSize,
                  get_docstring( "print_observation_sets_start_index_and_size" ).c_str( ) )
            .def( "remove_single_observation_sets",
                  py::overload_cast< std::shared_ptr< tom::ObservationCollectionParser > >(
                  &tom::ObservationCollection< double, TIME_TYPE >::removeSingleObservationSets ),
                  py::arg("observation_parser"),
                  get_docstring( "remove_single_observation_sets" ).c_str( ) )
            .def("set_reference_point",
                 py::overload_cast< tss::SystemOfBodies&, const Eigen::Vector3d&, const std::string&,
                 const std::string&, const tom::LinkEndType, const std::shared_ptr< tom::ObservationCollectionParser > >( &tom::ObservationCollection< double, TIME_TYPE>::setReferencePoint ),
                 py::arg( "bodies" ),
                 py::arg( "antenna_position" ),
                 py::arg( "antenna_name" ),
                 py::arg( "spacecraft_name" ),
                 py::arg( "link_end_type" ),
                 py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                 get_docstring("set_reference_point").c_str( ) )
            .def("set_reference_points",
                  py::overload_cast< tss::SystemOfBodies&, const std::map< double, Eigen::Vector3d >&, const std::string&,
                  const tom::LinkEndType, const std::shared_ptr< tom::ObservationCollectionParser > >(
                  &tom::ObservationCollection< double, TIME_TYPE >::setReferencePoints ),
                  py::arg("bodies"),
                  py::arg("antenna_switch_history"),
                  py::arg("spacecraft_name"),
                  py::arg("link_end_type"),
                  py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "set_reference_points" ).c_str( ) )
            .def("set_reference_point",
                 py::overload_cast< tss::SystemOfBodies&, const std::shared_ptr< te::Ephemeris >, const std::string&,
                 const std::string&, const tom::LinkEndType, const std::shared_ptr< tom::ObservationCollectionParser > >(
                 &tom::ObservationCollection< double, TIME_TYPE >::setReferencePoint ),
                 py::arg("bodies"),
                 py::arg("antenna_body_fixed_ephemeris"),
                 py::arg("antenna_name"),
                 py::arg("spacecraft_name"),
                 py::arg("link_end_type"),
                 py::arg("observation_parser") = std::make_shared< tom::ObservationCollectionParser >( ),
                 get_docstring( "set_reference_points" ).c_str( ) )
            .def( "remove_empty_observation_sets",
                  &tom::ObservationCollection< double, TIME_TYPE >::removeEmptySingleObservationSets,
                  get_docstring( "remove_empty_observation_sets" ).c_str( ) )
            .def( "add_dependent_variable",
                  &tom::ObservationCollection< double, TIME_TYPE >::addDependentVariable,
                  py::arg( "dependent_variable_settings" ),
                  py::arg( "bodies" ),
                  py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                  get_docstring( "add_dependent_variable" ).c_str( ) )
            .def( "dependent_variable",
                &tom::ObservationCollection< double, TIME_TYPE >::getDependentVariables,
                py::arg( "dependent_variable_settings" ),
                py::arg( "first_compatible_settings" ) = false,
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                get_docstring( "dependent_variable" ).c_str( ) )
            .def( "concatenated_dependent_variable",
                &tom::ObservationCollection< double, TIME_TYPE >::getConcatenatedDependentVariables,
                py::arg( "dependent_variable_settings" ),
                py::arg( "first_compatible_settings" ) = false,
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                get_docstring( "concatenated_dependent_variable" ).c_str( ) )
            .def( "compatible_dependent_variable_settings",
                &tom::ObservationCollection< double, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                py::arg( "dependent_variable_settings" ),
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                get_docstring( "compatible_dependent_variable_settings" ).c_str( ) )
            .def( "compatible_dependent_variables_list",
                &tom::ObservationCollection< double, TIME_TYPE >::getAllCompatibleDependentVariables,
                py::arg( "dependent_variable_settings" ),
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                get_docstring( "compatible_dependent_variables_list" ).c_str( ) )
            .def( "dependent_variable_history_per_set",
                &tom::ObservationCollection< double, TIME_TYPE >::getDependentVariableHistoryPerObservationSet,
                py::arg( "dependent_variable_settings" ),
                py::arg( "first_compatible_settings" ) = false,
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                get_docstring( "dependent_variable_history_per_set" ).c_str( ) )
            .def( "dependent_variable_history",
                &tom::ObservationCollection< double, TIME_TYPE >::getDependentVariableHistory,
                py::arg( "dependent_variable_settings" ),
                py::arg( "first_compatible_settings" ) = false,
                py::arg( "observation_parser" ) = std::make_shared< tom::ObservationCollectionParser >( ),
                get_docstring( "dependent_variable_history" ).c_str( ) );


    m.def( "create_single_observation_set",
           py::overload_cast<
                const tom::ObservableType, const tom::LinkEnds&,
                const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >&,
                const std::vector< TIME_TYPE >, const tom::LinkEndType,
                const std::shared_ptr< tom::ObservationAncilliarySimulationSettings > >(
                &tom::createSingleObservationSet< double, TIME_TYPE > ),
           py::arg( "observable_type" ),
           py::arg( "link_ends" ),
           py::arg( "observations" ),
           py::arg( "observation_times" ),
           py::arg( "reference_link_end" ),
           py::arg( "ancillary_settings" ),
           get_docstring( "create_single_observation_set" ).c_str( ) );

    m.def("filter_observations",
          py::overload_cast<
                  const std::shared_ptr< tom::SingleObservationSet< double, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationFilterBase >, const bool >(
                  &tom::filterObservations< double, TIME_TYPE > ),
            py::arg( "original_observation_set" ),
            py::arg( "observation_filter" ),
            py::arg( "save_filtered_observations" ) = false,
            get_docstring( "filter_observations" ).c_str( ) );

    m.def("split_observation_set",
          py::overload_cast<
                  const std::shared_ptr< tom::SingleObservationSet< double, TIME_TYPE > >,
                  const std::shared_ptr< tom::ObservationSetSplitterBase > >(
                  &tom::splitObservationSet< double, TIME_TYPE > ),
          py::arg( "original_observation_set" ),
          py::arg( "observation_splitter" ),
          get_docstring( "split_observation_set" ).c_str( ) );

    py::class_< tom::SingleObservationSet<double, TIME_TYPE>,
            std::shared_ptr<tom::SingleObservationSet<double, TIME_TYPE>>>(m, "SingleObservationSet",
                                                          get_docstring("SingleObservationSet").c_str() )
            .def("set_observations", py::overload_cast< const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >&>(
                    &tom::SingleObservationSet<double, TIME_TYPE>::setObservations),
                    py::arg("observations"), get_docstring("SingleObservationSet.set_observations").c_str() )
            .def("set_observations", py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                    &tom::SingleObservationSet<double, TIME_TYPE>::setObservations),
                    py::arg("observations"), get_docstring("SingleObservationSet.set_observations").c_str() )
            .def("set_residuals", py::overload_cast< const std::vector< Eigen::Matrix< double, Eigen::Dynamic, 1 > >&>(
                    &tom::SingleObservationSet<double, TIME_TYPE>::setResiduals),
                    py::arg("residuals"), get_docstring("SingleObservationSet.set_residuals").c_str() )
            .def("set_residuals", py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >(
                    &tom::SingleObservationSet<double, TIME_TYPE>::setResiduals),
                    py::arg("residuals"), get_docstring("SingleObservationSet.set_residuals").c_str() )
            .def("set_constant_weight", py::overload_cast< const double >( &tom::SingleObservationSet<double, TIME_TYPE>::setConstantWeight ),
                    py::arg("weight"), get_docstring("SingleObservationSet.set_constant_weight").c_str() )
            .def("set_constant_weight", py::overload_cast< const Eigen::Matrix< double, Eigen::Dynamic, 1 >& >( &tom::SingleObservationSet<double, TIME_TYPE>::setConstantWeight ),
                 py::arg("weight"), get_docstring("SingleObservationSet.set_constant_weight").c_str() )
            .def("set_tabulated_weights", py::overload_cast< const Eigen::VectorXd& >( &tom::SingleObservationSet<double, TIME_TYPE>::setTabulatedWeights ),
                 py::arg("weights"), get_docstring("SingleObservationSet.set_tabulated_weights").c_str() )
            .def("filter_observations", &tom::SingleObservationSet<double, TIME_TYPE>::filterObservations,
                 py::arg("filter"), py::arg("save_filtered_obs") = true,
                 get_docstring("SingleObservationSet.filter_observations").c_str() )
            .def_property_readonly("observable_type", &tom::SingleObservationSet<double, TIME_TYPE>::getObservableType,
                                   get_docstring("SingleObservationSet.observable_type").c_str() )
            .def_property("link_definition",
                          &tom::SingleObservationSet<double, TIME_TYPE>::getLinkEnds,
                          &tom::SingleObservationSet<double, TIME_TYPE>::setLinkEnds,
                          get_docstring("SingleObservationSet.link_definition").c_str() )
            .def_property_readonly("reference_link_end", &tom::SingleObservationSet<double, TIME_TYPE>::getReferenceLinkEnd,
                                   get_docstring("SingleObservationSet.reference_link_end").c_str() )
            .def_property_readonly("number_of_observables", &tom::SingleObservationSet<double, TIME_TYPE>::getNumberOfObservables,
                                   get_docstring("SingleObservationSet.number_of_observables").c_str())
            .def_property_readonly("single_observable_size", &tom::SingleObservationSet<double, TIME_TYPE>::getSingleObservableSize,
                                   get_docstring("SingleObservationSet.single_observaable_size").c_str())
            .def_property_readonly("total_observation_set_size", &tom::SingleObservationSet<double, TIME_TYPE>::getTotalObservationSetSize,
                                   get_docstring("SingleObservationSet.total_observation_set_size").c_str())
            .def_property_readonly("time_bounds", &tom::SingleObservationSet<double, TIME_TYPE>::getTimeBounds,
                                   get_docstring("SingleObservationSet.time_bounds").c_str())
            .def_property_readonly("list_of_observations", &tom::SingleObservationSet<double, TIME_TYPE>::getObservations,
                                   get_docstring("SingleObservationSet.list_of_observations").c_str() )
            .def_property_readonly("observation_times", &tom::SingleObservationSet<double, TIME_TYPE>::getObservationTimes,
                                   get_docstring("SingleObservationSet.observation_times").c_str() )
            .def_property_readonly("concatenated_observations", &tom::SingleObservationSet<double, TIME_TYPE>::getObservationsVector,
                                   get_docstring("SingleObservationSet.concatenated_observations").c_str() )
            .def_property_readonly("computed_observations", &tom::SingleObservationSet<double, TIME_TYPE>::getComputedObservations,
                                   get_docstring("SingleObservationSet.computed_observations").c_str() )
            .def_property_readonly("concatenated_computed_observations",
                                   &tom::SingleObservationSet<double, TIME_TYPE>::getComputedObservationsVector,
                                   get_docstring("SingleObservationSet.concatenated_computed_observations").c_str() )
            .def_property_readonly("residuals", &tom::SingleObservationSet<double, TIME_TYPE>::getResiduals,
                                   get_docstring("SingleObservationSet.residuals").c_str() )
            .def_property_readonly("concatenated_residuals", &tom::SingleObservationSet<double, TIME_TYPE>::getResidualsVector,
                                   get_docstring("SingleObservationSet.concatenated_residuals").c_str() )
            .def_property_readonly("weights", &tom::SingleObservationSet<double, TIME_TYPE>::getWeights,
                                   get_docstring("SingleObservationSet.weights").c_str())
            .def_property_readonly("concatenad_weights", &tom::SingleObservationSet<double, TIME_TYPE>::getWeightsVector,
                                   get_docstring("SingleObservationSet.concatenad_weights").c_str())
            .def_property("dependent_variables", &tom::SingleObservationSet< double, TIME_TYPE >::getObservationsDependentVariables,
                                                 &tom::SingleObservationSet< double, TIME_TYPE >::setObservationsDependentVariables,
                                   get_docstring("SingleObservationSet.dependent_variables").c_str())
            .def_property_readonly("dependent_variables_history", &tom::SingleObservationSet<double, TIME_TYPE>::getDependentVariableHistory,
                                   get_docstring("SingleObservationSet.dependent_variables_history").c_str())
            .def_property_readonly("observations_history", &tom::SingleObservationSet<double, TIME_TYPE>::getObservationsHistory,
                                   get_docstring("SingleObservationSet.observations_history").c_str() )
            .def_property_readonly("ancilliary_settings", &tom::SingleObservationSet<double, TIME_TYPE>::getAncilliarySettings,
                                   get_docstring("SingleObservationSet.ancilliary_settings").c_str() )
            .def_property("weights_vector", &tom::SingleObservationSet<double, TIME_TYPE>::getWeightsVector,
                                            &tom::SingleObservationSet<double, TIME_TYPE>::setTabulatedWeights,
                                   get_docstring("SingleObservationSet.weights_vector").c_str() )
            .def_property_readonly("filtered_observation_set", &tom::SingleObservationSet<double, TIME_TYPE>::getFilteredObservationSet,
                                    get_docstring("SingleObservationSet.filtered_observation_set").c_str() )
            .def_property_readonly("number_filtered_observations", &tom::SingleObservationSet<double, TIME_TYPE>::getNumberOfFilteredObservations,
                                   get_docstring("SingleObservationSet.number_filtered_observations").c_str() )
            .def( "single_dependent_variable",
                  py::overload_cast< std::shared_ptr< tss::ObservationDependentVariableSettings >, const bool >(
                          &tom::SingleObservationSet< double, TIME_TYPE >::getSingleDependentVariable ),
                          py::arg( "dependent_variable_settings" ),
                          py::arg( "return_first_compatible_settings" ) = false,
                          get_docstring("SingleObservationSet.single_dependent_variable").c_str() )
            .def( "compatible_dependent_variable_settings",
                  &tom::SingleObservationSet< double, TIME_TYPE >::getCompatibleDependentVariablesSettingsList,
                  get_docstring( "SingleObservationSet.compatible_dependent_variable_settings" ).c_str() )
            .def( "compatible_dependent_variables_list",
                  &tom::SingleObservationSet< double, TIME_TYPE >::getAllCompatibleDependentVariables,
                  get_docstring( "SingleObservationSet.compatible_dependent_variables_list" ).c_str() )
            .def( "single_dependent_variable_history",
                  &tom::SingleObservationSet< double, TIME_TYPE >::getSingleDependentVariableHistory,
                  get_docstring( "SingleObservationSet.single_dependent_variable_history" ).c_str() )
            .def_property_readonly( "dependent_variables_matrix",
                                    &tom::SingleObservationSet< double, TIME_TYPE >::getObservationsDependentVariablesMatrix,
                                    get_docstring( "SingleObservationSet.dependent_variables_matrix" ).c_str() );


    m.def("single_observation_set",
          &tss::singleObservationSetWithoutDependentVariables< double, TIME_TYPE >,
          py::arg("observable_type"),
          py::arg("link_definition" ),
          py::arg("observations" ),
          py::arg("observation_times"),
          py::arg("reference_link_end"),
          py::arg("ancilliary_settings") = nullptr,
          get_docstring("single_observation_set").c_str() );

    /*!
     *************** STATE TRANSITION INTERFACE ***************
     */

    py::class_<
            tp::CombinedStateTransitionAndSensitivityMatrixInterface,
            std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                m, "CombinedStateTransitionAndSensitivityMatrixInterface",
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface").c_str() )
            .def("state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"),
                 py::arg("add_central_body_dependency") = true,
                 py::arg("arc_defining_bodies" ) = std::vector< std::string >( ),
                 get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_sensitivity_at_epoch").c_str() )
            .def("full_state_transition_sensitivity_at_epoch",
                 &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                 getFullCombinedStateTransitionAndSensitivityMatrix,
                 py::arg("time"),
                 py::arg("add_central_body_dependency") = true,
                 py::arg("arc_defining_bodies" ) = std::vector< std::string >( ),
                 get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.full_state_transition_sensitivity_at_epoch").c_str() )
            .def_property_readonly(
                "state_transition_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getStateTransitionMatrixSize,
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.state_transition_size").c_str() )
            .def_property_readonly(
                "sensitivity_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getSensitivityMatrixSize,
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.sensitivity_size").c_str() )
            .def_property_readonly(
                "full_parameter_size",
                &tp::CombinedStateTransitionAndSensitivityMatrixInterface::getFullParameterVectorSize,
                get_docstring("CombinedStateTransitionAndSensitivityMatrixInterface.full_parameter_size").c_str() );

    /*!
     *************** COVARIANCE ***************
     */

    m.def("propagate_covariance_rsw_split_output",
          &tp::propagateCovarianceVectorsRsw,
          py::arg("covariance_output"),
          py::arg("estimator"),
          py::arg("output_times"),
          get_docstring("propagate_covariance_rsw_split_output").c_str( ) );


    m.def("propagate_formal_errors_rsw_split_output",
          &tp::propagateFormalErrorVectorsRsw,
          py::arg("covariance_output"),
          py::arg("estimator"),
          py::arg("output_times"),
          get_docstring("propagate_formal_errors_rsw_split_output").c_str( ) );

    m.def("propagate_covariance_split_output",
          py::overload_cast<
          const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >(
              &tp::propagateCovarianceVectors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_covariance_split_output").c_str() );

    m.def("propagate_covariance",
          py::overload_cast<
          const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >(
              &tp::propagateCovariance ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_covariance").c_str() );

    m.def("propagate_formal_errors_split_output",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateFormalErrorVectors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_formal_errors").c_str() );


    m.def("propagate_formal_errors",
          py::overload_cast< const Eigen::MatrixXd,
          const std::shared_ptr< tp::CombinedStateTransitionAndSensitivityMatrixInterface >,
          const std::vector< double > >( &tp::propagateFormalErrors ),
          py::arg("initial_covariance"),
          py::arg("state_transition_interface"),
          py::arg("output_times"),
          get_docstring("propagate_formal_errors").c_str() );



    /*!
     *************** ESTIMATION ***************
     */

    py::class_<
            tss::EstimationConvergenceChecker,
            std::shared_ptr<tss::EstimationConvergenceChecker>>(m, "EstimationConvergenceChecker",
                                                                get_docstring("EstimationConvergenceChecker").c_str() );

    m.def("estimation_convergence_checker",
          &tss::estimationConvergenceChecker,
          py::arg("maximum_iterations") = 5,
          py::arg("minimum_residual_change") = 0.0,
          py::arg("minimum_residual") = 0.0,
          py::arg("number_of_iterations_without_improvement") = 2,
          get_docstring("estimation_convergence_checker").c_str() );


    py::class_<
            tss::CovarianceAnalysisInput<double, TIME_TYPE>,
            std::shared_ptr<tss::CovarianceAnalysisInput<double, TIME_TYPE>>>(m, "CovarianceAnalysisInput",
                                                                           get_docstring("CovarianceAnalysisInput").c_str() )
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection<double, TIME_TYPE> >&,
                 const Eigen::MatrixXd,
                 const Eigen::MatrixXd  >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "consider_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 get_docstring("CovarianceAnalysisInput.ctor").c_str() )
            .def( "set_constant_weight",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantWeightsMatrix,
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_weight").c_str() )
            .def( "set_weights_from_observation_collection",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setWeightsFromObservationCollection,
                  get_docstring("CovarianceAnalysisInput.set_weights_from_observation_collection").c_str() )
            .def( "set_constant_single_observable_weight",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantSingleObservableWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_weight").c_str() )
            .def( "set_constant_single_observable_vector_weight",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantSingleObservableVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_vector_weight").c_str() )
            .def( "set_constant_single_observable_and_link_end_weight",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_and_link_end_weight").c_str() )
            .def( "set_constant_single_observable_and_link_end_vector_weight",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantSingleObservableAndLinkEndsVectorWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_single_observable_and_link_end_vector_weight").c_str() )
            .def( "set_total_single_observable_and_link_end_vector_weight",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setTabulatedSingleObservableAndLinkEndsWeights,
                  py::arg( "observable_type" ),
                  py::arg( "link_ends" ),
                  py::arg( "weight_vector" ),
                  get_docstring("CovarianceAnalysisInput.set_total_single_observable_and_link_end_vector_weight").c_str() )
            .def( "set_constant_weight_per_observable",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantPerObservableWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_weight_per_observable").c_str() )
            .def( "set_constant_vector_weight_per_observable",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setConstantPerObservableVectorWeightsMatrix,
                  py::arg( "weight_per_observable" ),
                  get_docstring("CovarianceAnalysisInput.set_constant_vector_weight_per_observable").c_str() )
            .def( "define_covariance_settings",
                  &tss::CovarianceAnalysisInput<double, TIME_TYPE>::defineCovarianceSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "limit_condition_number_for_warning" ) = 1.0E8,
                  get_docstring("CovarianceAnalysisInput.define_covariance_settings").c_str() )
            .def_property("weight_matrix_diagonal",
                          &tss::CovarianceAnalysisInput<double, TIME_TYPE>::getWeightsMatrixDiagonals,
                          &tss::CovarianceAnalysisInput<double, TIME_TYPE>::setWeightsMatrixDiagonals,
                          get_docstring("CovarianceAnalysisInput.weight_matrix_diagonal").c_str() );

    py::class_<
            tss::EstimationInput<double, TIME_TYPE>,
            std::shared_ptr<tss::EstimationInput<double, TIME_TYPE>>,
            tss::CovarianceAnalysisInput<double, TIME_TYPE>>(m, "EstimationInput",
                                                          get_docstring("EstimationInput").c_str() )
            .def(py::init<
                 const std::shared_ptr< tom::ObservationCollection<double, TIME_TYPE> >&,
                 const Eigen::MatrixXd,
                 std::shared_ptr< tss::EstimationConvergenceChecker >,
                 const Eigen::MatrixXd,
                 const Eigen::VectorXd,
                 const bool >( ),
                 py::arg( "observations_and_times" ),
                 py::arg( "inverse_apriori_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "convergence_checker" ) = std::make_shared< tss::EstimationConvergenceChecker >( ),
                 py::arg( "consider_covariance" ) = Eigen::MatrixXd::Zero( 0, 0 ),
                 py::arg( "consider_parameters_deviations" ) = Eigen::VectorXd::Zero( 0 ),
                 py::arg( "apply_final_parameter_correction" ) = true,
                 get_docstring("EstimationInput.ctor").c_str() )
            .def( "define_estimation_settings",
                  &tss::EstimationInput<double, TIME_TYPE>::defineEstimationSettings,
                  py::arg( "reintegrate_equations_on_first_iteration" ) = true,
                  py::arg( "reintegrate_variational_equations" ) = true,
                  py::arg( "save_design_matrix" ) = true,
                  py::arg( "print_output_to_terminal" ) = true,
                  py::arg( "save_residuals_and_parameters_per_iteration" ) = true,
                  py::arg( "save_state_history_per_iteration" ) = false,
                  py::arg( "limit_condition_number_for_warning" ) = 1.0E8,
                  py::arg( "condition_number_warning_each_iteration" ) = true,
                  get_docstring("EstimationInput.define_estimation_settings").c_str() );

    m.attr("PodInput") = m.attr("EstimationInput");


    py::class_<
            tss::CovarianceAnalysisOutput<double, TIME_TYPE>,
            std::shared_ptr<tss::CovarianceAnalysisOutput<double, TIME_TYPE>>>(m, "CovarianceAnalysisOutput",
                                                                            get_docstring("CovarianceAnalysisOutput").c_str() )
            .def_property_readonly("inverse_covariance",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getUnnormalizedInverseCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.inverse_covariance").c_str() )
            .def_property_readonly("covariance",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getUnnormalizedCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.covariance").c_str() )
            .def_property_readonly("inverse_normalized_covariance",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getNormalizedInverseCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.inverse_normalized_covariance").c_str() )
            .def_property_readonly("normalized_covariance",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getNormalizedCovarianceMatrix,
                                   get_docstring("CovarianceAnalysisOutput.normalized_covariance").c_str() )
            .def_property_readonly("formal_errors",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getFormalErrorVector,
                                   get_docstring("CovarianceAnalysisOutput.formal_errors").c_str() )
            .def_property_readonly("correlations",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getCorrelationMatrix,
                                   get_docstring("CovarianceAnalysisOutput.correlations").c_str() )
            .def_property_readonly("design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getUnnormalizedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.design_matrix").c_str() )
            .def_property_readonly("normalized_design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getNormalizedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.normalized_design_matrix").c_str() )
            .def_property_readonly("weighted_design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getUnnormalizedWeightedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.weighted_design_matrix").c_str() )
            .def_property_readonly("weighted_normalized_design_matrix",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getNormalizedWeightedDesignMatrix,
                                   get_docstring("CovarianceAnalysisOutput.weighted_normalized_design_matrix").c_str() )
            .def_property_readonly("consider_covariance_contribution",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getConsiderCovarianceContribution,
                                   get_docstring("CovarianceAnalysisOutput.consider_covariance_contribution").c_str() )
            .def_property_readonly("normalized_covariance_with_consider_parameters",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getNormalizedCovarianceWithConsiderParameters,
                                   get_docstring("CovarianceAnalysisOutput.normalized_covariance_with_consider_parameters").c_str() )
            .def_property_readonly("unnormalized_covariance_with_consider_parameters",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getUnnormalizedCovarianceWithConsiderParameters,
                                   get_docstring("CovarianceAnalysisOutput.unnormalized_covariance_with_consider_parameters").c_str() )
            .def_property_readonly("normalized_design_matrix_consider_parameters",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getNormalizedDesignMatrixConsiderParameters,
                                   get_docstring("CovarianceAnalysisOutput.normalized_design_matrix_consider_parameters").c_str() )
            .def_property_readonly("consider_normalization_factors",
                                   &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::getConsiderNormalizationFactors,
                                   get_docstring("CovarianceAnalysisOutput.consider_normalization_factors").c_str() )
            .def_readonly("normalization_terms",
                          &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::designMatrixTransformationDiagonal_,
                          get_docstring("CovarianceAnalysisOutput.normalization_terms").c_str() );


    py::class_<
            tss::EstimationOutput<double, TIME_TYPE>,
            std::shared_ptr<tss::EstimationOutput<double, TIME_TYPE>>,
            tss::CovarianceAnalysisOutput<double, TIME_TYPE>>(m, "EstimationOutput",
                                                           get_docstring("EstimationOutput").c_str() )
            .def_property_readonly("residual_history",
                                   &tss::EstimationOutput<double, TIME_TYPE>::getResidualHistoryMatrix,
                                   get_docstring("EstimationOutput.residual_history").c_str() )
            .def_property_readonly("parameter_history",
                                   &tss::EstimationOutput<double, TIME_TYPE>::getParameterHistoryMatrix,
                                   get_docstring("EstimationOutput.parameter_history").c_str() )
            .def_property_readonly("simulation_results_per_iteration",
                                   &tss::EstimationOutput<double, TIME_TYPE>::getSimulationResults,
                                   get_docstring("EstimationOutput.simulation_results_per_iteration").c_str() )
            .def_readonly("final_residuals",
                          &tss::EstimationOutput<double, TIME_TYPE>::residuals_,
                          get_docstring("EstimationOutput.final_residuals").c_str() )
            .def_readonly("final_parameters",
                          &tss::EstimationOutput<double, TIME_TYPE>::parameterEstimate_,
                          get_docstring("EstimationOutput.final_parameters").c_str() )
            .def_readonly("best_iteration",
                          &tss::EstimationOutput<double, TIME_TYPE>::bestIteration_,
                          get_docstring("EstimationOutput.best_iteration").c_str() );

    m.attr("PodOutput") = m.attr("EstimationOutput");



}

}
}
}// namespace tudatpy
