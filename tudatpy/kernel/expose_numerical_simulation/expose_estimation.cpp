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


namespace tudat {

    namespace propagators {

        std::map<double, Eigen::MatrixXd> propagateCovarianceRsw(
            const std::shared_ptr<
                tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>>
                orbitDeterminationManager,
            const std::vector<double> evaluationTimes) {
            std::map<double, Eigen::MatrixXd> propagatedCovariance;
            tp::propagateCovariance(
                propagatedCovariance,
                estimationOutput->getUnnormalizedCovarianceMatrix(),
                orbitDeterminationManager
                    ->getStateTransitionAndSensitivityMatrixInterface(),
                evaluationTimes);

            tss::SystemOfBodies bodies = orbitDeterminationManager->getBodies();

            std::shared_ptr<tep::EstimatableParameterSet<STATE_SCALAR_TYPE>>
                parameterSet =
                    orbitDeterminationManager->getParametersToEstimate();

            std::map<int,
                     std::shared_ptr<tep::EstimatableParameter<
                         Eigen::Matrix<STATE_SCALAR_TYPE, Eigen::Dynamic, 1>>>>
                initialStates = parameterSet->getInitialStateParameters();
            std::map<std::pair<std::string, std::string>, std::vector<int>>
                transformationList;
            for(auto it : initialStates) {
                if(std::dynamic_pointer_cast<
                       tep::InitialTranslationalStateParameter<
                           STATE_SCALAR_TYPE>>(it.second)) {
                    std::shared_ptr<tep::InitialTranslationalStateParameter<
                        STATE_SCALAR_TYPE>>
                        currentInitialState = std::dynamic_pointer_cast<
                            tep::InitialTranslationalStateParameter<
                                STATE_SCALAR_TYPE>>(it.second);
                    transformationList
                        [std::make_pair(currentInitialState->getParameterName()
                                            .second.first,
                                        currentInitialState->getCentralBody())]
                            .push_back(it.first);

                } else if(std::dynamic_pointer_cast<
                              tep::ArcWiseInitialTranslationalStateParameter<
                                  STATE_SCALAR_TYPE>>(it.second)) {
                    throw std::runtime_error(
                        "Error, multi-arc not yet supported in automatic "
                        "covariance conversion");
                }
            }

            Eigen::Matrix3d currentInertialToRswPosition;
            Eigen::Matrix6d currentInertialToRswState;
            Eigen::MatrixXd currentFullInertialToRswState =
                Eigen::MatrixXd::Zero(6, 6);

            std::map<double, Eigen::MatrixXd> propagatedRswCovariance;
            for(auto it : propagatedCovariance) {
                double currentTime = static_cast<double>(it.first);
                Eigen::MatrixXd currentCovariance = it.second;
                currentFullInertialToRswState.setZero();

                for(auto it_body : transformationList) {
                    Eigen::Vector6d relativeState =
                        bodies.getBody(it_body.first.first)
                            ->getStateInBaseFrameFromEphemeris(currentTime) -
                        bodies.getBody(it_body.first.second)
                            ->getStateInBaseFrameFromEphemeris(currentTime);
                    currentInertialToRswPosition = trf::
                        getInertialToRswSatelliteCenteredFrameRotationMatrix(
                            relativeState);
                    currentInertialToRswState.block(0, 0, 3, 3) =
                        currentInertialToRswPosition;
                    currentInertialToRswState.block(3, 3, 3, 3) =
                        currentInertialToRswPosition;
                    for(unsigned int j = 0; j < it_body.second.size(); j++) {
                        int currentStartIndex = it_body.second.at(j);
                        currentFullInertialToRswState.block(
                            currentStartIndex, currentStartIndex, 6, 6) =
                            currentInertialToRswState;
                    }
                }
                propagatedRswCovariance[currentTime] =
                    currentFullInertialToRswState * currentCovariance *
                    currentFullInertialToRswState.transpose();
            }
            return propagatedRswCovariance;
        }


        std::pair<std::vector<double>, std::vector<Eigen::MatrixXd>>
        propagateCovarianceVectorsRsw(
            const std::shared_ptr<
                tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>>
                orbitDeterminationManager,
            const std::vector<double> evaluationTimes) {
            std::map<double, Eigen::MatrixXd> propagatedRswCovariance =
                propagateCovarianceRsw(estimationOutput,
                                       orbitDeterminationManager,
                                       evaluationTimes);

            return std::make_pair(
                utilities::createVectorFromMapKeys(propagatedRswCovariance),
                utilities::createVectorFromMapValues(propagatedRswCovariance));
        }

        std::map<double, Eigen::VectorXd> propagateFormalErrorsRsw(
            const std::shared_ptr<
                tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>>
                orbitDeterminationManager,
            const std::vector<double> evaluationTimes) {
            std::map<double, Eigen::MatrixXd> propagatedCovariance;
            std::map<double, Eigen::VectorXd> propagatedFormalErrors;

            propagatedCovariance = propagateCovarianceRsw(
                estimationOutput, orbitDeterminationManager, evaluationTimes);
            tp::convertCovarianceHistoryToFormalErrorHistory(
                propagatedFormalErrors, propagatedCovariance);

            return propagatedFormalErrors;
        }

        std::pair<std::vector<double>, std::vector<Eigen::VectorXd>>
        propagateFormalErrorVectorsRsw(
            const std::shared_ptr<
                tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<STATE_SCALAR_TYPE, TIME_TYPE>>
                orbitDeterminationManager,
            const std::vector<double> evaluationTimes) {
            std::map<double, Eigen::VectorXd> propagatedFormalErrors =
                propagateFormalErrorsRsw(estimationOutput,
                                         orbitDeterminationManager,
                                         evaluationTimes);
            tp::propagateFormalErrorsRsw(
                estimationOutput, orbitDeterminationManager, evaluationTimes);
            return std::make_pair(
                utilities::createVectorFromMapKeys(propagatedFormalErrors),
                utilities::createVectorFromMapValues(propagatedFormalErrors));
        }


        std::pair<std::vector<double>, std::vector<Eigen::MatrixXd>>
        propagateCovarianceVectors(
            const Eigen::MatrixXd initialCovariance,
            const std::shared_ptr<
                tp::CombinedStateTransitionAndSensitivityMatrixInterface>
                stateTransitionInterface,
            const std::vector<double> evaluationTimes) {
            std::map<double, Eigen::MatrixXd> propagatedCovariance;
            tp::propagateCovariance(propagatedCovariance, initialCovariance,
                                    stateTransitionInterface, evaluationTimes);
            return std::make_pair(
                utilities::createVectorFromMapKeys(propagatedCovariance),
                utilities::createVectorFromMapValues(propagatedCovariance));
        }

        std::pair<std::vector<double>, std::vector<Eigen::VectorXd>>
        propagateFormalErrorVectors(
            const Eigen::MatrixXd initialCovariance,
            const std::shared_ptr<
                tp::CombinedStateTransitionAndSensitivityMatrixInterface>
                stateTransitionInterface,
            const std::vector<double> evaluationTimes) {
            std::map<double, Eigen::VectorXd> propagatedFormalErrors;
            tp::propagateFormalErrors(propagatedFormalErrors, initialCovariance,
                                      stateTransitionInterface,
                                      evaluationTimes);
            return std::make_pair(
                utilities::createVectorFromMapKeys(propagatedFormalErrors),
                utilities::createVectorFromMapValues(propagatedFormalErrors));
        }

    }  // namespace propagators

    namespace simulation_setup {

        std::pair<std::vector<double>, std::vector<Eigen::VectorXd>>
        getTargetAnglesAndRangeVector(
            const simulation_setup::SystemOfBodies& bodies,
            const std::pair<std::string, std::string> groundStationId,
            const std::string& targetBody, const std::vector<double> times,
            const bool transmittingToTarget) {
            std::map<double, Eigen::VectorXd> targetAnglesAndRange =
                getTargetAnglesAndRange(bodies, groundStationId, targetBody,
                                        times, transmittingToTarget);
            return std::make_pair(
                utilities::createVectorFromMapKeys(targetAnglesAndRange),
                utilities::createVectorFromMapValues(targetAnglesAndRange));
        }

        template <typename ObservationScalarType = double,
                  typename TimeType = double>
        std::shared_ptr<
            tom::SingleObservationSet<ObservationScalarType, TimeType>>
        singleObservationSetWithoutDependentVariables(
            const tom::ObservableType observableType,
            const tom::LinkDefinition& linkEnds,
            const std::vector<Eigen::Matrix<ObservationScalarType,
                                            Eigen::Dynamic, 1>>& observations,
            const std::vector<TimeType> observationTimes,
            const tom::LinkEndType referenceLinkEnd,
            const std::shared_ptr<
                observation_models::ObservationAncilliarySimulationSettings>
                ancilliarySettings = nullptr) {
            return std::make_shared<
                tom::SingleObservationSet<ObservationScalarType, TimeType>>(
                observableType, linkEnds, observations, observationTimes,
                referenceLinkEnd, std::vector<Eigen::VectorXd>(), nullptr,
                ancilliarySettings);
        }


    }  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace estimation {


            void expose_estimation(py::module& m) {
                /*!
                 *************** PARAMETERS ***************
                 */


                py::class_<tep::EstimatableParameterSet<STATE_SCALAR_TYPE>,
                           std::shared_ptr<tep::EstimatableParameterSet<
                               STATE_SCALAR_TYPE>>>(m,
                                                    "EstimatableParameterSet",
                                                    R"doc(

        Class containing a consolidated set of estimatable parameters.

        Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
        The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.create_parameters_to_estimate` function,
        whose output is indeed an *EstimatableParameterSet* object.





     )doc")
                    .def_property_readonly(
                        "parameter_set_size",
                        &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE>::getEstimatedParameterSetSize,
                        R"doc(

        **read-only**

        Size of the parameter set, i.e. amount of estimatable parameters contained in the set.

        :type: int
     )doc")
                    .def_property_readonly(
                        "initial_states_size",
                        &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::
                            getInitialDynamicalStateParameterSize,
                        R"doc(

        **read-only**

        Amount of initial state parameters contained in the set.

        :type: int
     )doc")
                    .def_property_readonly(
                        "initial_single_arc_states_size",
                        &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::
                            getInitialDynamicalSingleArcStateParameterSize,
                        R"doc(

        **read-only**

        Amount of initial state parameters in the set, which are treated in a single-arc fashion.

        :type: int
     )doc")
                    .def_property_readonly(
                        "initial_multi_arc_states_size",
                        &tep::EstimatableParameterSet<STATE_SCALAR_TYPE>::
                            getInitialDynamicalMultiArcStateParameterSize,
                        R"doc(

        **read-only**

        Amount of initial state parameters in the set, which are treated in a multi-arc fashion.

        :type: int
     )doc")
                    .def_property_readonly(
                        "constraints_size",
                        &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE>::getConstraintSize,
                        R"doc(

        **read-only**

        Total size of linear constraint that is to be applied during estimation.

        :type: int
     )doc")
                    .def_property(
                        "parameter_vector",
                        &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE>::getFullParameterValues<double>,
                        &tep::EstimatableParameterSet<
                            STATE_SCALAR_TYPE>::resetParameterValues<double>,
                        R"doc(

        Vector containing the parameter values of all parameters in the set.

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc")
                    .def("indices_for_parameter_type",
                         &tep::EstimatableParameterSet<
                             STATE_SCALAR_TYPE>::getIndicesForParameterType,
                         py::arg("parameter_type"),
                         R"doc(

        Function to retrieve the indices of a given type of parameter.

        Function to retrieve the index of all parameters of a given type from the parameter set.
        This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.


        Parameters
        ----------
        parameter_type : Tuple[ :class:`~tudatpy.numerical_simulation.estimation_setup.parameter.EstimatableParameterTypes`, Tuple[str, str] ]
            help
        Returns
        -------
        List[ Tuple[int, int] ]
            help





    )doc");

                /*!
                 *************** OBSERVATIONS ***************
                 */

                py::enum_<tom::ObservationFilterType>(
                    m, "ObservationFilterType",
                    R"doc(No documentation found.)doc")
                    .value("residual_filtering",
                           tom::ObservationFilterType::residual_filtering)
                    .value("absolute_value_filtering",
                           tom::ObservationFilterType::absolute_value_filtering)
                    .value("epochs_filtering",
                           tom::ObservationFilterType::epochs_filtering)
                    .value("time_bounds_filtering",
                           tom::ObservationFilterType::time_bounds_filtering)
                    .value("dependent_variable_filtering",
                           tom::ObservationFilterType::
                               dependent_variable_filtering)
                    .export_values();

                py::enum_<tom::ObservationSetSplitterType>(
                    m, "ObservationSetSplitterType",
                    R"doc(No documentation found.)doc")
                    .value("time_tags_splitter",
                           tom::ObservationSetSplitterType::time_tags_splitter)
                    .value(
                        "time_interval_splitter",
                        tom::ObservationSetSplitterType::time_interval_splitter)
                    .value("time_span_splitter",
                           tom::ObservationSetSplitterType::time_span_splitter)
                    .value("nb_observations_splitter",
                           tom::ObservationSetSplitterType::
                               nb_observations_splitter)
                    .export_values();

                py::enum_<tom::ObservationParserType>(
                    m, "ObservationParserType",
                    R"doc(No documentation found.)doc")
                    .value("empty_parser",
                           tom::ObservationParserType::empty_parser)
                    .value("observable_type_parser",
                           tom::ObservationParserType::observable_type_parser)
                    .value("link_ends_parser",
                           tom::ObservationParserType::link_ends_parser)
                    .value("link_end_str_parser",
                           tom::ObservationParserType::link_end_string_parser)
                    .value("link_end_id_parser",
                           tom::ObservationParserType::link_end_id_parser)
                    .value("link_end_type_parser",
                           tom::ObservationParserType::link_end_type_parser)
                    .value("single_link_end_parser",
                           tom::ObservationParserType::single_link_end_parser)
                    .value("time_bounds_parser",
                           tom::ObservationParserType::time_bounds_parser)
                    .value(
                        "ancillary_settings_parser",
                        tom::ObservationParserType::ancillary_settings_parser)
                    .value("multi_type_parser",
                           tom::ObservationParserType::multi_type_parser)
                    .export_values();

                py::class_<tom::ObservationFilterBase,
                           std::shared_ptr<tom::ObservationFilterBase>>(
                    m, "ObservationFilterBase",
                    R"doc(No documentation found.)doc");

                m.def("observation_filter",
                      py::overload_cast<tom::ObservationFilterType,
                                        const double, const bool, const bool>(
                          &tom::observationFilter),
                      py::arg("filter_type"), py::arg("filter_value"),
                      py::arg("filter_out") = true,
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_filter",
                      py::overload_cast<tom::ObservationFilterType,
                                        const std::vector<double>, const bool,
                                        const bool>(&tom::observationFilter),
                      py::arg("filter_type"), py::arg("filter_value"),
                      py::arg("filter_out") = true,
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

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
                    py::overload_cast<tom::ObservationFilterType, const double,
                                      const double, const bool, const bool>(
                        &tom::observationFilter),
                    py::arg("filter_type"), py::arg("first_filter_value"),
                    py::arg("second_filter_value"),
                    py::arg("filter_out") = true,
                    py::arg("use_opposite_condition") = false,
                    R"doc(No documentation found.)doc");

                m.def("observation_filter",
                      py::overload_cast<tom::ObservationFilterType,
                                        const Eigen::VectorXd&, const bool,
                                        const bool>(&tom::observationFilter),
                      py::arg("filter_type"), py::arg("filter_value"),
                      py::arg("filter_out") = true,
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_filter",
                      py::overload_cast<
                          std::shared_ptr<
                              tss::ObservationDependentVariableSettings>,
                          const Eigen::VectorXd&, const bool, const bool>(
                          &tom::observationFilter),
                      py::arg("dependent_variable_settings"),
                      py::arg("filter_value"), py::arg("filter_out") = true,
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                py::class_<tom::ObservationSetSplitterBase,
                           std::shared_ptr<tom::ObservationSetSplitterBase>>(
                    m, "ObservationSetSplitterBase",
                    R"doc(No documentation found.)doc");

                m.def("observation_set_splitter",
                      py::overload_cast<tom::ObservationSetSplitterType,
                                        const std::vector<double>, const int>(
                          &tom::observationSetSplitter),
                      py::arg("splitter_type"), py::arg("splitter_value"),
                      py::arg("min_number_observations") = 0,
                      R"doc(No documentation found.)doc");

                m.def("observation_set_splitter",
                      py::overload_cast<tom::ObservationSetSplitterType,
                                        const double, const int>(
                          &tom::observationSetSplitter),
                      py::arg("splitter_type"), py::arg("splitter_value"),
                      py::arg("min_number_observations") = 0,
                      R"doc(No documentation found.)doc");

                m.def("observation_set_splitter",
                      py::overload_cast<tom::ObservationSetSplitterType,
                                        const int, const int>(
                          &tom::observationSetSplitter),
                      py::arg("splitter_type"), py::arg("splitter_value"),
                      py::arg("min_number_observations") = 0,
                      R"doc(No documentation found.)doc");

                py::class_<tom::ObservationCollectionParser,
                           std::shared_ptr<tom::ObservationCollectionParser>>(
                    m, "ObservationCollectionParser",
                    R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<>(&tom::observationParser),
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const tom::ObservableType, const bool>(
                          &tom::observationParser),
                      py::arg("observable_type"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::vector<tom::ObservableType>&,
                                        const bool>(&tom::observationParser),
                      py::arg("observable_type_vector"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const tom::LinkEnds, const bool>(
                          &tom::observationParser),
                      py::arg("link_ends"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::vector<tom::LinkEnds>&,
                                        const bool>(&tom::observationParser),
                      py::arg("link_ends_vector"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::string, const bool,
                                        const bool>(&tom::observationParser),
                      py::arg("link_ends_str"),
                      py::arg("is_reference_point") = false,
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::vector<std::string>&,
                                        const bool, const bool>(
                          &tom::observationParser),
                      py::arg("link_ends_str_vector"),
                      py::arg("is_reference_point") = false,
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def(
                    "observation_parser",
                    py::overload_cast<
                        const std::pair<std::string, std::string>&, const bool>(
                        &tom::observationParser),
                    py::arg("link_end_id"),
                    py::arg("use_opposite_condition") = false,
                    R"doc(No documentation found.)doc");

                m.def(
                    "observation_parser",
                    py::overload_cast<
                        const std::vector<std::pair<std::string, std::string>>&,
                        const bool>(&tom::observationParser),
                    py::arg("link_end_ids_vector"),
                    py::arg("use_opposite_condition") = false,
                    R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const tom::LinkEndType&, const bool>(
                          &tom::observationParser),
                      py::arg("link_end_type"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::vector<tom::LinkEndType>&,
                                        const bool>(&tom::observationParser),
                      py::arg("link_end_types_vector"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<
                          const std::pair<tom::LinkEndType, tom::LinkEndId>&,
                          const bool>(&tom::observationParser),
                      py::arg("single_link_end"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::vector<std::pair<
                                            tom::LinkEndType, tom::LinkEndId>>&,
                                        const bool>(&tom::observationParser),
                      py::arg("single_link_ends_vector"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::pair<double, double>&,
                                        const bool>(&tom::observationParser),
                      py::arg("time_bounds"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<
                          const std::vector<std::pair<double, double>>&,
                          const bool>(&tom::observationParser),
                      py::arg("time_bounds_vector"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<
                          const std::shared_ptr<
                              tom::ObservationAncilliarySimulationSettings>,
                          const bool>(&tom::observationParser),
                      py::arg("ancillary_settings"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<
                          const std::vector<std::shared_ptr<
                              tom::ObservationAncilliarySimulationSettings>>&,
                          const bool>(&tom::observationParser),
                      py::arg("ancillary_settings_vector"),
                      py::arg("use_opposite_condition") = false,
                      R"doc(No documentation found.)doc");

                m.def("observation_parser",
                      py::overload_cast<const std::vector<std::shared_ptr<
                                            tom::ObservationCollectionParser>>&,
                                        const bool>(&tom::observationParser),
                      py::arg("observation_parsers"),
                      py::arg("combine_conditions") = false,
                      R"doc(No documentation found.)doc");


                py::class_<
                    tom::ObservationViabilityCalculator,
                    std::shared_ptr<tom::ObservationViabilityCalculator>>(
                    m, "ObservationViabilityCalculator",
                    R"doc(

        Template class for observation viability calculators.

        Template class for classes which conducts viability calculations on simulated observations.
        Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
        The user typically does not interact directly with this class.





     )doc")
                    .def("is_observation_viable",
                         &tom::ObservationViabilityCalculator::
                             isObservationViable,
                         py::arg("link_end_states"), py::arg("link_end_times"),
                         R"doc(

        Function to check whether an observation is viable.

        Function to check whether an observation is viable.
        The calculation is performed based on the given times and link end states.
        Note, that this function is called automatically during the simulation of observations.
        Direct calls to this function are generally not required.


        Parameters
        ----------
        link_end_states : List[ numpy.ndarray[numpy.float64[6, 1]] ]
            Vector of states of the link ends involved in the observation.
        link_end_times : List[float]
            Vector of times at the link ends involved in the observation.
        Returns
        -------
        bool
            True if observation is viable, false if not.





    )doc");

                py::class_<
                    tom::ObservationSimulatorBase<STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::ObservationSimulatorBase<
                        STATE_SCALAR_TYPE, TIME_TYPE>>>(m,
                                                        "ObservationSimulator",
                                                        R"doc(

        Class hosting the functionality for simulating observations.

        Class hosting the functionality for simulating a given observable over a defined link geometry.
        Instances of this class are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` objects upon instantiation of the :class:`~tudatpy.numerical_simulation.Estimator` class.





     )doc");

                py::class_<
                    tom::ObservationSimulator<1, STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::ObservationSimulator<
                        1, STATE_SCALAR_TYPE, TIME_TYPE>>,
                    tom::ObservationSimulatorBase<STATE_SCALAR_TYPE,
                                                  TIME_TYPE>>(
                    m, "ObservationSimulator_1",
                    R"doc(No documentation found.)doc");

                py::class_<
                    tom::ObservationSimulator<2, STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::ObservationSimulator<
                        2, STATE_SCALAR_TYPE, TIME_TYPE>>,
                    tom::ObservationSimulatorBase<STATE_SCALAR_TYPE,
                                                  TIME_TYPE>>(
                    m, "ObservationSimulator_2",
                    R"doc(No documentation found.)doc");

                py::class_<
                    tom::ObservationSimulator<3, STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::ObservationSimulator<
                        3, STATE_SCALAR_TYPE, TIME_TYPE>>,
                    tom::ObservationSimulatorBase<STATE_SCALAR_TYPE,
                                                  TIME_TYPE>>(
                    m, "ObservationSimulator_3",
                    R"doc(No documentation found.)doc");

                py::class_<
                    tom::ObservationSimulator<6, STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::ObservationSimulator<
                        6, STATE_SCALAR_TYPE, TIME_TYPE>>,
                    tom::ObservationSimulatorBase<STATE_SCALAR_TYPE,
                                                  TIME_TYPE>>(
                    m, "ObservationSimulator_6",
                    R"doc(No documentation found.)doc");

                m.def("simulate_observations",
                      &tss::simulateObservations<STATE_SCALAR_TYPE, TIME_TYPE>,
                      py::arg("simulation_settings"),
                      py::arg("observation_simulators"), py::arg("bodies"),
                      R"doc(

Function to simulate observations.

Function to simulate observations from set observation simulators and observation simulator settings.
Automatically iterates over all provided observation simulators, generating the full set of simulated observations.


Parameters
----------
observation_to_simulate : List[ :class:`ObservationSimulationSettings` ]
    List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.

observation_simulators : List[ :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` ]
    List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.

bodies : :class:`~tudatpy.numerical_simulation.environment.SystemOfBodies`
    Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.

Returns
-------
:class:`~tudatpy.numerical_simulation.estimation.ObservationCollection`
    Object collecting all products of the observation simulation.






    )doc");

                m.def("compute_residuals_and_dependent_variables",
                      &tss::computeResidualsAndDependentVariables<
                          STATE_SCALAR_TYPE, TIME_TYPE>,
                      py::arg("observation_collection"),
                      py::arg("observation_simulators"), py::arg("bodies"),
                      R"doc(No documentation found.)doc");


                m.def("create_pseudo_observations_and_models",
                      &tss::simulatePseudoObservations<TIME_TYPE,
                                                       STATE_SCALAR_TYPE>,
                      py::arg("bodies"), py::arg("observed_bodies"),
                      py::arg("central_bodies"), py::arg("initial_time"),
                      py::arg("final_time"), py::arg("time_step"),
                      R"doc(No documentation found.)doc");

                m.def("create_best_fit_to_ephemeris",
                      &tss::createBestFitToCurrentEphemeris<TIME_TYPE,
                                                            STATE_SCALAR_TYPE>,
                      py::arg("bodies"), py::arg("acceleration_models"),
                      py::arg("observed_bodies"), py::arg("central_bodies"),
                      py::arg("integrator_settings"), py::arg("initial_time"),
                      py::arg("final_time"), py::arg("data_point_interval"),
                      py::arg("additional_parameter_names") = std::vector<
                          std::shared_ptr<tep::EstimatableParameterSettings>>(),
                      py::arg("number_of_iterations") = 3,
                      py::arg("reintegrate_variational_equations") = true,
                      py::arg("results_print_frequency") = 0.0,
                      R"doc(No documentation found.)doc");

                m.def(
                    "set_existing_observations",
                    &tss::setExistingObservations<STATE_SCALAR_TYPE, TIME_TYPE>,
                    py::arg("observations"), py::arg("reference_link_end"),
                    py::arg("ancilliary_settings_per_observatble") = std::map<
                        tom::ObservableType,
                        std::shared_ptr<
                            tom::ObservationAncilliarySimulationSettings>>());

                m.def("compute_target_angles_and_range",
                      &tss::getTargetAnglesAndRange, py::arg("bodies"),
                      py::arg("station_id"), py::arg("target_body"),
                      py::arg("observation_times"),
                      py::arg("is_station_transmitting"),
                      R"doc(

Function to compute the azimuth angle, elevation angle and range at a ground station.

Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
as observed at that station


Parameters
----------
bodies : SystemOfBodies
    System of bodies that defines the full physical environment

station_id : tuple[ str, str]
    Identifier for the observing station, as a pair of strings: the body name and the station name.

target_body : str
    Name of body which is observed by ground station

observation_times : list[float]
    List of times at which the ground station observations are to be analyzed

is_station_transmitting : Bool
    Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
    This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

Returns
-------
dict[float,numpy.ndarray[numpy.float64[3, 1]]]
    Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






    )doc");

                m.def("compute_target_angles_and_range_vectors",
                      &tss::getTargetAnglesAndRangeVector, py::arg("bodies"),
                      py::arg("station_id"), py::arg("target_body"),
                      py::arg("observation_times"),
                      py::arg("is_station_transmitting"),
                      R"doc(

Function to compute the azimuth angle, elevation angle and range at a ground station.

Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
as observed at that station


Parameters
----------
bodies : SystemOfBodies
    System of bodies that defines the full physical environment

station_id : tuple[ str, str]
    Identifier for the observing station, as a pair of strings: the body name and the station name.

target_body : str
    Name of body which is observed by ground station

observation_times : list[float]
    List of times at which the ground station observations are to be analyzed

is_station_transmitting : Bool
    Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station.
    This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target

Returns
-------
dict[float,numpy.ndarray[numpy.float64[3, 1]]]
    Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range






    )doc");

                m.def(
                    "create_filtered_observation_collection",
                    py::overload_cast<
                        const std::shared_ptr<tom::ObservationCollection<
                            STATE_SCALAR_TYPE, TIME_TYPE>>,
                        const std::map<
                            std::shared_ptr<tom::ObservationCollectionParser>,
                            std::shared_ptr<tom::ObservationFilterBase>>&>(
                        &tom::filterObservations<STATE_SCALAR_TYPE, TIME_TYPE>),
                    py::arg("original_observation_collection"),
                    py::arg("observation_filters_map"),
                    R"doc(No documentation found.)doc");

                m.def(
                    "create_filtered_observation_collection",
                    py::overload_cast<
                        const std::shared_ptr<tom::ObservationCollection<
                            STATE_SCALAR_TYPE, TIME_TYPE>>,
                        const std::shared_ptr<tom::ObservationFilterBase>,
                        const std::shared_ptr<
                            tom::ObservationCollectionParser>>(
                        &tom::filterObservations<STATE_SCALAR_TYPE, TIME_TYPE>),
                    py::arg("original_observation_collection"),
                    py::arg("observation_filter"),
                    py::arg("observation_parser") =
                        std::make_shared<tom::ObservationCollectionParser>(),
                    R"doc(No documentation found.)doc");

                m.def("split_observation_collection",
                      &tom::splitObservationSets<STATE_SCALAR_TYPE, TIME_TYPE>,
                      py::arg("original_observation_collection"),
                      py::arg("observation_set_splitter"),
                      py::arg("observation_parser") =
                          std::make_shared<tom::ObservationCollectionParser>(),
                      R"doc(No documentation found.)doc");

                m.def("create_new_observation_collection",
                      &tom::createNewObservationCollection<STATE_SCALAR_TYPE,
                                                           TIME_TYPE>,
                      py::arg("original_observation_collection"),
                      py::arg("observation_parser") =
                          std::make_shared<tom::ObservationCollectionParser>(),
                      R"doc(No documentation found.)doc");


                py::class_<
                    tom::ObservationCollection<STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::ObservationCollection<
                        STATE_SCALAR_TYPE, TIME_TYPE>>>(m,
                                                        "ObservationCollection",
                                                        R"doc(

        Class collecting all observations and associated data for use in an estimation.

        Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
        this class is instantiated via a call to the :func:`~tudatpy.numerical_simulation.estimation.simulate_observations` function. More information is provided
        on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_





     )doc")
                    .def(py::init<std::vector<
                             std::shared_ptr<tom::SingleObservationSet<
                                 STATE_SCALAR_TYPE, TIME_TYPE>>>>(),
                         py::arg("observation_sets"))
                    .def_property_readonly(
                        "concatenated_times",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getConcatenatedTimeVector,
                        R"doc(

        **read-only**

        Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc")
                    .def_property_readonly(
                        "concatenated_float_times",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getConcatenatedDoubleTimeVector,
                        R"doc(

        **read-only**

        Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc")
                    .def_property_readonly(
                        "concatenated_weights",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getUnparsedConcatenatedWeights,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "concatenated_observations",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getObservationVector,
                        R"doc(

        **read-only**

        Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc")
                    .def_property_readonly(
                        "concatenated_link_definition_ids",
                        py::overload_cast<>(
                            &tom::ObservationCollection<
                                STATE_SCALAR_TYPE,
                                TIME_TYPE>::getConcatenatedLinkEndIds),
                        R"doc(

        **read-only**

        Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order of the present vector.

        :type: numpy.ndarray[ int ]
     )doc")
                    .def_property_readonly(
                        "link_definition_ids",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getInverseLinkEndIdentifierMap,
                        R"doc(

        **read-only**

        Dictionaty mapping a link end integer identifier to the specific link ends

        :type: dict[ int, dict[ LinkEndType, LinkEndId ] ]
     )doc")
                    .def_property_readonly(
                        "observable_type_start_index_and_size",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getObservationTypeStartAndSize,
                        R"doc(

        **read-only**

        Dictionary defining per obervable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type

        :type: dict[ ObservableType, [ int, int ] ]
     )doc")
                    .def_property_readonly(
                        "observation_set_start_index_and_size",
                        &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                    TIME_TYPE>::
                            getObservationSetStartAndSizePerLinkEndIndex,
                        R"doc(

        **read-only**

        The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the correspoding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per obervable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.

        :type: dict[ ObservableType, dict[ int, list[ int, int ] ] ]
     )doc")
                    .def_property_readonly(
                        "observation_vector_size",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getTotalObservableSize,
                        R"doc(

        **read-only**

        Length of the total vector of observations

        :type: int
     )doc")
                    .def_property_readonly(
                        "sorted_observation_sets",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSortedObservationSets,
                        R"doc(

        **read-only**

        The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_

        :type: dict[ ObservableType, dict[ int, list[ SingleObservationSet ] ] ]
     )doc")
                    .def_property_readonly(
                        "link_ends_per_observable_type",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getLinkEndsPerObservableType,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "link_definitions_per_observable",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getLinkDefinitionsPerObservable,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "time_bounds",
                        &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                    TIME_TYPE>::getTimeBounds,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "sorted_per_set_time_bounds",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSortedObservationSetsTimeBounds,
                        R"doc(No documentation found.)doc")
                    .def(
                        "set_observations",
                        py::overload_cast<const Eigen::Matrix<
                            STATE_SCALAR_TYPE, Eigen::Dynamic, 1>&>(
                            &tom::ObservationCollection<
                                STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                        py::arg("observations"),
                        R"doc(No documentation found.)doc")
                    .def(
                        "set_observations",
                        py::overload_cast<
                            const Eigen::Matrix<STATE_SCALAR_TYPE,
                                                Eigen::Dynamic, 1>&,
                            const std::shared_ptr<
                                tom::ObservationCollectionParser>>(
                            &tom::ObservationCollection<
                                STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                        py::arg("observations"), py::arg("observation_parser"),
                        R"doc(No documentation found.)doc")
                    .def(
                        "set_observations",
                        py::overload_cast<const std::map<
                            std::shared_ptr<tom::ObservationCollectionParser>,
                            Eigen::Matrix<STATE_SCALAR_TYPE, Eigen::Dynamic,
                                          1>>&>(
                            &tom::ObservationCollection<
                                STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                        py::arg("observations_per_parser"),
                        R"doc(No documentation found.)doc")
                    .def("set_residuals",
                         py::overload_cast<const Eigen::Matrix<
                             STATE_SCALAR_TYPE, Eigen::Dynamic, 1>&>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                         py::arg("residuals"),
                         R"doc(No documentation found.)doc")
                    .def("set_residuals",
                         py::overload_cast<
                             const Eigen::Matrix<STATE_SCALAR_TYPE,
                                                 Eigen::Dynamic, 1>&,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                         py::arg("residuals"), py::arg("observation_parser"),
                         R"doc(No documentation found.)doc")
                    .def("set_residuals",
                         py::overload_cast<const std::map<
                             std::shared_ptr<tom::ObservationCollectionParser>,
                             Eigen::Matrix<STATE_SCALAR_TYPE, Eigen::Dynamic,
                                           1>>&>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                         py::arg("residuals_per_parser"),
                         R"doc(No documentation found.)doc")
                    .def("get_link_definitions_for_observables",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getLinkDefinitionsForSingleObservable,
                         py::arg("observable_type"),
                         R"doc(No documentation found.)doc")
                    .def("get_single_link_and_type_observations",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getSingleLinkAndTypeObservationSets,
                         py::arg("observable_type"), py::arg("link_definition"),
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





    )doc")
                    .def("get_observable_types",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::getObservableTypes,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_bodies_in_link_ends",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::getBodiesInLinkEnds,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_reference_points_in_link_ends",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getReferencePointsInLinkEnds,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_time_bounds_list",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::getTimeBoundsList,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_time_bounds_per_set",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::getTimeBoundsPerSet,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def(
                        "get_observations",
                        &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                    TIME_TYPE>::getObservations,
                        py::arg("observation_parser") = std::make_shared<
                            tom::ObservationCollectionParser>(),
                        R"doc(No documentation found.)doc")
                    .def("get_concatenated_observations",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedObservations,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_observation_times",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::getObservationTimes,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_observation_times",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedObservationTimes,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_float_observation_times",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedDoubleObservationTimes,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_observations_and_times",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getObservationsAndTimes,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_observations_and_times",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedObservationsAndTimes,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_link_definition_ids",
                         py::overload_cast<
                             std::shared_ptr<tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::getConcatenatedLinkEndIds),
                         py::arg("observation_parser"),
                         R"doc(No documentation found.)doc")
                    .def("get_weights",
                         &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                     TIME_TYPE>::getWeights,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_weights",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedWeights,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_residuals",
                         &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                     TIME_TYPE>::getResiduals,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_residuals",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedResiduals,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def(
                        "get_rms_residuals",
                        &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                    TIME_TYPE>::getRmsResiduals,
                        py::arg("observation_parser") = std::make_shared<
                            tom::ObservationCollectionParser>(),
                        R"doc(No documentation found.)doc")
                    .def("get_mean_residuals",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::getMeanResiduals,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_computed_observations",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getComputedObservations,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_concatenated_computed_observations",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedComputedObservations,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_weight",
                         py::overload_cast<
                             const double,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setConstantWeight),
                         py::arg("weight"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_weight",
                         py::overload_cast<
                             const Eigen::VectorXd,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setConstantWeight),
                         py::arg("weight"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_weight_per_observation_parser",
                         py::overload_cast<std::map<
                             std::shared_ptr<tom::ObservationCollectionParser>,
                             double>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setConstantWeightPerObservable),
                         py::arg("weights_per_observation_parser"),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_weight_per_observation_parser",
                         py::overload_cast<std::map<
                             std::shared_ptr<tom::ObservationCollectionParser>,
                             Eigen::VectorXd>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setConstantWeightPerObservable),
                         py::arg("weights_per_observation_parser"),
                         R"doc(No documentation found.)doc")
                    .def("set_tabulated_weights",
                         py::overload_cast<
                             const Eigen::VectorXd,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setTabulatedWeights),
                         py::arg("tabulated_weights"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_tabulated_weights",
                         py::overload_cast<std::map<
                             std::shared_ptr<tom::ObservationCollectionParser>,
                             Eigen::VectorXd>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setTabulatedWeights),
                         py::arg("tabulated_weights"),
                         R"doc(No documentation found.)doc")
                    .def("append",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::appendObservationCollection,
                         py::arg("observation_collection_to_append"))
                    .def("filter_observations",
                         py::overload_cast<
                             const std::map<
                                 std::shared_ptr<
                                     tom::ObservationCollectionParser>,
                                 std::shared_ptr<tom::ObservationFilterBase>>&,
                             const bool>(&tom::ObservationCollection<
                                         STATE_SCALAR_TYPE,
                                         TIME_TYPE>::filterObservations),
                         py::arg("observation_filters"),
                         py::arg("save_filtered_observations") = true,
                         R"doc(No documentation found.)doc")
                    .def("filter_observations",
                         py::overload_cast<
                             std::shared_ptr<tom::ObservationFilterBase>,
                             std::shared_ptr<tom::ObservationCollectionParser>,
                             const bool>(&tom::ObservationCollection<
                                         STATE_SCALAR_TYPE,
                                         TIME_TYPE>::filterObservations),
                         py::arg("observation_filters"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         py::arg("save_filtered_observations") = true,
                         R"doc(No documentation found.)doc")
                    .def("split_observation_sets",
                         py::overload_cast<
                             std::shared_ptr<tom::ObservationSetSplitterBase>,
                             std::shared_ptr<tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::splitObservationSets),
                         py::arg("observation_set_splitter"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("get_single_observation_sets",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getSingleObservationSets,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("print_observation_sets_start_and_size",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::printObservationSetsStartAndSize,
                         R"doc(No documentation found.)doc")
                    .def("remove_single_observation_sets",
                         py::overload_cast<
                             std::shared_ptr<tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::removeSingleObservationSets),
                         py::arg("observation_parser"),
                         R"doc(No documentation found.)doc")
                    .def("set_reference_point",
                         py::overload_cast<
                             tss::SystemOfBodies&, const Eigen::Vector3d&,
                             const std::string&, const std::string&,
                             const tom::LinkEndType,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setReferencePoint),
                         py::arg("bodies"), py::arg("antenna_position"),
                         py::arg("antenna_name"), py::arg("spacecraft_name"),
                         py::arg("link_end_type"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_reference_points",
                         py::overload_cast<
                             tss::SystemOfBodies&,
                             const std::map<double, Eigen::Vector3d>&,
                             const std::string&, const tom::LinkEndType,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setReferencePoints),
                         py::arg("bodies"), py::arg("antenna_switch_history"),
                         py::arg("spacecraft_name"), py::arg("link_end_type"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_reference_point",
                         py::overload_cast<
                             tss::SystemOfBodies&,
                             const std::shared_ptr<te::Ephemeris>,
                             const std::string&, const std::string&,
                             const tom::LinkEndType,
                             const std::shared_ptr<
                                 tom::ObservationCollectionParser>>(
                             &tom::ObservationCollection<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setReferencePoint),
                         py::arg("bodies"),
                         py::arg("antenna_body_fixed_ephemeris"),
                         py::arg("antenna_name"), py::arg("spacecraft_name"),
                         py::arg("link_end_type"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("set_transponder_delay",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE, TIME_TYPE>::setTransponderDelay,
                         py::arg("spacecraft_name"),
                         py::arg("transponder_delay"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("remove_empty_observation_sets",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::removeEmptySingleObservationSets,
                         R"doc(No documentation found.)doc")
                    .def(
                        "add_dependent_variable",
                        &tom::ObservationCollection<
                            STATE_SCALAR_TYPE, TIME_TYPE>::addDependentVariable,
                        py::arg("dependent_variable_settings"),
                        py::arg("bodies"),
                        py::arg("observation_parser") = std::make_shared<
                            tom::ObservationCollectionParser>(),
                        R"doc(No documentation found.)doc")
                    .def("dependent_variable",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getDependentVariables,
                         py::arg("dependent_variable_settings"),
                         py::arg("first_compatible_settings") = false,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("concatenated_dependent_variable",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getConcatenatedDependentVariables,
                         py::arg("dependent_variable_settings"),
                         py::arg("first_compatible_settings") = false,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("compatible_dependent_variable_settings",
                         &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                     TIME_TYPE>::
                             getCompatibleDependentVariablesSettingsList,
                         py::arg("dependent_variable_settings"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("compatible_dependent_variables_list",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getAllCompatibleDependentVariables,
                         py::arg("dependent_variable_settings"),
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("dependent_variable_history_per_set",
                         &tom::ObservationCollection<STATE_SCALAR_TYPE,
                                                     TIME_TYPE>::
                             getDependentVariableHistoryPerObservationSet,
                         py::arg("dependent_variable_settings"),
                         py::arg("first_compatible_settings") = false,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc")
                    .def("dependent_variable_history",
                         &tom::ObservationCollection<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getDependentVariableHistory,
                         py::arg("dependent_variable_settings"),
                         py::arg("first_compatible_settings") = false,
                         py::arg("observation_parser") = std::make_shared<
                             tom::ObservationCollectionParser>(),
                         R"doc(No documentation found.)doc");

                m.def("merge_observation_collections",
                      &tss::mergeObservationCollections<STATE_SCALAR_TYPE,
                                                        TIME_TYPE>,
                      py::arg("observation_collection_list"));

                m.def("create_single_observation_set",
                      py::overload_cast<
                          const tom::ObservableType, const tom::LinkEnds&,
                          const std::vector<Eigen::Matrix<STATE_SCALAR_TYPE,
                                                          Eigen::Dynamic, 1>>&,
                          const std::vector<TIME_TYPE>, const tom::LinkEndType,
                          const std::shared_ptr<
                              tom::ObservationAncilliarySimulationSettings>>(
                          &tom::createSingleObservationSet<STATE_SCALAR_TYPE,
                                                           TIME_TYPE>),
                      py::arg("observable_type"), py::arg("link_ends"),
                      py::arg("observations"), py::arg("observation_times"),
                      py::arg("reference_link_end"),
                      py::arg("ancillary_settings"),
                      R"doc(No documentation found.)doc");

                m.def(
                    "filter_observations",
                    py::overload_cast<
                        const std::shared_ptr<tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>>,
                        const std::shared_ptr<tom::ObservationFilterBase>,
                        const bool>(
                        &tom::filterObservations<STATE_SCALAR_TYPE, TIME_TYPE>),
                    py::arg("original_observation_set"),
                    py::arg("observation_filter"),
                    py::arg("save_filtered_observations") = false,
                    R"doc(No documentation found.)doc");

                m.def(
                    "split_observation_set",
                    py::overload_cast<
                        const std::shared_ptr<tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>>,
                        const std::shared_ptr<tom::ObservationSetSplitterBase>,
                        const bool>(&tom::splitObservationSet<STATE_SCALAR_TYPE,
                                                              TIME_TYPE>),
                    py::arg("original_observation_set"),
                    py::arg("observation_splitter"),
                    py::arg("print_warning") = true,
                    R"doc(No documentation found.)doc");

                py::class_<
                    tom::SingleObservationSet<STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                              TIME_TYPE>>>(
                    m, "SingleObservationSet",
                    R"doc(

        Class collecting a single set of observations and associated data, of a given observable type, link ends, and ancilliary data.





     )doc")
                    .def(
                        "set_observations",
                        py::overload_cast<const std::vector<Eigen::Matrix<
                            STATE_SCALAR_TYPE, Eigen::Dynamic, 1>>&>(
                            &tom::SingleObservationSet<
                                STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                        py::arg("observations"),
                        R"doc(No documentation found.)doc")
                    .def(
                        "set_observations",
                        py::overload_cast<const Eigen::Matrix<
                            STATE_SCALAR_TYPE, Eigen::Dynamic, 1>&>(
                            &tom::SingleObservationSet<
                                STATE_SCALAR_TYPE, TIME_TYPE>::setObservations),
                        py::arg("observations"),
                        R"doc(No documentation found.)doc")
                    .def("set_residuals",
                         py::overload_cast<const std::vector<Eigen::Matrix<
                             STATE_SCALAR_TYPE, Eigen::Dynamic, 1>>&>(
                             &tom::SingleObservationSet<
                                 STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                         py::arg("residuals"),
                         R"doc(No documentation found.)doc")
                    .def("set_residuals",
                         py::overload_cast<const Eigen::Matrix<
                             STATE_SCALAR_TYPE, Eigen::Dynamic, 1>&>(
                             &tom::SingleObservationSet<
                                 STATE_SCALAR_TYPE, TIME_TYPE>::setResiduals),
                         py::arg("residuals"),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_weight",
                         py::overload_cast<const double>(
                             &tom::SingleObservationSet<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setConstantWeight),
                         py::arg("weight"), R"doc(No documentation found.)doc")
                    .def("set_constant_weight",
                         py::overload_cast<
                             const Eigen::Matrix<double, Eigen::Dynamic, 1>&>(
                             &tom::SingleObservationSet<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setConstantWeight),
                         py::arg("weight"), R"doc(No documentation found.)doc")
                    .def("set_tabulated_weights",
                         py::overload_cast<const Eigen::VectorXd&>(
                             &tom::SingleObservationSet<
                                 STATE_SCALAR_TYPE,
                                 TIME_TYPE>::setTabulatedWeights),
                         py::arg("weights"), R"doc(No documentation found.)doc")
                    .def("filter_observations",
                         &tom::SingleObservationSet<
                             STATE_SCALAR_TYPE, TIME_TYPE>::filterObservations,
                         py::arg("filter"), py::arg("save_filtered_obs") = true,
                         R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "observable_type",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getObservableType,
                        R"doc(

        **read-only**

        Type of observable for which the object stores observations

        :type: ObservableType
     )doc")
                    .def_property(
                        "link_definition",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getLinkEnds,
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::setLinkEnds,
                        R"doc(

        **read-only**

        Definition of the link ends for which the object stores observations

        :type: LinkDefinition
     )doc")
                    .def_property_readonly(
                        "reference_link_end",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getReferenceLinkEnd,
                        R"doc(

        **read-only**

        Reference link end for all stored observations

        :type: LinkEndType
     )doc")
                    .def_property_readonly(
                        "number_of_observables",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getNumberOfObservables,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "single_observable_size",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getSingleObservableSize,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "total_observation_set_size",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getTotalObservationSetSize,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "time_bounds",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getTimeBounds,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "list_of_observations",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getObservations,
                        R"doc(

        **read-only**

        List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1

        :type: list[ numpy.ndarray[numpy.float64[m, 1]] ]
     )doc")
                    .def_property_readonly(
                        "observation_times",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getObservationTimes,
                        R"doc(

        **read-only**

        Reference time for each of the observations in ``list_of_observations``

        :type: list[ float]
     )doc")
                    .def_property_readonly(
                        "concatenated_observations",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getObservationsVector,
                        R"doc(

        **read-only**

        Concatenated vector of all stored observations

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc")
                    .def_property_readonly(
                        "computed_observations",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getComputedObservations,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "concatenated_computed_observations",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getComputedObservationsVector,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "residuals",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getResiduals,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "concatenated_residuals",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getResidualsVector,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "rms_residuals",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getRmsResiduals,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "mean_residuals",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getMeanResiduals,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "weights",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getWeights,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "concatenad_weights",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getWeightsVector,
                        R"doc(No documentation found.)doc")
                    .def_property(
                        "dependent_variables",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getObservationsDependentVariables,
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::setObservationsDependentVariables,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "dependent_variables_history",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getDependentVariableHistory,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "observations_history",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getObservationsHistory,
                        R"doc(

        **read-only**

        Dictionary of observations sorted by time. Created by making a dictionaty with ``observation_times`` as keys and  ``list_of_observations`` as values

        :type: dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
     )doc")
                    .def_property_readonly(
                        "ancilliary_settings",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getAncilliarySettings,
                        R"doc(

        **read-only**

        Ancilliary settings all stored observations

        :type: ObservationAncilliarySimulationSettings
     )doc")
                    .def_property(
                        "weights_vector",
                        &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                   TIME_TYPE>::getWeightsVector,
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE, TIME_TYPE>::setTabulatedWeights,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "filtered_observation_set",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getFilteredObservationSet,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "number_filtered_observations",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getNumberOfFilteredObservations,
                        R"doc(No documentation found.)doc")
                    .def(
                        "single_dependent_variable",
                        py::overload_cast<
                            std::shared_ptr<
                                tss::ObservationDependentVariableSettings>,
                            const bool>(&tom::SingleObservationSet<
                                        STATE_SCALAR_TYPE,
                                        TIME_TYPE>::getSingleDependentVariable),
                        py::arg("dependent_variable_settings"),
                        py::arg("return_first_compatible_settings") = false,
                        R"doc(No documentation found.)doc")
                    .def("compatible_dependent_variable_settings",
                         &tom::SingleObservationSet<STATE_SCALAR_TYPE,
                                                    TIME_TYPE>::
                             getCompatibleDependentVariablesSettingsList,
                         R"doc(No documentation found.)doc")
                    .def("compatible_dependent_variables_list",
                         &tom::SingleObservationSet<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getAllCompatibleDependentVariables,
                         R"doc(No documentation found.)doc")
                    .def("single_dependent_variable_history",
                         &tom::SingleObservationSet<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::getSingleDependentVariableHistory,
                         R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "dependent_variables_matrix",
                        &tom::SingleObservationSet<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getObservationsDependentVariablesMatrix,
                        R"doc(No documentation found.)doc");


                m.def("single_observation_set",
                      &tss::singleObservationSetWithoutDependentVariables<
                          STATE_SCALAR_TYPE, TIME_TYPE>,
                      py::arg("observable_type"), py::arg("link_definition"),
                      py::arg("observations"), py::arg("observation_times"),
                      py::arg("reference_link_end"),
                      py::arg("ancilliary_settings") = nullptr,
                      R"doc(No documentation found.)doc");

                /*!
                 *************** STATE TRANSITION INTERFACE ***************
                 */

                py::class_<
                    tp::CombinedStateTransitionAndSensitivityMatrixInterface,
                    std::shared_ptr<
                        tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                    m, "CombinedStateTransitionAndSensitivityMatrixInterface",
                    R"doc(

        Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.

        Class establishing an interface to the State Transition and Sensitivity Matrices.
        Instances of this class are instantiated automatically upon creation of :class:`~tudatpy.numerical_simulation.Estimator` objects,
        using the simulation information in the observation, propagation and integration settings that the :class:`~tudatpy.numerical_simulation.Estimator` instance is linked to.





     )doc")
                    .def(
                        "state_transition_sensitivity_at_epoch",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getCombinedStateTransitionAndSensitivityMatrix,
                        py::arg("time"),
                        py::arg("add_central_body_dependency") = true,
                        py::arg("arc_defining_bodies") =
                            std::vector<std::string>(),
                        R"doc(

        Function to get the concatenated state transition and sensitivity matrix at a given time.

        Function to get the concatenated state transition and sensitivity matrix at a given time.
        Entries corresponding to parameters which are not active at the current arc are omitted.


        Parameters
        ----------
        time : float
            Time at which concatenated state transition and sensitivity matrix are to be retrieved.
        Returns
        -------
        numpy.ndarray[numpy.float64[m, n]]
            Concatenated state transition and sensitivity matrix at a given time.





    )doc")
                    .def(
                        "full_state_transition_sensitivity_at_epoch",
                        &tp::CombinedStateTransitionAndSensitivityMatrixInterface::
                            getFullCombinedStateTransitionAndSensitivityMatrix,
                        py::arg("time"),
                        py::arg("add_central_body_dependency") = true,
                        py::arg("arc_defining_bodies") =
                            std::vector<std::string>(),
                        R"doc(


        Parameters
        ----------
        time : float
            Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
        Returns
        -------
        numpy.ndarray[numpy.float64[m, n]]
            Full concatenated state transition and sensitivity matrix at a given time.





    )doc")
                    .def_property_readonly(
                        "state_transition_size",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getStateTransitionMatrixSize,
                        R"doc(

        **read-only**

        Size of the (square) state transition matrix.

        :type: int
     )doc")
                    .def_property_readonly(
                        "sensitivity_size",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getSensitivityMatrixSize,
                        R"doc(

        **read-only**

        Number of columns in the sensitivity matrix.

        :type: int
     )doc")
                    .def_property_readonly(
                        "full_parameter_size",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getFullParameterVectorSize,
                        R"doc(

        **read-only**

        Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.

        :type: int
     )doc");

                /*!
                 *************** COVARIANCE ***************
                 */

                m.def("propagate_covariance_rsw_split_output",
                      &tp::propagateCovarianceVectorsRsw,
                      py::arg("covariance_output"), py::arg("estimator"),
                      py::arg("output_times"),
                      R"doc(No documentation found.)doc");


                m.def("propagate_formal_errors_rsw_split_output",
                      &tp::propagateFormalErrorVectorsRsw,
                      py::arg("covariance_output"), py::arg("estimator"),
                      py::arg("output_times"),
                      R"doc(No documentation found.)doc");

                m.def(
                    "propagate_covariance_split_output",
                    py::overload_cast<
                        const Eigen::MatrixXd,
                        const std::shared_ptr<
                            tp::CombinedStateTransitionAndSensitivityMatrixInterface>,
                        const std::vector<double>>(
                        &tp::propagateCovarianceVectors),
                    py::arg("initial_covariance"),
                    py::arg("state_transition_interface"),
                    py::arg("output_times"),
                    R"doc(No documentation found.)doc");

                m.def(
                    "propagate_covariance",
                    py::overload_cast<
                        const Eigen::MatrixXd,
                        const std::shared_ptr<
                            tp::CombinedStateTransitionAndSensitivityMatrixInterface>,
                        const std::vector<double>>(&tp::propagateCovariance),
                    py::arg("initial_covariance"),
                    py::arg("state_transition_interface"),
                    py::arg("output_times"),
                    R"doc(

Function to propagate system covariance through time.

Function to propagate the covariance of a given system through time.
The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


Parameters
----------
initial_covariance : numpy.ndarray[numpy.float64[m, n]]
    System covariance matrix (symmetric and positive semi-definite) at initial time.
    Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
    Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

output_times : List[ float ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the covariance propagation,
    which always adheres to the integrator settings that the `state_transition_interface` links to.
    Output times which do not coincide with integration time steps are calculated via interpolation.

Returns
-------
Dict[ float, numpy.ndarray[numpy.float64[m, n]] ]
    Dictionary reporting the propagated covariances at each output time.






    )doc");

                m.def(
                    "propagate_formal_errors_split_output",
                    py::overload_cast<
                        const Eigen::MatrixXd,
                        const std::shared_ptr<
                            tp::CombinedStateTransitionAndSensitivityMatrixInterface>,
                        const std::vector<double>>(
                        &tp::propagateFormalErrorVectors),
                    py::arg("initial_covariance"),
                    py::arg("state_transition_interface"),
                    py::arg("output_times"),
                    R"doc(

Function to propagate system formal errors through time.

Function to propagate the formal errors of a given system through time.
Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


Parameters
----------
initial_covariance : numpy.ndarray[numpy.float64[m, n]]
    System covariance matrix (symmetric and positive semi-definite) at initial time.
    Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
    Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

output_times : List[ float ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the covariance propagation,
    which always adheres to the integrator settings that the `state_transition_interface` links to.
    Output times which do not coincide with integration time steps are calculated via interpolation.

Returns
-------
Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
    Dictionary reporting the propagated formal errors at each output time.






    )doc");


                m.def(
                    "propagate_formal_errors",
                    py::overload_cast<
                        const Eigen::MatrixXd,
                        const std::shared_ptr<
                            tp::CombinedStateTransitionAndSensitivityMatrixInterface>,
                        const std::vector<double>>(&tp::propagateFormalErrors),
                    py::arg("initial_covariance"),
                    py::arg("state_transition_interface"),
                    py::arg("output_times"),
                    R"doc(

Function to propagate system formal errors through time.

Function to propagate the formal errors of a given system through time.
Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.


Parameters
----------
initial_covariance : numpy.ndarray[numpy.float64[m, n]]
    System covariance matrix (symmetric and positive semi-definite) at initial time.
    Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)

state_transition_interface : :class:`~tudatpy.numerical_simulation.estimation.CombinedStateTransitionAndSensitivityMatrixInterface`
    Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.

output_times : List[ float ]
    Times at which the propagated covariance matrix shall be reported.
    Note that this argument has no impact on the integration time-steps of the covariance propagation,
    which always adheres to the integrator settings that the `state_transition_interface` links to.
    Output times which do not coincide with integration time steps are calculated via interpolation.

Returns
-------
Dict[ float, numpy.ndarray[numpy.float64[m, 1]] ]
    Dictionary reporting the propagated formal errors at each output time.






    )doc");


                /*!
                 *************** ESTIMATION ***************
                 */

                py::class_<tss::EstimationConvergenceChecker,
                           std::shared_ptr<tss::EstimationConvergenceChecker>>(
                    m, "EstimationConvergenceChecker",
                    R"doc(

        Class defining the convergence criteria for an estimation.

        Class defining the convergence criteria for an estimation.
        The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker` function.





     )doc");

                m.def("estimation_convergence_checker",
                      &tss::estimationConvergenceChecker,
                      py::arg("maximum_iterations") = 5,
                      py::arg("minimum_residual_change") = 0.0,
                      py::arg("minimum_residual") = 0.0,
                      py::arg("number_of_iterations_without_improvement") = 2,
                      R"doc(

Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object.

Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.


Parameters
----------
maximum_iterations : int, default = 5
    Maximum number of allowed iterations for estimation.
minimum_residual_change : float, default = 0.0
    Minimum required change in residual between two iterations.
minimum_residual : float, default = 0.0
    Minimum value of observation residual below which estimation is converged.
number_of_iterations_without_improvement : int, default = 2
    Number of iterations without reduction of residual.
Returns
-------
:class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker`
    Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation.






    )doc");


                py::class_<
                    tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tss::CovarianceAnalysisInput<
                        STATE_SCALAR_TYPE, TIME_TYPE>>>(
                    m, "CovarianceAnalysisInput",
                    R"doc(

        Class for defining all specific inputs to a covariance analysis.





     )doc")
                    .def(py::init<
                             const std::shared_ptr<tom::ObservationCollection<
                                 STATE_SCALAR_TYPE, TIME_TYPE>>&,
                             const Eigen::MatrixXd, const Eigen::MatrixXd>(),
                         py::arg("observations_and_times"),
                         py::arg("inverse_apriori_covariance") =
                             Eigen::MatrixXd::Zero(0, 0),
                         py::arg("consider_covariance") =
                             Eigen::MatrixXd::Zero(0, 0),
                         R"doc(

        Class constructor.

        Constructor through which the user can create instances of this class. Note that the weight are all initiated as 1.0, and the default settings of ``define_covariance_settings`` are used.


        Parameters
        ----------
        observations_and_times : ObservationCollection
            Total data structure of observations and associated times/link ends/type/etc.
        inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
            A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput`
            Instance of the :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.





    )doc")
                    .def("set_constant_weight",
                         &tss::CovarianceAnalysisInput<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::setConstantWeightsMatrix,
                         py::arg("weight"),
                         R"doc(

        Function to set a constant weight matrix for all observables.

        Function to set a constant weight matrix for all observables.
        The weights are applied to all observations managed by the given PodInput object.


        Parameters
        ----------
        constant_weight : float
            Constant weight factor that is to be applied to all observations.
        Returns
        -------
        None
            Function modifies the object in-place.





    )doc")
                    .def("set_weights_from_observation_collection",
                         &tss::CovarianceAnalysisInput<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::setWeightsFromObservationCollection,
                         R"doc(No documentation found.)doc")
                    .def("set_constant_single_observable_weight",
                         &tss::CovarianceAnalysisInput<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::setConstantSingleObservableWeights,
                         py::arg("observable_type"), py::arg("weight"),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_single_observable_vector_weight",
                         &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE,
                                                       TIME_TYPE>::
                             setConstantSingleObservableVectorWeights,
                         py::arg("observable_type"), py::arg("weight"),
                         R"doc(No documentation found.)doc")
                    .def("set_constant_single_observable_and_link_end_weight",
                         &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE,
                                                       TIME_TYPE>::
                             setConstantSingleObservableAndLinkEndsWeights,
                         py::arg("observable_type"), py::arg("link_ends"),
                         py::arg("weight"), R"doc(No documentation found.)doc")
                    .def(
                        "set_constant_single_observable_and_link_end_vector_"
                        "weight",
                        &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE,
                                                      TIME_TYPE>::
                            setConstantSingleObservableAndLinkEndsVectorWeights,
                        py::arg("observable_type"), py::arg("link_ends"),
                        py::arg("weight"), R"doc(No documentation found.)doc")
                    .def(
                        "set_total_single_observable_and_link_end_vector_"
                        "weight",
                        &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE,
                                                      TIME_TYPE>::
                            setTabulatedSingleObservableAndLinkEndsWeights,
                        py::arg("observable_type"), py::arg("link_ends"),
                        py::arg("weight_vector"),
                        R"doc(No documentation found.)doc")
                    .def("set_constant_weight_per_observable",
                         &tss::CovarianceAnalysisInput<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::setConstantPerObservableWeightsMatrix,
                         py::arg("weight_per_observable"),
                         R"doc(

        Function to set a constant weight matrix for a given type of observable.

        Function to set a constant weight matrix for a given type of observable.
        The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.


        Parameters
        ----------
        constant_weight : Dict[ :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservableType`, float ]
            Constant weight factor that is to be applied to all observations.
        Returns
        -------
        None
            Function modifies the object in-place.





    )doc")
                    .def("set_constant_vector_weight_per_observable",
                         &tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE,
                                                       TIME_TYPE>::
                             setConstantPerObservableVectorWeightsMatrix,
                         py::arg("weight_per_observable"),
                         R"doc(No documentation found.)doc")
                    .def("define_covariance_settings",
                         &tss::CovarianceAnalysisInput<
                             STATE_SCALAR_TYPE,
                             TIME_TYPE>::defineCovarianceSettings,
                         py::arg("reintegrate_equations_on_first_iteration") =
                             true,
                         py::arg("reintegrate_variational_equations") = true,
                         py::arg("save_design_matrix") = true,
                         py::arg("print_output_to_terminal") = true,
                         py::arg("limit_condition_number_for_warning") = 1.0E8,
                         R"doc(

        Function to define specific settings for covariance analysis process

        Function to define specific settings for covariance analysis process


        Parameters
        ----------
        reintegrate_equations : bool, default = True
            Boolean denoting whether the dynamics and variational equations are to be reintegrated
            or if existing values are to be used to perform first iteration.

        reintegrate_variational_equations : bool, default = True
            Boolean denoting whether the variational equations are to be reintegrated during estimation
            (if this is set to False, and ``reintegrate_equations`` to true, only the dynamics are re-integrated)

        save_design_matrix : bool, default = True
            Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
            :math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.

        print_output_to_terminal : bool, default = True
            Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.

        Returns
        -------
        None
            Function modifies the object in-place.





    )doc")
                    .def_property("weight_matrix_diagonal",
                                  &tss::CovarianceAnalysisInput<
                                      STATE_SCALAR_TYPE,
                                      TIME_TYPE>::getWeightsMatrixDiagonals,
                                  &tss::CovarianceAnalysisInput<
                                      STATE_SCALAR_TYPE,
                                      TIME_TYPE>::setWeightsMatrixDiagonals,
                                  R"doc(

        **read-only**

        Complete diagonal of the weights matrix that is to be used

        :type: numpy.ndarray[numpy.float64[n, 1]]
     )doc");

                py::class_<
                    tss::EstimationInput<STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<
                        tss::EstimationInput<STATE_SCALAR_TYPE, TIME_TYPE>>,
                    tss::CovarianceAnalysisInput<STATE_SCALAR_TYPE, TIME_TYPE>>(
                    m, "EstimationInput",
                    R"doc(

        Class for defining all inputs to the estimation.





     )doc")
                    .def(py::init<
                             const std::shared_ptr<tom::ObservationCollection<
                                 STATE_SCALAR_TYPE, TIME_TYPE>>&,
                             const Eigen::MatrixXd,
                             std::shared_ptr<tss::EstimationConvergenceChecker>,
                             const Eigen::MatrixXd, const Eigen::VectorXd,
                             const bool>(),
                         py::arg("observations_and_times"),
                         py::arg("inverse_apriori_covariance") =
                             Eigen::MatrixXd::Zero(0, 0),
                         py::arg("convergence_checker") = std::make_shared<
                             tss::EstimationConvergenceChecker>(),
                         py::arg("consider_covariance") =
                             Eigen::MatrixXd::Zero(0, 0),
                         py::arg("consider_parameters_deviations") =
                             Eigen::VectorXd::Zero(0),
                         py::arg("apply_final_parameter_correction") = true,
                         R"doc(

        Class constructor.

        Constructor through which the user can create instances of this class.


        Parameters
        ----------
        observations_and_times : ObservationCollection
            Total data structure of observations and associated times/link ends/type/etc.
        inverse_apriori_covariance : numpy.ndarray[numpy.float64[m, n]], default = [ ]
            A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
        convergence_checker : :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker`, default = :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker`
            Object defining when the estimation is converged.
        Returns
        -------
        :class:`~tudatpy.numerical_simulation.estimation.EstimationInput`
            Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationInput` class, defining the data and other settings to be used for the estimation.





    )doc")
                    .def(
                        "define_estimation_settings",
                        &tss::EstimationInput<STATE_SCALAR_TYPE, TIME_TYPE>::
                            defineEstimationSettings,
                        py::arg("reintegrate_equations_on_first_iteration") =
                            true,
                        py::arg("reintegrate_variational_equations") = true,
                        py::arg("save_design_matrix") = true,
                        py::arg("print_output_to_terminal") = true,
                        py::arg("save_residuals_and_parameters_per_iteration") =
                            true,
                        py::arg("save_state_history_per_iteration") = false,
                        py::arg("limit_condition_number_for_warning") = 1.0E8,
                        py::arg("condition_number_warning_each_iteration") =
                            true,
                        R"doc(

        Function to define specific settings for the estimation process

        Function to define specific settings for covariance analysis process


        Parameters
        ----------
        reintegrate_equations_on_first_iteration : bool, default = True
            Boolean denoting whether the dynamics and variational equations are to be reintegrated
            or if existing values are to be used to perform first iteration.

        reintegrate_variational_equations : bool, default = True
            Boolean denoting whether the variational equations are to be reintegrated during estimation
            (if this is set to False, and ``reintegrate_equations_on_first_iteration`` to true, only the dynamics are re-integrated)

        save_design_matrix : bool, default = True
            Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
            :math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.

        print_output_to_terminal : bool, default = True
            Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.

        save_residuals_and_parameters_per_iteration : bool, default = True
            Boolean denoting whether the residuals and parameters from the each iteration are to be saved.

        save_state_history_per_iteration : bool, default = False
            Boolean denoting whether the state history and dependent variables are to be saved on each iteration.

        Returns
        -------
        None
            Function modifies the object in-place.





    )doc");

                m.attr("PodInput") = m.attr("EstimationInput");


                py::class_<
                    tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE, TIME_TYPE>,
                    std::shared_ptr<tss::CovarianceAnalysisOutput<
                        STATE_SCALAR_TYPE, TIME_TYPE>>>(
                    m, "CovarianceAnalysisOutput",
                    R"doc(

        Class collecting all outputs from the covariance analysis process.





     )doc")
                    .def_property_readonly(
                        "inverse_covariance",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getUnnormalizedInverseCovarianceMatrix,
                        R"doc(

        **read-only**

        (Unnormalized) inverse estimation covariance matrix :math:`\mathbf{P}^{-1}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )doc")
                    .def_property_readonly(
                        "covariance",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getUnnormalizedCovarianceMatrix,
                        R"doc(

        **read-only**

        (Unnormalized) estimation covariance matrix :math:`\mathbf{P}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )doc")
                    .def_property_readonly(
                        "inverse_normalized_covariance",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getNormalizedInverseCovarianceMatrix,
                        R"doc(

        **read-only**

        Normalized inverse estimation covariance matrix :math:`\mathbf{\tilde{P}}^{-1}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )doc")
                    .def_property_readonly(
                        "normalized_covariance",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getNormalizedCovarianceMatrix,
                        R"doc(

        **read-only**

        Normalized estimation covariance matrix :math:`\mathbf{\tilde{P}}`.

        :type: numpy.ndarray[numpy.float64[m, m]]
     )doc")
                    .def_property_readonly(
                        "formal_errors",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getFormalErrorVector,
                        R"doc(

        **read-only**

        Formal error vector :math:`\boldsymbol{\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance)s

        :type: numpy.ndarray[numpy.float64[m, 1]]s
     )doc")
                    .def_property_readonly(
                        "correlations",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE, TIME_TYPE>::getCorrelationMatrix,
                        R"doc(

        **read-only**

        Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\sigma_{i}\sigma_{j})`

        :type: numpy.ndarray[numpy.float64[m, m]]
     )doc")
                    .def_property_readonly(
                        "design_matrix",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getUnnormalizedDesignMatrix,
                        R"doc(

        **read-only**

        Matrix of unnormalized partial derivatives :math:`\mathbf{H}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}}`.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )doc")
                    .def_property_readonly(
                        "normalized_design_matrix",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getNormalizedDesignMatrix,
                        R"doc(

        **read-only**

        Matrix of normalized partial derivatives :math:`\tilde{\mathbf{H}}`.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )doc")
                    .def_property_readonly(
                        "weighted_design_matrix",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getUnnormalizedWeightedDesignMatrix,
                        R"doc(

        **read-only**

        Matrix of weighted partial derivatives, equal to :math:`\mathbf{W}^{1/2}{\mathbf{H}}`

        :type: numpy.ndarray[numpy.float64[m, n]]
     )doc")
                    .def_property_readonly(
                        "weighted_normalized_design_matrix",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getNormalizedWeightedDesignMatrix,
                        R"doc(

        **read-only**

        Matrix of weighted, normalized partial derivatives, equal to :math:`\mathbf{W}^{1/2}\tilde{\mathbf{H}}`

        :type: numpy.ndarray[numpy.float64[m, n]]
     )doc")
                    .def_property_readonly(
                        "consider_covariance_contribution",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getConsiderCovarianceContribution,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "normalized_covariance_with_consider_parameters",
                        &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE,
                                                       TIME_TYPE>::
                            getNormalizedCovarianceWithConsiderParameters,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "unnormalized_covariance_with_consider_parameters",
                        &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE,
                                                       TIME_TYPE>::
                            getUnnormalizedCovarianceWithConsiderParameters,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "normalized_design_matrix_consider_parameters",
                        &tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE,
                                                       TIME_TYPE>::
                            getNormalizedDesignMatrixConsiderParameters,
                        R"doc(No documentation found.)doc")
                    .def_property_readonly(
                        "consider_normalization_factors",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::getConsiderNormalizationFactors,
                        R"doc(No documentation found.)doc")
                    .def_readonly(
                        "normalization_terms",
                        &tss::CovarianceAnalysisOutput<
                            STATE_SCALAR_TYPE,
                            TIME_TYPE>::designMatrixTransformationDiagonal_,
                        R"doc(

        **read-only**

        Vector of normalization terms used for covariance and design matrix

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc");


                py::class_<tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>,
                           std::shared_ptr<tss::EstimationOutput<
                               STATE_SCALAR_TYPE, TIME_TYPE>>,
                           tss::CovarianceAnalysisOutput<STATE_SCALAR_TYPE,
                                                         TIME_TYPE>>(
                    m, "EstimationOutput",
                    R"doc(

        Class collecting all outputs from the iterative estimation process.





     )doc")
                    .def_property_readonly(
                        "residual_history",
                        &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::
                            getResidualHistoryMatrix,
                        R"doc(

        **read-only**

        Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )doc")
                    .def_property_readonly(
                        "parameter_history",
                        &tss::EstimationOutput<STATE_SCALAR_TYPE, TIME_TYPE>::
                            getParameterHistoryMatrix,
                        R"doc(

        **read-only**

        Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.

        :type: numpy.ndarray[numpy.float64[m, n]]
     )doc")
                    .def_property_readonly(
                        "simulation_results_per_iteration",
                        &tss::EstimationOutput<STATE_SCALAR_TYPE,
                                               TIME_TYPE>::getSimulationResults,
                        R"doc(

        **read-only**

        List of complete numerical propagation results, with the :math:`i^{th}` entry of thee list thee results of the :math:`i^{th}` propagation

        :type: list[SimulationResults]
     )doc")
                    .def_readonly("final_residuals",
                                  &tss::EstimationOutput<STATE_SCALAR_TYPE,
                                                         TIME_TYPE>::residuals_,
                                  R"doc(

        **read-only**

        Vector of post-fit observation residuals, for the iteration with the lowest rms residuals.

        :type: numpy.ndarray[numpy.float64[m, 1]]
     )doc")
                    .def_readonly(
                        "final_parameters",
                        &tss::EstimationOutput<STATE_SCALAR_TYPE,
                                               TIME_TYPE>::parameterEstimate_,
                        R"doc(No documentation found.)doc")
                    .def_readonly(
                        "best_iteration",
                        &tss::EstimationOutput<STATE_SCALAR_TYPE,
                                               TIME_TYPE>::bestIteration_,
                        R"doc(No documentation found.)doc");

                m.attr("PodOutput") = m.attr("EstimationOutput");
            }

        }  // namespace estimation
    }  // namespace numerical_simulation
}  // namespace tudatpy
