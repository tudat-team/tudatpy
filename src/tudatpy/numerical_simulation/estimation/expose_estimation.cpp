/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <pybind11/chrono.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <tudat/astro/propagators/propagateCovariance.h>
#include <tudat/basics/utilities.h>
#include <tudat/simulation/estimation_setup/fitOrbitToEphemeris.h>

#include "tudatpy/docstrings.h"
#include "tudatpy/scalarTypes.h"

namespace py = pybind11;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace ts = tudat::statistics;
namespace tss = tudat::simulation_setup;
namespace tni = tudat::numerical_integrators;
namespace tom = tudat::observation_models;
namespace trf = tudat::reference_frames;


namespace tudat {

    namespace propagators {

        std::map<double, Eigen::MatrixXd> propagateCovarianceRsw(
            const std::shared_ptr<
                tss::CovarianceAnalysisOutput<double, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<double, TIME_TYPE>>
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

            std::shared_ptr<tep::EstimatableParameterSet<double>> parameterSet =
                orbitDeterminationManager->getParametersToEstimate();

            std::map<int, std::shared_ptr<
                              tep::EstimatableParameter<Eigen::VectorXd>>>
                initialStates = parameterSet->getInitialStateParameters();
            std::map<std::pair<std::string, std::string>, std::vector<int>>
                transformationList;
            for(auto it : initialStates) {
                if(std::dynamic_pointer_cast<
                       tep::InitialTranslationalStateParameter<double>>(
                       it.second)) {
                    std::shared_ptr<
                        tep::InitialTranslationalStateParameter<double>>
                        currentInitialState = std::dynamic_pointer_cast<
                            tep::InitialTranslationalStateParameter<double>>(
                            it.second);
                    transformationList
                        [std::make_pair(currentInitialState->getParameterName()
                                            .second.first,
                                        currentInitialState->getCentralBody())]
                            .push_back(it.first);

                } else if(std::dynamic_pointer_cast<
                              tep::ArcWiseInitialTranslationalStateParameter<
                                  double>>(it.second)) {
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
                double currentTime = it.first;
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
                tss::CovarianceAnalysisOutput<double, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<double, TIME_TYPE>>
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
                tss::CovarianceAnalysisOutput<double, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<double, TIME_TYPE>>
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
                tss::CovarianceAnalysisOutput<double, TIME_TYPE>>
                estimationOutput,
            const std::shared_ptr<
                tss::OrbitDeterminationManager<double, TIME_TYPE>>
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
                referenceLinkEnd,
                std::vector<
                    Eigen::Matrix<ObservationScalarType, Eigen::Dynamic, 1>>(),
                nullptr, ancilliarySettings);
        }


    }  // namespace simulation_setup

}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace estimation {


            PYBIND11_MODULE(expose_estimation, m) {
                /*!
                 *************** PARAMETERS ***************
                 */


                py::class_<
                    tep::EstimatableParameterSet<double>,
                    std::shared_ptr<tep::EstimatableParameterSet<double>>>(
                    m, "EstimatableParameterSet",
R"doc(Class containing a consolidated set of estimatable parameters.

	Class containing a consolidated set of estimatable parameters, linked to the environment and acceleration settings of the simulation.
	The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation_setup.create_parameters_to_estimate` factory function.
	
)doc")
                    .def_property_readonly(
                        "parameter_set_size",
                        &tep::EstimatableParameterSet<
                            double>::getEstimatedParameterSetSize,
R"doc(Size of the parameter set, i.e. amount of estimatable parameters contained in the set.
	)doc")
                    .def_property_readonly(
                        "initial_states_size",
                        &tep::EstimatableParameterSet<
                            double>::getInitialDynamicalStateParameterSize,
R"doc(Amount of initial state parameters contained in the set.
	)doc")
                    .def_property_readonly(
                        "initial_single_arc_states_size",
                        &tep::EstimatableParameterSet<double>::
                            getInitialDynamicalSingleArcStateParameterSize,
R"doc(Amount of initial state parameters in the set, which are treated in a single-arc fashion.
	)doc")
                    .def_property_readonly(
                        "initial_multi_arc_states_size",
                        &tep::EstimatableParameterSet<double>::
                            getInitialDynamicalMultiArcStateParameterSize,
R"doc(Amount of initial state parameters in the set, which are treated in a multi-arc fashion.
	)doc")
                    .def_property_readonly(
                        "constraints_size",
                        &tep::EstimatableParameterSet<
                            double>::getConstraintSize,
R"doc(Total size of linear constraint that is to be applied during estimation.
	)doc")
                    .def_property(
                        "parameter_vector",
                        &tep::EstimatableParameterSet<
                            double>::getFullParameterValues<double>,
                        &tep::EstimatableParameterSet<
                            double>::resetParameterValues<double>,
R"doc(Vector containing the parameter values of all parameters in the set.
	)doc")
                    .def("indices_for_parameter_type",
                         &tep::EstimatableParameterSet<
                             double>::getIndicesForParameterType,
                         py::arg("parameter_type"),
R"doc(Function to retrieve the indices of a given type of parameter.

	Function to retrieve the index of all parameters of a given type from the parameter set.
	This function can be very useful, since the order of parameters within the parameter set does not necessarily correspond to the order in which the elements were added to the set.
	

	:param parameter_type:
		help
	:return:
		help
)doc");

                /*!
                 *************** OBSERVATIONS ***************
                 */

                py::class_<
                    tom::ObservationViabilityCalculator,
                    std::shared_ptr<tom::ObservationViabilityCalculator>>(
                    m, "ObservationViabilityCalculator",
R"doc(Template class for observation viability calculators.

	Template class for classes which conducts viability calculations on simulated observations.
	Instances of the applicable ObservationViabilityCalculators are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSimulationSettings` objects during the simulation of observations (:func:`~tudatpy.numerical_simulation.estimation.simulate_observations`).
	The user typically does not interact directly with this class.
	
)doc")
                    .def("is_observation_viable",
                         &tom::ObservationViabilityCalculator::
                             isObservationViable,
                         py::arg("link_end_states"), py::arg("link_end_times"),
R"doc(Function to check whether an observation is viable.

	Function to check whether an observation is viable.
	The calculation is performed based on the given times and link end states.
	Note, that this function is called automatically during the simulation of observations.
	Direct calls to this function are generally not required.
	

	:param link_end_states:
		Vector of states of the link ends involved in the observation.
	:param link_end_times:
		Vector of times at the link ends involved in the observation.
	:return:
		True if observation is viable, false if not.
)doc");

                py::class_<tom::ObservationSimulatorBase<double, TIME_TYPE>,
                           std::shared_ptr<tom::ObservationSimulatorBase<
                               double, TIME_TYPE>>>(
                    m, "ObservationSimulator",
R"doc(Class hosting the functionality for simulating observations.

	Class hosting the functionality for simulating a given observable over a defined link geometry.
	Instances of this class are automatically created from the given :class:`~tudatpy.numerical_simulation.estimation_setup.observation.ObservationSettings` objects upon instantiation of the :class:`~tudatpy.numerical_simulation.Estimator` class.
	
)doc");

                py::class_<tom::ObservationSimulator<1, double, TIME_TYPE>,
                           std::shared_ptr<
                               tom::ObservationSimulator<1, double, TIME_TYPE>>,
                           tom::ObservationSimulatorBase<double, TIME_TYPE>>(
                    m, "ObservationSimulator_1",
get_docstring("ObservationSimulator_1").c_str());

                py::class_<tom::ObservationSimulator<2, double, TIME_TYPE>,
                           std::shared_ptr<
                               tom::ObservationSimulator<2, double, TIME_TYPE>>,
                           tom::ObservationSimulatorBase<double, TIME_TYPE>>(
                    m, "ObservationSimulator_2",
get_docstring("ObservationSimulator_2").c_str());

                py::class_<tom::ObservationSimulator<3, double, TIME_TYPE>,
                           std::shared_ptr<
                               tom::ObservationSimulator<3, double, TIME_TYPE>>,
                           tom::ObservationSimulatorBase<double, TIME_TYPE>>(
                    m, "ObservationSimulator_3",
get_docstring("ObservationSimulator_3").c_str());

                py::class_<tom::ObservationSimulator<6, double, TIME_TYPE>,
                           std::shared_ptr<
                               tom::ObservationSimulator<6, double, TIME_TYPE>>,
                           tom::ObservationSimulatorBase<double, TIME_TYPE>>(
                    m, "ObservationSimulator_6",
get_docstring("ObservationSimulator_6").c_str());

                m.def("simulate_observations",
                      &tss::simulateObservations<double, TIME_TYPE>,
                      py::arg("simulation_settings"),
                      py::arg("observation_simulators"), py::arg("bodies"),
R"doc(Function to simulate observations.

	Function to simulate observations from set observation simulators and observation simulator settings.
	Automatically iterates over all provided observation simulators, generating the full set of simulated observations.
	

	:param observation_to_simulate:
		List of settings objects, each object providing the observation time settings for simulating one type of observable and link end set.
		
	:param observation_simulators:
		List of :class:`~tudatpy.numerical_simulation.estimation.ObservationSimulator` objects, each object hosting the functionality for simulating one type of observable and link end set.
		
	:param bodies:
		Object consolidating all bodies and environment models, including ground station models, that constitute the physical environment.
		
	:return:
		Object collecting all products of the observation simulation.
)doc");


                m.def("create_pseudo_observations_and_models",
                      &tss::simulatePseudoObservations<TIME_TYPE, double>,
                      py::arg("bodies"), py::arg("observed_bodies"),
                      py::arg("central_bodies"), py::arg("initial_time"),
                      py::arg("final_time"), py::arg("time_step"),
get_docstring("create_pseudo_observations_and_models").c_str());

                m.def("set_existing_observations",
                      &tss::setExistingObservations<double, TIME_TYPE>,
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
R"doc(Function to compute the azimuth angle, elevation angle and range at a ground station.

	Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
	convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
	takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
	as observed at that station   
	

	:param bodies:
		System of bodies that defines the full physical environment
		
	:param station_id:
		Identifier for the observing station, as a pair of strings: the body name and the station name.
		
	:param target_body:
		Name of body which is observed by ground station
		
	:param observation_times:
		List of times at which the ground station observations are to be analyzed
		
	:param is_station_transmitting:
		Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station. 
		This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target
		
	:return:
		Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range
)doc");

                m.def("compute_target_angles_and_range_vectors",
                      &tss::getTargetAnglesAndRangeVector, py::arg("bodies"),
                      py::arg("station_id"), py::arg("target_body"),
                      py::arg("observation_times"),
                      py::arg("is_station_transmitting"),
R"doc(Function to compute the azimuth angle, elevation angle and range at a ground station.

	Function to compute the azimuth angle, elevation angle and range at a ground station. This functions is provided as a function of
	convenience, to prevent users having to manually define the relevant settings for this often-needed functionality. This function
	takes an observing station and a target body as input, and provides the observed angles and current range (without correction for aberrations, with correction for light time)
	as observed at that station   
	

	:param bodies:
		System of bodies that defines the full physical environment
		
	:param station_id:
		Identifier for the observing station, as a pair of strings: the body name and the station name.
		
	:param target_body:
		Name of body which is observed by ground station
		
	:param observation_times:
		List of times at which the ground station observations are to be analyzed
		
	:param is_station_transmitting:
		Boolean defining whether the observation times define times at which the station is transmitting to, or receiving from, the ground station. 
		This has an impact on the whether the light-time is computed forward or backward in time from the ground station to the target
		
	:return:
		Dictionary with the required output. Key defines the observation time, the value is an array of size three containing entry 0 - elevation angle, entry 1 - azimuth angle, entry 2 - range
)doc");

                py::class_<tom::ObservationCollection<double, TIME_TYPE>,
                           std::shared_ptr<
                               tom::ObservationCollection<double, TIME_TYPE>>>(
                    m, "ObservationCollection",
R"doc(Class collecting all observations and associated data for use in an estimation.

	Class containing the full set of observations and associated data, typically for input into the estimation. When using simulated data,
	this class is instantiated via a call to the :func:`~tudatpy.numerical_simulation.estimation.simulate_observations` function. More information is provided
	on the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_
	
)doc")
                    .def(py::init<std::vector<std::shared_ptr<
                             tom::SingleObservationSet<double, TIME_TYPE>>>>(),
                         py::arg("observation_sets"))
                    .def_property_readonly(
                        "concatenated_times",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getConcatenatedTimeVector,
R"doc(Vector containing concatenated observation times. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
	)doc")
                    .def_property_readonly(
                        "concatenated_observations",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getObservationVector,
R"doc(Vector containing concatenated observable values. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order
	)doc")
                    .def_property_readonly(
                        "concatenated_link_definition_ids",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getConcatenatedLinkEndIds,
R"doc(Vector containing concatenated indices identifying the link ends. Each set of link ends is assigned a unique integer identifier (for a given instance of this class). The definition of a given integer identifier with the link ends is given by this class' :func:`link_definition_ids` function. See `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_ for details on storage order of the present vector.
	)doc")
                    .def_property_readonly(
                        "link_definition_ids",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getInverseLinkEndIdentifierMap,
R"doc(Dictionaty mapping a link end integer identifier to the specific link ends
	)doc")
                    .def_property_readonly(
                        "observable_type_start_index_and_size",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getObservationTypeStartAndSize,
R"doc(Dictionary defining per obervable type (dict key), the index in the full observation vector (:func:`concatenated_observations`) where the given observable type starts, and the number of subsequent entries in this vector containing a value of an observable of this type
	)doc")
                    .def_property_readonly(
                        "observation_set_start_index_and_size",
                        &tom::ObservationCollection<double, TIME_TYPE>::
                            getObservationSetStartAndSizePerLinkEndIndex,
R"doc(The nested dictionary/list returned by this property mirrors the structure of the :func:`sorted_observation_sets` property of this class. The present function provides the start index and size of the observables in the full observation vector that come from the correspoding `SingleObservationSet` in the :func:`sorted_observation_sets` Consequently, the present property returns a nested dictionary defining per obervable type, link end identifier, and `SingleObservationSet` index (for the given observable type and link end identifier), where the observables in the given `SingleObservationSet` starts, and the number of subsequent entries in this vector containing data from it.
	)doc")
                    .def_property_readonly(
                        "observation_vector_size",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getTotalObservableSize,
R"doc(Length of the total vector of observations
	)doc")
                    .def_property_readonly(
                        "sorted_observation_sets",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getSortedObservationSets,
R"doc(The nested dictionary/list contains the list of `SingleObservationSet` objects, in the same method as they are stored internally in the present class. Specifics on the storage order are given in the `user guide <https://docs.tudat.space/en/stable/_src_user_guide/state_estimation/observation_simulation.html#accessing-and-analyzing-the-observations>`_
	)doc")
                    .def_property_readonly(
                        "link_ends_per_observable_type",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getLinkEndsPerObservableType,
get_docstring("ObservationCollection.link_ends_per_observable_type").c_str())
                    .def_property_readonly(
                        "link_definitions_per_observable",
                        &tom::ObservationCollection<
                            double, TIME_TYPE>::getLinkDefinitionsPerObservable,
get_docstring("ObservationCollection.link_definitions_per_observable").c_str())
                    .def("get_link_definitions_for_observables",
                         &tom::ObservationCollection<double, TIME_TYPE>::
                             getLinkDefinitionsForSingleObservable,
                         py::arg("observable_type"),
get_docstring("ObservationCollection.get_link_definitions_for_observables").c_str())
                    .def("get_single_link_and_type_observations",
                         &tom::ObservationCollection<double, TIME_TYPE>::
                             getSingleLinkAndTypeObservationSets,
                         py::arg("observable_type"), py::arg("link_definition"),
R"doc(Function to get all observation sets for a given observable type and link definition.

	:param observable_type:
		Observable type of which observations are to be simulated.
	:param link_ends:
		Link ends for which observations are to be simulated.
	:return:
		List of observation sets for given observable type and link definition.
)doc");

                py::class_<tom::SingleObservationSet<double, TIME_TYPE>,
                           std::shared_ptr<
                               tom::SingleObservationSet<double, TIME_TYPE>>>(
                    m, "SingleObservationSet",
R"doc(Class collecting a single set of observations and associated data, of a given observable type, link ends, and ancilliary data.

)doc")
                    .def_property_readonly(
                        "observable_type",
                        &tom::SingleObservationSet<
                            double, TIME_TYPE>::getObservableType,
R"doc(Type of observable for which the object stores observations
	)doc")
                    .def_property_readonly(
                        "link_definition",
                        &tom::SingleObservationSet<double,
                                                   TIME_TYPE>::getLinkEnds,
R"doc(Definition of the link ends for which the object stores observations
	)doc")
                    .def_property_readonly(
                        "reference_link_end",
                        &tom::SingleObservationSet<
                            double, TIME_TYPE>::getReferenceLinkEnd,
R"doc(Reference link end for all stored observations
	)doc")
                    .def_property_readonly(
                        "list_of_observations",
                        &tom::SingleObservationSet<double,
                                                   TIME_TYPE>::getObservations,
R"doc(List of separate stored observations. Each entry of this list is a vector containing a single observation. In cases where the observation is single-valued (range, Doppler), the vector is size 1, but for multi-valued observations such as angular position, each vector in the list will have size >1
	)doc")
                    .def_property_readonly(
                        "observation_times",
                        &tom::SingleObservationSet<
                            double, TIME_TYPE>::getObservationTimes,
R"doc(Reference time for each of the observations in ``list_of_observations``
	)doc")
                    .def_property_readonly(
                        "concatenated_observations",
                        &tom::SingleObservationSet<
                            double, TIME_TYPE>::getObservationsVector,
R"doc(Concatenated vector of all stored observations
	)doc")
                    .def_property_readonly(
                        "observations_history",
                        &tom::SingleObservationSet<
                            double, TIME_TYPE>::getObservationsHistory,
R"doc(Dictionary of observations sorted by time. Created by making a dictionaty with ``observation_times`` as keys and  ``list_of_observations`` as values
	)doc")
                    .def_property_readonly(
                        "ancilliary_settings",
                        &tom::SingleObservationSet<
                            double, TIME_TYPE>::getAncilliarySettings,
R"doc(Ancilliary settings all stored observations
	)doc")
                    .def_property(
                        "weights_vector",
                        &tom::SingleObservationSet<double,
                                                   TIME_TYPE>::getWeightsVector,
                        &tom::SingleObservationSet<double,
                                                   TIME_TYPE>::setWeightsVector,
get_docstring("SingleObservationSet.weights_vector").c_str());


                m.def("single_observation_set",
                      &tss::singleObservationSetWithoutDependentVariables<
                          double, TIME_TYPE>,
                      py::arg("observable_type"), py::arg("link_definition"),
                      py::arg("observations"), py::arg("observation_times"),
                      py::arg("reference_link_end"),
                      py::arg("ancilliary_settings") = nullptr,
get_docstring("single_observation_set").c_str());

                /*!
                 *************** STATE TRANSITION INTERFACE ***************
                 */

                py::class_<
                    tp::CombinedStateTransitionAndSensitivityMatrixInterface,
                    std::shared_ptr<
                        tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
                    m, "CombinedStateTransitionAndSensitivityMatrixInterface",
R"doc(Class establishing an interface with the simulation's State Transition and Sensitivity Matrices.

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
R"doc(Function to get the concatenated state transition and sensitivity matrix at a given time.

	Function to get the concatenated state transition and sensitivity matrix at a given time.
	Entries corresponding to parameters which are not active at the current arc are omitted.
	

	:param time:
		Time at which concatenated state transition and sensitivity matrix are to be retrieved.
	:return:
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
	:param time:
		Time at which full concatenated state transition and sensitivity matrix are to be retrieved.
	:return:
		Full concatenated state transition and sensitivity matrix at a given time.
)doc")
                    .def_property_readonly(
                        "state_transition_size",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getStateTransitionMatrixSize,
R"doc(Size of the (square) state transition matrix.
	)doc")
                    .def_property_readonly(
                        "sensitivity_size",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getSensitivityMatrixSize,
R"doc(Number of columns in the sensitivity matrix.
	)doc")
                    .def_property_readonly(
                        "full_parameter_size",
                        &tp::
                            CombinedStateTransitionAndSensitivityMatrixInterface::
                                getFullParameterVectorSize,
R"doc(Full amount of parameters w.r.t. which partials have been set up via State Transition and Sensitivity Matrices.
	)doc");

                /*!
                 *************** COVARIANCE ***************
                 */

                m.def("propagate_covariance_rsw_split_output",
                      &tp::propagateCovarianceVectorsRsw,
                      py::arg("covariance_output"), py::arg("estimator"),
                      py::arg("output_times"),
get_docstring("propagate_covariance_rsw_split_output").c_str());


                m.def("propagate_formal_errors_rsw_split_output",
                      &tp::propagateFormalErrorVectorsRsw,
                      py::arg("covariance_output"), py::arg("estimator"),
                      py::arg("output_times"),
get_docstring("propagate_formal_errors_rsw_split_output").c_str());

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
get_docstring("propagate_covariance_split_output").c_str());

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
R"doc(Function to propagate system covariance through time.

	Function to propagate the covariance of a given system through time.
	The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
	

	:param initial_covariance:
		System covariance matrix (symmetric and positive semi-definite) at initial time.
		Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
		
	:param state_transition_interface:
		Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
		
	:param output_times:
		Times at which the propagated covariance matrix shall be reported.
		Note that this argument has no impact on the integration time-steps of the covariance propagation,
		which always adheres to the integrator settings that the `state_transition_interface` links to.
		Output times which do not coincide with integration time steps are calculated via interpolation.
		
	:return:
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
R"doc(Function to propagate system formal errors through time.

	Function to propagate the formal errors of a given system through time.
	Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
	The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
	

	:param initial_covariance:
		System covariance matrix (symmetric and positive semi-definite) at initial time.
		Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
		
	:param state_transition_interface:
		Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
		
	:param output_times:
		Times at which the propagated covariance matrix shall be reported.
		Note that this argument has no impact on the integration time-steps of the covariance propagation,
		which always adheres to the integrator settings that the `state_transition_interface` links to.
		Output times which do not coincide with integration time steps are calculated via interpolation.
		
	:return:
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
R"doc(Function to propagate system formal errors through time.

	Function to propagate the formal errors of a given system through time.
	Note that in practice the entire covariance matrix is propagated, but only the formal errors (variances) are reported at the output times.
	The system dynamics and numerical settings of the propagation are prescribed by the `state_transition_interface` parameter.
	

	:param initial_covariance:
		System covariance matrix (symmetric and positive semi-definite) at initial time.
		Dimensions have to be consistent with estimatable parameters in the system (specified by `state_transition_interface`)
		
	:param state_transition_interface:
		Interface to the variational equations of the system dynamics, handling the propagation of the covariance matrix through time.
		
	:param output_times:
		Times at which the propagated covariance matrix shall be reported.
		Note that this argument has no impact on the integration time-steps of the covariance propagation,
		which always adheres to the integrator settings that the `state_transition_interface` links to.
		Output times which do not coincide with integration time steps are calculated via interpolation.
		
	:return:
		Dictionary reporting the propagated formal errors at each output time.
)doc");


                /*!
                 *************** ESTIMATION ***************
                 */

                py::class_<tss::EstimationConvergenceChecker,
                           std::shared_ptr<tss::EstimationConvergenceChecker>>(
                    m, "EstimationConvergenceChecker",
R"doc(Class defining the convergence criteria for an estimation.

	Class defining the convergence criteria for an estimation.
	The user typically creates instances of this class via the :func:`~tudatpy.numerical_simulation.estimation.estimation_convergence_checker` factory function.
	
)doc");

                m.def("estimation_convergence_checker",
                      &tss::estimationConvergenceChecker,
                      py::arg("maximum_iterations") = 5,
                      py::arg("minimum_residual_change") = 0.0,
                      py::arg("minimum_residual") = 0.0,
                      py::arg("number_of_iterations_without_improvement") = 2,
R"doc(Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object.

	Function for creating an :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` object, which is required for defining the convergence criteria of an estimation.
	

	:param maximum_iterations:
		Maximum number of allowed iterations for estimation.
	:param minimum_residual_change:
		Minimum required change in residual between two iterations.
	:param minimum_residual:
		Minimum value of observation residual below which estimation is converged.
	:param number_of_iterations_without_improvement:
		Number of iterations without reduction of residual.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationConvergenceChecker` class, defining the convergence criteria for an estimation.
)doc");


                py::class_<tss::CovarianceAnalysisInput<double, TIME_TYPE>,
                           std::shared_ptr<tss::CovarianceAnalysisInput<
                               double, TIME_TYPE>>>(
                    m, "CovarianceAnalysisInput",
R"doc(Class for defining all specific inputs to a covariance analysis.

)doc")
                    .def(
                        py::init<
                            const std::shared_ptr<
                                tom::ObservationCollection<double, TIME_TYPE>>&,
                            const Eigen::MatrixXd, const Eigen::MatrixXd>(),
                        py::arg("observations_and_times"),
                        py::arg("inverse_apriori_covariance") =
                            Eigen::MatrixXd::Zero(0, 0),
                        py::arg("consider_covariance") =
                            Eigen::MatrixXd::Zero(0, 0),
R"doc(Class constructor.

	Constructor through which the user can create instances of this class. Note that the weight are all initiated as 1.0, and the default settings of ``define_covariance_settings`` are used.
	

	:param observations_and_times:
		Total data structure of observations and associated times/link ends/type/etc.
	:param inverse_apriori_covariance:
		A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation.CovarianceAnalysisInput` class, defining the data and other settings to be used for the covariance analysis.
)doc")
                    .def("set_constant_weight",
                         &tss::CovarianceAnalysisInput<
                             double, TIME_TYPE>::setConstantWeightsMatrix,
                         py::arg("weight"),
R"doc(Function to set a constant weight matrix for all observables.

	Function to set a constant weight matrix for all observables.
	The weights are applied to all observations managed by the given PodInput object.
	

	:param constant_weight:
		Constant weight factor that is to be applied to all observations.
	:return:
		Function modifies the object in-place.
)doc")
                    .def("set_constant_single_observable_weight",
                         &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                             setConstantSingleObservableWeights,
                         py::arg("observable_type"), py::arg("weight"),
get_docstring("CovarianceAnalysisInput.set_constant_single_observable_weight").c_str())
                    .def("set_constant_single_observable_vector_weight",
                         &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                             setConstantSingleObservableVectorWeights,
                         py::arg("observable_type"), py::arg("weight"),
get_docstring("CovarianceAnalysisInput.set_constant_single_observable_vector_weight").c_str())
                    .def("set_constant_single_observable_and_link_end_weight",
                         &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                             setConstantSingleObservableAndLinkEndsWeights,
                         py::arg("observable_type"), py::arg("link_ends"),
                         py::arg("weight"),
get_docstring("CovarianceAnalysisInput.set_constant_single_observable_and_link_end_weight").c_str())
                    .def(
                        "set_constant_single_observable_and_link_end_vector_"
                        "weight",
                        &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                            setConstantSingleObservableAndLinkEndsVectorWeights,
                        py::arg("observable_type"), py::arg("link_ends"),
                        py::arg("weight"),
get_docstring("CovarianceAnalysisInput.set_constant_single_observable_and_link_end_vector_weight").c_str())
                    .def(
                        "set_total_single_observable_and_link_end_vector_"
                        "weight",
                        &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                            setTabulatedSingleObservableAndLinkEndsWeights,
                        py::arg("observable_type"), py::arg("link_ends"),
                        py::arg("weight_vector"),
get_docstring("CovarianceAnalysisInput.set_total_single_observable_and_link_end_vector_weight").c_str())
                    .def("set_constant_weight_per_observable",
                         &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                             setConstantPerObservableWeightsMatrix,
                         py::arg("weight_per_observable"),
R"doc(Function to set a constant weight matrix for a given type of observable.

	Function to set a constant weight matrix for a given type of observable.
	The weights are applied to all observations of the observable type specified by the `weight_per_observable` parameter.
	

	:param constant_weight:
		Constant weight factor that is to be applied to all observations.
	:return:
		Function modifies the object in-place.
)doc")
                    .def("set_constant_vector_weight_per_observable",
                         &tss::CovarianceAnalysisInput<double, TIME_TYPE>::
                             setConstantPerObservableVectorWeightsMatrix,
                         py::arg("weight_per_observable"),
get_docstring("CovarianceAnalysisInput.set_constant_vector_weight_per_observable").c_str())
                    .def("define_covariance_settings",
                         &tss::CovarianceAnalysisInput<
                             double, TIME_TYPE>::defineCovarianceSettings,
                         py::arg("reintegrate_equations_on_first_iteration") =
                             true,
                         py::arg("reintegrate_variational_equations") = true,
                         py::arg("save_design_matrix") = true,
                         py::arg("print_output_to_terminal") = true,
                         py::arg("limit_condition_number_for_warning") = 1.0E8,
R"doc(Function to define specific settings for covariance analysis process

	Function to define specific settings for covariance analysis process
	

	:param reintegrate_equations:
		Boolean denoting whether the dynamics and variational equations are to be reintegrated
		or if existing values are to be used to perform first iteration.
		
	:param reintegrate_variational_equations:
		Boolean denoting whether the variational equations are to be reintegrated during estimation 
		(if this is set to False, and ``reintegrate_equations`` to true, only the dynamics are re-integrated)
		
	:param save_design_matrix:
		Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
		:math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.
		
	:param print_output_to_terminal:
		Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.
		
	:return:
		Function modifies the object in-place.
)doc")
                    .def_property(
                        "weight_matrix_diagonal",
                        &tss::CovarianceAnalysisInput<
                            double, TIME_TYPE>::getWeightsMatrixDiagonals,
                        &tss::CovarianceAnalysisInput<
                            double, TIME_TYPE>::setWeightsMatrixDiagonals,
R"doc(Complete diagonal of the weights matrix that is to be used
	)doc");

                py::class_<
                    tss::EstimationInput<double, TIME_TYPE>,
                    std::shared_ptr<tss::EstimationInput<double, TIME_TYPE>>,
                    tss::CovarianceAnalysisInput<double, TIME_TYPE>>(
                    m, "EstimationInput",
R"doc(Class for defining all inputs to the estimation.

)doc")
                    .def(
                        py::init<
                            const std::shared_ptr<
                                tom::ObservationCollection<double, TIME_TYPE>>&,
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
R"doc(Class constructor.

	Constructor through which the user can create instances of this class.
	

	:param observations_and_times:
		Total data structure of observations and associated times/link ends/type/etc.
	:param inverse_apriori_covariance:
		A priori covariance matrix (unnormalized) of estimated parameters. This should be either a size 0x0 matrix (no a priori information), or a square matrix with the same size as the number of parameters that are considered
	:param convergence_checker:
		Object defining when the estimation is converged.
	:return:
		Instance of the :class:`~tudatpy.numerical_simulation.estimation.EstimationInput` class, defining the data and other settings to be used for the estimation.
)doc")
                    .def(
                        "define_estimation_settings",
                        &tss::EstimationInput<
                            double, TIME_TYPE>::defineEstimationSettings,
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
R"doc(Function to define specific settings for the estimation process

	Function to define specific settings for covariance analysis process
	

	:param reintegrate_equations_on_first_iteration:
		Boolean denoting whether the dynamics and variational equations are to be reintegrated
		or if existing values are to be used to perform first iteration.
		
	:param reintegrate_variational_equations:
		Boolean denoting whether the variational equations are to be reintegrated during estimation 
		(if this is set to False, and ``reintegrate_equations_on_first_iteration`` to true, only the dynamics are re-integrated)
		
	:param save_design_matrix:
		Boolean denoting whether to save the partials matrix (also called design matrix) :math:`\mathbf{H}` in the output. Setting this to false makes the
		:math:`\mathbf{H}` matrix unavailable to the user, with the advantage of lower RAM usage.
		
	:param print_output_to_terminal:
		Boolean denoting whether to print covariance-analysis-specific output to the terminal when running the estimation.
		
	:param save_residuals_and_parameters_per_iteration:
		Boolean denoting whether the residuals and parameters from the each iteration are to be saved.
		
	:param save_state_history_per_iteration:
		Boolean denoting whether the state history and dependent variables are to be saved on each iteration.
		
	:return:
		Function modifies the object in-place.
)doc");

                m.attr("PodInput") = m.attr("EstimationInput");


                py::class_<tss::CovarianceAnalysisOutput<double, TIME_TYPE>,
                           std::shared_ptr<tss::CovarianceAnalysisOutput<
                               double, TIME_TYPE>>>(
                    m, "CovarianceAnalysisOutput",
R"doc(Class collecting all outputs from the covariance analysis process.

)doc")
                    .def_property_readonly(
                        "inverse_covariance",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getUnnormalizedInverseCovarianceMatrix,
R"doc((Unnormalized) inverse estimation covariance matrix :math:`\mathbf{P}^{-1}`.
	)doc")
                    .def_property_readonly(
                        "covariance",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getUnnormalizedCovarianceMatrix,
R"doc((Unnormalized) estimation covariance matrix :math:`\mathbf{P}`.
	)doc")
                    .def_property_readonly(
                        "inverse_normalized_covariance",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getNormalizedInverseCovarianceMatrix,
R"doc(Normalized inverse estimation covariance matrix :math:`\mathbf{\tilde{P}}^{-1}`.
	)doc")
                    .def_property_readonly(
                        "normalized_covariance",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getNormalizedCovarianceMatrix,
R"doc(Normalized estimation covariance matrix :math:`\mathbf{\tilde{P}}`.
	)doc")
                    .def_property_readonly(
                        "formal_errors",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getFormalErrorVector,
R"doc(Formal error vector :math:`\boldsymbol{\sigma}` of the estimation result (e.g. square root of diagonal entries of covariance)s
	)doc")
                    .def_property_readonly(
                        "correlations",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getCorrelationMatrix,
R"doc(Correlation matrix of the estimation result. Entry :math:`i,j` is equal to :math:`P_{i,j}/(\sigma_{i}\sigma_{j})`
	)doc")
                    .def_property_readonly(
                        "design_matrix",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getUnnormalizedDesignMatrix,
R"doc(Matrix of unnormalized partial derivatives :math:`\mathbf{H}=\frac{\partial\mathbf{h}}{\partial\mathbf{p}}`.
	)doc")
                    .def_property_readonly(
                        "normalized_design_matrix",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getNormalizedDesignMatrix,
R"doc(Matrix of normalized partial derivatives :math:`\tilde{\mathbf{H}}`.
	)doc")
                    .def_property_readonly(
                        "weighted_design_matrix",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getUnnormalizedWeightedDesignMatrix,
R"doc(Matrix of weighted partial derivatives, equal to :math:`\mathbf{W}^{1/2}{\mathbf{H}}`
	)doc")
                    .def_property_readonly(
                        "weighted_normalized_design_matrix",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getNormalizedWeightedDesignMatrix,
R"doc(Matrix of weighted, normalized partial derivatives, equal to :math:`\mathbf{W}^{1/2}\tilde{\mathbf{H}}`
	)doc")
                    .def_property_readonly(
                        "consider_covariance_contribution",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getConsiderCovarianceContribution,
get_docstring("CovarianceAnalysisOutput.consider_covariance_contribution").c_str())
                    .def_property_readonly(
                        "normalized_covariance_with_consider_parameters",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getNormalizedCovarianceWithConsiderParameters,
get_docstring("CovarianceAnalysisOutput.normalized_covariance_with_consider_parameters").c_str())
                    .def_property_readonly(
                        "unnormalized_covariance_with_consider_parameters",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getUnnormalizedCovarianceWithConsiderParameters,
get_docstring("CovarianceAnalysisOutput.unnormalized_covariance_with_consider_parameters").c_str())
                    .def_property_readonly(
                        "normalized_design_matrix_consider_parameters",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            getNormalizedDesignMatrixConsiderParameters,
get_docstring("CovarianceAnalysisOutput.normalized_design_matrix_consider_parameters").c_str())
                    .def_property_readonly(
                        "consider_normalization_factors",
                        &tss::CovarianceAnalysisOutput<
                            double, TIME_TYPE>::getConsiderNormalizationFactors,
get_docstring("CovarianceAnalysisOutput.consider_normalization_factors").c_str())
                    .def_readonly(
                        "normalization_terms",
                        &tss::CovarianceAnalysisOutput<double, TIME_TYPE>::
                            designMatrixTransformationDiagonal_,
R"doc(Vector of normalization terms used for covariance and design matrix
	)doc");


                py::class_<
                    tss::EstimationOutput<double, TIME_TYPE>,
                    std::shared_ptr<tss::EstimationOutput<double, TIME_TYPE>>,
                    tss::CovarianceAnalysisOutput<double, TIME_TYPE>>(
                    m, "EstimationOutput",
R"doc(Class collecting all outputs from the iterative estimation process.

)doc")
                    .def_property_readonly(
                        "residual_history",
                        &tss::EstimationOutput<
                            double, TIME_TYPE>::getResidualHistoryMatrix,
R"doc(Residual vectors, concatenated per iteration into a matrix; the :math:`i^{th}` column has the residuals from the :math:`i^{th}` iteration.
	)doc")
                    .def_property_readonly(
                        "parameter_history",
                        &tss::EstimationOutput<
                            double, TIME_TYPE>::getParameterHistoryMatrix,
R"doc(Parameter vectors, concatenated per iteration into a matrix. Column 0 contains pre-estimation values. The :math:`(i+1)^{th}` column has the residuals from the :math:`i^{th}` iteration.
	)doc")
                    .def_property_readonly(
                        "simulation_results_per_iteration",
                        &tss::EstimationOutput<double,
                                               TIME_TYPE>::getSimulationResults,
R"doc(List of complete numerical propagation results, with the :math:`i^{th}` entry of thee list thee results of the :math:`i^{th}` propagation
	)doc")
                    .def_readonly(
                        "final_residuals",
                        &tss::EstimationOutput<double, TIME_TYPE>::residuals_,
R"doc(Vector of post-fit observation residuals.
	)doc")
                    .def_readonly(
                        "final_parameters",
                        &tss::EstimationOutput<double,
                                               TIME_TYPE>::parameterEstimate_,
get_docstring("EstimationOutput.final_parameters").c_str())
                    .def_readonly(
                        "best_iteration",
                        &tss::EstimationOutput<double,
                                               TIME_TYPE>::bestIteration_,
get_docstring("EstimationOutput.best_iteration").c_str());

                m.attr("PodOutput") = m.attr("EstimationOutput");
            }

        }  // namespace estimation
    }  // namespace numerical_simulation
}  // namespace tudatpy
