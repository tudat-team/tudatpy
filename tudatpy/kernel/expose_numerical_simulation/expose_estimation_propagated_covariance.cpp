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


}  // namespace tudat

namespace tudatpy {
    namespace numerical_simulation {
        namespace estimation {


            void expose_estimation_propagated_covariance(py::module& m) {

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
    }      // namespace numerical_simulation
}  // namespace tudatpy
