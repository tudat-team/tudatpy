/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "../../wasm_module.h"
#include "../../eigen_wasm.h"
#include "../../stl_wasm.h"
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/estimation_setup/orbitDeterminationManager.h>
#include <tudat/astro/orbit_determination/podInputOutputTypes.h>
#include <tudat/astro/propagators/propagateCovariance.h>
#include <tudat/basics/utilities.h>
#include <tudat/astro/reference_frames/referenceFrameTransformations.h>

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace tep = tudat::estimatable_parameters;
namespace tp = tudat::propagators;
namespace trf = tudat::reference_frames;

// ============================================================================
// Covariance Propagation Helper Functions (matching Python bindings)
// ============================================================================

namespace tudat {
namespace propagators {

// Propagate covariance and convert to RSW frame
std::map<double, Eigen::MatrixXd> propagateCovarianceRsw(
    const Eigen::MatrixXd initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>> orbitDeterminationManager,
    const std::vector<double> evaluationTimes)
{
    std::map<double, Eigen::MatrixXd> propagatedCovariance;
    tp::propagateCovariance(
        propagatedCovariance,
        initialCovariance,
        orbitDeterminationManager->getStateTransitionAndSensitivityMatrixInterface(),
        evaluationTimes);

    tss::SystemOfBodies bodies = orbitDeterminationManager->getBodies();

    std::shared_ptr<tep::EstimatableParameterSet<double>> parameterSet =
        orbitDeterminationManager->getParametersToEstimate();

    std::map<int, std::shared_ptr<tep::EstimatableParameter<Eigen::Matrix<double, Eigen::Dynamic, 1>>>>
        initialStates = parameterSet->getInitialStateParameters();

    std::map<std::pair<std::string, std::string>, std::vector<int>> transformationList;
    for (auto it : initialStates)
    {
        if (std::dynamic_pointer_cast<tep::InitialTranslationalStateParameter<double>>(it.second))
        {
            std::shared_ptr<tep::InitialTranslationalStateParameter<double>> currentInitialState =
                std::dynamic_pointer_cast<tep::InitialTranslationalStateParameter<double>>(it.second);
            transformationList[std::make_pair(
                currentInitialState->getParameterName().second.first,
                currentInitialState->getCentralBody())].push_back(it.first);
        }
        else if (std::dynamic_pointer_cast<tep::ArcWiseInitialTranslationalStateParameter<double>>(it.second))
        {
            throw std::runtime_error(
                "Error, multi-arc not yet supported in automatic covariance conversion");
        }
    }

    Eigen::Matrix3d currentInertialToRswPosition = Eigen::Matrix3d::Zero();
    Eigen::Matrix6d currentInertialToRswState = Eigen::Matrix6d::Zero();
    Eigen::MatrixXd currentFullInertialToRswState = Eigen::MatrixXd::Zero(6, 6);

    std::map<double, Eigen::MatrixXd> propagatedRswCovariance;
    for (auto it : propagatedCovariance)
    {
        double currentTime = static_cast<double>(it.first);
        Eigen::MatrixXd currentCovariance = it.second;
        currentFullInertialToRswState.setZero();

        for (auto it_body : transformationList)
        {
            Eigen::Vector6d relativeState =
                bodies.getBody(it_body.first.first)->getStateInBaseFrameFromEphemeris(currentTime) -
                bodies.getBody(it_body.first.second)->getStateInBaseFrameFromEphemeris(currentTime);
            currentInertialToRswPosition =
                trf::getInertialToRswSatelliteCenteredFrameRotationMatrix(relativeState);
            currentInertialToRswState.block(0, 0, 3, 3) = currentInertialToRswPosition;
            currentInertialToRswState.block(3, 3, 3, 3) = currentInertialToRswPosition;
            for (unsigned int j = 0; j < it_body.second.size(); j++)
            {
                int currentStartIndex = it_body.second.at(j);
                currentFullInertialToRswState.block(currentStartIndex, currentStartIndex, 6, 6) =
                    currentInertialToRswState;
            }
        }
        propagatedRswCovariance[currentTime] = currentFullInertialToRswState * currentCovariance *
            currentFullInertialToRswState.transpose();
    }
    return propagatedRswCovariance;
}

// Split output version of RSW covariance propagation
std::pair<std::vector<double>, std::vector<Eigen::MatrixXd>> propagateCovarianceVectorsRsw(
    const Eigen::MatrixXd initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>> orbitDeterminationManager,
    const std::vector<double> evaluationTimes)
{
    std::map<double, Eigen::MatrixXd> propagatedRswCovariance =
        propagateCovarianceRsw(initialCovariance, orbitDeterminationManager, evaluationTimes);

    return std::make_pair(
        utilities::createVectorFromMapKeys(propagatedRswCovariance),
        utilities::createVectorFromMapValues(propagatedRswCovariance));
}

// Propagate formal errors in RSW frame
std::map<double, Eigen::VectorXd> propagateFormalErrorsRsw(
    const Eigen::MatrixXd initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>> orbitDeterminationManager,
    const std::vector<double> evaluationTimes)
{
    std::map<double, Eigen::MatrixXd> propagatedCovariance;
    std::map<double, Eigen::VectorXd> propagatedFormalErrors;

    propagatedCovariance =
        propagateCovarianceRsw(initialCovariance, orbitDeterminationManager, evaluationTimes);
    tp::convertCovarianceHistoryToFormalErrorHistory(propagatedFormalErrors, propagatedCovariance);

    return propagatedFormalErrors;
}

// Split output version of RSW formal errors propagation
std::pair<std::vector<double>, std::vector<Eigen::VectorXd>> propagateFormalErrorVectorsRsw(
    const Eigen::MatrixXd initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>> orbitDeterminationManager,
    const std::vector<double> evaluationTimes)
{
    std::map<double, Eigen::VectorXd> propagatedFormalErrors = propagateFormalErrorsRsw(
        initialCovariance, orbitDeterminationManager, evaluationTimes);
    return std::make_pair(
        utilities::createVectorFromMapKeys(propagatedFormalErrors),
        utilities::createVectorFromMapValues(propagatedFormalErrors));
}

// Propagate covariance (split output)
std::pair<std::vector<double>, std::vector<Eigen::MatrixXd>> propagateCovarianceVectors(
    const Eigen::MatrixXd initialCovariance,
    const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface> stateTransitionInterface,
    const std::vector<double> evaluationTimes)
{
    std::map<double, Eigen::MatrixXd> propagatedCovariance;
    tp::propagateCovariance(
        propagatedCovariance, initialCovariance, stateTransitionInterface, evaluationTimes);
    return std::make_pair(
        utilities::createVectorFromMapKeys(propagatedCovariance),
        utilities::createVectorFromMapValues(propagatedCovariance));
}

// Propagate formal errors (split output)
std::pair<std::vector<double>, std::vector<Eigen::VectorXd>> propagateFormalErrorVectors(
    const Eigen::MatrixXd initialCovariance,
    const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface> stateTransitionInterface,
    const std::vector<double> evaluationTimes)
{
    std::map<double, Eigen::VectorXd> propagatedFormalErrors;
    tp::propagateFormalErrors(
        propagatedFormalErrors, initialCovariance, stateTransitionInterface, evaluationTimes);
    return std::make_pair(
        utilities::createVectorFromMapKeys(propagatedFormalErrors),
        utilities::createVectorFromMapValues(propagatedFormalErrors));
}

} // namespace propagators
} // namespace tudat

// ============================================================================
// WASM Wrapper Functions (convert Eigen types to/from JS-friendly wrappers)
// ============================================================================

namespace tudatpy_wasm {

using namespace tudatpy_wasm;

// Wrapper for propagate_covariance_rsw returning JS-friendly types
std::map<double, MatrixXdWrapper> propagateCovarianceRswWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>>& orbitDeterminationManager,
    const std::vector<double>& evaluationTimes)
{
    auto result = tp::propagateCovarianceRsw(
        initialCovariance.data, orbitDeterminationManager, evaluationTimes);

    std::map<double, MatrixXdWrapper> wrappedResult;
    for (const auto& entry : result) {
        wrappedResult[entry.first] = MatrixXdWrapper(entry.second);
    }
    return wrappedResult;
}

// Wrapper for propagate_covariance_rsw_split_output
std::pair<std::vector<double>, std::vector<MatrixXdWrapper>> propagateCovarianceRswSplitOutputWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>>& orbitDeterminationManager,
    const std::vector<double>& evaluationTimes)
{
    auto result = tp::propagateCovarianceVectorsRsw(
        initialCovariance.data, orbitDeterminationManager, evaluationTimes);

    std::vector<MatrixXdWrapper> wrappedMatrices;
    for (const auto& mat : result.second) {
        wrappedMatrices.push_back(MatrixXdWrapper(mat));
    }
    return std::make_pair(result.first, wrappedMatrices);
}

// Wrapper for propagate_formal_errors_rsw
std::map<double, VectorXdWrapper> propagateFormalErrorsRswWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>>& orbitDeterminationManager,
    const std::vector<double>& evaluationTimes)
{
    auto result = tp::propagateFormalErrorsRsw(
        initialCovariance.data, orbitDeterminationManager, evaluationTimes);

    std::map<double, VectorXdWrapper> wrappedResult;
    for (const auto& entry : result) {
        wrappedResult[entry.first] = VectorXdWrapper(entry.second);
    }
    return wrappedResult;
}

// Wrapper for propagate_formal_errors_rsw_split_output
std::pair<std::vector<double>, std::vector<VectorXdWrapper>> propagateFormalErrorsRswSplitOutputWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tss::OrbitDeterminationManager<double, double>>& orbitDeterminationManager,
    const std::vector<double>& evaluationTimes)
{
    auto result = tp::propagateFormalErrorVectorsRsw(
        initialCovariance.data, orbitDeterminationManager, evaluationTimes);

    std::vector<VectorXdWrapper> wrappedVectors;
    for (const auto& vec : result.second) {
        wrappedVectors.push_back(VectorXdWrapper(vec));
    }
    return std::make_pair(result.first, wrappedVectors);
}

// Wrapper for propagate_covariance (using state transition interface)
std::map<double, MatrixXdWrapper> propagateCovarianceWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>& stateTransitionInterface,
    const std::vector<double>& evaluationTimes)
{
    std::map<double, Eigen::MatrixXd> result;
    tp::propagateCovariance(result, initialCovariance.data, stateTransitionInterface, evaluationTimes);

    std::map<double, MatrixXdWrapper> wrappedResult;
    for (const auto& entry : result) {
        wrappedResult[entry.first] = MatrixXdWrapper(entry.second);
    }
    return wrappedResult;
}

// Wrapper for propagate_covariance_split_output
std::pair<std::vector<double>, std::vector<MatrixXdWrapper>> propagateCovarianceSplitOutputWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>& stateTransitionInterface,
    const std::vector<double>& evaluationTimes)
{
    auto result = tp::propagateCovarianceVectors(
        initialCovariance.data, stateTransitionInterface, evaluationTimes);

    std::vector<MatrixXdWrapper> wrappedMatrices;
    for (const auto& mat : result.second) {
        wrappedMatrices.push_back(MatrixXdWrapper(mat));
    }
    return std::make_pair(result.first, wrappedMatrices);
}

// Wrapper for propagate_formal_errors
std::map<double, VectorXdWrapper> propagateFormalErrorsWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>& stateTransitionInterface,
    const std::vector<double>& evaluationTimes)
{
    std::map<double, Eigen::VectorXd> result;
    tp::propagateFormalErrors(result, initialCovariance.data, stateTransitionInterface, evaluationTimes);

    std::map<double, VectorXdWrapper> wrappedResult;
    for (const auto& entry : result) {
        wrappedResult[entry.first] = VectorXdWrapper(entry.second);
    }
    return wrappedResult;
}

// Wrapper for propagate_formal_errors_split_output
std::pair<std::vector<double>, std::vector<VectorXdWrapper>> propagateFormalErrorsSplitOutputWrapper(
    const MatrixXdWrapper& initialCovariance,
    const std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>& stateTransitionInterface,
    const std::vector<double>& evaluationTimes)
{
    auto result = tp::propagateFormalErrorVectors(
        initialCovariance.data, stateTransitionInterface, evaluationTimes);

    std::vector<VectorXdWrapper> wrappedVectors;
    for (const auto& vec : result.second) {
        wrappedVectors.push_back(VectorXdWrapper(vec));
    }
    return std::make_pair(result.first, wrappedVectors);
}

} // namespace tudatpy_wasm

WASM_MODULE_PATH("estimation_estimation_analysis")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_estimation_analysis) {
    using namespace emscripten;

    // EstimationConvergenceChecker class
    class_<tss::EstimationConvergenceChecker>("estimation_estimation_analysis_EstimationConvergenceChecker")
        .smart_ptr<std::shared_ptr<tss::EstimationConvergenceChecker>>(
            "shared_ptr_EstimationConvergenceChecker")
        .constructor<const unsigned int, const double, const double, const int>()
        .function("isEstimationConverged", &tss::EstimationConvergenceChecker::isEstimationConverged);

    // Factory function for convergence checker
    function("estimation_estimation_analysis_estimation_convergence_checker",
        &tss::estimationConvergenceChecker);

    // CovarianceAnalysisInput class
    class_<tss::CovarianceAnalysisInput<double, double>>(
        "estimation_estimation_analysis_CovarianceAnalysisInput")
        .smart_ptr<std::shared_ptr<tss::CovarianceAnalysisInput<double, double>>>(
            "shared_ptr_CovarianceAnalysisInput")
        .function("setConstantWeightsMatrix",
            &tss::CovarianceAnalysisInput<double, double>::setConstantWeightsMatrix)
        .function("defineCovarianceSettings",
            &tss::CovarianceAnalysisInput<double, double>::defineCovarianceSettings)
        .function("getInverseOfAprioriCovariance",
            select_overload<Eigen::MatrixXd()>(
                &tss::CovarianceAnalysisInput<double, double>::getInverseOfAprioriCovariance))
        .function("getReintegrateEquationsOnFirstIteration",
            &tss::CovarianceAnalysisInput<double, double>::getReintegrateEquationsOnFirstIteration)
        .function("getSaveDesignMatrix",
            &tss::CovarianceAnalysisInput<double, double>::getSaveDesignMatrix);

    // EstimationInput class
    class_<tss::EstimationInput<double, double>, base<tss::CovarianceAnalysisInput<double, double>>>(
        "estimation_estimation_analysis_EstimationInput")
        .smart_ptr<std::shared_ptr<tss::EstimationInput<double, double>>>(
            "shared_ptr_EstimationInput")
        .function("defineEstimationSettings",
            &tss::EstimationInput<double, double>::defineEstimationSettings)
        .function("getSaveResidualsAndParametersFromEachIteration",
            &tss::EstimationInput<double, double>::getSaveResidualsAndParametersFromEachIteration)
        .function("getConvergenceChecker",
            &tss::EstimationInput<double, double>::getConvergenceChecker)
        .function("setConvergenceChecker",
            &tss::EstimationInput<double, double>::setConvergenceChecker);

    // CovarianceAnalysisOutput struct
    class_<tss::CovarianceAnalysisOutput<double, double>>(
        "estimation_estimation_analysis_CovarianceAnalysisOutput")
        .smart_ptr<std::shared_ptr<tss::CovarianceAnalysisOutput<double, double>>>(
            "shared_ptr_CovarianceAnalysisOutput")
        .function("getNormalizationTerms",
            &tss::CovarianceAnalysisOutput<double, double>::getNormalizationTerms)
        .function("getNormalizedCovarianceMatrix",
            &tss::CovarianceAnalysisOutput<double, double>::getNormalizedCovarianceMatrix)
        .function("getUnnormalizedCovarianceMatrix",
            &tss::CovarianceAnalysisOutput<double, double>::getUnnormalizedCovarianceMatrix)
        .function("getFormalErrorVector",
            &tss::CovarianceAnalysisOutput<double, double>::getFormalErrorVector)
        .function("getCorrelationMatrix",
            &tss::CovarianceAnalysisOutput<double, double>::getCorrelationMatrix)
        .function("getUnnormalizedDesignMatrix",
            &tss::CovarianceAnalysisOutput<double, double>::getUnnormalizedDesignMatrix)
        .function("getNormalizedDesignMatrix",
            &tss::CovarianceAnalysisOutput<double, double>::getNormalizedDesignMatrix);

    // EstimationOutput struct
    class_<tss::EstimationOutput<double, double>, base<tss::CovarianceAnalysisOutput<double, double>>>(
        "estimation_estimation_analysis_EstimationOutput")
        .smart_ptr<std::shared_ptr<tss::EstimationOutput<double, double>>>(
            "shared_ptr_EstimationOutput")
        .function("getResidualHistoryMatrix",
            &tss::EstimationOutput<double, double>::getResidualHistoryMatrix)
        .function("getParameterHistoryMatrix",
            &tss::EstimationOutput<double, double>::getParameterHistoryMatrix)
        .property("parameterEstimate", &tss::EstimationOutput<double, double>::parameterEstimate_)
        .property("residuals", &tss::EstimationOutput<double, double>::residuals_)
        .property("bestIteration", &tss::EstimationOutput<double, double>::bestIteration_)
        .property("residualStandardDeviation", &tss::EstimationOutput<double, double>::residualStandardDeviation_);

    // OrbitDeterminationManager (Estimator) class
    // This is the main class for orbit determination / parameter estimation
    class_<tss::OrbitDeterminationManager<double, double>>(
        "estimation_estimation_analysis_Estimator")
        .smart_ptr<std::shared_ptr<tss::OrbitDeterminationManager<double, double>>>(
            "shared_ptr_Estimator")
        .constructor<const tss::SystemOfBodies&,
                     const std::shared_ptr<tep::EstimatableParameterSet<double>>,
                     const std::vector<std::shared_ptr<tss::ObservationModelSettings>>&,
                     const std::shared_ptr<tudat::propagators::PropagatorSettings<double>>,
                     const bool>()
        // Methods matching Python bindings
        .function("perform_estimation",
            &tss::OrbitDeterminationManager<double, double>::estimateParameters)
        .function("compute_covariance",
            &tss::OrbitDeterminationManager<double, double>::computeCovariance)
        // Properties matching Python bindings
        .function("getObservationSimulators",
            &tss::OrbitDeterminationManager<double, double>::getObservationSimulators)
        .function("getObservationManagers",
            &tss::OrbitDeterminationManager<double, double>::getObservationManagers)
        .function("getVariationalEquationsSolver",
            &tss::OrbitDeterminationManager<double, double>::getVariationalEquationsSolver)
        .function("getStateTransitionAndSensitivityMatrixInterface",
            &tss::OrbitDeterminationManager<double, double>::getStateTransitionAndSensitivityMatrixInterface)
        // Additional utility methods
        .function("getParametersToEstimate",
            &tss::OrbitDeterminationManager<double, double>::getParametersToEstimate)
        .function("getBodies",
            &tss::OrbitDeterminationManager<double, double>::getBodies);

    // CombinedStateTransitionAndSensitivityMatrixInterface class
    class_<tp::CombinedStateTransitionAndSensitivityMatrixInterface>(
        "estimation_estimation_analysis_CombinedStateTransitionAndSensitivityMatrixInterface")
        .smart_ptr<std::shared_ptr<tp::CombinedStateTransitionAndSensitivityMatrixInterface>>(
            "shared_ptr_CombinedStateTransitionAndSensitivityMatrixInterface");

    // ============================================================================
    // Covariance Propagation Functions (using WASM-friendly wrapper types)
    // ============================================================================

    // propagate_covariance_rsw - Map output (estimator-based, converts to RSW frame)
    function("estimation_estimation_analysis_propagate_covariance_rsw",
        &tudatpy_wasm::propagateCovarianceRswWrapper);

    // propagate_covariance_rsw_split_output - Tuple output
    function("estimation_estimation_analysis_propagate_covariance_rsw_split_output",
        &tudatpy_wasm::propagateCovarianceRswSplitOutputWrapper);

    // propagate_formal_errors_rsw - Map output
    function("estimation_estimation_analysis_propagate_formal_errors_rsw",
        &tudatpy_wasm::propagateFormalErrorsRswWrapper);

    // propagate_formal_errors_rsw_split_output - Tuple output
    function("estimation_estimation_analysis_propagate_formal_errors_rsw_split_output",
        &tudatpy_wasm::propagateFormalErrorsRswSplitOutputWrapper);

    // propagate_covariance - Map output (uses state transition interface directly)
    function("estimation_estimation_analysis_propagate_covariance",
        &tudatpy_wasm::propagateCovarianceWrapper);

    // propagate_covariance_split_output - Tuple output
    function("estimation_estimation_analysis_propagate_covariance_split_output",
        &tudatpy_wasm::propagateCovarianceSplitOutputWrapper);

    // propagate_formal_errors - Map output
    function("estimation_estimation_analysis_propagate_formal_errors",
        &tudatpy_wasm::propagateFormalErrorsWrapper);

    // propagate_formal_errors_split_output - Tuple output
    function("estimation_estimation_analysis_propagate_formal_errors_split_output",
        &tudatpy_wasm::propagateFormalErrorsSplitOutputWrapper);
}

#endif
