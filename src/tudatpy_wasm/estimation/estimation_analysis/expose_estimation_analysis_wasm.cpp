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

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;
namespace tep = tudat::estimatable_parameters;

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
    class_<tss::OrbitDeterminationManager<double, double>>(
        "estimation_estimation_analysis_Estimator")
        .smart_ptr<std::shared_ptr<tss::OrbitDeterminationManager<double, double>>>(
            "shared_ptr_Estimator")
        .constructor<const tss::SystemOfBodies&,
                     const std::shared_ptr<tep::EstimatableParameterSet<double>>,
                     const std::vector<std::shared_ptr<tss::ObservationModelSettings>>&,
                     const std::shared_ptr<tudat::propagators::PropagatorSettings<double>>,
                     const bool>()
        .function("estimateParameters",
            &tss::OrbitDeterminationManager<double, double>::estimateParameters)
        .function("computeCovariance",
            &tss::OrbitDeterminationManager<double, double>::computeCovariance)
        .function("getParametersToEstimate",
            &tss::OrbitDeterminationManager<double, double>::getParametersToEstimate)
        .function("getBodies",
            &tss::OrbitDeterminationManager<double, double>::getBodies)
        .function("getStateTransitionAndSensitivityMatrixInterface",
            &tss::OrbitDeterminationManager<double, double>::getStateTransitionAndSensitivityMatrixInterface);
}

#endif
