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

#include <tudat/simulation/propagation_setup/propagationResults.h>
#include <tudat/simulation/propagation_setup/dependentVariablesInterface.h>

namespace tp = tudat::propagators;

WASM_MODULE_PATH("dynamics_propagation")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation) {
    using namespace emscripten;

    // DependentVariablesInterface base class
    class_<tp::DependentVariablesInterface<double>>("dynamics_propagation_DependentVariablesInterface")
        .smart_ptr<std::shared_ptr<tp::DependentVariablesInterface<double>>>(
            "shared_ptr_DependentVariablesInterface");

    // SingleArcDependentVariablesInterface
    class_<tp::SingleArcDependentVariablesInterface<double>,
           base<tp::DependentVariablesInterface<double>>>(
        "dynamics_propagation_SingleArcDependentVariablesInterface")
        .smart_ptr<std::shared_ptr<tp::SingleArcDependentVariablesInterface<double>>>(
            "shared_ptr_SingleArcDependentVariablesInterface")
        .function("getDependentVariables", &tp::SingleArcDependentVariablesInterface<double>::getDependentVariables)
        .function("getDependentVariableIds", &tp::SingleArcDependentVariablesInterface<double>::getDependentVariableIds)
        .function("getDependentVariablesize", &tp::SingleArcDependentVariablesInterface<double>::getDependentVariablesize);

    // SimulationResults base class
    class_<tp::SimulationResults<double, double>>("dynamics_propagation_SimulationResults")
        .smart_ptr<std::shared_ptr<tp::SimulationResults<double, double>>>(
            "shared_ptr_SimulationResults")
        .function("getDependentVariablesInterface", &tp::SimulationResults<double, double>::getDependentVariablesInterface);

    // SingleArcSimulationResults
    class_<tp::SingleArcSimulationResults<double, double>,
           base<tp::SimulationResults<double, double>>>(
        "dynamics_propagation_SingleArcSimulationResults")
        .smart_ptr<std::shared_ptr<tp::SingleArcSimulationResults<double, double>>>(
            "shared_ptr_SingleArcSimulationResults")
        .function("getEquationsOfMotionNumericalSolution",
            &tp::SingleArcSimulationResults<double, double>::getEquationsOfMotionNumericalSolution)
        .function("getEquationsOfMotionNumericalSolutionRaw",
            &tp::SingleArcSimulationResults<double, double>::getEquationsOfMotionNumericalSolutionRaw)
        .function("getDependentVariableHistory",
            &tp::SingleArcSimulationResults<double, double>::getDependentVariableHistory)
        .function("getCumulativeComputationTimeHistory",
            &tp::SingleArcSimulationResults<double, double>::getCumulativeComputationTimeHistory)
        .function("getCumulativeNumberOfFunctionEvaluations",
            &tp::SingleArcSimulationResults<double, double>::getCumulativeNumberOfFunctionEvaluations)
        .function("getTotalComputationRuntime",
            &tp::SingleArcSimulationResults<double, double>::getTotalComputationRuntime)
        .function("getTotalNumberOfFunctionEvaluations",
            &tp::SingleArcSimulationResults<double, double>::getTotalNumberOfFunctionEvaluations)
        .function("getPropagationTerminationReason",
            &tp::SingleArcSimulationResults<double, double>::getPropagationTerminationReason)
        .function("integrationCompletedSuccessfully",
            &tp::SingleArcSimulationResults<double, double>::integrationCompletedSuccessfully)
        .function("getDependentVariableId",
            &tp::SingleArcSimulationResults<double, double>::getDependentVariableId)
        .function("getProcessedStateIds",
            &tp::SingleArcSimulationResults<double, double>::getProcessedStateIds)
        .function("getPropagatedStateIds",
            &tp::SingleArcSimulationResults<double, double>::getPropagatedStateIds)
        .function("getPropagatedStateSize",
            &tp::SingleArcSimulationResults<double, double>::getPropagatedStateSize)
        .function("getPropagationIsPerformed",
            &tp::SingleArcSimulationResults<double, double>::getPropagationIsPerformed)
        .function("getSolutionIsCleared",
            &tp::SingleArcSimulationResults<double, double>::getSolutionIsCleared)
        .function("getArcInitialAndFinalTime",
            &tp::SingleArcSimulationResults<double, double>::getArcInitialAndFinalTime)
        .function("getSingleArcDependentVariablesInterface",
            &tp::SingleArcSimulationResults<double, double>::getSingleArcDependentVariablesInterface);

    // SingleArcVariationalSimulationResults
    class_<tp::SingleArcVariationalSimulationResults<double, double>,
           base<tp::SimulationResults<double, double>>>(
        "dynamics_propagation_SingleArcVariationalSimulationResults")
        .smart_ptr<std::shared_ptr<tp::SingleArcVariationalSimulationResults<double, double>>>(
            "shared_ptr_SingleArcVariationalSimulationResults")
        .function("getStateTransitionSolution",
            &tp::SingleArcVariationalSimulationResults<double, double>::getStateTransitionSolution)
        .function("getSensitivitySolution",
            &tp::SingleArcVariationalSimulationResults<double, double>::getSensitivitySolution)
        .function("getDynamicsResults",
            &tp::SingleArcVariationalSimulationResults<double, double>::getDynamicsResults)
        .function("getStateTransitionMatrixSize",
            &tp::SingleArcVariationalSimulationResults<double, double>::getStateTransitionMatrixSize)
        .function("getSensitivityMatrixSize",
            &tp::SingleArcVariationalSimulationResults<double, double>::getSensitivityMatrixSize);

    // MultiArcSimulationResults
    class_<tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>,
           base<tp::SimulationResults<double, double>>>(
        "dynamics_propagation_MultiArcSimulationResults")
        .smart_ptr<std::shared_ptr<tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>>>(
            "shared_ptr_MultiArcSimulationResults")
        .function("getSingleArcResults",
            &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getSingleArcResults)
        .function("getArcStartTimes",
            &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getArcStartTimes)
        .function("getArcEndTimes",
            &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getArcEndTimes)
        .function("getPropagationIsPerformed",
            &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getPropagationIsPerformed)
        .function("getSolutionIsCleared",
            &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getSolutionIsCleared)
        .function("integrationCompletedSuccessfully",
            &tp::MultiArcSimulationResults<tp::SingleArcSimulationResults, double, double>::integrationCompletedSuccessfully);

    // HybridArcSimulationResults
    class_<tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, double>,
           base<tp::SimulationResults<double, double>>>(
        "dynamics_propagation_HybridArcSimulationResults")
        .smart_ptr<std::shared_ptr<tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, double>>>(
            "shared_ptr_HybridArcSimulationResults")
        .function("getSingleArcResults",
            &tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getSingleArcResults)
        .function("getMultiArcResults",
            &tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, double>::getMultiArcResults)
        .function("integrationCompletedSuccessfully",
            &tp::HybridArcSimulationResults<tp::SingleArcSimulationResults, double, double>::integrationCompletedSuccessfully);
}

#endif
