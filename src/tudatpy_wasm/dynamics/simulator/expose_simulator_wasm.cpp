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

#include <tudat/simulation/propagation_setup/dynamicsSimulator.h>
#include <tudat/simulation/estimation_setup/createNumericalSimulator.h>

namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;

// Wrapper function to create dynamics simulator (template function with default arg)
std::shared_ptr<tp::DynamicsSimulator<double, double>> createDynamicsSimulatorWrapper(
    const tss::SystemOfBodies& bodies,
    const std::shared_ptr<tp::PropagatorSettings<double>> propagatorSettings,
    const bool areEquationsOfMotionToBeIntegrated)
{
    return tss::createDynamicsSimulator<double, double>(bodies, propagatorSettings, areEquationsOfMotionToBeIntegrated);
}

WASM_MODULE_PATH("dynamics_simulator")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_simulator) {
    using namespace emscripten;

    // DynamicsSimulator base class
    class_<tp::DynamicsSimulator<double, double>>("dynamics_simulator_DynamicsSimulator")
        .smart_ptr<std::shared_ptr<tp::DynamicsSimulator<double, double>>>(
            "shared_ptr_DynamicsSimulator")
        .function("getSystemOfBodies", &tp::DynamicsSimulator<double, double>::getSystemOfBodies)
        .function("getSetIntegratedResult", &tp::DynamicsSimulator<double, double>::getSetIntegratedResult)
        .function("getPropagationResults", &tp::DynamicsSimulator<double, double>::getPropagationResults);

    // SingleArcDynamicsSimulator
    class_<tp::SingleArcDynamicsSimulator<double, double>,
           base<tp::DynamicsSimulator<double, double>>>(
        "dynamics_simulator_SingleArcDynamicsSimulator")
        .smart_ptr<std::shared_ptr<tp::SingleArcDynamicsSimulator<double, double>>>(
            "shared_ptr_SingleArcDynamicsSimulator")
        .constructor<const tss::SystemOfBodies&,
                     const std::shared_ptr<tp::SingleArcPropagatorSettings<double, double>>,
                     const bool>()
        .function("integrateEquationsOfMotion",
            static_cast<void(tp::SingleArcDynamicsSimulator<double, double>::*)(
                const Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>&)>(
                &tp::SingleArcDynamicsSimulator<double, double>::integrateEquationsOfMotion))
        .function("getEquationsOfMotionNumericalSolution",
            &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolution)
        .function("getEquationsOfMotionNumericalSolutionRaw",
            &tp::SingleArcDynamicsSimulator<double, double>::getEquationsOfMotionNumericalSolutionRaw)
        .function("getDependentVariableHistory",
            &tp::SingleArcDynamicsSimulator<double, double>::getDependentVariableHistory)
        .function("getCumulativeComputationTimeHistory",
            &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeComputationTimeHistory)
        .function("getCumulativeNumberOfFunctionEvaluations",
            &tp::SingleArcDynamicsSimulator<double, double>::getCumulativeNumberOfFunctionEvaluations)
        .function("getIntegratorSettings",
            &tp::SingleArcDynamicsSimulator<double, double>::getIntegratorSettings)
        .function("getSingleArcPropagationResults",
            &tp::SingleArcDynamicsSimulator<double, double>::getSingleArcPropagationResults)
        .function("getPropagatorSettings",
            &tp::SingleArcDynamicsSimulator<double, double>::getPropagatorSettings)
        .function("getInitialPropagationTime",
            &tp::SingleArcDynamicsSimulator<double, double>::getInitialPropagationTime)
        .function("integrationCompletedSuccessfully",
            &tp::SingleArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully);

    // MultiArcDynamicsSimulator
    class_<tp::MultiArcDynamicsSimulator<double, double>,
           base<tp::DynamicsSimulator<double, double>>>(
        "dynamics_simulator_MultiArcDynamicsSimulator")
        .smart_ptr<std::shared_ptr<tp::MultiArcDynamicsSimulator<double, double>>>(
            "shared_ptr_MultiArcDynamicsSimulator")
        .constructor<const tss::SystemOfBodies&,
                     const std::shared_ptr<tp::MultiArcPropagatorSettings<double, double>>,
                     const bool>()
        .function("getSingleArcDynamicsSimulators",
            &tp::MultiArcDynamicsSimulator<double, double>::getSingleArcDynamicsSimulators)
        .function("getMultiArcPropagationResults",
            &tp::MultiArcDynamicsSimulator<double, double>::getMultiArcPropagationResults)
        .function("integrationCompletedSuccessfully",
            &tp::MultiArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully);

    // HybridArcDynamicsSimulator
    class_<tp::HybridArcDynamicsSimulator<double, double>,
           base<tp::DynamicsSimulator<double, double>>>(
        "dynamics_simulator_HybridArcDynamicsSimulator")
        .smart_ptr<std::shared_ptr<tp::HybridArcDynamicsSimulator<double, double>>>(
            "shared_ptr_HybridArcDynamicsSimulator")
        .constructor<const tss::SystemOfBodies&,
                     const std::shared_ptr<tp::HybridArcPropagatorSettings<double, double>>,
                     const bool>()
        .function("getSingleArcDynamicsSimulator",
            &tp::HybridArcDynamicsSimulator<double, double>::getSingleArcDynamicsSimulator)
        .function("getMultiArcDynamicsSimulator",
            &tp::HybridArcDynamicsSimulator<double, double>::getMultiArcDynamicsSimulator)
        .function("getHybridArcPropagationResults",
            &tp::HybridArcDynamicsSimulator<double, double>::getHybridArcPropagationResults)
        .function("integrationCompletedSuccessfully",
            &tp::HybridArcDynamicsSimulator<double, double>::integrationCompletedSuccessfully);

    // Factory function for creating dynamics simulator
    function("dynamics_simulator_create_dynamics_simulator",
        &createDynamicsSimulatorWrapper);
}

#endif
