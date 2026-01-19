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
#include "../../../wasm_module.h"
#include "../../../eigen_wasm.h"
#include "../../../stl_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/propagation_setup/propagationSettings.h>
#include <tudat/simulation/propagation_setup/propagationTerminationSettings.h>
#include <tudat/simulation/propagation_setup/propagationTermination.h>
#include <tudat/simulation/propagation_setup/propagationProcessingSettings.h>
#include <tudat/astro/propagators/nBodyStateDerivative.h>
#include <tudat/astro/propagators/rotationalMotionStateDerivative.h>
#include <tudat/astro/propagators/singleStateTypeDerivative.h>

namespace tp = tudat::propagators;
namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;
namespace tni = tudat::numerical_integrators;
namespace tmrf = tudat::root_finders;

WASM_MODULE_PATH("dynamics_propagation_setup_propagator")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_propagator) {
    using namespace emscripten;

    // ========================================================================
    // TranslationalPropagatorType enum
    // ========================================================================
    enum_<tp::TranslationalPropagatorType>("dynamics_propagation_setup_propagator_TranslationalPropagatorType")
        .value("undefined_translational_propagator", tp::undefined_translational_propagator)
        .value("cowell", tp::cowell)
        .value("encke", tp::encke)
        .value("gauss_keplerian", tp::gauss_keplerian)
        .value("gauss_modified_equinoctial", tp::gauss_modified_equinoctial)
        .value("unified_state_model_quaternions", tp::unified_state_model_quaternions)
        .value("unified_state_model_modified_rodrigues_parameters", tp::unified_state_model_modified_rodrigues_parameters)
        .value("unified_state_model_exponential_map", tp::unified_state_model_exponential_map);

    // ========================================================================
    // RotationalPropagatorType enum
    // ========================================================================
    enum_<tp::RotationalPropagatorType>("dynamics_propagation_setup_propagator_RotationalPropagatorType")
        .value("undefined_rotational_propagator", tp::undefined_rotational_propagator)
        .value("quaternions", tp::quaternions)
        .value("modified_rodrigues_parameters", tp::modified_rodrigues_parameters)
        .value("exponential_map", tp::exponential_map);

    // ========================================================================
    // PropagationTerminationTypes enum
    // ========================================================================
    enum_<tp::PropagationTerminationTypes>("dynamics_propagation_setup_propagator_PropagationTerminationTypes")
        .value("time_stopping_condition_type", tp::time_stopping_condition)
        .value("cpu_time_stopping_condition_type", tp::cpu_time_stopping_condition)
        .value("dependent_variable_stopping_condition_type", tp::dependent_variable_stopping_condition)
        .value("hybrid_stopping_condition_type", tp::hybrid_stopping_condition)
        .value("custom_stopping_condition_type", tp::custom_stopping_condition);

    // ========================================================================
    // IntegratedStateType enum (exposed as StateType in Python)
    // ========================================================================
    enum_<tp::IntegratedStateType>("dynamics_propagation_setup_propagator_IntegratedStateType")
        .value("hybrid_type", tp::hybrid)
        .value("translational_type", tp::translational_state)
        .value("rotational_type", tp::rotational_state)
        .value("mass_type", tp::body_mass_state)
        .value("custom_type", tp::custom_state);

    // ========================================================================
    // PropagationTerminationReason enum
    // ========================================================================
    enum_<tp::PropagationTerminationReason>("dynamics_propagation_setup_propagator_PropagationTerminationReason")
        .value("propagation_never_run", tp::propagation_never_run)
        .value("unknown_propagation_termination_reason", tp::unknown_propagation_termination_reason)
        .value("termination_condition_reached", tp::termination_condition_reached)
        .value("runtime_error_caught_in_propagation", tp::runtime_error_caught_in_propagation)
        .value("nan_or_inf_detected_in_state", tp::nan_or_inf_detected_in_state);

    // ========================================================================
    // PropagationTerminationDetails class
    // ========================================================================
    class_<tp::PropagationTerminationDetails>("dynamics_propagation_setup_propagator_PropagationTerminationDetails")
        .smart_ptr<std::shared_ptr<tp::PropagationTerminationDetails>>("shared_ptr_PropagationTerminationDetails")
        .function("getPropagationTerminationReason", &tp::PropagationTerminationDetails::getPropagationTerminationReason)
        .function("getTerminationOnExactCondition", &tp::PropagationTerminationDetails::getTerminationOnExactCondition);

    // ========================================================================
    // PropagationPrintSettings class
    // ========================================================================
    class_<tp::PropagationPrintSettings>("dynamics_propagation_setup_propagator_PropagationPrintSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationPrintSettings>>("shared_ptr_PropagationPrintSettings")
        .constructor<>()
        // Boolean properties with getter/setter
        .function("getPrintDependentVariableData", &tp::PropagationPrintSettings::getPrintDependentVariableData)
        .function("setPrintDependentVariableData", &tp::PropagationPrintSettings::setPrintDependentVariableData)
        .function("getPrintPropagatedStateData", &tp::PropagationPrintSettings::getPrintPropagatedStateData)
        .function("setPrintPropagatedStateData", &tp::PropagationPrintSettings::setPrintPropagatedStateData)
        .function("getPrintProcessedStateData", &tp::PropagationPrintSettings::getPrintProcessedStateData)
        .function("setPrintProcessedStateData", &tp::PropagationPrintSettings::setPrintProcessedStateData)
        .function("getPrintNumberOfFunctionEvaluations", &tp::PropagationPrintSettings::getPrintNumberOfFunctionEvaluations)
        .function("setPrintNumberOfFunctionEvaluations", &tp::PropagationPrintSettings::setPrintNumberOfFunctionEvaluations)
        .function("getPrintPropagationTime", &tp::PropagationPrintSettings::getPrintPropagationTime)
        .function("setPrintPropagationTime", &tp::PropagationPrintSettings::setPrintPropagationTime)
        .function("getPrintTerminationReason", &tp::PropagationPrintSettings::getPrintTerminationReason)
        .function("setPrintTerminationReason", &tp::PropagationPrintSettings::setPrintTerminationReason)
        .function("getPrintInitialAndFinalConditions", &tp::PropagationPrintSettings::getPrintInitialAndFinalConditions)
        .function("setPrintInitialAndFinalConditions", &tp::PropagationPrintSettings::setPrintInitialAndFinalConditions)
        // Frequency properties
        .function("getResultsPrintFrequencyInSeconds", &tp::PropagationPrintSettings::getResultsPrintFrequencyInSeconds)
        .function("setResultsPrintFrequencyInSeconds", &tp::PropagationPrintSettings::setResultsPrintFrequencyInSeconds)
        .function("getResultsPrintFrequencyInSteps", &tp::PropagationPrintSettings::getResultsPrintFrequencyInSteps)
        .function("setResultsPrintFrequencyInSteps", &tp::PropagationPrintSettings::setResultsPrintFrequencyInSteps)
        .function("getPrintDependentVariableDuringPropagation", &tp::PropagationPrintSettings::getPrintDependentVariableDuringPropagation)
        .function("setPrintDependentVariableDuringPropagation", &tp::PropagationPrintSettings::setPrintDependentVariableDuringPropagation)
        // Enable/disable methods
        .function("enableAllPrinting", select_overload<void()>(&tp::PropagationPrintSettings::enableAllPrinting))
        .function("enableAllPrintingWithFrequency", select_overload<void(const double, const int)>(&tp::PropagationPrintSettings::enableAllPrinting))
        .function("disableAllPrinting", &tp::PropagationPrintSettings::disableAllPrinting);

    // ========================================================================
    // PropagatorProcessingSettings base class
    // ========================================================================
    class_<tp::PropagatorProcessingSettings>("dynamics_propagation_setup_propagator_PropagatorProcessingSettings")
        .smart_ptr<std::shared_ptr<tp::PropagatorProcessingSettings>>("shared_ptr_PropagatorProcessingSettings")
        .function("getClearNumericalSolutions", &tp::PropagatorProcessingSettings::getClearNumericalSolutions)
        .function("setClearNumericalSolutions", &tp::PropagatorProcessingSettings::setClearNumericalSolutions)
        .function("getSetIntegratedResult", &tp::PropagatorProcessingSettings::getSetIntegratedResult)
        .function("setIntegratedResult", &tp::PropagatorProcessingSettings::setIntegratedResult)
        .function("getUpdateDependentVariableInterpolator", &tp::PropagatorProcessingSettings::getUpdateDependentVariableInterpolator)
        .function("setUpdateDependentVariableInterpolator", &tp::PropagatorProcessingSettings::setUpdateDependentVariableInterpolator);

    // ========================================================================
    // SingleArcPropagatorProcessingSettings
    // ========================================================================
    class_<tp::SingleArcPropagatorProcessingSettings, base<tp::PropagatorProcessingSettings>>(
        "dynamics_propagation_setup_propagator_SingleArcPropagatorProcessingSettings")
        .smart_ptr<std::shared_ptr<tp::SingleArcPropagatorProcessingSettings>>(
            "shared_ptr_SingleArcPropagatorProcessingSettings")
        .function("getPrintSettings", &tp::SingleArcPropagatorProcessingSettings::getPrintSettings)
        .function("getResultsSaveFrequencyInSteps", &tp::SingleArcPropagatorProcessingSettings::getResultsSaveFrequencyInSteps)
        .function("setResultsSaveFrequencyInSteps", &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSteps)
        .function("getResultsSaveFrequencyInSeconds", &tp::SingleArcPropagatorProcessingSettings::getResultsSaveFrequencyInSeconds)
        .function("setResultsSaveFrequencyInSeconds", &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSeconds);

    // ========================================================================
    // MultiArcPropagatorProcessingSettings
    // ========================================================================
    class_<tp::MultiArcPropagatorProcessingSettings, base<tp::PropagatorProcessingSettings>>(
        "dynamics_propagation_setup_propagator_MultiArcPropagatorProcessingSettings")
        .smart_ptr<std::shared_ptr<tp::MultiArcPropagatorProcessingSettings>>(
            "shared_ptr_MultiArcPropagatorProcessingSettings")
        .function("resetAndApplyConsistentSingleArcPrintSettings",
            &tp::MultiArcPropagatorProcessingSettings::resetAndApplyConsistentSingleArcPrintSettings)
        .function("getPrintFirstArcOnly", &tp::MultiArcPropagatorProcessingSettings::getPrintFirstArcOnly)
        .function("resetPrintFirstArcOnly", &tp::MultiArcPropagatorProcessingSettings::resetPrintFirstArcOnly)
        .function("useIdenticalSettings", &tp::MultiArcPropagatorProcessingSettings::useIdenticalSettings)
        .function("getSingleArcSettings", &tp::MultiArcPropagatorProcessingSettings::getSingleArcSettings);

    // ========================================================================
    // HybridArcPropagatorProcessingSettings
    // ========================================================================
    class_<tp::HybridArcPropagatorProcessingSettings, base<tp::PropagatorProcessingSettings>>(
        "dynamics_propagation_setup_propagator_HybridArcPropagatorProcessingSettings")
        .smart_ptr<std::shared_ptr<tp::HybridArcPropagatorProcessingSettings>>(
            "shared_ptr_HybridArcPropagatorProcessingSettings")
        .function("resetAndApplyConsistentPrintSettings",
            &tp::HybridArcPropagatorProcessingSettings::resetAndApplyConsistentPrintSettings)
        .function("getSingleArcSettings", &tp::HybridArcPropagatorProcessingSettings::getSingleArcSettings)
        .function("getMultiArcSettings", &tp::HybridArcPropagatorProcessingSettings::getMultiArcSettings);

    // ========================================================================
    // PropagationTerminationSettings base class
    // ========================================================================
    class_<tp::PropagationTerminationSettings>("dynamics_propagation_setup_propagator_PropagationTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationTerminationSettings>>("shared_ptr_PropagationTerminationSettings");

    // ========================================================================
    // Derived termination settings classes
    // ========================================================================
    class_<tp::PropagationDependentVariableTerminationSettings, base<tp::PropagationTerminationSettings>>(
        "dynamics_propagation_setup_propagator_PropagationDependentVariableTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationDependentVariableTerminationSettings>>(
            "shared_ptr_PropagationDependentVariableTerminationSettings");

    class_<tp::PropagationTimeTerminationSettings, base<tp::PropagationTerminationSettings>>(
        "dynamics_propagation_setup_propagator_PropagationTimeTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationTimeTerminationSettings>>(
            "shared_ptr_PropagationTimeTerminationSettings");

    class_<tp::PropagationCPUTimeTerminationSettings, base<tp::PropagationTerminationSettings>>(
        "dynamics_propagation_setup_propagator_PropagationCPUTimeTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationCPUTimeTerminationSettings>>(
            "shared_ptr_PropagationCPUTimeTerminationSettings");

    class_<tp::PropagationCustomTerminationSettings, base<tp::PropagationTerminationSettings>>(
        "dynamics_propagation_setup_propagator_PropagationCustomTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationCustomTerminationSettings>>(
            "shared_ptr_PropagationCustomTerminationSettings");

    class_<tp::PropagationHybridTerminationSettings, base<tp::PropagationTerminationSettings>>(
        "dynamics_propagation_setup_propagator_PropagationHybridTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationHybridTerminationSettings>>(
            "shared_ptr_PropagationHybridTerminationSettings");

    class_<tp::NonSequentialPropagationTerminationSettings, base<tp::PropagationTerminationSettings>>(
        "dynamics_propagation_setup_propagator_NonSequentialPropagationTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::NonSequentialPropagationTerminationSettings>>(
            "shared_ptr_NonSequentialPropagationTerminationSettings");

    // ========================================================================
    // PropagatorSettings base class (template instantiation for double)
    // ========================================================================
    class_<tp::PropagatorSettings<double>>("dynamics_propagation_setup_propagator_PropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::PropagatorSettings<double>>>("shared_ptr_PropagatorSettings")
        .function("getInitialStates", &tp::PropagatorSettings<double>::getInitialStates)
        .function("resetInitialStates", &tp::PropagatorSettings<double>::resetInitialStates)
        .function("getPropagatedStateSize", &tp::PropagatorSettings<double>::getPropagatedStateSize)
        .function("getIsMultiArc", &tp::PropagatorSettings<double>::getIsMultiArc);

    // ========================================================================
    // SingleArcPropagatorSettings
    // ========================================================================
    class_<tp::SingleArcPropagatorSettings<double, double>, base<tp::PropagatorSettings<double>>>(
        "dynamics_propagation_setup_propagator_SingleArcPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::SingleArcPropagatorSettings<double, double>>>(
            "shared_ptr_SingleArcPropagatorSettings")
        .function("getStateType", &tp::SingleArcPropagatorSettings<double, double>::getStateType)
        .function("getTerminationSettings", &tp::SingleArcPropagatorSettings<double, double>::getTerminationSettings)
        .function("resetTerminationSettings", &tp::SingleArcPropagatorSettings<double, double>::resetTerminationSettings)
        .function("getInitialTime", &tp::SingleArcPropagatorSettings<double, double>::getInitialTime)
        .function("resetInitialTime", &tp::SingleArcPropagatorSettings<double, double>::resetInitialTime)
        .function("getIntegratorSettings", &tp::SingleArcPropagatorSettings<double, double>::getIntegratorSettings)
        .function("setIntegratorSettings", &tp::SingleArcPropagatorSettings<double, double>::setIntegratorSettings)
        .function("getOutputSettings", &tp::SingleArcPropagatorSettings<double, double>::getOutputSettings)
        .function("getPrintSettings", &tp::SingleArcPropagatorSettings<double, double>::getPrintSettings);

    // ========================================================================
    // TranslationalStatePropagatorSettings
    // ========================================================================
    class_<tp::TranslationalStatePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_TranslationalStatePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::TranslationalStatePropagatorSettings<double, double>>>(
            "shared_ptr_TranslationalStatePropagatorSettings")
        .function("getPropagatedStateSize", &tp::TranslationalStatePropagatorSettings<double, double>::getPropagatedStateSize)
        .function("resetAccelerationModelsMap", &tp::TranslationalStatePropagatorSettings<double, double>::resetAccelerationModelsMap);

    // ========================================================================
    // RotationalStatePropagatorSettings
    // ========================================================================
    class_<tp::RotationalStatePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_RotationalStatePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::RotationalStatePropagatorSettings<double, double>>>(
            "shared_ptr_RotationalStatePropagatorSettings");

    // ========================================================================
    // MassPropagatorSettings
    // ========================================================================
    class_<tp::MassPropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_MassPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::MassPropagatorSettings<double, double>>>(
            "shared_ptr_MassPropagatorSettings");

    // ========================================================================
    // CustomStatePropagatorSettings
    // ========================================================================
    class_<tp::CustomStatePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_CustomStatePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::CustomStatePropagatorSettings<double, double>>>(
            "shared_ptr_CustomStatePropagatorSettings");

    // ========================================================================
    // MultiTypePropagatorSettings
    // ========================================================================
    class_<tp::MultiTypePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_MultiTypePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::MultiTypePropagatorSettings<double, double>>>(
            "shared_ptr_MultiTypePropagatorSettings")
        .function("resetInitialStates", &tp::MultiTypePropagatorSettings<double, double>::resetInitialStates)
        .function("resetIntegratedStateModels", &tp::MultiTypePropagatorSettings<double, double>::resetIntegratedStateModels)
        .function("getSingleTypePropagatorSettings", &tp::MultiTypePropagatorSettings<double, double>::getSingleTypePropagatorSettings)
        .function("getPropagatorSettingsMap", &tp::MultiTypePropagatorSettings<double, double>::getPropagatorSettingsMap);

    // ========================================================================
    // MultiArcPropagatorSettings
    // ========================================================================
    class_<tp::MultiArcPropagatorSettings<double, double>, base<tp::PropagatorSettings<double>>>(
        "dynamics_propagation_setup_propagator_MultiArcPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::MultiArcPropagatorSettings<double, double>>>(
            "shared_ptr_MultiArcPropagatorSettings")
        .function("getSingleArcSettings", &tp::MultiArcPropagatorSettings<double, double>::getSingleArcSettings)
        .function("getArcStartTimes", &tp::MultiArcPropagatorSettings<double, double>::getArcStartTimes)
        .function("getOutputSettings", &tp::MultiArcPropagatorSettings<double, double>::getOutputSettings)
        .function("getInitialStateList", &tp::MultiArcPropagatorSettings<double, double>::getInitialStateList);

    // ========================================================================
    // HybridArcPropagatorSettings
    // ========================================================================
    class_<tp::HybridArcPropagatorSettings<double, double>, base<tp::PropagatorSettings<double>>>(
        "dynamics_propagation_setup_propagator_HybridArcPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::HybridArcPropagatorSettings<double, double>>>(
            "shared_ptr_HybridArcPropagatorSettings")
        .function("getOutputSettings", &tp::HybridArcPropagatorSettings<double, double>::getOutputSettings);

    // ========================================================================
    // Termination settings factory functions
    // ========================================================================
    function("dynamics_propagation_setup_propagator_time_termination",
        &tp::propagationTimeTerminationSettings);

    function("dynamics_propagation_setup_propagator_cpu_time_termination",
        &tp::propagationCPUTimeTerminationSettings);

    function("dynamics_propagation_setup_propagator_dependent_variable_termination",
        &tp::propagationDependentVariableTerminationSettings);

    function("dynamics_propagation_setup_propagator_hybrid_termination",
        &tp::propagationHybridTerminationSettings);

    function("dynamics_propagation_setup_propagator_non_sequential_termination",
        &tp::nonSequentialPropagationTerminationSettings);

    // ========================================================================
    // Propagator settings factory functions
    // ========================================================================

    // Translational propagator
    function("dynamics_propagation_setup_propagator_translational",
        &tp::translationalStatePropagatorSettings<double, double>);

    // Rotational propagator
    function("dynamics_propagation_setup_propagator_rotational",
        &tp::rotationalStatePropagatorSettings<double, double>);

    // Mass propagator
    function("dynamics_propagation_setup_propagator_mass",
        &tp::massPropagatorSettings<double, double>);

    // Multitype propagator
    function("dynamics_propagation_setup_propagator_multitype",
        &tp::multiTypePropagatorSettings<double, double>);

    // Multi-arc propagator
    function("dynamics_propagation_setup_propagator_multi_arc",
        &tp::multiArcPropagatorSettings<double, double>);

    // Hybrid-arc propagator
    function("dynamics_propagation_setup_propagator_hybrid_arc",
        &tp::hybridArcPropagatorSettings<double, double>);

    // ========================================================================
    // Utility functions
    // ========================================================================
    function("dynamics_propagation_setup_propagator_get_single_integration_size",
        &tp::getSingleIntegrationSize);

    function("dynamics_propagation_setup_propagator_get_integrated_type_and_body_list",
        &tp::getIntegratedTypeAndBodyList<double, double>);

    function("dynamics_propagation_setup_propagator_add_dependent_variable_settings",
        &tp::addDepedentVariableSettings<double>);
}

#endif
