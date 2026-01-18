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

WASM_MODULE_PATH("dynamics_propagation_setup_propagator")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_propagator) {
    using namespace emscripten;

    // TranslationalPropagatorType enum
    enum_<tp::TranslationalPropagatorType>("dynamics_propagation_setup_propagator_TranslationalPropagatorType")
        .value("undefined_translational_propagator", tp::undefined_translational_propagator)
        .value("cowell", tp::cowell)
        .value("encke", tp::encke)
        .value("gauss_keplerian", tp::gauss_keplerian)
        .value("gauss_modified_equinoctial", tp::gauss_modified_equinoctial)
        .value("unified_state_model_quaternions", tp::unified_state_model_quaternions)
        .value("unified_state_model_modified_rodrigues_parameters", tp::unified_state_model_modified_rodrigues_parameters)
        .value("unified_state_model_exponential_map", tp::unified_state_model_exponential_map);

    // RotationalPropagatorType enum
    enum_<tp::RotationalPropagatorType>("dynamics_propagation_setup_propagator_RotationalPropagatorType")
        .value("undefined_rotational_propagator", tp::undefined_rotational_propagator)
        .value("quaternions", tp::quaternions)
        .value("modified_rodrigues_parameters", tp::modified_rodrigues_parameters)
        .value("exponential_map", tp::exponential_map);

    // IntegratedStateType enum
    enum_<tp::IntegratedStateType>("dynamics_propagation_setup_propagator_IntegratedStateType")
        .value("hybrid", tp::hybrid)
        .value("translational_state", tp::translational_state)
        .value("rotational_state", tp::rotational_state)
        .value("body_mass_state", tp::body_mass_state)
        .value("custom_state", tp::custom_state);

    // PropagationTerminationReason enum
    enum_<tp::PropagationTerminationReason>("dynamics_propagation_setup_propagator_PropagationTerminationReason")
        .value("propagation_never_run", tp::propagation_never_run)
        .value("unknown_propagation_termination_reason", tp::unknown_propagation_termination_reason)
        .value("termination_condition_reached", tp::termination_condition_reached)
        .value("runtime_error_caught_in_propagation", tp::runtime_error_caught_in_propagation)
        .value("nan_or_inf_detected_in_state", tp::nan_or_inf_detected_in_state);

    // PropagationTerminationDetails class
    class_<tp::PropagationTerminationDetails>("dynamics_propagation_setup_propagator_PropagationTerminationDetails")
        .smart_ptr<std::shared_ptr<tp::PropagationTerminationDetails>>("shared_ptr_PropagationTerminationDetails")
        .function("getPropagationTerminationReason", &tp::PropagationTerminationDetails::getPropagationTerminationReason)
        .function("getTerminationOnExactCondition", &tp::PropagationTerminationDetails::getTerminationOnExactCondition);

    // PropagationTerminationSettings base class
    class_<tp::PropagationTerminationSettings>("dynamics_propagation_setup_propagator_PropagationTerminationSettings")
        .smart_ptr<std::shared_ptr<tp::PropagationTerminationSettings>>("shared_ptr_PropagationTerminationSettings");

    // PropagatorSettings base class (template instantiation for double)
    class_<tp::PropagatorSettings<double>>("dynamics_propagation_setup_propagator_PropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::PropagatorSettings<double>>>("shared_ptr_PropagatorSettings")
        .function("getInitialStates", &tp::PropagatorSettings<double>::getInitialStates)
        .function("resetInitialStates", &tp::PropagatorSettings<double>::resetInitialStates)
        .function("getPropagatedStateSize", &tp::PropagatorSettings<double>::getPropagatedStateSize)
        .function("getIsMultiArc", &tp::PropagatorSettings<double>::getIsMultiArc);

    // SingleArcPropagatorSettings
    class_<tp::SingleArcPropagatorSettings<double, double>, base<tp::PropagatorSettings<double>>>(
        "dynamics_propagation_setup_propagator_SingleArcPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::SingleArcPropagatorSettings<double, double>>>(
            "shared_ptr_SingleArcPropagatorSettings")
        .function("getStateType", &tp::SingleArcPropagatorSettings<double, double>::getStateType)
        .function("getTerminationSettings", &tp::SingleArcPropagatorSettings<double, double>::getTerminationSettings)
        .function("resetTerminationSettings", &tp::SingleArcPropagatorSettings<double, double>::resetTerminationSettings)
        .function("getInitialTime", &tp::SingleArcPropagatorSettings<double, double>::getInitialTime)
        .function("resetInitialTime", &tp::SingleArcPropagatorSettings<double, double>::resetInitialTime);

    // TranslationalStatePropagatorSettings
    class_<tp::TranslationalStatePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_TranslationalStatePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::TranslationalStatePropagatorSettings<double, double>>>(
            "shared_ptr_TranslationalStatePropagatorSettings");

    // RotationalStatePropagatorSettings
    class_<tp::RotationalStatePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_RotationalStatePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::RotationalStatePropagatorSettings<double, double>>>(
            "shared_ptr_RotationalStatePropagatorSettings");

    // MassPropagatorSettings
    class_<tp::MassPropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_MassPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::MassPropagatorSettings<double, double>>>(
            "shared_ptr_MassPropagatorSettings");

    // MultiTypePropagatorSettings
    class_<tp::MultiTypePropagatorSettings<double, double>,
           base<tp::SingleArcPropagatorSettings<double, double>>>(
        "dynamics_propagation_setup_propagator_MultiTypePropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::MultiTypePropagatorSettings<double, double>>>(
            "shared_ptr_MultiTypePropagatorSettings");

    // MultiArcPropagatorSettings
    class_<tp::MultiArcPropagatorSettings<double, double>, base<tp::PropagatorSettings<double>>>(
        "dynamics_propagation_setup_propagator_MultiArcPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::MultiArcPropagatorSettings<double, double>>>(
            "shared_ptr_MultiArcPropagatorSettings")
        .function("getSingleArcSettings", &tp::MultiArcPropagatorSettings<double, double>::getSingleArcSettings)
        .function("getArcStartTimes", &tp::MultiArcPropagatorSettings<double, double>::getArcStartTimes);

    // HybridArcPropagatorSettings
    class_<tp::HybridArcPropagatorSettings<double, double>, base<tp::PropagatorSettings<double>>>(
        "dynamics_propagation_setup_propagator_HybridArcPropagatorSettings")
        .smart_ptr<std::shared_ptr<tp::HybridArcPropagatorSettings<double, double>>>(
            "shared_ptr_HybridArcPropagatorSettings");

    // PropagatorProcessingSettings base class
    class_<tp::PropagatorProcessingSettings>("dynamics_propagation_setup_propagator_PropagatorProcessingSettings")
        .smart_ptr<std::shared_ptr<tp::PropagatorProcessingSettings>>("shared_ptr_PropagatorProcessingSettings")
        .function("getClearNumericalSolutions", &tp::PropagatorProcessingSettings::getClearNumericalSolutions)
        .function("getSetIntegratedResult", &tp::PropagatorProcessingSettings::getSetIntegratedResult);

    // SingleArcPropagatorProcessingSettings
    class_<tp::SingleArcPropagatorProcessingSettings, base<tp::PropagatorProcessingSettings>>(
        "dynamics_propagation_setup_propagator_SingleArcPropagatorProcessingSettings")
        .smart_ptr<std::shared_ptr<tp::SingleArcPropagatorProcessingSettings>>(
            "shared_ptr_SingleArcPropagatorProcessingSettings")
        .function("setResultsSaveFrequencyInSteps", &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSteps)
        .function("setResultsSaveFrequencyInSeconds", &tp::SingleArcPropagatorProcessingSettings::setResultsSaveFrequencyInSeconds);

    // Termination settings factory functions
    function("dynamics_propagation_setup_propagator_time_termination",
        select_overload<std::shared_ptr<tp::PropagationTerminationSettings>(
            const double, const bool)>(&tp::propagationTimeTerminationSettings));

    function("dynamics_propagation_setup_propagator_cpu_time_termination",
        &tp::propagationCPUTimeTerminationSettings);

    function("dynamics_propagation_setup_propagator_dependent_variable_termination",
        select_overload<std::shared_ptr<tp::PropagationTerminationSettings>(
            const std::shared_ptr<tp::SingleDependentVariableSaveSettings>,
            const double, const bool, const bool,
            const std::shared_ptr<tudat::root_finders::RootFinderSettings>)>(
            &tp::propagationDependentVariableTerminationSettings));

    function("dynamics_propagation_setup_propagator_hybrid_termination",
        select_overload<std::shared_ptr<tp::PropagationTerminationSettings>(
            const std::vector<std::shared_ptr<tp::PropagationTerminationSettings>>,
            const bool)>(&tp::propagationHybridTerminationSettings));
}

#endif
