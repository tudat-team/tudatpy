/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Integrator settings module bindings for WASM.
 *    Mirrors: src/tudatpy/dynamics/propagation_setup/integrator/expose_integrator.cpp
 *
 *    Provides settings for numerical integrators used in orbit propagation:
 *    - Runge-Kutta methods (RK4, RKF45, RKF78, etc.)
 *    - Bulirsch-Stoer
 *    - Adams-Bashforth-Moulton
 */

#ifdef __EMSCRIPTEN__

#include <emscripten/bind.h>
#include "expose_integrator_wasm.h"
#include "../../../wasm_module.h"
#include "../../../shared_ptr_wasm.h"

// Tudat headers
#include <tudat/simulation/propagation_setup.h>

namespace tni = tudat::numerical_integrators;

// Use double for both state and time (standard precision)
using TIME_TYPE = double;

// Register this module path
WASM_MODULE_PATH("dynamics_propagation_setup_integrator")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_integrator) {
    using namespace emscripten;

    // ========================================================================
    // Enumerations
    // ========================================================================

    // Minimum step handling
    enum_<tni::MinimumIntegrationTimeStepHandling>("dynamics_propagation_setup_integrator_MinimumIntegrationTimeStepHandling")
        .value("throw_exception_below_minimum", tni::MinimumIntegrationTimeStepHandling::throw_exception_below_minimum)
        .value("set_to_minimum_step_silently", tni::MinimumIntegrationTimeStepHandling::set_to_minimum_step_silently)
        .value("set_to_minimum_step_single_warning", tni::MinimumIntegrationTimeStepHandling::set_to_minimum_step_single_warning)
        .value("set_to_minimum_step_every_time_warning", tni::MinimumIntegrationTimeStepHandling::set_to_minimum_step_every_time_warning);

    // Available integrator types
    enum_<tni::AvailableIntegrators>("dynamics_propagation_setup_integrator_AvailableIntegrators")
        .value("runge_kutta_fixed_step_size_type", tni::AvailableIntegrators::rungeKuttaFixedStepSize)
        .value("runge_kutta_variable_step_size_type", tni::AvailableIntegrators::rungeKuttaVariableStepSize)
        .value("bulirsch_stoer_type", tni::AvailableIntegrators::bulirschStoer)
        .value("adams_bashforth_moulton_type", tni::AvailableIntegrators::adamsBashforthMoulton);

    // Runge-Kutta coefficient sets
    enum_<tni::CoefficientSets>("dynamics_propagation_setup_integrator_CoefficientSets")
        .value("euler_forward", tni::forwardEuler)
        .value("rk_4", tni::rungeKutta4Classic)
        .value("explicit_mid_point", tni::explicitMidPoint)
        .value("explicit_trapezoid_rule", tni::explicitTrapezoidRule)
        .value("ralston", tni::ralston)
        .value("rk_3", tni::rungeKutta3)
        .value("ralston_3", tni::ralston3)
        .value("SSPRK3", tni::SSPRK3)
        .value("ralston_4", tni::ralston4)
        .value("three_eight_rule_rk_4", tni::threeEighthRuleRK4)
        .value("heun_euler", tni::heunEuler)
        .value("rkf_12", tni::rungeKuttaFehlberg12)
        .value("rkf_45", tni::rungeKuttaFehlberg45)
        .value("rkf_56", tni::rungeKuttaFehlberg56)
        .value("rkf_78", tni::rungeKuttaFehlberg78)
        .value("rkdp_87", tni::rungeKutta87DormandPrince)
        .value("rkf_89", tni::rungeKuttaFehlberg89)
        .value("rkv_89", tni::rungeKuttaVerner89)
        .value("rkf_108", tni::rungeKuttaFeagin108)
        .value("rkf_1210", tni::rungeKuttaFeagin1210)
        .value("rkf_1412", tni::rungeKuttaFeagin1412);

    // Order to integrate
    enum_<tni::RungeKuttaCoefficients::OrderEstimateToIntegrate>("dynamics_propagation_setup_integrator_OrderToIntegrate")
        .value("lower", tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::lower)
        .value("higher", tni::RungeKuttaCoefficients::OrderEstimateToIntegrate::higher);

    // Extrapolation method step sequences
    enum_<tni::ExtrapolationMethodStepSequences>("dynamics_propagation_setup_integrator_ExtrapolationMethodStepSequences")
        .value("bulirsch_stoer_sequence", tni::ExtrapolationMethodStepSequences::bulirsch_stoer_sequence)
        .value("deufelhard_sequence", tni::ExtrapolationMethodStepSequences::deufelhard_sequence);

    // ========================================================================
    // Settings Classes
    // ========================================================================

    // Base IntegratorSettings class
    class_<tni::IntegratorSettings<TIME_TYPE>>("dynamics_propagation_setup_integrator_IntegratorSettings")
        .smart_ptr<std::shared_ptr<tni::IntegratorSettings<TIME_TYPE>>>("shared_ptr_IntegratorSettings");

    // RungeKuttaFixedStepSizeSettings
    class_<tni::RungeKuttaFixedStepSizeSettings<TIME_TYPE>,
           base<tni::IntegratorSettings<TIME_TYPE>>>("dynamics_propagation_setup_integrator_RungeKuttaFixedStepSizeSettings")
        .smart_ptr<std::shared_ptr<tni::RungeKuttaFixedStepSizeSettings<TIME_TYPE>>>("shared_ptr_RungeKuttaFixedStepSizeSettings");

    // Step size control settings
    class_<tni::IntegratorStepSizeControlSettings>("dynamics_propagation_setup_integrator_IntegratorStepSizeControlSettings")
        .smart_ptr<std::shared_ptr<tni::IntegratorStepSizeControlSettings>>("shared_ptr_IntegratorStepSizeControlSettings");

    // Step size validation settings
    class_<tni::IntegratorStepSizeValidationSettings>("dynamics_propagation_setup_integrator_IntegratorStepSizeValidationSettings")
        .smart_ptr<std::shared_ptr<tni::IntegratorStepSizeValidationSettings>>("shared_ptr_IntegratorStepSizeValidationSettings");

    // BulirschStoerIntegratorSettings
    class_<tni::BulirschStoerIntegratorSettings<TIME_TYPE>,
           base<tni::IntegratorSettings<TIME_TYPE>>>("dynamics_propagation_setup_integrator_BulirschStoerIntegratorSettings")
        .smart_ptr<std::shared_ptr<tni::BulirschStoerIntegratorSettings<TIME_TYPE>>>("shared_ptr_BulirschStoerIntegratorSettings");

    // AdamsBashforthMoultonSettings
    class_<tni::AdamsBashforthMoultonSettings<TIME_TYPE>,
           base<tni::IntegratorSettings<TIME_TYPE>>>("dynamics_propagation_setup_integrator_AdamsBashforthMoultonSettings")
        .smart_ptr<std::shared_ptr<tni::AdamsBashforthMoultonSettings<TIME_TYPE>>>("shared_ptr_AdamsBashforthMoultonSettings");

    // ========================================================================
    // Factory Functions
    // ========================================================================

    // Step size validation settings factory
    function("dynamics_propagation_setup_integrator_step_size_validation",
        &tni::stepSizeValidationSettings);

    // Step size control settings factories
    function("dynamics_propagation_setup_integrator_step_size_control_elementwise_scalar_tolerance",
        &tni::perElementIntegratorStepSizeControlSettings<double>);

    // Runge-Kutta fixed step
    function("dynamics_propagation_setup_integrator_runge_kutta_fixed_step",
        &tni::rungeKuttaFixedStepSettings<TIME_TYPE>);

    // Runge-Kutta variable step
    function("dynamics_propagation_setup_integrator_runge_kutta_variable_step",
        &tni::multiStageVariableStepSizeSettings<TIME_TYPE>);

    // Bulirsch-Stoer variable step
    function("dynamics_propagation_setup_integrator_bulirsch_stoer_variable_step",
        &tni::bulirschStoerVariableStepIntegratorSettings<TIME_TYPE>);

    // Bulirsch-Stoer fixed step
    function("dynamics_propagation_setup_integrator_bulirsch_stoer_fixed_step",
        &tni::bulirschStoerFixedStepIntegratorSettings<TIME_TYPE>);

    // Adams-Bashforth-Moulton
    function("dynamics_propagation_setup_integrator_adams_bashforth_moulton",
        &tni::adamsBashforthMoultonSettings<TIME_TYPE>);

    function("dynamics_propagation_setup_integrator_adams_bashforth_moulton_fixed_order",
        &tni::adamsBashforthMoultonSettingsFixedOrder<TIME_TYPE>);

    function("dynamics_propagation_setup_integrator_adams_bashforth_moulton_fixed_step",
        &tni::adamsBashforthMoultonSettingsFixedStep<TIME_TYPE>);

    function("dynamics_propagation_setup_integrator_adams_bashforth_moulton_fixed_step_fixed_order",
        &tni::adamsBashforthMoultonSettingsFixedStepFixedOrder<TIME_TYPE>);

    // ========================================================================
    // Convenience functions for simple integrators
    // ========================================================================

    // Euler
    function("dynamics_propagation_setup_integrator_euler",
        &tni::eulerSettings<TIME_TYPE>);

    // RK4
    function("dynamics_propagation_setup_integrator_runge_kutta_4",
        &tni::rungeKutta4Settings<TIME_TYPE>);
}

#endif // __EMSCRIPTEN__
