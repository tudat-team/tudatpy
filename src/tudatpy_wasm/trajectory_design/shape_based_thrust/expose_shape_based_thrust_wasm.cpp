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

#include <tudat/astro/low_thrust/shape_based/baseFunctionsHodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/getRecommendedBaseFunctionsHodographicShaping.h>
#include <tudat/astro/low_thrust/shape_based/hodographicShapingLeg.h>
#include <tudat/astro/low_thrust/shape_based/sphericalShapingLeg.h>

namespace tsbm = tudat::shape_based_methods;

WASM_MODULE_PATH("trajectory_design_shape_based_thrust")

EMSCRIPTEN_BINDINGS(tudatpy_trajectory_design_shape_based_thrust) {
    using namespace emscripten;

    // baseFunctionHodographicShapingType enum
    enum_<tsbm::baseFunctionHodographicShapingType>(
        "trajectory_design_shape_based_thrust_BaseFunctionHodographicShapingType")
        .value("constant", tsbm::constant)
        .value("sine", tsbm::sine)
        .value("cosine", tsbm::cosine)
        .value("exponential", tsbm::exponential)
        .value("scaledExponential", tsbm::scaledExponential)
        .value("exponentialSine", tsbm::exponentialSine)
        .value("scaledExponentialSine", tsbm::scaledExponentialSine)
        .value("exponentialCosine", tsbm::exponentialCosine)
        .value("scaledExponentialCosine", tsbm::scaledExponentialCosine)
        .value("power", tsbm::power)
        .value("scaledPower", tsbm::scaledPower)
        .value("powerCosine", tsbm::powerCosine)
        .value("scaledPowerCosine", tsbm::scaledPowerCosine)
        .value("powerSine", tsbm::powerSine)
        .value("scaledPowerSine", tsbm::scaledPowerSine);

    // BaseFunctionHodographicShaping base class
    class_<tsbm::BaseFunctionHodographicShaping>(
        "trajectory_design_shape_based_thrust_BaseFunctionHodographicShaping")
        .smart_ptr<std::shared_ptr<tsbm::BaseFunctionHodographicShaping>>(
            "shared_ptr_BaseFunctionHodographicShaping")
        .function("evaluateFunction", &tsbm::BaseFunctionHodographicShaping::evaluateFunction)
        .function("evaluateDerivative", &tsbm::BaseFunctionHodographicShaping::evaluateDerivative)
        .function("evaluateIntegral", &tsbm::BaseFunctionHodographicShaping::evaluateIntegral);

    // Factory functions for hodographic shaping base functions
    function("trajectory_design_shape_based_thrust_hodograph_constant",
        &tsbm::hodographConstant);

    function("trajectory_design_shape_based_thrust_hodograph_sine",
        &tsbm::hodographSine);

    function("trajectory_design_shape_based_thrust_hodograph_cosine",
        &tsbm::hodographCosine);

    function("trajectory_design_shape_based_thrust_hodograph_exponential",
        &tsbm::hodographExponential);

    function("trajectory_design_shape_based_thrust_hodograph_scaled_exponential",
        &tsbm::hodographScaledExponential);

    function("trajectory_design_shape_based_thrust_hodograph_exponential_sine",
        &tsbm::hodographExponentialSine);

    function("trajectory_design_shape_based_thrust_hodograph_scaled_exponential_sine",
        &tsbm::hodographScaledExponentialSine);

    function("trajectory_design_shape_based_thrust_hodograph_exponential_cosine",
        &tsbm::hodographExponentialCosine);

    function("trajectory_design_shape_based_thrust_hodograph_scaled_exponential_cosine",
        &tsbm::hodographScaledExponentialCosine);

    function("trajectory_design_shape_based_thrust_hodograph_power",
        &tsbm::hodographPower);

    function("trajectory_design_shape_based_thrust_hodograph_scaled_power",
        &tsbm::hodographScaledPower);

    function("trajectory_design_shape_based_thrust_hodograph_power_sine",
        &tsbm::hodographPowerSine);

    function("trajectory_design_shape_based_thrust_hodograph_scaled_power_sine",
        &tsbm::hodographScaledPowerSine);

    function("trajectory_design_shape_based_thrust_hodograph_power_cosine",
        &tsbm::hodographPowerCosine);

    function("trajectory_design_shape_based_thrust_hodograph_scaled_power_cosine",
        &tsbm::hodographScaledPowerCosine);

    // Recommended base functions
    function("trajectory_design_shape_based_thrust_recommended_radial_hodograph_functions",
        select_overload<std::vector<std::shared_ptr<tsbm::BaseFunctionHodographicShaping>>(
            const double)>(&tsbm::getRecommendedRadialVelocityBaseFunctions));

    function("trajectory_design_shape_based_thrust_recommended_normal_hodograph_functions",
        select_overload<std::vector<std::shared_ptr<tsbm::BaseFunctionHodographicShaping>>(
            const double)>(&tsbm::getRecommendedNormalBaseFunctions));

    function("trajectory_design_shape_based_thrust_recommended_axial_hodograph_functions",
        select_overload<std::vector<std::shared_ptr<tsbm::BaseFunctionHodographicShaping>>(
            const double, const int)>(&tsbm::getRecommendedAxialVelocityBaseFunctions));
}

#endif
