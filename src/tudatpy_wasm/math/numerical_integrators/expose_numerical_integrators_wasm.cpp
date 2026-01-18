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
#include "../../shared_ptr_wasm.h"

#include <tudat/math/integrators/createNumericalIntegrator.h>
#include <tudat/math/integrators/rungeKuttaCoefficients.h>

namespace tni = tudat::numerical_integrators;

WASM_MODULE_PATH("math_numerical_integrators")

EMSCRIPTEN_BINDINGS(tudatpy_math_numerical_integrators) {
    using namespace emscripten;

    // AvailableIntegrators enum
    enum_<tni::AvailableIntegrators>("math_numerical_integrators_AvailableIntegrators")
        .value("rungeKutta4", tni::rungeKutta4)
        .value("euler", tni::euler)
        .value("rungeKuttaVariableStepSize", tni::rungeKuttaVariableStepSize)
        .value("adamsBashforthMoulton", tni::adamsBashforthMoulton)
        .value("bulirschStoer", tni::bulirschStoer);

    // CoefficientSets enum (at namespace level, not inside RungeKuttaCoefficients)
    enum_<tni::CoefficientSets>("math_numerical_integrators_CoefficientSets")
        .value("rungeKutta4Classic", tni::rungeKutta4Classic)
        .value("rungeKuttaFehlberg45", tni::rungeKuttaFehlberg45)
        .value("rungeKuttaFehlberg56", tni::rungeKuttaFehlberg56)
        .value("rungeKuttaFehlberg78", tni::rungeKuttaFehlberg78)
        .value("rungeKutta87DormandPrince", tni::rungeKutta87DormandPrince)
        .value("forwardEuler", tni::forwardEuler);

    // IntegratorSettings base class
    class_<tni::IntegratorSettings<double>>(
        "math_numerical_integrators_IntegratorSettings")
        .smart_ptr<std::shared_ptr<tni::IntegratorSettings<double>>>(
            "shared_ptr_IntegratorSettings");

    // RungeKuttaVariableStepSizeSettings derived class
    class_<tni::RungeKuttaVariableStepSizeSettings<double>,
           base<tni::IntegratorSettings<double>>>(
        "math_numerical_integrators_RungeKuttaVariableStepSizeSettings")
        .smart_ptr<std::shared_ptr<tni::RungeKuttaVariableStepSizeSettings<double>>>(
            "shared_ptr_RungeKuttaVariableStepSizeSettings");

    // BulirschStoerIntegratorSettings derived class
    class_<tni::BulirschStoerIntegratorSettings<double>,
           base<tni::IntegratorSettings<double>>>(
        "math_numerical_integrators_BulirschStoerIntegratorSettings")
        .smart_ptr<std::shared_ptr<tni::BulirschStoerIntegratorSettings<double>>>(
            "shared_ptr_BulirschStoerIntegratorSettings");

    // AdamsBashforthMoultonSettings derived class
    class_<tni::AdamsBashforthMoultonSettings<double>,
           base<tni::IntegratorSettings<double>>>(
        "math_numerical_integrators_AdamsBashforthMoultonSettings")
        .smart_ptr<std::shared_ptr<tni::AdamsBashforthMoultonSettings<double>>>(
            "shared_ptr_AdamsBashforthMoultonSettings");
}

#endif
