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

#include <tudat/math/root_finders.h>

namespace trf = tudat::root_finders;

WASM_MODULE_PATH("math_root_finders")

EMSCRIPTEN_BINDINGS(tudatpy_math_root_finders) {
    using namespace emscripten;

    // Enums
    enum_<trf::MaximumIterationHandling>("math_root_finders_MaximumIterationHandling")
        .value("accept_result", trf::accept_result)
        .value("accept_result_with_warning", trf::accept_result_with_warning)
        .value("throw_exception", trf::throw_exception);

    // RootFinderSettings base class
    class_<trf::RootFinderSettings>("math_root_finders_RootFinderSettings")
        .smart_ptr<std::shared_ptr<trf::RootFinderSettings>>("shared_ptr_RootFinderSettings");

    // Factory functions
    function("math_root_finders_bisection",
        &trf::bisectionRootFinderSettings);

    function("math_root_finders_newton_raphson",
        &trf::newtonRaphsonRootFinderSettings);

    function("math_root_finders_halley",
        &trf::halleyRootFinderSettings);

    function("math_root_finders_secant",
        &trf::secantRootFinderSettings);
}

#endif
