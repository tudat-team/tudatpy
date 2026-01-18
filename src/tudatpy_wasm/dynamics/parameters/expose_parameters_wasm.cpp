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

#include <tudat/simulation/estimation_setup/createEstimatableParameters.h>

namespace tep = tudat::estimatable_parameters;

WASM_MODULE_PATH("dynamics_parameters")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_parameters) {
    using namespace emscripten;

    // EstimatableParameterSet class
    class_<tep::EstimatableParameterSet<double>>(
        "dynamics_parameters_EstimatableParameterSet")
        .smart_ptr<std::shared_ptr<tep::EstimatableParameterSet<double>>>(
            "shared_ptr_EstimatableParameterSet")
        .function("parameter_set_size",
            &tep::EstimatableParameterSet<double>::getEstimatedParameterSetSize)
        .function("initial_states_size",
            &tep::EstimatableParameterSet<double>::getInitialDynamicalStateParameterSize)
        .function("initial_single_arc_states_size",
            &tep::EstimatableParameterSet<double>::getInitialDynamicalSingleArcStateParameterSize)
        .function("initial_multi_arc_states_size",
            &tep::EstimatableParameterSet<double>::getInitialDynamicalMultiArcStateParameterSize)
        .function("constraints_size",
            &tep::EstimatableParameterSet<double>::getConstraintSize)
        .function("get_parameter_descriptions",
            &tep::EstimatableParameterSet<double>::getParametersDescriptions);

    // Helper function to print parameter names
    function("dynamics_parameters_print_parameter_names",
        &tep::printEstimatableParameterEntries<double>);
}

#endif
