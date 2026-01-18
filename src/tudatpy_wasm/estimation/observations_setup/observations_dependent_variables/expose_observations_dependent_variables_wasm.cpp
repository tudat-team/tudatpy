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
#include "../../../stl_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/estimation_setup/observationOutputSettings.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("estimation_observations_setup_observations_dependent_variables")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_setup_observations_dependent_variables) {
    using namespace emscripten;

    // ObservationDependentVariableSettings base class
    class_<tss::ObservationDependentVariableSettings>(
        "estimation_observations_setup_observations_dependent_variables_ObservationDependentVariableSettings")
        .smart_ptr<std::shared_ptr<tss::ObservationDependentVariableSettings>>(
            "shared_ptr_ObservationDependentVariableSettings");

    // Factory functions for observation dependent variables
    function("estimation_observations_setup_observations_dependent_variables_elevation_angle",
        &tss::elevationAngleDependentVariable);

    function("estimation_observations_setup_observations_dependent_variables_azimuth_angle",
        &tss::azimuthAngleDependentVariable);

    function("estimation_observations_setup_observations_dependent_variables_target_range",
        &tss::targetRangeBetweenLinkEndsDependentVariable);
}

#endif
