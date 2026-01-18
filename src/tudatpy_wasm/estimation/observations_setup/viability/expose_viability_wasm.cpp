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

#include <tudat/simulation/estimation_setup/createObservationViability.h>

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observations_setup_viability")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_setup_viability) {
    using namespace emscripten;

    // ObservationViabilityType enum
    enum_<tom::ObservationViabilityType>(
        "estimation_observations_setup_viability_ObservationViabilityType")
        .value("minimum_elevation_angle", tom::minimum_elevation_angle)
        .value("body_avoidance_angle", tom::body_avoidance_angle)
        .value("body_occultation", tom::body_occultation);

    // ObservationViabilitySettings base class (in observation_models namespace)
    class_<tom::ObservationViabilitySettings>(
        "estimation_observations_setup_viability_ObservationViabilitySettings")
        .smart_ptr<std::shared_ptr<tom::ObservationViabilitySettings>>(
            "shared_ptr_ObservationViabilitySettings");

    // Factory functions (in observation_models namespace)
    // Single link end versions take pair<string, string>
    function("estimation_observations_setup_viability_elevation_angle_viability",
        select_overload<std::shared_ptr<tom::ObservationViabilitySettings>(
            const std::pair<std::string, std::string>, const double)>(
            &tom::elevationAngleViabilitySettings));

    function("estimation_observations_setup_viability_body_avoidance_viability",
        select_overload<std::shared_ptr<tom::ObservationViabilitySettings>(
            const std::pair<std::string, std::string>, const std::string, const double)>(
            &tom::bodyAvoidanceAngleViabilitySettings));

    function("estimation_observations_setup_viability_body_occultation_viability",
        select_overload<std::shared_ptr<tom::ObservationViabilitySettings>(
            const std::pair<std::string, std::string>, const std::string)>(
            &tom::bodyOccultationViabilitySettings));
}

#endif
