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

#include <tudat/simulation/estimation_setup/simulateObservations.h>

namespace tss = tudat::simulation_setup;
namespace tom = tudat::observation_models;

// Wrapper function for simulateObservations (template function)
std::shared_ptr<tom::ObservationCollection<double, double>> simulateObservationsWrapper(
    const std::vector<std::shared_ptr<tss::ObservationSimulationSettings<double>>>& observationsToSimulate,
    const std::vector<std::shared_ptr<tom::ObservationSimulatorBase<double, double>>>& observationSimulators,
    const tss::SystemOfBodies& bodies)
{
    return tss::simulateObservations<double, double>(observationsToSimulate, observationSimulators, bodies);
}

WASM_MODULE_PATH("estimation_observations_setup_observations_wrapper")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations_setup_observations_wrapper) {
    using namespace emscripten;

    // Factory function for simulating observations
    function("estimation_observations_setup_observations_wrapper_simulate_observations",
        &simulateObservationsWrapper);
}

#endif
