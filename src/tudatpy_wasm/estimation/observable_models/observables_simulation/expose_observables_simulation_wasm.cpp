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

namespace tom = tudat::observation_models;

WASM_MODULE_PATH("estimation_observable_models_observables_simulation")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observable_models_observables_simulation) {
    using namespace emscripten;

    // ObservationViabilityCalculator class
    class_<tom::ObservationViabilityCalculator>(
        "estimation_observable_models_observables_simulation_ObservationViabilityCalculator")
        .smart_ptr<std::shared_ptr<tom::ObservationViabilityCalculator>>(
            "shared_ptr_ObservationViabilityCalculator")
        .function("isObservationViable",
            &tom::ObservationViabilityCalculator::isObservationViable);

    // ObservationSimulatorBase class
    class_<tom::ObservationSimulatorBase<double, double>>(
        "estimation_observable_models_observables_simulation_ObservationSimulatorBase")
        .smart_ptr<std::shared_ptr<tom::ObservationSimulatorBase<double, double>>>(
            "shared_ptr_ObservationSimulatorBase");

    // ObservationSimulator class
    class_<tom::ObservationSimulator<1, double, double>,
           base<tom::ObservationSimulatorBase<double, double>>>(
        "estimation_observable_models_observables_simulation_ObservationSimulator1")
        .smart_ptr<std::shared_ptr<tom::ObservationSimulator<1, double, double>>>(
            "shared_ptr_ObservationSimulator1");

    class_<tom::ObservationSimulator<2, double, double>,
           base<tom::ObservationSimulatorBase<double, double>>>(
        "estimation_observable_models_observables_simulation_ObservationSimulator2")
        .smart_ptr<std::shared_ptr<tom::ObservationSimulator<2, double, double>>>(
            "shared_ptr_ObservationSimulator2");

    class_<tom::ObservationSimulator<3, double, double>,
           base<tom::ObservationSimulatorBase<double, double>>>(
        "estimation_observable_models_observables_simulation_ObservationSimulator3")
        .smart_ptr<std::shared_ptr<tom::ObservationSimulator<3, double, double>>>(
            "shared_ptr_ObservationSimulator3");

    class_<tom::ObservationSimulator<6, double, double>,
           base<tom::ObservationSimulatorBase<double, double>>>(
        "estimation_observable_models_observables_simulation_ObservationSimulator6")
        .smart_ptr<std::shared_ptr<tom::ObservationSimulator<6, double, double>>>(
            "shared_ptr_ObservationSimulator6");
}

#endif
