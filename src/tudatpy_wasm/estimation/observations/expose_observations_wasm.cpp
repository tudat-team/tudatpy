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

#include <tudat/simulation/estimation_setup/singleObservationSet.h>
#include <tudat/simulation/estimation_setup/observationCollection.h>

namespace tom = tudat::observation_models;
namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("estimation_observations")

EMSCRIPTEN_BINDINGS(tudatpy_estimation_observations) {
    using namespace emscripten;

    // SingleObservationSet class
    class_<tom::SingleObservationSet<double, double>>("estimation_observations_SingleObservationSet")
        .smart_ptr<std::shared_ptr<tom::SingleObservationSet<double, double>>>(
            "shared_ptr_SingleObservationSet")
        .function("getObservableType", &tom::SingleObservationSet<double, double>::getObservableType)
        .function("getLinkEnds", &tom::SingleObservationSet<double, double>::getLinkEnds)
        .function("getObservationTimes", &tom::SingleObservationSet<double, double>::getObservationTimes)
        .function("getObservations", &tom::SingleObservationSet<double, double>::getObservations)
        .function("getNumberOfObservables", &tom::SingleObservationSet<double, double>::getNumberOfObservables);

    // ObservationCollection class
    class_<tom::ObservationCollection<double, double>>("estimation_observations_ObservationCollection")
        .smart_ptr<std::shared_ptr<tom::ObservationCollection<double, double>>>(
            "shared_ptr_ObservationCollection")
        .function("getTotalObservableSize", &tom::ObservationCollection<double, double>::getTotalObservableSize)
        .function("getConcatenatedTimeVector", &tom::ObservationCollection<double, double>::getConcatenatedTimeVector)
        .function("getConcatenatedObservations", &tom::ObservationCollection<double, double>::getConcatenatedObservations);
}

#endif
