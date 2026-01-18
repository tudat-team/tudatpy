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

#include <tudat/astro/basic_astro/missionGeometry.h>

namespace tmg = tudat::mission_geometry;

WASM_MODULE_PATH("astro_fundamentals")

namespace {
using namespace tudatpy_wasm;

double computeShadowFunction(const Vector3dWrapper& occultedBodyPosition,
                              double occultedBodyRadius,
                              const Vector3dWrapper& occultingBodyPosition,
                              double occultingBodyRadius,
                              const Vector3dWrapper& satellitePosition) {
    return tmg::computeShadowFunction(
        occultedBodyPosition.eigen(),
        occultedBodyRadius,
        occultingBodyPosition.eigen(),
        occultingBodyRadius,
        satellitePosition.eigen());
}

}

EMSCRIPTEN_BINDINGS(tudatpy_astro_fundamentals) {
    using namespace emscripten;

    function("astro_fundamentals_compute_shadow_function", &computeShadowFunction);
}

#endif
