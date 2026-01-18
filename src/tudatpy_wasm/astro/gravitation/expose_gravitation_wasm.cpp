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

#include <tudat/math/basic/legendrePolynomials.h>

namespace tbm = tudat::basic_mathematics;

WASM_MODULE_PATH("astro_gravitation")

EMSCRIPTEN_BINDINGS(tudatpy_astro_gravitation) {
    using namespace emscripten;

    function("astro_gravitation_legendre_normalization_factor",
        &tbm::calculateLegendreGeodesyNormalizationFactor);
}

#endif
