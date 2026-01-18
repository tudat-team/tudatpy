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

#include <tudat/math/geometric.h>

namespace tgs = tudat::geometric_shapes;

WASM_MODULE_PATH("math_geometry")

EMSCRIPTEN_BINDINGS(tudatpy_math_geometry) {
    using namespace emscripten;

    // SurfaceGeometry base class (note: in tudat:: namespace, not geometric_shapes)
    class_<tudat::SurfaceGeometry>("math_geometry_SurfaceGeometry")
        .smart_ptr<std::shared_ptr<tudat::SurfaceGeometry>>("shared_ptr_SurfaceGeometry");

    // CompositeSurfaceGeometry class
    class_<tgs::CompositeSurfaceGeometry, base<tudat::SurfaceGeometry>>("math_geometry_CompositeSurfaceGeometry")
        .smart_ptr<std::shared_ptr<tgs::CompositeSurfaceGeometry>>("shared_ptr_CompositeSurfaceGeometry");

    // Capsule class
    class_<tgs::Capsule, base<tgs::CompositeSurfaceGeometry>>("math_geometry_Capsule")
        .smart_ptr<std::shared_ptr<tgs::Capsule>>("shared_ptr_Capsule")
        .constructor<double, double, double, double, double>()
        .function("getMiddleRadius", &tgs::Capsule::getMiddleRadius)
        .function("getLength", &tgs::Capsule::getLength);
}

#endif
