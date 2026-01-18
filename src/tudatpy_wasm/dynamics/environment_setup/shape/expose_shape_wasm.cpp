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
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/createBodyShapeModel.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup_shape")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_shape) {
    using namespace emscripten;

    // BodyShapeSettings base class
    class_<tss::BodyShapeSettings>("dynamics_environment_setup_shape_BodyShapeSettings")
        .smart_ptr<std::shared_ptr<tss::BodyShapeSettings>>("shared_ptr_BodyShapeSettings");

    // SphericalBodyShapeSettings derived class
    class_<tss::SphericalBodyShapeSettings, base<tss::BodyShapeSettings>>(
        "dynamics_environment_setup_shape_SphericalBodyShapeSettings")
        .smart_ptr<std::shared_ptr<tss::SphericalBodyShapeSettings>>(
            "shared_ptr_SphericalBodyShapeSettings")
        .function("getRadius", &tss::SphericalBodyShapeSettings::getRadius);

    // OblateSphericalBodyShapeSettings derived class
    class_<tss::OblateSphericalBodyShapeSettings, base<tss::BodyShapeSettings>>(
        "dynamics_environment_setup_shape_OblateSphericalBodyShapeSettings")
        .smart_ptr<std::shared_ptr<tss::OblateSphericalBodyShapeSettings>>(
            "shared_ptr_OblateSphericalBodyShapeSettings")
        .function("getEquatorialRadius", &tss::OblateSphericalBodyShapeSettings::getEquatorialRadius)
        .function("getFlattening", &tss::OblateSphericalBodyShapeSettings::getFlattening);

    // PolyhedronBodyShapeSettings derived class
    class_<tss::PolyhedronBodyShapeSettings, base<tss::BodyShapeSettings>>(
        "dynamics_environment_setup_shape_PolyhedronBodyShapeSettings")
        .smart_ptr<std::shared_ptr<tss::PolyhedronBodyShapeSettings>>(
            "shared_ptr_PolyhedronBodyShapeSettings");

    // HybridBodyShapeSettings derived class
    class_<tss::HybridBodyShapeSettings, base<tss::BodyShapeSettings>>(
        "dynamics_environment_setup_shape_HybridBodyShapeSettings")
        .smart_ptr<std::shared_ptr<tss::HybridBodyShapeSettings>>(
            "shared_ptr_HybridBodyShapeSettings")
        .function("getSwitchoverAltitude", &tss::HybridBodyShapeSettings::getSwitchoverAltitude);

    // Factory functions
    function("dynamics_environment_setup_shape_spherical",
        &tss::sphericalBodyShapeSettings);

    function("dynamics_environment_setup_shape_spherical_spice",
        &tss::fromSpiceSphericalBodyShapeSettings);

    function("dynamics_environment_setup_shape_oblate_spherical",
        &tss::oblateSphericalBodyShapeSettings);

    function("dynamics_environment_setup_shape_polyhedron",
        &tss::polyhedronBodyShapeSettings);

    function("dynamics_environment_setup_shape_hybrid",
        &tss::hybridBodyShapeSettings);
}

#endif
