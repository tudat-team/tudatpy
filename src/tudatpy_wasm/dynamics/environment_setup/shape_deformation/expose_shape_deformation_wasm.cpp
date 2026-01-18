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

#include <tudat/simulation/environment_setup/createBodyDeformationModel.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup_shape_deformation")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_shape_deformation) {
    using namespace emscripten;

    // BodyDeformationSettings base class
    class_<tss::BodyDeformationSettings>(
        "dynamics_environment_setup_shape_deformation_BodyDeformationSettings")
        .smart_ptr<std::shared_ptr<tss::BodyDeformationSettings>>(
            "shared_ptr_BodyDeformationSettings");

    // BasicSolidBodyDeformationSettings derived class
    class_<tss::BasicSolidBodyDeformationSettings, base<tss::BodyDeformationSettings>>(
        "dynamics_environment_setup_shape_deformation_BasicSolidBodyDeformationSettings")
        .smart_ptr<std::shared_ptr<tss::BasicSolidBodyDeformationSettings>>(
            "shared_ptr_BasicSolidBodyDeformationSettings");

    // Factory functions
    function("dynamics_environment_setup_shape_deformation_basic_solid_body_tidal",
        &tss::basicTidalBodyShapeDeformation);

    function("dynamics_environment_setup_shape_deformation_degree_two_basic_solid_body_tidal",
        &tss::degreeTwoBasicTidalBodyShapeDeformation);

    function("dynamics_environment_setup_shape_deformation_iers_2010_solid_body_tidal",
        &tss::iers2010TidalBodyShapeDeformation);

    function("dynamics_environment_setup_shape_deformation_pole_tidal",
        &tss::poleTideBodyShapeDeformation);

    function("dynamics_environment_setup_shape_deformation_ocean_tidal",
        &tss::oceanTideBodyShapeDeformation);
}

#endif
