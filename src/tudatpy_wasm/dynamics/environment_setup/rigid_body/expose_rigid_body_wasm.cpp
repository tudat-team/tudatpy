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
#include "../../../eigen_wasm.h"
#include "../../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/createBodies.h>

namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_environment_setup_rigid_body")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment_setup_rigid_body) {
    using namespace emscripten;

    // RigidBodyPropertiesSettings base class
    class_<tss::RigidBodyPropertiesSettings>(
        "dynamics_environment_setup_rigid_body_RigidBodyPropertiesSettings")
        .smart_ptr<std::shared_ptr<tss::RigidBodyPropertiesSettings>>(
            "shared_ptr_RigidBodyPropertiesSettings")
        .function("getRigidBodyPropertiesType",
            &tss::RigidBodyPropertiesSettings::getRigidBodyPropertiesType);

    // ========================================================================
    // Factory functions
    // ========================================================================

    // Constant rigid body properties
    function("dynamics_environment_setup_rigid_body_constant_rigid_body_properties",
        &tss::constantRigidBodyPropertiesSettings);

    // Custom time-dependent rigid body properties
    function("dynamics_environment_setup_rigid_body_custom_time_dependent_rigid_body_properties",
        &tss::fromFunctionRigidBodyPropertiesSettings);

    // Custom mass-dependent rigid body properties
    function("dynamics_environment_setup_rigid_body_custom_mass_dependent_rigid_body_properties",
        &tss::massDependentMassDistributionSettings);
}

#endif
