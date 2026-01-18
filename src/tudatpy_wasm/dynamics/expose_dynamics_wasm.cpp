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
#include "../wasm_module.h"

// This file serves as the main entry point for the dynamics module bindings.
// All submodules are compiled separately and linked together.
// The module hierarchy is managed through prefixed function/class names.

WASM_MODULE_PATH("dynamics")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics) {
    using namespace emscripten;

    // The dynamics module contains the following submodules:
    // - environment: Body, SystemOfBodies, Ephemeris, etc.
    // - environment_setup: BodySettings, SystemOfBodySettings, and submodules:
    //   - atmosphere
    //   - ephemeris
    //   - gravity_field
    //   - gravity_field_variation
    //   - aerodynamic_coefficients
    //   - radiation_pressure
    //   - rotation_model
    //   - shape
    //   - shape_deformation
    //   - rigid_body
    //   - ground_station
    //   - vehicle_systems
    // - propagation_setup: Model creation functions, and submodules:
    //   - acceleration
    //   - dependent_variable
    //   - integrator
    //   - mass_rate
    //   - propagator
    //   - thrust
    //   - torque
    // - propagation: State propagation results and helper functions
    // - simulator: DynamicsSimulator, VariationalSimulator, etc.
    // - parameters: EstimatableParameterSet
    // - parameters_setup: EstimatableParameterSettings and factory functions

    // All submodule bindings are registered through their individual
    // EMSCRIPTEN_BINDINGS blocks with prefixed names.
}

#endif
