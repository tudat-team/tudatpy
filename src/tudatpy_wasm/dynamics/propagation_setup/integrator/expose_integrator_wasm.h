/*    Copyright (c) 2010-2025, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Integrator settings module bindings for WASM.
 *    Mirrors: src/tudatpy/dynamics/propagation_setup/integrator/expose_integrator.cpp
 */

#ifndef TUDATPY_WASM_INTEGRATOR_H
#define TUDATPY_WASM_INTEGRATOR_H

#ifdef __EMSCRIPTEN__

namespace tudatpy_wasm {
namespace dynamics {
namespace propagation_setup {
namespace integrator {

// Embind bindings are registered via EMSCRIPTEN_BINDINGS macro

} // namespace integrator
} // namespace propagation_setup
} // namespace dynamics
} // namespace tudatpy_wasm

#endif // __EMSCRIPTEN__

#endif // TUDATPY_WASM_INTEGRATOR_H
