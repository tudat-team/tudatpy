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
#include "../../shared_ptr_wasm.h"

#include <tudat/simulation/environment_setup/body.h>
#include <tudat/astro/ephemerides/ephemeris.h>
#include <tudat/astro/gravitation/gravityFieldModel.h>
#include <tudat/astro/aerodynamics/atmosphereModel.h>

namespace tss = tudat::simulation_setup;
namespace te = tudat::ephemerides;
namespace tg = tudat::gravitation;
namespace ta = tudat::aerodynamics;

WASM_MODULE_PATH("dynamics_environment")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_environment) {
    using namespace emscripten;

    // Ephemeris base class
    class_<te::Ephemeris>("dynamics_environment_Ephemeris")
        .smart_ptr<std::shared_ptr<te::Ephemeris>>("shared_ptr_Ephemeris")
        .function("getCartesianState", &te::Ephemeris::getCartesianState)
        .function("getReferenceFrameOrigin", &te::Ephemeris::getReferenceFrameOrigin)
        .function("getReferenceFrameOrientation", &te::Ephemeris::getReferenceFrameOrientation);

    // GravityFieldModel base class
    class_<tg::GravityFieldModel>("dynamics_environment_GravityFieldModel")
        .smart_ptr<std::shared_ptr<tg::GravityFieldModel>>("shared_ptr_GravityFieldModel")
        .function("getGravitationalParameter", &tg::GravityFieldModel::getGravitationalParameter);

    // AtmosphereModel base class
    class_<ta::AtmosphereModel>("dynamics_environment_AtmosphereModel")
        .smart_ptr<std::shared_ptr<ta::AtmosphereModel>>("shared_ptr_AtmosphereModel");

    // Body class
    class_<tss::Body>("dynamics_environment_Body")
        .smart_ptr<std::shared_ptr<tss::Body>>("shared_ptr_Body")
        .function("getState", &tss::Body::getState)
        .function("getPosition", &tss::Body::getPosition)
        .function("getVelocity", &tss::Body::getVelocity)
        .function("getEphemeris", &tss::Body::getEphemeris)
        .function("getGravityFieldModel", &tss::Body::getGravityFieldModel)
        .function("getAtmosphereModel", &tss::Body::getAtmosphereModel)
        .function("getBodyMass", &tss::Body::getBodyMass)
        .function("setConstantMass", &tss::Body::setConstantBodyMass);

    // SystemOfBodies class
    class_<tss::SystemOfBodies>("dynamics_environment_SystemOfBodies")
        .smart_ptr<std::shared_ptr<tss::SystemOfBodies>>("shared_ptr_SystemOfBodies")
        .function("getBody", &tss::SystemOfBodies::getBody)
        .function("at", &tss::SystemOfBodies::at)
        .function("getFrameOrigin", &tss::SystemOfBodies::getFrameOrigin)
        .function("getFrameOrientation", &tss::SystemOfBodies::getFrameOrientation);
}

#endif
