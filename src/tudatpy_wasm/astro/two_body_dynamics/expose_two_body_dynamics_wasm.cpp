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

#include <tudat/astro/mission_segments/lambertTargeterIzzo.h>
#include <tudat/astro/basic_astro/keplerPropagator.h>
#include <tudat/astro/basic_astro/orbitalElementConversions.h>
#include <tudat/astro/basic_astro/astrodynamicsFunctions.h>

namespace tms = tudat::mission_segments;
namespace toec = tudat::orbital_element_conversions;
namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("astro_two_body_dynamics")

namespace {
using namespace tudatpy_wasm;

Vector6dWrapper propagateKeplerOrbit(const Vector6dWrapper& initialKeplerElements,
                                      double propagationTime,
                                      double gravitationalParameter) {
    return Vector6dWrapper(toec::propagateKeplerOrbit<double>(
        initialKeplerElements.eigen(), propagationTime, gravitationalParameter));
}

}

EMSCRIPTEN_BINDINGS(tudatpy_astro_two_body_dynamics) {
    using namespace emscripten;

    // LambertTargeter base class
    class_<tms::LambertTargeter>("astro_two_body_dynamics_LambertTargeter")
        .smart_ptr<std::shared_ptr<tms::LambertTargeter>>("shared_ptr_LambertTargeter");

    // LambertTargeterIzzo
    class_<tms::LambertTargeterIzzo, base<tms::LambertTargeter>>("astro_two_body_dynamics_LambertTargeterIzzo")
        .smart_ptr<std::shared_ptr<tms::LambertTargeterIzzo>>("shared_ptr_LambertTargeterIzzo")
        .constructor<const Eigen::Vector3d&, const Eigen::Vector3d&, double, double, bool, double, int>()
        .function("getInertialVelocityAtDeparture", &tms::LambertTargeterIzzo::getInertialVelocityAtDeparture)
        .function("getInertialVelocityAtArrival", &tms::LambertTargeterIzzo::getInertialVelocityAtArrival);

    // Free functions
    function("astro_two_body_dynamics_propagate_kepler_orbit", &propagateKeplerOrbit);
    function("astro_two_body_dynamics_compute_kepler_orbit_period", &tba::computeKeplerOrbitalPeriod);
    function("astro_two_body_dynamics_compute_kepler_mean_motion", &tba::computeKeplerMeanMotion);
}

#endif
