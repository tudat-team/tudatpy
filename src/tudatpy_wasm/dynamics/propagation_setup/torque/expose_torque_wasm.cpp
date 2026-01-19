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

#include <tudat/simulation/propagation_setup/torqueSettings.h>

namespace tss = tudat::simulation_setup;
namespace tba = tudat::basic_astrodynamics;

WASM_MODULE_PATH("dynamics_propagation_setup_torque")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_torque) {
    using namespace emscripten;

    // AvailableTorque enum - Full list matching Python bindings
    enum_<tba::AvailableTorque>("dynamics_propagation_setup_torque_AvailableTorque")
        .value("torque_free", tba::torque_free)
        .value("underfined_torque", tba::underfined_torque)
        .value("inertial_torque", tba::inertial_torque)
        .value("second_order_gravitational_torque", tba::second_order_gravitational_torque)
        .value("spherical_harmonic_gravitational_torque", tba::spherical_harmonic_gravitational_torque)
        .value("aerodynamic_torque", tba::aerodynamic_torque)
        .value("radiation_pressure_torque", tba::radiation_pressure_torque)
        .value("dissipative_torque", tba::dissipative_torque)
        .value("custom_torque", tba::custom_torque);

    // TorqueSettings base class
    class_<tss::TorqueSettings>("dynamics_propagation_setup_torque_TorqueSettings")
        .smart_ptr<std::shared_ptr<tss::TorqueSettings>>("shared_ptr_TorqueSettings")
        .property("torqueType", &tss::TorqueSettings::torqueType_);

    // SphericalHarmonicTorqueSettings
    class_<tss::SphericalHarmonicTorqueSettings, base<tss::TorqueSettings>>(
        "dynamics_propagation_setup_torque_SphericalHarmonicTorqueSettings")
        .smart_ptr<std::shared_ptr<tss::SphericalHarmonicTorqueSettings>>(
            "shared_ptr_SphericalHarmonicTorqueSettings");

    // CustomTorqueSettings derived class
    class_<tss::CustomTorqueSettings, base<tss::TorqueSettings>>(
        "dynamics_propagation_setup_torque_CustomTorqueSettings")
        .smart_ptr<std::shared_ptr<tss::CustomTorqueSettings>>(
            "shared_ptr_CustomTorqueSettings");

    // Factory functions
    function("dynamics_propagation_setup_torque_aerodynamic",
        &tss::aerodynamicTorque);

    function("dynamics_propagation_setup_torque_second_degree_gravitational",
        &tss::secondDegreeGravitationalTorque);

    function("dynamics_propagation_setup_torque_spherical_harmonic_gravitational",
        &tss::sphericalHarmonicGravitationalTorque);

    function("dynamics_propagation_setup_torque_radiation_pressure",
        &tss::radiationPressureTorque);

    function("dynamics_propagation_setup_torque_dissipative",
        &tss::dissipativeTorque);

    // Note: custom_torque is not exposed because it requires std::function which is
    // complex to bind in Emscripten. Users can create custom torques on the JS side.
}

#endif
