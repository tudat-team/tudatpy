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

#include <tudat/simulation/propagation_setup/accelerationSettings.h>
#include <tudat/astro/basic_astro/accelerationModelTypes.h>

namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_propagation_setup_acceleration")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_acceleration) {
    using namespace emscripten;

    // AvailableAcceleration enum
    enum_<tba::AvailableAcceleration>(
        "dynamics_propagation_setup_acceleration_AvailableAcceleration")
        .value("undefined_acceleration_type", tba::undefined_acceleration)
        .value("point_mass_gravity_type", tba::point_mass_gravity)
        .value("aerodynamic_type", tba::aerodynamic)
        .value("cannonball_radiation_pressure_type", tba::cannon_ball_radiation_pressure)
        .value("spherical_harmonic_gravity_type", tba::spherical_harmonic_gravity)
        .value("mutual_spherical_harmonic_gravity_type", tba::mutual_spherical_harmonic_gravity)
        .value("polyhedron_gravity_type", tba::polyhedron_gravity)
        .value("thrust_acceleration_type", tba::thrust_acceleration)
        .value("relativistic_correction_acceleration_type", tba::relativistic_correction_acceleration)
        .value("empirical_acceleration_type", tba::empirical_acceleration)
        .value("custom_acceleration_type", tba::custom_acceleration)
        .value("radiation_pressure_type", tba::radiation_pressure);

    // AccelerationSettings base class
    class_<tss::AccelerationSettings>(
        "dynamics_propagation_setup_acceleration_AccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::AccelerationSettings>>(
            "shared_ptr_AccelerationSettings");

    // SphericalHarmonicAccelerationSettings derived class
    class_<tss::SphericalHarmonicAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_SphericalHarmonicAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::SphericalHarmonicAccelerationSettings>>(
            "shared_ptr_SphericalHarmonicAccelerationSettings");

    // MutualSphericalHarmonicAccelerationSettings derived class
    class_<tss::MutualSphericalHarmonicAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_MutualSphericalHarmonicAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::MutualSphericalHarmonicAccelerationSettings>>(
            "shared_ptr_MutualSphericalHarmonicAccelerationSettings");

    // EmpiricalAccelerationSettings derived class
    class_<tss::EmpiricalAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_EmpiricalAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::EmpiricalAccelerationSettings>>(
            "shared_ptr_EmpiricalAccelerationSettings");

    // ThrustAccelerationSettings derived class
    class_<tss::ThrustAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_ThrustAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustAccelerationSettings>>(
            "shared_ptr_ThrustAccelerationSettings");

    // Factory functions
    function("dynamics_propagation_setup_acceleration_point_mass_gravity",
        &tss::pointMassGravityAcceleration);

    function("dynamics_propagation_setup_acceleration_spherical_harmonic_gravity",
        &tss::sphericalHarmonicAcceleration);

    function("dynamics_propagation_setup_acceleration_mutual_spherical_harmonic_gravity",
        &tss::mutualSphericalHarmonicAcceleration);

    function("dynamics_propagation_setup_acceleration_aerodynamic",
        &tss::aerodynamicAcceleration);

    function("dynamics_propagation_setup_acceleration_radiation_pressure",
        &tss::radiationPressureAcceleration);

    function("dynamics_propagation_setup_acceleration_relativistic_correction",
        &tss::relativisticAccelerationCorrection);

    function("dynamics_propagation_setup_acceleration_empirical",
        &tss::empiricalAcceleration);

    function("dynamics_propagation_setup_acceleration_thrust_from_engines",
        &tss::thrustAcceleration);

    function("dynamics_propagation_setup_acceleration_polyhedron_gravity",
        &tss::polyhedronAcceleration);
}

#endif
