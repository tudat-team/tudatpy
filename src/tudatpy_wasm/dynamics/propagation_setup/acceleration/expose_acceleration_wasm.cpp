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
#include "../../../eigen_wasm.h"

#include <tudat/simulation/propagation_setup/accelerationSettings.h>
#include <tudat/astro/basic_astro/accelerationModelTypes.h>

namespace tba = tudat::basic_astrodynamics;
namespace tss = tudat::simulation_setup;

WASM_MODULE_PATH("dynamics_propagation_setup_acceleration")

EMSCRIPTEN_BINDINGS(tudatpy_dynamics_propagation_setup_acceleration) {
    using namespace emscripten;

    // ========================================================================
    // AvailableAcceleration enum - Full list matching Python bindings
    // ========================================================================
    enum_<tba::AvailableAcceleration>(
        "dynamics_propagation_setup_acceleration_AvailableAcceleration")
        .value("undefined_acceleration_type", tba::undefined_acceleration)
        .value("point_mass_gravity_type", tba::point_mass_gravity)
        .value("aerodynamic_type", tba::aerodynamic)
        .value("cannonball_radiation_pressure_type", tba::cannon_ball_radiation_pressure)
        .value("spherical_harmonic_gravity_type", tba::spherical_harmonic_gravity)
        .value("mutual_spherical_harmonic_gravity_type", tba::mutual_spherical_harmonic_gravity)
        .value("polyhedron_gravity_type", tba::polyhedron_gravity)
        .value("ring_gravity_type", tba::ring_gravity)
        .value("thrust_acceleration_type", tba::thrust_acceleration)
        .value("relativistic_correction_acceleration_type", tba::relativistic_correction_acceleration)
        .value("empirical_acceleration_type", tba::empirical_acceleration)
        .value("rtg_acceleration_type", tba::rtg_acceleration)
        .value("direct_tidal_dissipation_in_central_body_acceleration_type", tba::direct_tidal_dissipation_in_central_body_acceleration)
        .value("direct_tidal_dissipation_in_orbiting_body_acceleration_type", tba::direct_tidal_dissipation_in_orbiting_body_acceleration)
        .value("quasi_impulsive_shots_acceleration_type", tba::momentum_wheel_desaturation_acceleration)
        .value("custom_acceleration_type", tba::custom_acceleration)
        .value("radiation_pressure_type", tba::radiation_pressure)
        .value("einstein_infeld_hoffmann_acceleration_type", tba::einstein_infeld_hoffmann_acceleration)
        .value("yarkovsky_acceleration_type", tba::yarkovsky_acceleration);

    // ========================================================================
    // AccelerationSettings base class and derived classes
    // ========================================================================

    // AccelerationSettings base class
    class_<tss::AccelerationSettings>(
        "dynamics_propagation_setup_acceleration_AccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::AccelerationSettings>>(
            "shared_ptr_AccelerationSettings");

    // SphericalHarmonicAccelerationSettings derived class
    class_<tss::SphericalHarmonicAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_SphericalHarmonicAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::SphericalHarmonicAccelerationSettings>>(
            "shared_ptr_SphericalHarmonicAccelerationSettings")
        .property("maximumDegree", &tss::SphericalHarmonicAccelerationSettings::maximumDegree_)
        .property("maximumOrder", &tss::SphericalHarmonicAccelerationSettings::maximumOrder_);

    // MutualSphericalHarmonicAccelerationSettings derived class
    class_<tss::MutualSphericalHarmonicAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_MutualSphericalHarmonicAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::MutualSphericalHarmonicAccelerationSettings>>(
            "shared_ptr_MutualSphericalHarmonicAccelerationSettings")
        .property("maximumDegreeOfBodyExertingAcceleration",
            &tss::MutualSphericalHarmonicAccelerationSettings::maximumDegreeOfBodyExertingAcceleration_)
        .property("maximumOrderOfBodyExertingAcceleration",
            &tss::MutualSphericalHarmonicAccelerationSettings::maximumOrderOfBodyExertingAcceleration_)
        .property("maximumDegreeOfBodyUndergoingAcceleration",
            &tss::MutualSphericalHarmonicAccelerationSettings::maximumDegreeOfBodyUndergoingAcceleration_)
        .property("maximumOrderOfBodyUndergoingAcceleration",
            &tss::MutualSphericalHarmonicAccelerationSettings::maximumOrderOfBodyUndergoingAcceleration_)
        .property("maximumDegreeOfCentralBody",
            &tss::MutualSphericalHarmonicAccelerationSettings::maximumDegreeOfCentralBody_)
        .property("maximumOrderOfCentralBody",
            &tss::MutualSphericalHarmonicAccelerationSettings::maximumOrderOfCentralBody_);

    // EmpiricalAccelerationSettings derived class
    class_<tss::EmpiricalAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_EmpiricalAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::EmpiricalAccelerationSettings>>(
            "shared_ptr_EmpiricalAccelerationSettings")
        .property("constantAcceleration", &tss::EmpiricalAccelerationSettings::constantAcceleration_)
        .property("sineAcceleration", &tss::EmpiricalAccelerationSettings::sineAcceleration_)
        .property("cosineAcceleration", &tss::EmpiricalAccelerationSettings::cosineAcceleration_);

    // RTGAccelerationSettings derived class
    class_<tss::RTGAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_RTGAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::RTGAccelerationSettings>>(
            "shared_ptr_RTGAccelerationSettings")
        .property("bodyFixedForceVectorAtReferenceEpoch",
            &tss::RTGAccelerationSettings::bodyFixedForceVectorAtReferenceEpoch_)
        .property("decayScaleFactor",
            &tss::RTGAccelerationSettings::decayScaleFactor_)
        .property("referenceEpoch",
            &tss::RTGAccelerationSettings::referenceEpoch_);

    // RelativisticAccelerationCorrectionSettings derived class
    class_<tss::RelativisticAccelerationCorrectionSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_RelativisticAccelerationCorrectionSettings")
        .smart_ptr<std::shared_ptr<tss::RelativisticAccelerationCorrectionSettings>>(
            "shared_ptr_RelativisticAccelerationCorrectionSettings")
        .property("calculateSchwarzschildCorrection",
            &tss::RelativisticAccelerationCorrectionSettings::calculateSchwarzschildCorrection_)
        .property("calculateLenseThirringCorrection",
            &tss::RelativisticAccelerationCorrectionSettings::calculateLenseThirringCorrection_)
        .property("calculateDeSitterCorrection",
            &tss::RelativisticAccelerationCorrectionSettings::calculateDeSitterCorrection_);

    // CustomAccelerationSettings derived class
    class_<tss::CustomAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_CustomAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::CustomAccelerationSettings>>(
            "shared_ptr_CustomAccelerationSettings");

    // DirectTidalDissipationAccelerationSettings derived class
    class_<tss::DirectTidalDissipationAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_DirectTidalDissipationAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::DirectTidalDissipationAccelerationSettings>>(
            "shared_ptr_DirectTidalDissipationAccelerationSettings")
        .property("k2LoveNumber", &tss::DirectTidalDissipationAccelerationSettings::k2LoveNumber_)
        .property("timeLag", &tss::DirectTidalDissipationAccelerationSettings::timeLag_)
        .property("inverseTidalQualityFactor", &tss::DirectTidalDissipationAccelerationSettings::inverseTidalQualityFactor_)
        .property("tidalPeriod", &tss::DirectTidalDissipationAccelerationSettings::tidalPeriod_)
        .property("includeDirectRadialComponent", &tss::DirectTidalDissipationAccelerationSettings::includeDirectRadialComponent_)
        .property("useTideRaisedOnPlanet", &tss::DirectTidalDissipationAccelerationSettings::useTideRaisedOnPlanet_)
        .property("explicitLibraionalTideOnSatellite", &tss::DirectTidalDissipationAccelerationSettings::explicitLibraionalTideOnSatellite_);

    // MomentumWheelDesaturationAccelerationSettings derived class (quasi-impulsive shots)
    class_<tss::MomentumWheelDesaturationAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_MomentumWheelDesaturationAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::MomentumWheelDesaturationAccelerationSettings>>(
            "shared_ptr_MomentumWheelDesaturationAccelerationSettings")
        .property("thrustMidTimes", &tss::MomentumWheelDesaturationAccelerationSettings::thrustMidTimes_)
        .property("totalManeuverTime", &tss::MomentumWheelDesaturationAccelerationSettings::totalManeuverTime_)
        .property("maneuverRiseTime", &tss::MomentumWheelDesaturationAccelerationSettings::maneuverRiseTime_);

    // ThrustAccelerationSettings derived class
    class_<tss::ThrustAccelerationSettings, base<tss::AccelerationSettings>>(
        "dynamics_propagation_setup_acceleration_ThrustAccelerationSettings")
        .smart_ptr<std::shared_ptr<tss::ThrustAccelerationSettings>>(
            "shared_ptr_ThrustAccelerationSettings");

    // ========================================================================
    // Factory functions - matching Python bindings exactly
    // ========================================================================

    // Point mass gravity
    function("dynamics_propagation_setup_acceleration_point_mass_gravity",
        &tss::pointMassGravityAcceleration);

    // Einstein Infeld Hofmann (EIH) relativistic gravity
    function("dynamics_propagation_setup_acceleration_einstein_infeld_hofmann",
        &tss::einsteinInfledHoffmannGravityAcceleration);

    // Aerodynamic
    function("dynamics_propagation_setup_acceleration_aerodynamic",
        &tss::aerodynamicAcceleration);

    // Radiation pressure (general)
    function("dynamics_propagation_setup_acceleration_radiation_pressure",
        &tss::radiationPressureAcceleration);

    // Cannonball radiation pressure (legacy interface)
    function("dynamics_propagation_setup_acceleration_cannonball_radiation_pressure",
        &tss::cannonBallRadiationPressureAcceleration);

    // Spherical harmonic gravity
    function("dynamics_propagation_setup_acceleration_spherical_harmonic_gravity",
        &tss::sphericalHarmonicAcceleration);

    // Mutual spherical harmonic gravity
    function("dynamics_propagation_setup_acceleration_mutual_spherical_harmonic_gravity",
        &tss::mutualSphericalHarmonicAcceleration);

    // Polyhedron gravity
    function("dynamics_propagation_setup_acceleration_polyhedron_gravity",
        &tss::polyhedronAcceleration);

    // Ring gravity
    function("dynamics_propagation_setup_acceleration_ring_gravity",
        &tss::ringAcceleration);

    // Relativistic correction
    function("dynamics_propagation_setup_acceleration_relativistic_correction",
        &tss::relativisticAccelerationCorrection);

    // Empirical
    function("dynamics_propagation_setup_acceleration_empirical",
        &tss::empiricalAcceleration);

    // RTG (Radioisotope Thermoelectric Generator)
    function("dynamics_propagation_setup_acceleration_rtg",
        &tss::rtgAcceleration);

    // Yarkovsky effect
    function("dynamics_propagation_setup_acceleration_yarkovsky",
        &tss::yarkovskyAcceleration);

    // Direct tidal dissipation (time lag variant)
    function("dynamics_propagation_setup_acceleration_direct_tidal_dissipation_acceleration",
        &tss::directTidalDissipationAcceleration);

    // Direct tidal dissipation (inverse Q variant)
    function("dynamics_propagation_setup_acceleration_direct_tidal_dissipation_acceleration_from_inv_q",
        &tss::directTidalDissipationAccelerationFromInvQ);

    // Quasi-impulsive shots (momentum wheel desaturation)
    function("dynamics_propagation_setup_acceleration_quasi_impulsive_shots_acceleration",
        &tss::momentumWheelDesaturationAcceleration);

    // Thrust from engines (list of engine names)
    function("dynamics_propagation_setup_acceleration_thrust_from_engines",
        &tss::thrustAcceleration);

    // Thrust from single engine
    function("dynamics_propagation_setup_acceleration_thrust_from_engine",
        &tss::thrustAccelerationFromSingleEngine);

    // Thrust from all engines
    function("dynamics_propagation_setup_acceleration_thrust_from_all_engines",
        &tss::thrustAccelerationFromAllEngines);
}

#endif
