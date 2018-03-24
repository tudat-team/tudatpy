/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_ACCELERATION_H
#define TUDAT_JSONINTERFACE_ACCELERATION_H

#include "Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Map of `AvailableAcceleration`s string representations.
/*!
 * Map of `AvailableAcceleration`s string representations.
 *
 * There is a naming inconsistency between "pointMassGravity" and `central_gravity`,
 * and between "thirdBodyPointMassGravity" and `third_body_central_gravity`.
 *
 * The values "thirdBodyPointMassGravity", "thirdBodySphericalHarmonicGravity" and
 * "thirdBodyMutualSphericalHarmonicGravity" are supported by the JSON interface, and are equivalent to
 * "pointMassGravity", "sphericalHarmonicGravity" and "mutualSphericalHarmonicGravity", respectively.
 * However, a warning will be printed indicating that the body causing the acceleration may not be treated as a
 * third body, requesting the user to remove the part "thridBody" in order to silence the warning.
 */
static std::map< AvailableAcceleration, std::string > accelerationTypes =
{
    { undefined_acceleration, "undefined" },
    { point_mass_gravity, "pointMassGravity" },
    { aerodynamic, "aerodynamic" },
    { cannon_ball_radiation_pressure, "cannonBallRadiationPressure" },
    { spherical_harmonic_gravity, "sphericalHarmonicGravity" },
    { mutual_spherical_harmonic_gravity, "mutualSphericalHarmonicGravity" },
    { third_body_point_mass_gravity, "thirdBodyPointMassGravity" },
    { third_body_spherical_harmonic_gravity, "thirdBodySphericalHarmonicGravity" },
    { third_body_mutual_spherical_harmonic_gravity, "thirdBodyMutualSphericalHarmonicGravity" },
    { thrust_acceleration, "thrust" },
    { relativistic_correction_acceleration, "relativisticCorrection" },
    { empirical_acceleration, "empirical" }
};

//! `AvailableAcceleration`s not supported by `json_interface`.
static std::vector< AvailableAcceleration > unsupportedAccelerationTypes =
{
    third_body_point_mass_gravity,
    third_body_spherical_harmonic_gravity,
    third_body_mutual_spherical_harmonic_gravity
};

//! Convert `AvailableAcceleration` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AvailableAcceleration& accelerationType )
{
    jsonObject = json_interface::stringFromEnum( accelerationType, accelerationTypes );
}

//! Convert `json` to `AvailableAcceleration`.
inline void from_json( const nlohmann::json& jsonObject, AvailableAcceleration& accelerationType )
{
    accelerationType = json_interface::enumFromString( jsonObject, accelerationTypes );
}

} // namespace basic_astrodynamics


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< AccelerationSettings >& accelerationSettings );

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< AccelerationSettings >& accelerationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ACCELERATION_H
