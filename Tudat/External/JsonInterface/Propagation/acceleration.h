/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Map of `AvailableAcceleration`s supported by `json_interface`.
/*!
 * Map of `AvailableAcceleration`s supported by `json_interface`.
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
static std::map< std::string, AvailableAcceleration > availableAccelerations =
{
    { "undefined",                               undefined_acceleration },
    { "pointMassGravity",                        point_mass_gravity },
    { "aerodynamic",                             aerodynamic },
    { "cannonBallRadiationPressure",             cannon_ball_radiation_pressure },
    { "sphericalHarmonicGravity",                spherical_harmonic_gravity },
    { "mutualSphericalHarmonicGravity",          mutual_spherical_harmonic_gravity },
    { "thirdBodyPointMassGravity",               third_body_point_mass_gravity },
    { "thirdBodySphericalHarmonicGravity",       third_body_spherical_harmonic_gravity },
    { "thirdBodyMutualSphericalHarmonicGravity", third_body_mutual_spherical_harmonic_gravity },
    // FIXME: { "thrust",                                  thrust_acceleration },
    { "relativisticCorrection",                  relativistic_correction_acceleration },
    { "empirical",                               empirical_acceleration }
};

//! Convert `AvailableAcceleration` to `json`.
void to_json( json& jsonObject, const AvailableAcceleration& availableAcceleration );

//! Convert `json` to `AvailableAcceleration`.
void from_json( const json& jsonObject, AvailableAcceleration& availableAcceleration );

} // namespace basic_astrodynamics


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `AccelerationSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< AccelerationSettings >& accelerationSettings );

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< AccelerationSettings >& accelerationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ACCELERATION_H
