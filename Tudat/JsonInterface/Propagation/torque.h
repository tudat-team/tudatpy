/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_TORQUE_H
#define TUDAT_JSONINTERFACE_TORQUE_H

#include "Tudat/SimulationSetup/PropagationSetup/torqueSettings.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Map of `AvailableTorque`s string representations.
static std::map< AvailableTorque, std::string > torqueTypes =
{
    { underfined_torque, "undefined" },
    { second_order_gravitational_torque, "secondOrderGravitational" },
    { aerodynamic_torque, "aerodynamic" }
};

//! `AvailableTorque`s not supported by `json_interface`.
static std::vector< AvailableTorque > unsupportedTorqueTypes = { };

//! Convert `AvailableTorque` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AvailableTorque& torqueType )
{
    jsonObject = json_interface::stringFromEnum( torqueType, torqueTypes );
}

//! Convert `json` to `AvailableTorque`.
inline void from_json( const nlohmann::json& jsonObject, AvailableTorque& torqueType )
{
    torqueType = json_interface::enumFromString( jsonObject, torqueTypes );
}

} // namespace basic_astrodynamics


namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `TorqueSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< TorqueSettings >& torqueSettings );

//! Create a shared pointer to a `TorqueSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< TorqueSettings >& torqueSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_TORQUE_H
