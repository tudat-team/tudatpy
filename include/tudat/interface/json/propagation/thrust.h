/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_THRUST_H
#define TUDAT_JSONINTERFACE_THRUST_H

#include "tudat/simulation/propagation_setup/accelerationSettings.h"
#include "tudat/simulation/propagation_setup/thrustSettings.h"
#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

// ThrustDirectionGuidanceTypes

//! Map of `ThrustDirectionGuidanceTypes` string representations.
static std::map< ThrustDirectionGuidanceTypes, std::string > thrustDirectionTypes =
{
    { colinear_with_state_segment_thrust_direction, "colinearWithStateSegment" },
    { thrust_direction_from_existing_body_orientation, "fromExistingBodyOrientation" },
    { custom_thrust_direction, "customDirection" },
    { custom_thrust_orientation, "customOrientation" }
};

//! `ThrustDirectionGuidanceTypes` not supported by `json_interface`.
static std::vector< ThrustDirectionGuidanceTypes > unsupportedThrustDirectionTypes =
{
    custom_thrust_direction,
    custom_thrust_orientation
};

//! Convert `ThrustDirectionGuidanceTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const ThrustDirectionGuidanceTypes& thrustDirectionType )
{
    jsonObject = json_interface::stringFromEnum( thrustDirectionType, thrustDirectionTypes );
}

//! Convert `json` to `ThrustDirectionGuidanceTypes`.
inline void from_json( const nlohmann::json& jsonObject, ThrustDirectionGuidanceTypes& thrustDirectionType )
{
    thrustDirectionType = json_interface::enumFromString( jsonObject, thrustDirectionTypes );
}


// ThrustDirectionGuidanceSettings

//! Create a `json` object from a shared pointer to a `ThrustDirectionGuidanceSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings );

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ThrustDirectionGuidanceSettings >& directionSettings );


// ThrustMagnitudeTypes

//! Map of `ThrustMagnitudeTypes` string representations.
static std::map< ThrustMagnitudeTypes, std::string > thrustMagnitudeTypes =
{
    { constant_thrust_magnitude, "constant" },
    { from_engine_properties_thrust_magnitude, "fromEngineProperties" },
    { thrust_magnitude_from_time_function, "timeDependent" },
    { thrust_magnitude_from_dependent_variables, "variableDependent" }
};

//! `ThrustMagnitudeTypes` not supported by `json_interface`.
static std::vector< ThrustMagnitudeTypes > unsupportedThrustMagnitudeTypes =
{
    thrust_magnitude_from_time_function,
    thrust_magnitude_from_dependent_variables
};

//! Convert `ThrustMagnitudeTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const ThrustMagnitudeTypes& thrustMagnitudeType )
{
    jsonObject = json_interface::stringFromEnum( thrustMagnitudeType, thrustMagnitudeTypes );
}

//! Convert `json` to `ThrustMagnitudeTypes`.
inline void from_json( const nlohmann::json& jsonObject, ThrustMagnitudeTypes& thrustMagnitudeType )
{
    thrustMagnitudeType = json_interface::enumFromString( jsonObject, thrustMagnitudeTypes );
}


// ThrustMagnitudeSettings

//! Create a `json` object from a shared pointer to a `ThrustMagnitudeSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ThrustMagnitudeSettings >& magnitudeSettings );

//! Create a shared pointer to a `AccelerationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ThrustMagnitudeSettings >& magnitudeSettings );



// Thrust

//! Create a `json` object from a shared pointer to a `ThrustAccelerationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< ThrustAccelerationSettings >& thrustAccelerationSettings );

//! Create a shared pointer to a `ThrustAccelerationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< ThrustAccelerationSettings >& thrustAccelerationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_THRUST_H
