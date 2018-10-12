/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_RADIATIONPRESSURE_H
#define TUDAT_JSONINTERFACE_RADIATIONPRESSURE_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `RadiationPressureType`s string representations.
static std::map< RadiationPressureType, std::string > radiationPressureTypes =
{
    { cannon_ball, "cannonBall" }
};

//! `RadiationPressureType`s not supported by `json_interface`.
static std::vector< RadiationPressureType > unsupportedRadiationPressureTypes = { };

//! Convert `RadiationPressureType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const RadiationPressureType& radiationPressureType )
{
    jsonObject = json_interface::stringFromEnum( radiationPressureType, radiationPressureTypes );
}

//! Convert `json` to `RadiationPressureType`.
inline void from_json( const nlohmann::json& jsonObject, RadiationPressureType& radiationPressureType )
{
    radiationPressureType = json_interface::enumFromString( jsonObject, radiationPressureTypes );
}

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void to_json( nlohmann::json& jsonObject,
              const std::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings );

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
void from_json( const nlohmann::json& jsonObject,
                std::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_RADIATIONPRESSURE_H
