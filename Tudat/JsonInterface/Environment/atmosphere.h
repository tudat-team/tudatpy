/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_ATMOSPHERE_H
#define TUDAT_JSONINTERFACE_ATMOSPHERE_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `AtmosphereTypes` string representations.
static std::map< AtmosphereTypes, std::string > atmosphereTypes =
{
    { exponential_atmosphere, "exponential" },
    { tabulated_atmosphere, "tabulated" },
    { nrlmsise00, "nrlmsise00" }
};

//! `AtmosphereTypes` not supported by `json_interface`.
static std::vector< AtmosphereTypes > unsupportedAtmosphereTypes = { };

//! Convert `AtmosphereTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AtmosphereTypes& atmosphereType )
{
    jsonObject = json_interface::stringFromEnum( atmosphereType, atmosphereTypes );
}

//! Convert `json` to `AtmosphereTypes`.
inline void from_json( const nlohmann::json& jsonObject, AtmosphereTypes& atmosphereType )
{
    atmosphereType = json_interface::enumFromString( jsonObject, atmosphereTypes );
}

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< AtmosphereSettings >& atmosphereSettings );

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< AtmosphereSettings >& atmosphereSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ATMOSPHERE_H
