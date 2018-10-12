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

namespace aerodynamics
{

// AtmosphereIndependentVariables

//! Map of `AtmosphereIndependentVariables` string representations.
static std::map< AtmosphereIndependentVariables, std::string > atmosphereIndependentVariables =
{
    { altitude_dependent_atmosphere, "altitude" },
    { longitude_dependent_atmosphere, "longitude" },
    { latitude_dependent_atmosphere, "latitude" },
    { time_dependent_atmosphere, "time" }
};

//! `AtmosphereIndependentVariables` not supported by `json_interface`.
static std::vector< AtmosphereIndependentVariables > unsupportedAtmosphereIndependentVariables = { };

//! Convert `AtmosphereIndependentVariables` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AtmosphereIndependentVariables& atmosphereIndependentVariable )
{
    jsonObject = json_interface::stringFromEnum( atmosphereIndependentVariable, atmosphereIndependentVariables );
}

//! Convert `json` to `AtmosphereIndependentVariables`.
inline void from_json( const nlohmann::json& jsonObject, AtmosphereIndependentVariables& atmosphereIndependentVariable )
{
    atmosphereIndependentVariable = json_interface::enumFromString( jsonObject, atmosphereIndependentVariables );
}


// AtmosphereDependentVariables

//! Map of `AtmosphereDependentVariables` string representations.
static std::map< AtmosphereDependentVariables, std::string > atmosphereDependentVariables =
{
    { density_dependent_atmosphere, "density" },
    { pressure_dependent_atmosphere, "pressure" },
    { temperature_dependent_atmosphere, "temperature" },
    { gas_constant_dependent_atmosphere, "gasConstant" },
    { specific_heat_ratio_dependent_atmosphere, "specificHeatRatio" },
    { molar_mass_dependent_atmosphere, "molarMass" }
};

//! `AtmosphereDependentVariables` not supported by `json_interface`.
static std::vector< AtmosphereDependentVariables > unsupportedAtmosphereDependentVariables = { };

//! Convert `AtmosphereDependentVariables` to `json`.
inline void to_json( nlohmann::json& jsonObject, const AtmosphereDependentVariables& atmosphereDependentVariable )
{
    jsonObject = json_interface::stringFromEnum( atmosphereDependentVariable, atmosphereDependentVariables );
}

//! Convert `json` to `AtmosphereDependentVariables`.
inline void from_json( const nlohmann::json& jsonObject, AtmosphereDependentVariables& atmosphereDependentVariable )
{
    atmosphereDependentVariable = json_interface::enumFromString( jsonObject, atmosphereDependentVariables );
}

}


namespace simulation_setup
{

// AtmosphereTypes

//! Map of `AtmosphereTypes` string representations.
static std::map< AtmosphereTypes, std::string > atmosphereTypes =
{
    { exponential_atmosphere, "exponential" },
    { tabulated_atmosphere, "tabulated" },
    { nrlmsise00, "nrlmsise00" }
};

//! `AtmosphereTypes` not supported by `json_interface`.
static std::vector< AtmosphereTypes > unsupportedAtmosphereTypes = { custom_constant_temperature_atmosphere };

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


// AtmosphereSettings

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< AtmosphereSettings >& atmosphereSettings );

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< AtmosphereSettings >& atmosphereSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ATMOSPHERE_H
