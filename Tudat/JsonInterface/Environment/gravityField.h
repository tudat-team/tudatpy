/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_GRAVITYFIELD_H
#define TUDAT_JSONINTERFACE_GRAVITYFIELD_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

// GravityFieldType

//! Map of `GravityFieldType`s string representations.
static std::map< GravityFieldType, std::string > gravityFieldTypes =
{
    { central, "pointMass" },
    { central_spice, "pointMassSpice" },
    { spherical_harmonic, "sphericalHarmonic" }
};

//! `GravityFieldType` not supported by `json_interface`.
static std::vector< GravityFieldType > unsupportedGravityFieldTypes = { };

//! Convert `GravityFieldType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const GravityFieldType& gravityFieldType )
{
    jsonObject = json_interface::stringFromEnum( gravityFieldType, gravityFieldTypes );
}

//! Convert `json` to `GravityFieldType`.
inline void from_json( const nlohmann::json& jsonObject, GravityFieldType& gravityFieldType )
{
    gravityFieldType = json_interface::enumFromString( jsonObject, gravityFieldTypes );
}


// SphericalHarmonicsModel

//! Map of `SphericalHarmonicsModel`s string representations.
static std::map< SphericalHarmonicsModel, std::string > sphericalHarmonicsModels =
{
    { customModel, "custom" },
    { egm96, "egm96" },
    { ggm02c, "ggm02c" },
    { ggm02s, "ggm02s" },
    { glgm3150, "glgm3150" },
    { lpe200, "lpe200" },
    { jgmro120d, "jgmro120d" }
};

//! `SphericalHarmonicsModel` not supported by `json_interface`.
static std::vector< SphericalHarmonicsModel > unsupportedSphericalHarmonicsModels = { };

//! Convert `SphericalHarmonicsModel` to `json`.
inline void to_json( nlohmann::json& jsonObject, const SphericalHarmonicsModel& sphericalHarmonicsModel )
{
    jsonObject = json_interface::stringFromEnum( sphericalHarmonicsModel, sphericalHarmonicsModels );
}

//! Convert `json` to `SphericalHarmonicsModel`.
inline void from_json( const nlohmann::json& jsonObject, SphericalHarmonicsModel& sphericalHarmonicsModel )
{
    sphericalHarmonicsModel =
            json_interface::enumFromString( jsonObject, sphericalHarmonicsModels );
}


// GravityFieldSettings

//! Create a `json` object from a shared pointer to a `GravityFieldSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< GravityFieldSettings >& gravityFieldSettings );

//! Create a shared pointer to a `GravityFieldSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< GravityFieldSettings >& gravityFieldSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_GRAVITYFIELD_H
