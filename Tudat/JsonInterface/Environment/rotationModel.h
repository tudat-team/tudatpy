/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_ROTATIONMODEL_H
#define TUDAT_JSONINTERFACE_ROTATIONMODEL_H

#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"


namespace tudat
{

namespace simulation_setup
{

//! Map of `RotationModelType`s string representations.
static std::map< RotationModelType, std::string > rotationModelTypes =
{
    { simple_rotation_model, "simple" },
    { spice_rotation_model, "spice" },
    { gcrs_to_itrs_rotation_model, "gcrsToItrs" }
};

//! `RotationModelType`s not supported by `json_interface`.
static std::vector< RotationModelType > unsupportedRotationModelTypes = { };

//! Convert `RotationModelType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const RotationModelType& rotationModelType )
{
    jsonObject = json_interface::stringFromEnum( rotationModelType, rotationModelTypes );
}

//! Convert `json` to `RotationModelType`.
inline void from_json( const nlohmann::json& jsonObject, RotationModelType& rotationModelType )
{
    rotationModelType = json_interface::enumFromString( jsonObject, rotationModelTypes );
}

//! Create a `json` object from a shared pointer to a `RotationModelSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< RotationModelSettings >& rotationModelSettings );

//! Create a shared pointer to a `RotationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< RotationModelSettings >& rotationModelSettings );

} // namespace simulation_setup

namespace basic_astrodynamics
{

//! Map of `RotationModelType`s string representations.
static std::map< IAUConventions, std::string > precessionNutationConventions =
{
    { iau_2000_a, "IAU2000a" },
    { iau_2000_b, "IAU2000b" },
    { iau_2006, "IAU2006" }
};

//! Convert `IAUConventions` to `json`.
inline void to_json( nlohmann::json& jsonObject, const IAUConventions& conventionType )
{
    jsonObject = json_interface::stringFromEnum( conventionType, precessionNutationConventions );
}

//! Convert `json` to `IAUConventions`.
inline void from_json( const nlohmann::json& jsonObject, IAUConventions& conventionType )
{
    conventionType = json_interface::enumFromString( jsonObject, precessionNutationConventions );
}

}

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ROTATIONMODEL_H
