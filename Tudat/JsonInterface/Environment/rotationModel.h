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

#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `RotationModelType`s string representations.
static std::map< RotationModelType, std::string > rotationModelTypes =
{
    { simple_rotation_model, "simple" },
    { spice_rotation_model, "spice" }
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
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< RotationModelSettings >& rotationModelSettings );

//! Create a shared pointer to a `RotationModelSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< RotationModelSettings >& rotationModelSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ROTATIONMODEL_H
