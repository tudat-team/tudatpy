/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_SHAPEMODEL_H
#define TUDAT_JSONINTERFACE_SHAPEMODEL_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `BodyShapeTypes` string representations.
static std::map< BodyShapeTypes, std::string > bodyShapeTypes =
{
    { spherical, "spherical" },
    { spherical_spice, "sphericalSpice" },
    { oblate_spheroid, "oblateSpheroid" }
};

//! `BodyShapeTypes` not supported by `json_interface`.
static std::vector< BodyShapeTypes > unsupportedBodyShapeTypes = { };

//! Convert `BodyShapeTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const BodyShapeTypes& bodyShapeType )
{
    jsonObject = json_interface::stringFromEnum( bodyShapeType, bodyShapeTypes );
}

//! Convert `json` to `BodyShapeTypes`.
inline void from_json( const nlohmann::json& jsonObject, BodyShapeTypes& bodyShapeType )
{
    bodyShapeType = json_interface::enumFromString( jsonObject, bodyShapeTypes );
}

//! Create a `json` object from a shared pointer to a `BodyShapeSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< BodyShapeSettings >& bodyShapeSettings );

//! Create a shared pointer to a `BodyShapeSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< BodyShapeSettings >& bodyShapeSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SHAPEMODEL_H
