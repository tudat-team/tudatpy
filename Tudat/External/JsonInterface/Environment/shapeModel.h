/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/SimulationSetup/EnvironmentSetup/createBodyShapeModel.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `BodyShapeTypes` supported by `json_interface`.
static std::map< std::string, BodyShapeTypes > bodyShapeTypes =
{
    { "spherical",      spherical },
    { "sphericalSpice", spherical_spice },
    { "oblateSpheroid", oblate_spheroid }
};

//! Convert `BodyShapeTypes` to `json`.
void to_json( json& jsonObject, const BodyShapeTypes& bodyShapeTypes );

//! Convert `json` to `BodyShapeTypes`.
void from_json( const json& jsonObject, BodyShapeTypes& bodyShapeTypes );

//! Create a `json` object from a shared pointer to a `BodyShapeSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( bodyShapeSettings )`.
void to_json( json& jsonObject, const boost::shared_ptr< BodyShapeSettings >& bodyShapeSettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `BodyShapeSettings` object from a `json` object.
/*!
 * Create a shared pointer to a `BodyShapeSettings` object from a `json` object.
 * \param settings `json` object containing the settings for one shape model.
 * \param keyTree Key tree at which the object containing the shape model settings can be accessed.
 * Empty if `settings` contains ONLY the shape model settings.
 * \return Shared pointer to a `BodyShapeSettings` object.
 */
boost::shared_ptr< simulation_setup::BodyShapeSettings > createShapeModelSettings(
        const json& settings, const KeyTree& keyTree = { } );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SHAPEMODEL_H
