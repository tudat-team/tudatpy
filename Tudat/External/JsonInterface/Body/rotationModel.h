/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h>

#include "Tudat/External/JsonInterface/jsonInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `RotationModelType`s supported by `json_interface`.
static std::map< std::string, RotationModelType > rotationModelTypes =
{
    { "simple", simple_rotation_model },
    { "spice",  spice_rotation_model }
};

//! Convert `RotationModelType` to `json`.
void to_json( json& jsonObject, const RotationModelType& rotationModelType );

//! Convert `json` to `RotationModelType`.
void from_json( const json& jsonObject, RotationModelType& rotationModelType );

//! Create a `json` object from a shared pointer to a `RotationModelSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( rotationModelSettings )`.
void to_json( json& jsonObject, const boost::shared_ptr< RotationModelSettings >& rotationModelSettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `RotationModelSettings` object from a `json` object.
/*!
 * Create a shared pointer to a `RotationModelSettings` object from a `json` object.
 * \param settings `json` object containing the settings for one rotational model.
 * \param keyTree Key tree at which the object containing the rotational model settings can be accessed.
 * Empty if `settings` contains ONLY the rotational model settings.
 * \return Shared pointer to a `RotationModelSettings` object.
 */
boost::shared_ptr< simulation_setup::RotationModelSettings > createRotationModelSettings(
        const json& settings, const KeyTree& keyTree = { } );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ROTATIONMODEL_H
