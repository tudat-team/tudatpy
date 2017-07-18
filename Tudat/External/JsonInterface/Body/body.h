/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_BODY_H
#define TUDAT_JSONINTERFACE_BODY_H

#include <Tudat/SimulationSetup/EnvironmentSetup/body.h>
#include <Tudat/SimulationSetup/EnvironmentSetup/createBodies.h>

#include "Tudat/External/JsonInterface/jsonInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `BodySettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( bodySettings )`.
void to_json( json& jsonObject, const boost::shared_ptr< BodySettings >& bodySettings );

} // namespace simulation_setup


namespace json_interface
{

//! Update a `BodySettings` object with the settings from a `json` object.
/*!
 * Update a `BodySettings` object with the settings from a `json` object.
 * Does not change the values already defined that are not provided in `settings`.
 * \param bodyMap Map with the bodies created so far.
 * If the body to be updated interfaces with celestial bodies, those must already be defined in `bodyMap`.
 * If the body to be updated does not exist, it will created with empty constructor.
 * \param bodyName The name of the body to be updated.
 * \param settings `json` object containing only the settings for one body.
 */
void updateBodySettings( std::map< std::string, boost::shared_ptr< simulation_setup::BodySettings > >& bodySettingsMap,
                         const std::string& bodyName, const json &settings );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_BODY_H
