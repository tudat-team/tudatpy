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

#include <Tudat/SimulationSetup/EnvironmentSetup/createBodies.h>

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

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

//! Create a `BodySettings` object with the settings from a `json` object.
/*!
 * Create a `BodySettings` object with the settings from a `json` object.
 * \param settings `json` object containing the settings for one body.
 * \param keyTree Key tree at which the object containing the body settings can be accessed.
 * Empty if `settings` contains ONLY the body settings.
 * \return Body settings object.
 */
boost::shared_ptr< simulation_setup::BodySettings > createBodySettings(
        const json& settings, const KeyTree& keyTree = { } );

//! Update a `BodySettings` object with the settings from a `json` object.
/*!
 * Update a `BodySettings` object with the settings from a `json` object.
 * Does not change the values already defined in `bodySettings` that are not provided in `settings`.
 * \param bodySettings Body settings object to be updated.
 * \param settings `json` object containing only the settings for one body.
 * \param keyTree Key tree at which the object containing the body settings can be accessed.
 * Empty if `settings` contains ONLY the body settings.
 */
void updateBodySettings( boost::shared_ptr< simulation_setup::BodySettings >& bodySettings,
                         const json& settings, const KeyTree& keyTree = { } );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_BODY_H
