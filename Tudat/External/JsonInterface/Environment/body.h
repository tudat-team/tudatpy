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
void to_json( json& jsonObject, const boost::shared_ptr< BodySettings >& bodySettings );

/*
//! Create a shared pointer to a `BodySettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< BodySettings >& bodySettings );
*/

} // namespace simulation_setup


namespace json_interface
{

//! Create a simulation_setup::BodySettings object with the settings from \p jsonObject.
/*!
 * @copybrief createBodySettings
 * \param jsonObject `json` object containing the settings for one body.
 * \return Body settings object.
 */
boost::shared_ptr< simulation_setup::BodySettings > createBodySettings( const json& jsonObject );

//! Update \p bodySettings with the settings from \p jsonObject.
/*!
 * @copybrief updateBodySettings
 * Does not change the values already defined in \p bodySettings that are not specified in \p jsonObject.
 * \param bodySettings Body settings object to be updated.
 * \param jsonObject `json` object containing only the settings for one body.
 */
void updateBodySettings( boost::shared_ptr< simulation_setup::BodySettings >& bodySettings, const json& jsonObject );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_BODY_H
