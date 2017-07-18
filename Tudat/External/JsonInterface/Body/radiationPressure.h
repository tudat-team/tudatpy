/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_RADIATIONPRESSURE_H
#define TUDAT_JSONINTERFACE_RADIATIONPRESSURE_H

#include <Tudat/SimulationSetup/EnvironmentSetup/createRadiationPressureInterface.h>

#include "Tudat/External/JsonInterface/jsonInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `RadiationPressureType`s supported by `json_interface`.
static std::map< std::string, RadiationPressureType > radiationPressureTypes =
{
    { "cannonBall", cannon_ball }
};
void to_json( json& jsonObject, const RadiationPressureType& radiationPressureType );
void from_json( const json& jsonObject, RadiationPressureType& radiationPressureType );

//! Create a `json` object from a shared pointer to a `RadiationPressureInterfaceSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as
//! `json( radiationPressureInterfaceSettings )`.
void to_json( json& jsonObject,
              const boost::shared_ptr< RadiationPressureInterfaceSettings >& radiationPressureInterfaceSettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `RadiationPressureInterfaceSettings` object from a `json` object.
/*!
 * Create a shared pointer to a `RadiationPressureInterfaceSettings` object from a `json` object.
 * \param settings `json` object containing only the settings for one radiation pressure interface.
 * \param sourceBodyName The name of the radiating body.
 * \param fallbackArea Fallback reference area to be used when no reference area is speciefied in `settings`.
 * \return Shared pointer to a `RadiationPressureInterfaceSettings` object.
 */
boost::shared_ptr< simulation_setup::RadiationPressureInterfaceSettings > createRadiationPressureInterfaceSettings(
    const json &settings, const std::string& sourceBodyName, const double& fallbackArea = TUDAT_NAN );

} // namespace json_interface


} // namespace tudat

#endif // TUDAT_JSONINTERFACE_RADIATIONPRESSURE_H
