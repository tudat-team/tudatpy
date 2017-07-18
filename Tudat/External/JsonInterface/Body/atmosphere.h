/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_ATMOSPHERE_H
#define TUDAT_JSONINTERFACE_ATMOSPHERE_H

#include <Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.h>

#include "Tudat/External/JsonInterface/jsonInterface.h"

namespace tudat
{

namespace simulation_setup
{

//! Map of `AtmosphereTypes` supported by `json_interface`.
static std::map< std::string, AtmosphereTypes > atmosphereTypes =
{
    { "exponential", exponential_atmosphere },
    { "tabulated",   tabulated_atmosphere },
    { "nrlmsise00",  nrlmsise00 }
};
void to_json( json& jsonObject, const AtmosphereTypes& atmosphereType );
void from_json( const json& jsonObject, AtmosphereTypes& atmosphereType );

//! Create a `json` object from a shared pointer to a `AtmosphereSettings` object.
//! Called automatically by `nlohmann::json` when using a constructor such as `json( atmosphereSettings )`.
void to_json( json& jsonObject, const boost::shared_ptr< AtmosphereSettings >& atmosphereSettings );

} // namespace simulation_setup


namespace json_interface
{

//! Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
/*!
 * Create a shared pointer to a `AtmosphereSettings` object from a `json` object.
 * \param settings `json` object containing only the settings for one atmosphere model.
 * \return Shared pointer to a `AtmosphereSettings` object.
 */
boost::shared_ptr< simulation_setup::AtmosphereSettings > createAtmosphereSettings( const json &settings );

} // namespace json_interface


} // namespace tudat

#endif // TUDAT_JSONINTERFACE_ATMOSPHERE_H
