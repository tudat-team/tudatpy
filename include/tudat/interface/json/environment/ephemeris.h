/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_EPHEMERIS_H
#define TUDAT_JSONINTERFACE_EPHEMERIS_H

#include "tudat/simulation/environment_setup/createEphemeris.h"
#include "tudat/interface/json/support/valueAccess.h"
#include "tudat/interface/json/support/valueConversions.h"

namespace tudat
{

namespace ephemerides
{

namespace simulation_setup
{

//! Map of `EphemerisType`s string representations.
static std::map< EphemerisType, std::string > ephemerisTypes =
{
    { approximate_planet_positions, "approximatePlanetPositions" },
    { direct_spice_ephemeris, "directSpice" },
    { tabulated_ephemeris, "tabulated" },
    { interpolated_spice, "interpolatedSpice" },
    { constant_ephemeris, "constant" },
    { kepler_ephemeris, "kepler" },
    { custom_ephemeris, "custom" }
};

//! `EphemerisType` not supported by `json_interface`.
static std::vector< EphemerisType > unsupportedEphemerisTypes =
{
    custom_ephemeris
};

//! Convert `EphemerisType` to `json`.
inline void to_json( nlohmann::json& jsonObject, const EphemerisType& ephemerisType )
{
    jsonObject = json_interface::stringFromEnum( ephemerisType, ephemerisTypes );
}

//! Convert `json` to `EphemerisType`.
inline void from_json( const nlohmann::json& jsonObject, EphemerisType& ephemerisType )
{
    ephemerisType = json_interface::enumFromString( jsonObject, ephemerisTypes );
}

//! Create a `json` object from a shared pointer to a `EphemerisSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< EphemerisSettings >& ephemerisSettings );

//! Create a shared pointer to a `EphemerisSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< EphemerisSettings >& ephemerisSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_EPHEMERIS_H
