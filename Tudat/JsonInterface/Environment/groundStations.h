/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_GROUNDSTATION_H
#define TUDAT_JSONINTERFACE_GROUNDSTATION_H

#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace ground_stations
{

//! Map of `PositionElementTypes` string representations.
static std::map< coordinate_conversions::PositionElementTypes, std::string > positionElementTypes =
{
    { coordinate_conversions::cartesian_position, "cartesianPosition" },
    { coordinate_conversions::spherical_position, "sphericalPosition" },
    { coordinate_conversions::geodetic_position, "geodeticPosition" }
};

//! `PositionElementTypes` not supported by `json_interface`.
static std::vector< coordinate_conversions::PositionElementTypes > unsupportedPositionElementTypes =
{
};


//! Convert `PositionElementTypes` to `json`.
inline void to_json( nlohmann::json& jsonObject, const coordinate_conversions::PositionElementTypes& positionElementType )
{
    jsonObject = json_interface::stringFromEnum( positionElementType, positionElementTypes );
}

//! Convert `json` to `PositionElementTypes`.
inline void from_json( const nlohmann::json& jsonObject, coordinate_conversions::PositionElementTypes& positionElementType )
{
    positionElementType = json_interface::enumFromString( jsonObject, positionElementTypes );
}

}

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GroundStationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< GroundStationSettings >& groundStationSettings );

//! Create a shared pointer to a `GroundStationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< GroundStationSettings >& groundStationSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_GROUNDSTATION_H
