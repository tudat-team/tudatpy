/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/JsonInterface/Environment/groundStations.h"

namespace tudat
{

namespace simulation_setup
{

//! Create a `json` object from a shared pointer to a `GroundStationSettings` object.
void to_json( nlohmann::json& jsonObject, const std::shared_ptr< GroundStationSettings >& groundStationSettings )
{
    if ( ! groundStationSettings )
    {
        return;
    }
    using namespace coordinate_conversions;
    using namespace json_interface;
    using K = Keys::Body::GroundStation;

    jsonObject[ K::positionElementType ] = groundStationSettings->getPositionElementType( );
    jsonObject[ K::stationName ] = groundStationSettings->getStationName( );
    jsonObject[ K::stationPosition ] = groundStationSettings->getGroundStationPosition( );
}

//! Create a shared pointer to a `GroundStationSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, std::shared_ptr< GroundStationSettings >& groundStationSettings )
{
    using namespace coordinate_conversions;
    using namespace json_interface;
    using K = Keys::Body::GroundStation;

    PositionElementTypes positionElementType =
            enumFromString( getValue< std::string >( jsonObject, K::positionElementType ), ground_stations::positionElementTypes );
    Eigen::Vector3d stationPosition =
            getValue< Eigen::Vector3d >( jsonObject, K::stationPosition );
    std::string stationName =
            getValue< std::string >( jsonObject, K::stationName );

    groundStationSettings = std::make_shared< GroundStationSettings >(
                stationName, stationPosition, positionElementType );

}

} // namespace simulation_setup

} // namespace tudat
