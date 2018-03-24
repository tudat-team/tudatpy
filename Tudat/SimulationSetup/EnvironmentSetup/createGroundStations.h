/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEGROUNDSTATIONS_H
#define TUDAT_CREATEGROUNDSTATIONS_H

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"

namespace tudat
{

namespace simulation_setup
{

class GroundStationSettings
{
public:
    GroundStationSettings(
            const std::string& stationName,
            const Eigen::Vector3d& groundStationPosition,
            const coordinate_conversions::PositionElementTypes positionElementType =
            coordinate_conversions::cartesian_position ):
        stationName_( stationName ),
        groundStationPosition_( groundStationPosition ),
        positionElementType_( positionElementType ){ }

    std::string getStationName( )
    {
        return stationName_;
    }

    Eigen::Vector3d getGroundStationPosition( )
    {
        return groundStationPosition_;
    }

    coordinate_conversions::PositionElementTypes getPositionElementType( )
    {
        return positionElementType_;
    }

protected:

    std::string stationName_;

    Eigen::Vector3d groundStationPosition_;

    coordinate_conversions::PositionElementTypes positionElementType_;
};

//! Function to create a ground station from pre-defined station state object, and add it to a Body object
/*!
 * Function to create a ground station from pre-defined station state object, and add it to a Body object
 * \param body Body object in which the newly created ground station is to be added.
 * \param groundStationName Name of ground station that is to be created
 * \param groundStationState Object defining the state of the ground-station in a body-fixed frame
 */
void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const boost::shared_ptr< ground_stations::GroundStationState > groundStationState );

//! Function to create a ground station and add it to a Body object
/*!
 * Function to create a ground station and add it to a Body object
 * \param body Body object in which the newly created ground station is to be added.
 * \param groundStationName Name of ground station that is to be created
 * \param groundStationPosition Position of ground station in body-fixed frame
 * \param positionElementType Element type (e.g. Cartesian, spherical, etc.) of groundStationPosition.
 */
void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType =
        coordinate_conversions::cartesian_position );

//! Function to create a set of ground stations and add them to the corresponding Body objects
/*!
 * Function to create a set of ground stations and add them to the corresponding Body objects
 * \param bodyMap List of body objects to which the ground stations are to be added
 * \param groundStationsWithPosition List of ground station positions, key is first: associated body; second: ground station
 * name
 * \param positionElementType Element type (e.g. Cartesian, spherical, etc.) of Vector3d in groundStationsWithPosition.
 */
void createGroundStations(
        const NamedBodyMap& bodyMap,
        const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
        const coordinate_conversions::PositionElementTypes positionElementType =
        coordinate_conversions::cartesian_position );

void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string& bodyName,
        const boost::shared_ptr< GroundStationSettings > groundStationSettings );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEGROUNDSTATIONS_H
