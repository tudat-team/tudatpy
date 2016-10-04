#ifndef TUDAT_CREATEGROUNDSTATIONS_H
#define TUDAT_CREATEGROUNDSTATIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace simulation_setup
{

void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const boost::shared_ptr< ground_stations::GroundStationState > groundStationState );

void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType = coordinate_conversions::cartesian_position );

void createGroundStations(
        const NamedBodyMap& bodyMap,
        const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
        const coordinate_conversions::PositionElementTypes positionElementType = coordinate_conversions::cartesian_position );


}

}

#endif // TUDAT_CREATEGROUNDSTATIONS_H
