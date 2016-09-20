#ifndef CREATEGROUNDSTATIONS_H
#define CREATEGROUNDSTATIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/GroundStations/groundStation.h"
#include "Tudat/InputOutput/basicInputOutput.h"

namespace tudat
{

namespace simulation_setup
{


void createGroundStations( const NamedBodyMap& bodyMap,
                           const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition );


void createGroundStations( const NamedBodyMap& bodyMap,
                           std::vector< std::pair< std::string, std::string > > groundStations );
}

}

#endif // CREATEGROUNDSTATIONS_H
