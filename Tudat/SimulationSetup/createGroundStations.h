#ifndef CREATEGROUNDSTATIONS_H
#define CREATEGROUNDSTATIONS_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/SimulationSetup/body.h"
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

void setSingleBodyGroundStationPositionVariationFunctions( boost::shared_ptr< Body > body,
                                                           gravitation::BodyDeformationTypes variationType,
                                                           const double initialTime,
                                                           const double finalTime,
                                                           const double timeStep = 600.0 );

void setSingleBodyGroundStationPositionVariationFunctions( boost::shared_ptr< Body > body,
                                                           std::vector< gravitation::BodyDeformationTypes > variationTypes,
                                                           const double initialTime,
                                                           const double finalTime,
                                                           const double timeStep = 600.0 );

boost::function< Eigen::Vector3d( const double ) > getGroundStationPositionVariationFunction(
        boost::shared_ptr< ground_stations::GroundStation > groundStation,
        boost::shared_ptr< Body > bodyWithStation,
        gravitation::BodyDeformationTypes variationType,
        const double initialTime,
        const double finalTime,
        const double timeStep );


void setGroundStationPositionVariationFunction( boost::shared_ptr< ground_stations::GroundStation > groundStation,
                                                boost::shared_ptr< Body > bodyWithStation,
                                                gravitation::BodyDeformationTypes variationType,
                                                const double initialTime,
                                                const double finalTime,
                                                const double timeStep = 600.0 );

void setGroundStationPositionVariationFunction( boost::shared_ptr< ground_stations::GroundStation > groundStation,
                                                boost::shared_ptr< Body > bodyWithStation,
                                                std::vector< gravitation::BodyDeformationTypes > variationTypes,
                                                const double initialTime,
                                                const double finalTime,
                                                const double timeStep = 600.0 );

void updateGroundStationStates( std::vector< std::pair< std::string, std::string > >,
                                NamedBodyMap& bodyMap );
}

}

#endif // CREATEGROUNDSTATIONS_H
