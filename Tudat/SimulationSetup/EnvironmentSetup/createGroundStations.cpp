#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"

namespace tudat
{

namespace simulation_setup
{

void createGroundStations( const NamedBodyMap& bodyMap,
                           const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition )
{
    using namespace tudat::ephemerides;
    using namespace tudat::coordinate_conversions;

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > currentBodyShapeModel;

    for( std::map< std::pair< std::string, std::string >, Eigen::Vector3d >::const_iterator stationIterator = groundStationsWithPosition.begin( );
         stationIterator != groundStationsWithPosition.end( ); stationIterator++ )
    {
        if( bodyMap.count( stationIterator->first.first ) > 0 )
        {
            currentBodyShapeModel =  bodyMap.at( stationIterator->first.first )->getShapeModel( );

            if( currentBodyShapeModel == NULL )
            {
                std::cerr<<"Error when making ground station "<<stationIterator->first.second<<" of "<<stationIterator->first.first<<", body "<<
                           " has no shape model."<<std::endl;
            }
            else
            {

                boost::shared_ptr< ground_stations::NominalGroundStationState > stationState;

                    stationState = boost::make_shared< ground_stations::NominalGroundStationState >(
                                stationIterator->second, currentBodyShapeModel  );

                bodyMap.at( stationIterator->first.first )->addGroundStation(
                            stationIterator->first.second, boost::make_shared< ground_stations::GroundStation >(
                                stationState, stationIterator->first.second ) );

            }
        }
    }
}

}

}
