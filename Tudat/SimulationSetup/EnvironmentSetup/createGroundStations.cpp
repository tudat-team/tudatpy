#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalBodyShapeModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"

namespace tudat
{

namespace simulation_setup
{

void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const boost::shared_ptr< ground_stations::GroundStationState > groundStationState )
{

    body->addGroundStation( groundStationName, boost::make_shared< ground_stations::GroundStation >(
                                groundStationState, groundStationName ) );
}

void createGroundStation(
        const boost::shared_ptr< Body >& body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType )
{

    boost::shared_ptr< basic_astrodynamics::BodyShapeModel > currentBodyShapeModel =  body->getShapeModel( );

    boost::shared_ptr< ground_stations::GroundStationState > stationState = boost::make_shared< ground_stations::GroundStationState >(
                groundStationPosition, positionElementType, currentBodyShapeModel );

    createGroundStation( body, groundStationName, stationState );

}

void createGroundStations( const NamedBodyMap& bodyMap,
                           const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
                           const coordinate_conversions::PositionElementTypes positionElementType )
{
    using namespace tudat::ephemerides;
    using namespace tudat::coordinate_conversions;

    for( std::map< std::pair< std::string, std::string >, Eigen::Vector3d >::const_iterator stationIterator = groundStationsWithPosition.begin( );
         stationIterator != groundStationsWithPosition.end( ); stationIterator++ )
    {
        if( bodyMap.count( stationIterator->first.first ) > 0 )
        {
            createGroundStation( bodyMap.at( stationIterator->first.first ), stationIterator->first.second,
                                 stationIterator->second, positionElementType );
        }
    }
}

}

}
