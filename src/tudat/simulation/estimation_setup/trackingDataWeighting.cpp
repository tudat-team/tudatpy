#include "tudat/simulation/estimation_setup/trackingDataWeighting.h"

#include "tudat/astro/ground_stations/groundStation.h"

namespace tudat
{

namespace observation_models
{

Eigen::Vector3d getGroundStationPositionForTrackingDataLinkEnd( const simulation_setup::SystemOfBodies& bodies, const LinkEndId& linkEndId )
{
    if( !bodies.doesBodyExist( linkEndId.bodyName_ ) )
    {
        throw std::runtime_error( "Error when retrieving ground station position for TrackingData weighing, body '" + linkEndId.bodyName_ +
                                  "' does not exist." );
    }

    if( linkEndId.bodyName_ != "Earth" )
    {
        return Eigen::Vector3d::Zero( );
    }

    const std::shared_ptr< simulation_setup::Body > body = bodies.getBody( linkEndId.bodyName_ );
    const std::map< std::string, std::shared_ptr< ground_stations::GroundStation > > groundStations = body->getGroundStationMap( );
    const auto stationIterator = groundStations.find( linkEndId.getReferencePointName( ) );
    if( stationIterator == groundStations.end( ) || stationIterator->second == nullptr )
    {
        throw std::runtime_error( "Error when retrieving ground station position for TrackingData weighing, station '" +
                                  linkEndId.getReferencePointName( ) + "' does not exist on Earth." );
    }

    return stationIterator->second->getNominalStationState( )->getNominalCartesianPosition( );
}

}  // namespace observation_models

}  // namespace tudat
