/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace simulation_setup
{

//! Function to create a ground station from pre-defined station state object, and add it to a Body object
void createGroundStation(
        const std::shared_ptr< Body > body,
        const std::string groundStationName,
        const std::shared_ptr< ground_stations::GroundStationState > groundStationState )
{
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            std::make_shared< ground_stations::PointingAnglesCalculator >(
                std::bind( &ephemerides::RotationalEphemeris::getRotationToTargetFrame, body->getRotationalEphemeris( ), std::placeholders::_1 ),
                std::bind( &ground_stations::GroundStationState::getRotationFromBodyFixedToTopocentricFrame, groundStationState, std::placeholders::_1 ) );
    body->addGroundStation( groundStationName, std::make_shared< ground_stations::GroundStation >(
                                groundStationState, pointingAnglesCalculator, groundStationName ) );
}

//! Function to create a ground station and add it to a Body object
void createGroundStation(
        const std::shared_ptr< Body > body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType )
{
    createGroundStation( body, groundStationName, std::make_shared< ground_stations::GroundStationState >(
                             groundStationPosition, positionElementType, body->getShapeModel( ) ) );

}

//! Function to create a set of ground stations and add them to the corresponding Body objects
void createGroundStations(
        const SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
        const coordinate_conversions::PositionElementTypes positionElementType )
{
    for( std::map< std::pair< std::string, std::string >, Eigen::Vector3d >::const_iterator
         stationIterator = groundStationsWithPosition.begin( );
         stationIterator != groundStationsWithPosition.end( ); stationIterator++ )
    {
        if( bodies.count( stationIterator->first.first ) > 0 )
        {
            createGroundStation( bodies.at( stationIterator->first.first ), stationIterator->first.second,
                                 stationIterator->second, positionElementType );
        }
    }
}

void createGroundStation(
        const std::shared_ptr< Body > body,
        const std::string& bodyName,
        const std::shared_ptr< GroundStationSettings > groundStationSettings )
{

    if( body->getGroundStationMap( ).count( groundStationSettings->getStationName( ) ) != 0 )
    {
        throw std::runtime_error(
                    "Error when creating ground station " + groundStationSettings->getStationName( ) +
                    " on body " + bodyName + ", station already exists." );
    }
    else
    {
        createGroundStation( body, groundStationSettings->getStationName( ),
                             groundStationSettings->getGroundStationPosition( ),
                             groundStationSettings->getPositionElementType( ) );
    }
}

std::vector< std::pair< std::string, std::string > > getGroundStationsLinkEndList(
        const std::shared_ptr< Body > body )
{
    std::vector< std::pair< std::string, std::string > > stationList;

    std::map<std::string, std::shared_ptr<ground_stations::GroundStation> > groundStationsMap = body->getGroundStationMap();
    for( auto stationIterator : groundStationsMap )
    {
        stationList.push_back( std::make_pair( body->getBodyName( ), stationIterator.first ) );
    }
    return stationList;
}

std::vector< double >  getTargetElevationAngles(
        const std::shared_ptr< Body > observingBody,
        const std::shared_ptr< Body > targetBody,
        const std::string groundStationName,
        const std::vector< double > times )
{
    if( observingBody->getGroundStationMap( ).count( groundStationName ) == 0 )
    {
        throw std::runtime_error( "Error when computing elevating angle, station " + groundStationName +
                                  " not found on body " + observingBody->getBodyName( ) );
    }


    std::function< Eigen::Vector6d( const double& ) > groundStationStateFunction =
            getLinkEndCompleteEphemerisFunction(
                observingBody, std::make_pair( observingBody->getBodyName( ), groundStationName ) );

    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
            observingBody->getGroundStationMap( ).at( groundStationName )->getPointingAnglesCalculator( );
    Eigen::Vector3d relativePosition;
    std::vector< double > elevationAngles;
    for( unsigned int i = 0; i < times.size( ); i++ )
    {
        targetBody->getStateInBaseFrameFromEphemeris( times.at( i ) );
        groundStationStateFunction( times.at( i ) );
        relativePosition = ( targetBody->getStateInBaseFrameFromEphemeris( times.at( i ) ) -
                groundStationStateFunction( times.at( i ) ) ).segment( 0, 3 );
        elevationAngles.push_back( pointingAnglesCalculator->calculateElevationAngle( relativePosition, times.at( i ) ) );
    }
    return elevationAngles;
}

}

}
