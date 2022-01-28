/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/ground_stations/groundStation.h"

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
        const std::shared_ptr< Body > body,
        const std::string groundStationName,
        const std::shared_ptr< ground_stations::GroundStationState > groundStationState );

//! Function to create a ground station and add it to a Body object
/*!
 * Function to create a ground station and add it to a Body object
 * \param body Body object in which the newly created ground station is to be added.
 * \param groundStationName Name of ground station that is to be created
 * \param groundStationPosition Position of ground station in body-fixed frame
 * \param positionElementType Element type (e.g. Cartesian, spherical, etc.) of groundStationPosition.
 */
void createGroundStation(
        const std::shared_ptr< Body > body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType =
        coordinate_conversions::cartesian_position );

//! Function to create a set of ground stations and add them to the corresponding Body objects
/*!
 * Function to create a set of ground stations and add them to the corresponding Body objects
 * \param bodies List of body objects to which the ground stations are to be added
 * \param groundStationsWithPosition List of ground station positions, key is first: associated body; second: ground station
 * name
 * \param positionElementType Element type (e.g. Cartesian, spherical, etc.) of Vector3d in groundStationsWithPosition.
 */
void createGroundStations(
        const SystemOfBodies& bodies,
        const std::map< std::pair< std::string, std::string >, Eigen::Vector3d >& groundStationsWithPosition,
        const coordinate_conversions::PositionElementTypes positionElementType =
        coordinate_conversions::cartesian_position );

void createGroundStation(
        const std::shared_ptr< Body > body,
        const std::string& bodyName,
        const std::shared_ptr< GroundStationSettings > groundStationSettings );

std::vector< std::pair< std::string, std::string > > getGroundStationsLinkEndList(
        const std::shared_ptr< Body > body );


//! Function to create an ephemeris for a reference point on a body
/*!
 *  Function to create an ephemeris for a reference point on a body, taking into account the time-variable rotation
 *  of the body and its global ephemeris
 *  \param bodyWithReferencePoint Body on which reference point is located
 *  \param bodyRotationModel Rotation model that is to be used for going from body-fixed to inertial frame
 *  \param referencePointStateFunction Function returning the state of the reference point on the body (in a body-fixed
 *  frame).
 *  \return Reference point ephemeris in global coordinates.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< ephemerides::Ephemeris > createReferencePointEphemeris(
        const std::shared_ptr< simulation_setup::Body > bodyWithReferencePoint,
        const std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > referencePointStateFunction )
{
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    if( bodyWithReferencePoint == nullptr )
    {
        throw std::runtime_error( "Error when creating reference point ephemeris, body is not provided" );
    }

    std::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel =
            bodyWithReferencePoint->getRotationalEphemeris( );

    if( bodyRotationModel == nullptr )
    {
        throw std::runtime_error( "Error when creating reference point ephemeris, no body rotation model is provided" );
    }

    if( referencePointStateFunction == nullptr )
    {
        throw std::runtime_error( "Error when creating reference point ephemeris, no reference point state function is provided" );
    }

    // Create list of state/rotation functions that are to be used
    std::map< int, std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > > stationEphemerisVector;
    stationEphemerisVector[ 2 ] = std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris
                                               < StateScalarType, TimeType >, bodyWithReferencePoint, std::placeholders::_1 );
    stationEphemerisVector[ 0 ] = referencePointStateFunction;

    std::map< int, std::function< StateType( const TimeType, const StateType& ) > > stationRotationVector;
    stationRotationVector[ 1 ] =  std::bind( &ephemerides::transformStateToInertialOrientation
                                               < StateScalarType, TimeType >, std::placeholders::_2, std::placeholders::_1, bodyRotationModel );

    // Create and return ephemeris
    return std::make_shared< ephemerides::CompositeEphemeris< TimeType, StateScalarType > >(
                stationEphemerisVector, stationRotationVector, "SSB", "ECLIPJ2000" );
}

//! Function to retrieve a state function for a link end (either a body center of mass or ground station).
/*!
 *  Function to retrieve a state function for a link end (either a body center of mass or ground station).
 *  \param bodyWithLinkEnd Body on/in which link end is situated.
 *  \param linkEndId Id of link end for which state function is to be created. First: name of body, second: name of
 *  reference point (empty if center of mass is to be used
 *  \return Requested state function
 */
template< typename TimeType = double, typename StateScalarType = double >
std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > getLinkEndCompleteEphemerisFunction(
        const std::shared_ptr< simulation_setup::Body > bodyWithLinkEnd,
        const std::pair< std::string, std::string >& linkEndId )
{
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    std::function< StateType( const TimeType& ) > linkEndCompleteEphemerisFunction;

    // Checking transmitter if a reference point is to be used
    if( linkEndId.second != "" )
    {
        if( bodyWithLinkEnd->getGroundStationMap( ).count( linkEndId.second ) == 0 )
        {
            std::string errorMessage = "Error when making ephemeris function for " + linkEndId.first + ", " +
                    linkEndId.second + ", station not found.";
            throw std::runtime_error( errorMessage );
        }

        // Retrieve function to calculate state of transmitter S/C
        linkEndCompleteEphemerisFunction =
                std::bind( &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< StateScalarType,TimeType >,
                             createReferencePointEphemeris< TimeType, StateScalarType >(
                                 bodyWithLinkEnd,
                                 std::bind( &ground_stations::GroundStation::getStateInPlanetFixedFrame
                                              < StateScalarType, TimeType >,
                                              bodyWithLinkEnd->getGroundStation( linkEndId.second ), std::placeholders::_1 ) ), std::placeholders::_1 );

    }
    // Else, create state function for center of mass
    else
    {
        // Create function to calculate state of transmitting ground station.
        linkEndCompleteEphemerisFunction =
                std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
                                                        bodyWithLinkEnd, std::placeholders::_1 );
    }
    return linkEndCompleteEphemerisFunction;
}

//! Function to create a state function of a link end, expressed in base frame.
/*!
 *  Function to create a state function of a link end, expressed in base frame.
 *  \param linkEndId Name of the reference point for which state function is to be created.
 *  \param bodies List of body objects that comprises the environment
 *  \return Requested state function.
 */
template< typename TimeType = double, typename StateScalarType = double >
std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > getLinkEndCompleteEphemerisFunction(
        const std::pair< std::string, std::string > linkEndId, const simulation_setup::SystemOfBodies& bodies )
{
    if( bodies.count( linkEndId.first ) == 0  )
    {
        std::string errorMessage = "Error when making ephemeris function for " + linkEndId.first + ", " +
                linkEndId.second + ", body not found.";
        throw std::runtime_error( errorMessage );
    }
    return getLinkEndCompleteEphemerisFunction< TimeType, StateScalarType >( bodies.at( linkEndId.first ), linkEndId );
}

std::vector< double >  getTargetElevationAngles(
        const std::shared_ptr< Body > observingBody,
        const std::shared_ptr< Body > targetBody,
        const std::string groundStationName,
        const std::vector< double > times );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEGROUNDSTATIONS_H
