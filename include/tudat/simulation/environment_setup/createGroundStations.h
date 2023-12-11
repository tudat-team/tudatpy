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
#include "tudat/astro/observation_models/linkTypeDefs.h"


namespace tudat
{

namespace simulation_setup
{

enum StationMotionModelTypes
{
    linear_station_motion,
    piecewise_constant_station_motion,
    custom_station_motion,
    body_deformation_station_motion
};

class GroundStationMotionSettings
{
public:
    GroundStationMotionSettings(
            const StationMotionModelTypes& modelType ):modelType_( modelType ){ }

    virtual ~GroundStationMotionSettings( ){ }

    StationMotionModelTypes getModelType( )
    {
        return modelType_;
    }

private:
    StationMotionModelTypes modelType_;
};

class LinearGroundStationMotionSettings: public GroundStationMotionSettings
{
public:
    LinearGroundStationMotionSettings(
            const Eigen::Vector3d& linearVelocity,
            const double referenceEpoch = 0.0 ):
        GroundStationMotionSettings( linear_station_motion ),
    linearVelocity_( linearVelocity ),
    referenceEpoch_( referenceEpoch ){ }

    virtual ~LinearGroundStationMotionSettings( ){ }

    Eigen::Vector3d linearVelocity_;

    double referenceEpoch_;
};

class PiecewiseConstantGroundStationMotionSettings: public GroundStationMotionSettings
{
public:
    PiecewiseConstantGroundStationMotionSettings(
            const std::map< double, Eigen::Vector3d >& displacementList ):
        GroundStationMotionSettings( piecewise_constant_station_motion ),
    displacementList_( displacementList ){ }

    virtual ~PiecewiseConstantGroundStationMotionSettings( ){ }

    std::map< double, Eigen::Vector3d > displacementList_;
};

class CustomGroundStationMotionSettings: public GroundStationMotionSettings
{
public:
//    CustomGroundStationMotionSettings(
//            const std::function< Eigen::Vector6d( const double ) > customDisplacementModel ):
//        GroundStationMotionSettings( custom_station_motion ),
//    customDisplacementModel_( customDisplacementModel ){ }

    CustomGroundStationMotionSettings(
            const std::function< Eigen::Vector3d( const double ) > customDisplacementModel ):
        GroundStationMotionSettings( custom_station_motion ),
        customDisplacementModel_( [=](const double time ){
        return ( Eigen::Vector6d( ) << customDisplacementModel( time ), Eigen::Vector3d::Zero( ) ).finished( ); } ){ }

    virtual ~CustomGroundStationMotionSettings( ){ }

    const std::function< Eigen::Vector6d( const double ) > customDisplacementModel_;
};


inline std::shared_ptr< GroundStationMotionSettings > linearGroundStationMotionSettings(
        const Eigen::Vector3d& linearVelocity,
        const double referenceEpoch =  0.0 )
{
    return std::make_shared< LinearGroundStationMotionSettings >(
                linearVelocity, referenceEpoch );
}

inline std::shared_ptr< GroundStationMotionSettings > piecewiseConstantGroundStationMotionSettings(
        const std::map< double, Eigen::Vector3d >& displacementList)
{
    return std::make_shared< PiecewiseConstantGroundStationMotionSettings >(
                displacementList );
}

inline std::shared_ptr< GroundStationMotionSettings > customGroundStationMotionSettings(
        const std::function< Eigen::Vector3d( const double ) > customDisplacementModel )
{
    return std::make_shared< CustomGroundStationMotionSettings >(
                customDisplacementModel );
}


class GroundStationSettings
{
public:
    GroundStationSettings(
            const std::string& stationName,
            const Eigen::Vector3d& groundStationPosition,
            const coordinate_conversions::PositionElementTypes positionElementType =
            coordinate_conversions::cartesian_position,
            const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings =
            std::vector< std::shared_ptr< GroundStationMotionSettings > >( ) ):
        stationName_( stationName ),
        groundStationPosition_( groundStationPosition ),
        positionElementType_( positionElementType ),
    stationMotionSettings_( stationMotionSettings ){ }

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

    std::vector< std::shared_ptr< GroundStationMotionSettings > > getStationMotionSettings( )
    {
        return stationMotionSettings_;
    }

    void addStationMotionSettings( const std::shared_ptr< GroundStationMotionSettings > stationMotionSetting )
    {
        stationMotionSettings_.push_back( stationMotionSetting );
    }

protected:

    std::string stationName_;

    Eigen::Vector3d groundStationPosition_;

    coordinate_conversions::PositionElementTypes positionElementType_;

    std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings_;
};

inline std::shared_ptr< GroundStationSettings > groundStationSettings(
        const std::string& stationName,
        const Eigen::Vector3d& groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType =
        coordinate_conversions::cartesian_position,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings =
        std::vector< std::shared_ptr< GroundStationMotionSettings > >( ) )
{
    return std::make_shared< GroundStationSettings >(
                stationName, groundStationPosition, positionElementType, stationMotionSettings );
}

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

std::shared_ptr< ground_stations::GroundStationState > createGroundStationState(
        const std::shared_ptr< Body > body,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType );

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
        coordinate_conversions::cartesian_position,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings =
        std::vector< std::shared_ptr< GroundStationMotionSettings > >( ) );

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
    stationEphemerisVector[ 2 ] = std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< StateScalarType, TimeType >,
        bodyWithReferencePoint, std::placeholders::_1 );
    stationEphemerisVector[ 0 ] = referencePointStateFunction;

    std::map< int, std::function< StateType( const TimeType, const StateType& ) > > stationRotationVector;
    stationRotationVector[ 1 ] =  std::bind( &ephemerides::transformStateToInertialOrientation< StateScalarType, TimeType >,
        std::placeholders::_2, std::placeholders::_1, bodyRotationModel );

    // Create and return ephemeris
    return std::make_shared< ephemerides::CompositeEphemeris< TimeType, StateScalarType > >(
                stationEphemerisVector, stationRotationVector, "SSB", "ECLIPJ2000" );
}

template< typename TimeType = double, typename StateScalarType = double >
std::shared_ptr< ephemerides::Ephemeris > createReferencePointEphemeris(
    const std::shared_ptr< simulation_setup::Body > bodyWithLinkEnd,
    const std::string& referencePointName )
{
    std::shared_ptr< ephemerides::Ephemeris > stationEphemeris;
    if( referencePointName != "" )
    {
        std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType& ) > referencePointStateFunction;

        bool isPointGroundStation = simulation_setup::isReferencePointGroundStation( bodyWithLinkEnd, referencePointName );

        if( isPointGroundStation )
        {
            if ( bodyWithLinkEnd->getGroundStationMap( ).count( referencePointName ) == 0 )
            {
                std::string errorMessage = "Error when making ephemeris for station " + bodyWithLinkEnd->getBodyName() + ", " +
                    referencePointName + ", station not found.";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                referencePointStateFunction = std::bind( &ground_stations::GroundStation::getStateInPlanetFixedFrame<StateScalarType, TimeType>,
                                                         bodyWithLinkEnd->getGroundStation( referencePointName ), std::placeholders::_1 );
            }
        }
        else
        {
            if ( bodyWithLinkEnd->getVehicleSystems( ) == nullptr )
            {
                std::string errorMessage = "Error when making ephemeris for reference point " + bodyWithLinkEnd->getBodyName() + ", " +
                                           referencePointName + ", no vehicle systems found ";
                throw std::runtime_error( errorMessage );
            }
            else if( bodyWithLinkEnd->getVehicleSystems( )->doesReferencePointExist( referencePointName ) == false )
            {
                std::string errorMessage = "Error when making ephemeris for reference point " + bodyWithLinkEnd->getBodyName() + ", " +
                                           referencePointName + ", reference point not found in vehicle systems. ";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                referencePointStateFunction = std::bind(
                    &system_models::VehicleSystems::getReferencePointStateInBodyFixedFrame< StateScalarType, TimeType >,
                    bodyWithLinkEnd->getVehicleSystems( ), referencePointName, std::placeholders::_1 );
            }
        }

        // Retrieve function to calculate state of transmitter S/C
        stationEphemeris = createReferencePointEphemeris<TimeType, StateScalarType>(
            bodyWithLinkEnd, referencePointStateFunction );
    }
    else
    {
        throw std::runtime_error( "Error when making ground station ephemeris, no station ID specified" );
    }
    return stationEphemeris;
}


template< typename StateScalarType = double >
Eigen::Matrix< StateScalarType, 3, 1 > getGroundStationPositionDuringPropagation(
        const std::shared_ptr< simulation_setup::Body > bodyWithLinkEnd,
        const std::string& stationName )
{
    if( bodyWithLinkEnd->getGroundStationMap( ).count( stationName ) == 0 )
    {
        std::string errorMessage = "Error when getting grounf station position for " + bodyWithLinkEnd->getBodyName( ) + ", " +
                stationName + ", station not found.";
        throw std::runtime_error( errorMessage );
    }

    return bodyWithLinkEnd->getPosition( ) + bodyWithLinkEnd->getCurrentRotationToGlobalFrame( ) *
            bodyWithLinkEnd->getGroundStation( stationName )->getNominalStationState( )->getCartesianPositionInTime(
                bodyWithLinkEnd->getDoubleTimeOfCurrentState( ) );

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
        const observation_models::LinkEndId& linkEndId )
{
    typedef Eigen::Matrix< StateScalarType, 6, 1 > StateType;

    std::function< StateType( const TimeType& ) > linkEndCompleteEphemerisFunction;

    if( linkEndId.bodyName_ != bodyWithLinkEnd->getBodyName( ) )
    {
        throw std::runtime_error( "Error when making ground station ephemeris function, input body names are inconsistent" );
    }

    // Checking transmitter if a reference point is to be used
    if( linkEndId.stationName_ != "" )
    {

        // Retrieve function to calculate state of transmitter S/C
        linkEndCompleteEphemerisFunction =
                std::bind( &ephemerides::Ephemeris::getTemplatedStateFromEphemeris< StateScalarType,TimeType >,
                             createReferencePointEphemeris< TimeType, StateScalarType >(
                                 bodyWithLinkEnd, linkEndId.stationName_ ), std::placeholders::_1 );

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
        const observation_models::LinkEndId linkEndId, const simulation_setup::SystemOfBodies& bodies )
{
    if( bodies.count( linkEndId.bodyName_ ) == 0  )
    {
        std::string errorMessage;
        if ( linkEndId.stationName_ != "" )
        {
            errorMessage = "Error when making ephemeris function for body " + linkEndId.bodyName_ + ", station " +
                    linkEndId.stationName_ + ": body not found.";
        }
        else
        {
            errorMessage = "Error when making ephemeris function for body " + linkEndId.bodyName_ + ": body not found.";
        }


        throw std::runtime_error( errorMessage );
    }
    return getLinkEndCompleteEphemerisFunction< TimeType, StateScalarType >( bodies.at( linkEndId.bodyName_ ), linkEndId );
}

std::vector< double >  getTargetElevationAngles(
        const std::shared_ptr< Body > observingBody,
        const std::shared_ptr< Body > targetBody,
        const std::string groundStationName,
        const std::vector< double > times );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEGROUNDSTATIONS_H
