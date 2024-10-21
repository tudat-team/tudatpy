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

std::shared_ptr< ground_stations::StationMotionModel > createGroundStationMotionModels(
        const std::shared_ptr< Body > body,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings =
        std::vector< std::shared_ptr< GroundStationMotionSettings > >( ),
        const simulation_setup::SystemOfBodies& bodies = simulation_setup::SystemOfBodies( ) )
{

    std::vector< std::shared_ptr< ground_stations::StationMotionModel > > stationMotionModelList;

    for( unsigned int i = 0; i < stationMotionSettings.size( ); i++ )
    {
        std::shared_ptr< ground_stations::StationMotionModel > currentStationMotionModel;
        switch( stationMotionSettings.at( i )->getModelType( ) )
        {
        case body_deformation_station_motion:
        {
            std::shared_ptr< BodyDeformationStationMotionSettings > bodyDeformationModelSettings =
                    std::dynamic_pointer_cast< BodyDeformationStationMotionSettings >( stationMotionSettings.at( i ) );
            if( bodyDeformationModelSettings == nullptr )
            {
                throw std::runtime_error( "Error when making body deformation station motion model, settings type is incompatible" );
            }
            currentStationMotionModel =
                std::make_shared< ground_stations::BodyDeformationStationMotionModel >(
                    std::bind( &Body::getBodyDeformationModelsReference, body ), bodyDeformationModelSettings->throwExceptionWhenNotAvailable_ );
            break;
        }
        case linear_station_motion:
        {
            std::shared_ptr< LinearGroundStationMotionSettings > linearStationMotionSettings =
                    std::dynamic_pointer_cast< LinearGroundStationMotionSettings >( stationMotionSettings.at( i ) );
            if( linearStationMotionSettings == nullptr )
            {
                throw std::runtime_error( "Error when making linear ground station motion model, settings type is incompatible" );
            }
            currentStationMotionModel = std::make_shared< ground_stations::LinearStationMotionModel >(
                        linearStationMotionSettings->linearVelocity_, linearStationMotionSettings->referenceEpoch_ );
            break;
        }
        case piecewise_constant_station_motion:
        {
            std::shared_ptr< PiecewiseConstantGroundStationMotionSettings > piecewiseConstantStationMotionSettings =
                    std::dynamic_pointer_cast< PiecewiseConstantGroundStationMotionSettings >( stationMotionSettings.at( i ) );
            if( piecewiseConstantStationMotionSettings == nullptr )
            {
                throw std::runtime_error( "Error when making piecewise constant ground station motion model, settings type is incompatible" );
            }
            currentStationMotionModel = std::make_shared< ground_stations::PiecewiseConstantStationMotionModel >(
                        piecewiseConstantStationMotionSettings->displacementList_ );
            break;
        }
        case bodycentric_to_barycentric_station_position_motion:
        {
            std::shared_ptr< BodyCentricToBarycentricGroundStationMotionSettings > relativisticStationMotionSettings =
                std::dynamic_pointer_cast< BodyCentricToBarycentricGroundStationMotionSettings >( stationMotionSettings.at( i ) );
            if( relativisticStationMotionSettings == nullptr )
            {
                throw std::runtime_error( "Error when making relativistic (body-centered to barycentered) ground station motion model, settings type is incompatible" );
            }


            std::function< Eigen::Vector6d( const double ) > bodyBarycentricStateFunction;
            if( bodies.getFrameOrigin( ) == "SSB" )
            {
                bodyBarycentricStateFunction = std::bind( &Body::getStateInBaseFrameFromEphemeris<double, double>, body, std::placeholders::_1 );
            }
            else
            {
                throw std::runtime_error( "Error when getting body state function for (body-centered to barycentered) ground station motion model, only SSB global frame origin is currently supported" );
            }

            std::function< Eigen::Quaterniond( const double ) > inertialToBodyFixedRotationFunction = std::bind(
                &ephemerides::RotationalEphemeris::getRotationToTargetFrame, body->getRotationalEphemeris( ), std::placeholders::_1 );

            std::function< Eigen::Vector3d( const double ) > centralBodyBarycentricPositionFunction = nullptr;
            std::function< double( ) > centralBodyGravitationalParameterFunction = nullptr;
            if( bodies.count( relativisticStationMotionSettings->centralBodyName_ ) == 0  && relativisticStationMotionSettings->useGeneralRelativisticCorrection_ )
            {
                throw std::runtime_error( "Error when making bodycentric to barycentric station position correction, body " + relativisticStationMotionSettings->centralBodyName_  + " not found" );
            }
            else if( relativisticStationMotionSettings->useGeneralRelativisticCorrection_ )
            {
                centralBodyBarycentricPositionFunction = std::bind( &Body::getPositionInBaseFrameFromEphemeris<double, double>, body, std::placeholders::_1 );
                if( bodies.at( relativisticStationMotionSettings->centralBodyName_ )->getGravityFieldModel( ) == nullptr )
                {
                    throw std::runtime_error( "Error when making bodycentric to barycentric station position correction, body " + relativisticStationMotionSettings->centralBodyName_  + " has no grvaity field" );

                }
                centralBodyGravitationalParameterFunction = std::bind(
                    &gravitation::GravityFieldModel::getGravitationalParameter, bodies.at( relativisticStationMotionSettings->centralBodyName_ )->getGravityFieldModel( ) );
            }


            currentStationMotionModel = std::make_shared< ground_stations::BodyCentricToBarycentricRelativisticStationMotion >(
            bodyBarycentricStateFunction, centralBodyBarycentricPositionFunction, inertialToBodyFixedRotationFunction, centralBodyGravitationalParameterFunction,
                relativisticStationMotionSettings->useGeneralRelativisticCorrection_ );
            break;
        }
        case custom_station_motion:
        {
            std::shared_ptr< CustomGroundStationMotionSettings > customStationMotionSettings =
                    std::dynamic_pointer_cast< CustomGroundStationMotionSettings >( stationMotionSettings.at( i ) );
            if( customStationMotionSettings == nullptr )
            {
                throw std::runtime_error( "Error when making custom ground station motion model, settings type is incompatible" );
            }
            currentStationMotionModel = std::make_shared< ground_stations::CustomStationMotionModel >(
                        customStationMotionSettings->customDisplacementModel_ );
            break;
        }
        default:
            throw std::runtime_error( "Error when making ground station motion model, settings type not recognized" );

        }
        stationMotionModelList.push_back( currentStationMotionModel );
    }
    return std::make_shared< ground_stations::CombinedStationMotionModel >( stationMotionModelList );
}

std::shared_ptr< ground_stations::GroundStationState > createGroundStationState(
        const std::shared_ptr< Body > body,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings =
        std::vector< std::shared_ptr< GroundStationMotionSettings > >( ),
        const simulation_setup::SystemOfBodies& bodies = simulation_setup::SystemOfBodies( ) )
{
    std::shared_ptr< ground_stations::StationMotionModel > stationMotionModel =
            createGroundStationMotionModels( body,  stationMotionSettings, bodies );
    return std::make_shared< ground_stations::GroundStationState >(
                                 groundStationPosition, positionElementType, body->getShapeModel( ), stationMotionModel );
}

//! Function to create a ground station and add it to a Body object
void createGroundStation(
        const std::shared_ptr< Body > body,
        const std::string groundStationName,
        const Eigen::Vector3d groundStationPosition,
        const coordinate_conversions::PositionElementTypes positionElementType,
        const std::vector< std::shared_ptr< GroundStationMotionSettings > > stationMotionSettings )
{
    createGroundStation( body, groundStationName, createGroundStationState(
                             body, groundStationPosition, positionElementType, stationMotionSettings ) );

}

//! Function to create a set of ground stations and add them to the corresponding Body objects
void    createGroundStations(
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
        const std::shared_ptr< GroundStationSettings > groundStationSettings )
{

    if( body->getGroundStationMap( ).count( groundStationSettings->getStationName( ) ) != 0 )
    {
        throw std::runtime_error(
                    "Error when creating ground station " + groundStationSettings->getStationName( ) +
                    " on body " + body->getBodyName( ) + ", station already exists." );
    }
    else
    {
        createGroundStation( body, groundStationSettings->getStationName( ),
                             groundStationSettings->getGroundStationPosition( ),
                             groundStationSettings->getPositionElementType( ),
                             groundStationSettings->getStationMotionSettings( ) );
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
//
//std::vector< double >  getTargetElevationAngles(
//        const std::shared_ptr< Body > observingBody,
//        const std::shared_ptr< Body > targetBody,
//        const std::string groundStationName,
//        const std::vector< double > times )
//{
//    if( observingBody->getGroundStationMap( ).count( groundStationName ) == 0 )
//    {
//        throw std::runtime_error( "Error when computing elevating angle, station " + groundStationName +
//                                  " not found on body " + observingBody->getBodyName( ) );
//    }
//
//
//    std::function< Eigen::Vector6d( const double& ) > groundStationStateFunction =
//            getLinkEndCompleteEphemerisFunction(
//                observingBody, std::make_pair( observingBody->getBodyName( ), groundStationName ) );
//
//    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAnglesCalculator =
//            observingBody->getGroundStationMap( ).at( groundStationName )->getPointingAnglesCalculator( );
//    Eigen::Vector3d relativePosition;
//    std::vector< double > elevationAngles;
//    for( unsigned int i = 0; i < times.size( ); i++ )
//    {
//        targetBody->getStateInBaseFrameFromEphemeris( times.at( i ) );
//        groundStationStateFunction( times.at( i ) );
//        relativePosition = ( targetBody->getStateInBaseFrameFromEphemeris( times.at( i ) ) -
//                groundStationStateFunction( times.at( i ) ) ).segment( 0, 3 );
//        elevationAngles.push_back( pointingAnglesCalculator->calculateElevationAngleFromInertialVector( relativePosition, times.at( i ) ) );
//    }
//    return elevationAngles;
//}

}

}
