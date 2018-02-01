/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/directTidalDissipationAccelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/directTidalTimeLag.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial derivative of direct tidal acceleration due to tide on planet w.r.t. position of satellite
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnPlanetWrtPosition(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const Eigen::Vector3d planetAngularVelocityVector,
        const double currentTidalAccelerationMultiplier, const double timeLag, const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativePositionUnitVector = relativePosition.normalized( );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double positionVelocityInnerProduct = relativePosition.dot( relativeVelocity );
    double distance = relativePosition.norm( );
    double distanceSquared = distance * distance;
    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return currentTidalAccelerationMultiplier * (
                radialComponentMultiplier * ( Eigen::Matrix3d::Identity( )
                                              - 8.0 * relativePositionUnitVector * relativePositionUnitVector.transpose( ) ) +
                timeLag * ( -20.0 * positionVelocityInnerProduct *
                            relativePositionUnitVector * relativePositionUnitVector.transpose( ) /
                            distanceSquared + 2.0 / distanceSquared * (
                                Eigen::Matrix3d::Identity( ) * positionVelocityInnerProduct +
                                relativePosition * relativeVelocity.transpose( ) ) -
                            8.0 / distance * ( relativePosition.cross( planetAngularVelocityVector) + relativeVelocity ) *
                            relativePositionUnitVector.transpose( ) -
                            linear_algebra::getCrossProductMatrix( planetAngularVelocityVector ) ) );
}

//! Function to compute partial derivative of direct tidal acceleration due to tide on planet w.r.t. velocity of satellite/
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnPlanetWrtVelocity(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const double currentTidalAccelerationMultiplier,
        const double timeLag )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativePositionUnitVector = relativePosition.normalized( );

    return currentTidalAccelerationMultiplier * timeLag *
            ( 2.0 * relativePositionUnitVector * relativePositionUnitVector.transpose( ) +
              Eigen::Matrix3d::Identity( ) );
}

//! Function to compute partial derivative of direct tidal acceleration due to tide on satellite w.r.t. position of satellite
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnSatelliteWrtPosition(
        const Eigen::Vector6d relativeStateOfBodyExertingTide,
        const double currentTidalAccelerationMultiplier, const double timeLag, const bool includeDirectRadialComponent )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativePositionUnitVector = relativePosition.normalized( );
    Eigen::Vector3d relativeVelocity = relativeStateOfBodyExertingTide.segment( 3, 3 );

    double positionVelocityInnerProduct = relativePosition.dot( relativeVelocity );
    double distance = relativePosition.norm( );
    double distanceSquared = distance * distance;
    double radialComponentMultiplier = ( includeDirectRadialComponent == true ) ? 1.0 : 0.0;

    return currentTidalAccelerationMultiplier * (
                2.0 * radialComponentMultiplier * (
                    Eigen::Matrix3d::Identity( ) - 8.0 * relativePositionUnitVector * relativePositionUnitVector.transpose( ) ) +
                timeLag * 3.5 * ( -20.0 * positionVelocityInnerProduct *
                                  relativePositionUnitVector * relativePositionUnitVector.transpose( ) / distanceSquared +
                                  2.0 / distanceSquared * ( Eigen::Matrix3d::Identity( ) * positionVelocityInnerProduct +
                                                            relativePosition * relativeVelocity.transpose( ) ) ) );
}

//! Function to compute partial derivative of direct tidal acceleration due to tide on satellite w.r.t. velocity of satellite
Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnSatelliteWrtVelocity(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const double currentTidalAccelerationMultiplier,
        const double timeLag )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativePositionUnitVector = relativePosition.normalized( );
    return currentTidalAccelerationMultiplier * timeLag * 3.5 *
            ( 2.0 * relativePositionUnitVector * relativePositionUnitVector.transpose( ) );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
DirectTidalDissipationAccelerationPartial::getParameterPartialFunction(
        boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )

{
    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies.
    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        // If parameter is gravitational parameter, check and create dependency function .
        partialFunctionPair = this->getGravitationalParameterPartialFunction( parameter->getParameterName( ) );
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::direct_dissipation_tidal_time_lag )
    {
        if( ( parameter->getParameterName( ).second.first == acceleratingBody_ ) &&
                tidalAcceleration_->getModelTideOnPlanet( ) )
        {
            boost::shared_ptr< estimatable_parameters::DirectTidalTimeLag > timeLagParameter =
                    boost::dynamic_pointer_cast< estimatable_parameters::DirectTidalTimeLag >( parameter );
            if( timeLagParameter == NULL )
            {
                throw std::runtime_error( "Error when getting partial of DirectTidalDissipationAcceleration w.r.t. DirectTidalTimeLag, models are inconsistent" );
            }
            else
            {
                std::vector< std::string > bodiesCausingDeformation = timeLagParameter->getBodiesCausingDeformation( );
                if( bodiesCausingDeformation.size( ) == 0 || (
                            std::find( bodiesCausingDeformation.begin( ), bodiesCausingDeformation.end( ), acceleratedBody_ ) !=
                            bodiesCausingDeformation.end( ) ) )
                {
                    partialFunctionPair = std::make_pair(
                                boost::bind( &DirectTidalDissipationAccelerationPartial::wrtTidalTimeLag, this, _1 ), 1 );
                }
            }

        }
        else if( ( parameter->getParameterName( ).second.first == acceleratedBody_ ) &&
                 !tidalAcceleration_->getModelTideOnPlanet( ) )
        {
            boost::shared_ptr< estimatable_parameters::DirectTidalTimeLag > timeLagParameter =
                    boost::dynamic_pointer_cast< estimatable_parameters::DirectTidalTimeLag >( parameter );
            if( timeLagParameter == NULL )
            {
                throw std::runtime_error( "Error when getting partial of DirectTidalDissipationAcceleration w.r.t. DirectTidalTimeLag, models are inconsistent" );
            }
            else
            {
                std::vector< std::string > bodiesCausingDeformation = timeLagParameter->getBodiesCausingDeformation( );
                if( bodiesCausingDeformation.size( ) == 0 || (
                            std::find( bodiesCausingDeformation.begin( ), bodiesCausingDeformation.end( ), acceleratingBody_ ) !=
                            bodiesCausingDeformation.end( ) ) )
                {
                    partialFunctionPair = std::make_pair(
                                boost::bind( &DirectTidalDissipationAccelerationPartial::wrtTidalTimeLag, this, _1 ), 1 );
                }
            }


        }
    }
    else
    {
        partialFunctionPair = std::make_pair( boost::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }


    return partialFunctionPair;
}

void DirectTidalDissipationAccelerationPartial::update( const double currentTime )
{
    tidalAcceleration_->updateMembers( currentTime );

    if( !( currentTime_ == currentTime ) )
    {
        currentRelativeBodyState_ = tidalAcceleration_->getCurrentRelativeState( );

        if( tidalAcceleration_->getModelTideOnPlanet( ) )
        {
            currentPartialWrtPosition_ = computeDirectTidalAccelerationDueToTideOnPlanetWrtPosition(
                        currentRelativeBodyState_,
                        tidalAcceleration_->getCurrentAngularVelocityVectorOfBodyUndergoingTide( ),
                        tidalAcceleration_->getCurrentTidalAccelerationMultiplier( ),
                        tidalAcceleration_->getTimeLag( ), tidalAcceleration_->getIncludeDirectRadialComponent( ) );
            currentPartialWrtVelocity_ = computeDirectTidalAccelerationDueToTideOnPlanetWrtVelocity(
                        currentRelativeBodyState_, tidalAcceleration_->getCurrentTidalAccelerationMultiplier( ),
                        tidalAcceleration_->getTimeLag( ) );
        }
        else
        {
            currentPartialWrtPosition_ = computeDirectTidalAccelerationDueToTideOnSatelliteWrtPosition(
                        currentRelativeBodyState_, tidalAcceleration_->getCurrentTidalAccelerationMultiplier( ),
                        tidalAcceleration_->getTimeLag( ), tidalAcceleration_->getIncludeDirectRadialComponent( ) );
            currentPartialWrtVelocity_ = computeDirectTidalAccelerationDueToTideOnSatelliteWrtVelocity(
                        currentRelativeBodyState_, tidalAcceleration_->getCurrentTidalAccelerationMultiplier( ),
                        tidalAcceleration_->getTimeLag( ) );
        }

        currentTime_ = currentTime;
    }
}

//! Function to create a function returning the current partial w.r.t. a gravitational parameter.
std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >
DirectTidalDissipationAccelerationPartial::getGravitationalParameterPartialFunction(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterId )
{
    boost::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    // Check if parameter is gravitational parameter.
    if( parameterId.first ==  estimatable_parameters::gravitational_parameter )
    {
        // Check if parameter body is central body.
        if( parameterId.second.first == acceleratingBody_ )
        {
            if( !tidalAcceleration_->getModelTideOnPlanet( ) )
            {
                partialFunction = boost::bind( &DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfPlanet,
                                               this, _1 );
                numberOfColumns = 1;
            }

        }

        // Check if parameter body is accelerated body, and if the mutual acceleration is used.
        if( parameterId.second.first == acceleratedBody_ )
        {

            partialFunction = boost::bind( &DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfSatellite,
                                           this, _1 );
            numberOfColumns = 1;
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfPlanet( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = tidalAcceleration_->getAcceleration( ) *
            2.0 / tidalAcceleration_->getGravitationalParameterFunctionOfBodyExertingTide( )( );
}

//! Function to compute derivative w.r.t. gravitational parameter of satellite
void DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfSatellite( Eigen::MatrixXd& gravitationalParameterPartial )
{
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = -tidalAcceleration_->getAcceleration( ) /
            tidalAcceleration_->getGravitationalParameterFunctionOfBodyUndergoingTide( )( );

    if( tidalAcceleration_->getModelTideOnPlanet( ) )
    {
        gravitationalParameterPartial *= -1.0;
    }
}

//! Function to compute derivative w.r.t. tidal time lag parameter.
void DirectTidalDissipationAccelerationPartial::wrtTidalTimeLag( Eigen::MatrixXd& gravitationalParameterPartial )
{
    Eigen::Vector3d tidalAccelerationWithoutScaling = tidalAcceleration_->getAcceleration( ) /
            tidalAcceleration_->getCurrentTidalAccelerationMultiplier( );
    if( tidalAcceleration_->getIncludeDirectRadialComponent( ) )
    {
        tidalAccelerationWithoutScaling -= ( 1.0 + static_cast< double >( !tidalAcceleration_->getModelTideOnPlanet( ) ) ) *
                tidalAcceleration_->getCurrentRelativeState( ).segment( 0, 3 );
    }
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) =  tidalAcceleration_->getCurrentTidalAccelerationMultiplier( ) *
            tidalAccelerationWithoutScaling / tidalAcceleration_->getTimeLag( );
}


}

}
