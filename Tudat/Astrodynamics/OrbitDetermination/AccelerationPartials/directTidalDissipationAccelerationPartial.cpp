/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/directTidalDissipationAccelerationPartial.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{


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
                radialComponentMultiplier * ( Eigen::Matrix3d::Identity( ) - 8.0 * relativePositionUnitVector * relativePositionUnitVector.transpose( ) ) +
                timeLag * ( -20.0 * positionVelocityInnerProduct *
                            relativePositionUnitVector * relativePositionUnitVector.transpose( ) / distanceSquared + 2.0 / distanceSquared * (
                                Eigen::Matrix3d::Identity( ) * positionVelocityInnerProduct +
                                relativePosition * relativeVelocity.transpose( ) ) -
                            8.0 / distance * ( relativePosition.cross( planetAngularVelocityVector) + relativeVelocity ) *
                            relativePositionUnitVector.transpose( ) -
                            linear_algebra::getCrossProductMatrix( planetAngularVelocityVector ) ) );
}

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

Eigen::Matrix3d computeDirectTidalAccelerationDueToTideOnSatelliteWrtVelocity(
        const Eigen::Vector6d relativeStateOfBodyExertingTide, const double currentTidalAccelerationMultiplier,
        const double timeLag )
{
    Eigen::Vector3d relativePosition = relativeStateOfBodyExertingTide.segment( 0, 3 );
    Eigen::Vector3d relativePositionUnitVector = relativePosition.normalized( );
    return currentTidalAccelerationMultiplier * timeLag * 3.5 *
                 ( 2.0 * relativePositionUnitVector * relativePositionUnitVector.transpose( ) );
}

//! Constructor
DirectTidalDissipationAccelerationPartial::DirectTidalDissipationAccelerationPartial(
        const boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > tidalAcceleration,
        const std::string acceleratedBody,
        const std::string acceleratingBody ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::direct_tidal_dissipation_acceleration ),
    tidalAcceleration_( tidalAcceleration )
{

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
            partialFunctionPair = std::make_pair(
                        boost::bind( &DirectTidalDissipationAccelerationPartial::wrtTidalTimeLag, this, _1 ), 1 );
        }
        else if( ( parameter->getParameterName( ).second.first == acceleratedBody_ ) &&
                !tidalAcceleration_->getModelTideOnPlanet( ) )
        {
            partialFunctionPair = std::make_pair(
                        boost::bind( &DirectTidalDissipationAccelerationPartial::wrtTidalTimeLag, this, _1 ), 1 );
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
            partialFunction = boost::bind( &DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfPlanet,
                                           this, _1 );
            numberOfColumns = 1;

        }

        // Check if parameter body is accelerated body, and if the mutual acceleration is used.
        if( parameterId.second.first == acceleratedBody_ )
        {
            if( !tidalAcceleration_->getModelTideOnPlanet( ) )
            {
                partialFunction = boost::bind( &DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfSatellite,
                                               this, _1 );
                numberOfColumns = 1;
            }
        }
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

//! Function to calculate central gravity partial w.r.t. central body gravitational parameter
void DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfPlanet( Eigen::MatrixXd& gravitationalParameterPartial )
{
    std::cout<<"Exerting: "<<tidalAcceleration_->getMassFunctionOfBodyExertingTide( )( )<<std::endl;
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = tidalAcceleration_->getAcceleration( ) *
            2.0 / tidalAcceleration_->getMassFunctionOfBodyExertingTide( )( );
}

void DirectTidalDissipationAccelerationPartial::wrtGravitationalParameterOfSatellite( Eigen::MatrixXd& gravitationalParameterPartial )
{
    std::cout<<"Undergoing: "<<tidalAcceleration_->getMassFunctionOfBodyUndergoingTide( )( )<<std::endl;
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = -tidalAcceleration_->getAcceleration( ) /
            tidalAcceleration_->getMassFunctionOfBodyUndergoingTide( )( );
}

void DirectTidalDissipationAccelerationPartial::wrtTidalTimeLag( Eigen::MatrixXd& gravitationalParameterPartial )
{
    Eigen::Vector3d tidalAccelerationWithoutScaling = tidalAcceleration_->getAcceleration( ) /
            tidalAcceleration_->getCurrentTidalAccelerationMultiplier( );
    if( tidalAcceleration_->getIncludeDirectRadialComponent( ) )
    {
        tidalAccelerationWithoutScaling -= tidalAcceleration_->getCurrentRelativeState( ).segment( 0, 3 );
    }
    gravitationalParameterPartial.block( 0, 0, 3, 1 ) = tidalAccelerationWithoutScaling / tidalAcceleration_->getTimeLag( );
}

void DirectTidalDissipationAccelerationPartial::wrtTidalLoveNumber( Eigen::MatrixXd& gravitationalParameterPartial )
{

}




}

}
