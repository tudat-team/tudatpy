/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/relativisticAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
void computePartialOfSchwarzschildAccelerationCorrectionWrtPosition(
        const Eigen::Vector6d& relativeState, Eigen::Vector3d& currentAcceleration, Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter, const double ppnParameterGamma, const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = 2.0 * ( ppnParameterGamma + ppnParameterBeta ) * gravitationalParameter / distance * (
                Eigen::Matrix3d::Identity( ) - position * position.transpose( ) / ( distance * distance ) );
    partialMatrix -= ppnParameterGamma * velocity.dot( velocity ) * Eigen::Matrix3d::Identity( );
    partialMatrix += 2.0 * ( 1.0 + ppnParameterGamma ) * velocity * velocity.transpose( );
    partialMatrix *= gravitationalParameter *  physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance );
    partialMatrix -= 3.0 * currentAcceleration * position.transpose( ) / ( distance * distance );
}

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
void computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity(
        const Eigen::Vector6d& relativeState, Eigen::Matrix3d& partialMatrix,
        const double gravitationalParameter, const double ppnParameterGamma )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = gravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( - 2.0 * ppnParameterGamma * position * velocity.transpose( ) +
              2.0 * ( 1.0 + ppnParameterGamma ) * (
                  position.dot( velocity ) * Eigen::Matrix3d::Identity( ) + velocity * position.transpose( ) ) );

}

//! Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
void computePartialOfSchwarzschildAccelerationCorrectionWrtGravitationalParameter(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix, const double ppnParameterGamma, const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( -ppnParameterGamma * ( velocity.dot( velocity ) ) * position +
              2.0 * ( 1.0 + ppnParameterGamma ) *
              ( position.dot( velocity ) ) * velocity
              + 4.0 * gravitationalParameter * ( ppnParameterGamma + ppnParameterBeta ) *
              position / distance );
}

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterGamma(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) * (
                ( 2.0 * gravitationalParameter / distance - velocity.dot( velocity ) ) * position + 2.0 * position.dot( velocity ) * velocity );
}

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterBeta(
        const Eigen::Vector6d& relativeState,
        const double gravitationalParameter,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) * (
                ( 2.0 * gravitationalParameter / distance ) * position );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
RelativisticAccelerationPartial::getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Create partial function if parameter is central body gravitational parameter
    if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::gravitational_parameter:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtGravitationalParameterOfCentralBody, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }
    // Create partial function if parameter is PPN parameter beta or gamma
    else if( parameter->getParameterName( ).second.first == "global_metric"  )
    {
        switch( parameter->getParameterName( ).first )
        {
        case estimatable_parameters::ppn_parameter_gamma:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtPpnParameterGamma, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        case estimatable_parameters::ppn_parameter_beta:
            partialFunction = std::bind( &RelativisticAccelerationPartial::wrtPpnParameterBeta, this, std::placeholders::_1 );
            numberOfRows = 1;
            break;
        default:
            break;
        }
    }
    return std::make_pair( partialFunction, numberOfRows );
}

//! Function for updating partial w.r.t. the bodies' states
void RelativisticAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRelativeState_ = ( acceleratedBodyState_( ) - centralBodyState_( ) );
        currentAcceleration_ = currentAccelerationFunction_( );

        computePartialOfSchwarzschildAccelerationCorrectionWrtPosition(
                    currentRelativeState_, currentAcceleration_, currentPartialWrtPosition_,
                    centralBodyGravitationalParameterFunction_( ),
                    ppnGammaParameterFunction_( ), ppnBetaParameterFunction_( ) );
        computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity(
                    currentRelativeState_, currentPartialWrtVelocity_,
                    centralBodyGravitationalParameterFunction_( ),
                    ppnGammaParameterFunction_( ) );
        currentTime_ = currentTime;
    }
}

}

}


