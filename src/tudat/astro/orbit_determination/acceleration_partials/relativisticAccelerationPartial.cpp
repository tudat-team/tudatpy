/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cmath>

#include "tudat/astro/orbit_determination/acceleration_partials/relativisticAccelerationPartial.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. position of body undergoing acceleration
void computePartialOfSchwarzschildAccelerationCorrectionWrtPosition( const Eigen::Vector6d& relativeState,
                                                                     Eigen::Vector3d& currentAcceleration,
                                                                     Eigen::Matrix3d& partialMatrix,
                                                                     const double gravitationalParameter,
                                                                     const double ppnParameterGamma,
                                                                     const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = 2.0 * ( ppnParameterGamma + ppnParameterBeta ) * gravitationalParameter / distance *
            ( Eigen::Matrix3d::Identity( ) - position * position.transpose( ) / ( distance * distance ) );
    partialMatrix -= ppnParameterGamma * velocity.dot( velocity ) * Eigen::Matrix3d::Identity( );
    partialMatrix += 2.0 * ( 1.0 + ppnParameterGamma ) * velocity * velocity.transpose( );
    partialMatrix *= gravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance );
    partialMatrix -= 3.0 * currentAcceleration * position.transpose( ) / ( distance * distance );
}

//! Function to compute partial of Schwarzschild acceleration correction w.r.t. velocity of body undergoing acceleration
void computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity( const Eigen::Vector6d& relativeState,
                                                                     Eigen::Matrix3d& partialMatrix,
                                                                     const double gravitationalParameter,
                                                                     const double ppnParameterGamma )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = gravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( -2.0 * ppnParameterGamma * position * velocity.transpose( ) +
              2.0 * ( 1.0 + ppnParameterGamma ) *
                      ( position.dot( velocity ) * Eigen::Matrix3d::Identity( ) + velocity * position.transpose( ) ) );
}

//! Function to compute partial derivative of Schwarzschild acceleration correction w.r.t. central body gravitational patameter.
void computePartialOfSchwarzschildAccelerationCorrectionWrtGravitationalParameter( const Eigen::Vector6d& relativeState,
                                                                                   const double gravitationalParameter,
                                                                                   Eigen::MatrixXd& partialMatrix,
                                                                                   const double ppnParameterGamma,
                                                                                   const double ppnParameterBeta )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / ( distance * distance * distance ) *
            ( -ppnParameterGamma * ( velocity.dot( velocity ) ) * position +
              2.0 * ( 1.0 + ppnParameterGamma ) * ( position.dot( velocity ) ) * velocity +
              4.0 * gravitationalParameter * ( ppnParameterGamma + ppnParameterBeta ) * position / distance );
}

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter gamma
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterGamma( const Eigen::Vector6d& relativeState,
                                                                              const double gravitationalParameter,
                                                                              Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    Eigen::Vector3d velocity = relativeState.segment( 3, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) *
            ( ( 2.0 * gravitationalParameter / distance - velocity.dot( velocity ) ) * position +
              2.0 * position.dot( velocity ) * velocity );
}

//! Function to compute the partial derivative of Schwarzschild acceleration correction w.r.t. PPN parameter beta
void computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterBeta( const Eigen::Vector6d& relativeState,
                                                                             const double gravitationalParameter,
                                                                             Eigen::MatrixXd& partialMatrix )
{
    Eigen::Vector3d position = relativeState.segment( 0, 3 );
    double distance = position.norm( );

    partialMatrix = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT * gravitationalParameter / ( distance * distance * distance ) *
            ( ( 2.0 * gravitationalParameter / distance ) * position );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > RelativisticAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    // Create partial function if parameter is primary body gravitational parameter (de Sitter term)
    if( calculateDeSitterCorrection_ && parameter->getParameterName( ).second.first == primaryBodyName_ &&
        parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter && primaryBodyName_ != acceleratingBody_ )
    {
        partialFunction =
                std::bind( &RelativisticAccelerationPartial::wrtGravitationalParameterOfPrimaryBody, this, std::placeholders::_1 );
        numberOfRows = 1;
    }
    // Create partial function if parameter is central body gravitational parameter
    else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
    {
        switch( parameter->getParameterName( ).first )
        {
            case estimatable_parameters::gravitational_parameter:
                partialFunction =
                        std::bind( &RelativisticAccelerationPartial::wrtGravitationalParameterOfCentralBody, this, std::placeholders::_1 );
                numberOfRows = 1;
                break;
            default:
                break;
        }
    }
    // Create partial function if parameter is PPN parameter beta or gamma
    else if( parameter->getParameterName( ).second.first == "global_metric" )
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

//! Function to compute partial derivative of relativistic acceleration correction w.r.t. central body gravitational parameter.
void RelativisticAccelerationPartial::wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, 1 );
    if( calculateSchwarzschildCorrection_ )
    {
        computePartialOfSchwarzschildAccelerationCorrectionWrtGravitationalParameter( currentRelativeState_,
                                                                                      centralBodyGravitationalParameterFunction_( ),
                                                                                      partialMatrix,
                                                                                      ppnGammaParameterFunction_( ),
                                                                                      ppnBetaParameterFunction_( ) );
    }
    // The Lense-Thirring term is linear in the central body gravitational parameter
    partialMatrix += currentLenseThirringAcceleration_ / centralBodyGravitationalParameterFunction_( );
}

//! Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter gamma
void RelativisticAccelerationPartial::wrtPpnParameterGamma( Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, 1 );
    if( calculateSchwarzschildCorrection_ )
    {
        computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterGamma(
                currentRelativeState_, centralBodyGravitationalParameterFunction_( ), partialMatrix );
    }
    // The Lense-Thirring and de Sitter terms are proportional to ( 1 + gamma ) and ( 1 + 2 gamma ), respectively
    const double ppnParameterGamma = ppnGammaParameterFunction_( );
    partialMatrix += currentLenseThirringAcceleration_ / ( 1.0 + ppnParameterGamma ) +
            2.0 * currentDeSitterAcceleration_ / ( 1.0 + 2.0 * ppnParameterGamma );
}

//! Function to compute partial derivative of relativistic acceleration correction w.r.t. PPN parameter beta
void RelativisticAccelerationPartial::wrtPpnParameterBeta( Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, 1 );
    if( calculateSchwarzschildCorrection_ )
    {
        computePartialOfSchwarzschildAccelerationCorrectionWrtPpnParameterBeta(
                currentRelativeState_, centralBodyGravitationalParameterFunction_( ), partialMatrix );
    }
}

//! Function for updating partial w.r.t. the bodies' states
void RelativisticAccelerationPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        currentRelativeState_ = ( acceleratedBodyState_( ) - centralBodyState_( ) );
        const double gravitationalParameter = centralBodyGravitationalParameterFunction_( );
        const double ppnParameterGamma = ppnGammaParameterFunction_( );
        currentPartialWrtPosition_.setZero( );
        currentPartialWrtVelocity_.setZero( );

        if( calculateSchwarzschildCorrection_ )
        {
            currentAcceleration_ = relativity::calculateScharzschildGravitationalAccelerationCorrection(
                    gravitationalParameter, currentRelativeState_, ppnParameterGamma, ppnBetaParameterFunction_( ) );
            computePartialOfSchwarzschildAccelerationCorrectionWrtPosition( currentRelativeState_,
                                                                            currentAcceleration_,
                                                                            currentPartialWrtPosition_,
                                                                            gravitationalParameter,
                                                                            ppnParameterGamma,
                                                                            ppnBetaParameterFunction_( ) );
            computePartialOfSchwarzschildAccelerationCorrectionWrtVelocity(
                    currentRelativeState_, currentPartialWrtVelocity_, gravitationalParameter, ppnParameterGamma );
        }

        // Lense-Thirring term a = k ( 3 / r^5 ( r x v ) ( r . J ) + 1 / r^3 ( v x J ) ), with k = ( 1 + gamma ) mu / c^2 and J the
        // (state-independent) angular momentum per unit mass of the central body
        currentLenseThirringAcceleration_.setZero( );
        if( calculateLenseThirringCorrection_ )
        {
            const Eigen::Vector3d position = currentRelativeState_.segment( 0, 3 );
            const Eigen::Vector3d velocity = currentRelativeState_.segment( 3, 3 );
            const Eigen::Vector3d angularMomentum = centralBodyAngularMomentumFunction_( );
            const double distance = position.norm( );
            const double coefficient =
                    ( 1.0 + ppnParameterGamma ) * gravitationalParameter * physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT;
            const double inverseDistanceFifth = 1.0 / std::pow( distance, 5 );
            const double positionDotAngularMomentum = position.dot( angularMomentum );
            const Eigen::Vector3d positionCrossVelocity = position.cross( velocity );

            currentLenseThirringAcceleration_ = relativity::calculateLenseThirringCorrectionAcceleration(
                    gravitationalParameter, currentRelativeState_, angularMomentum, ppnParameterGamma );
            currentPartialWrtPosition_ += coefficient *
                    ( 3.0 * inverseDistanceFifth *
                              ( -positionDotAngularMomentum * linear_algebra::getCrossProductMatrix( velocity ) +
                                positionCrossVelocity * angularMomentum.transpose( ) -
                                velocity.cross( angularMomentum ) * position.transpose( ) ) -
                      15.0 * inverseDistanceFifth / ( distance * distance ) * positionDotAngularMomentum * positionCrossVelocity *
                              position.transpose( ) );
            currentPartialWrtVelocity_ += coefficient *
                    ( 3.0 * inverseDistanceFifth * positionDotAngularMomentum * linear_algebra::getCrossProductMatrix( position ) -
                      linear_algebra::getCrossProductMatrix( angularMomentum ) / ( distance * distance * distance ) );
        }

        // de Sitter term a = -K ( ( V x R ) x v ), with K = ( 1 + 2 gamma ) mu_P / ( c^2 R^3 ), R and V the state of the central body
        // w.r.t. the primary body P, and v the relative velocity of the accelerated body
        currentDeSitterAcceleration_.setZero( );
        currentDeSitterPartialWrtCentralBodyPosition_.setZero( );
        currentDeSitterPartialWrtCentralBodyVelocity_.setZero( );
        if( calculateDeSitterCorrection_ )
        {
            const Eigen::Vector6d centralBodyStateWrtPrimary = centralBodyState_( ) - primaryBodyState_( );
            const Eigen::Vector3d primaryPosition = centralBodyStateWrtPrimary.segment( 0, 3 );
            const Eigen::Vector3d primaryVelocity = centralBodyStateWrtPrimary.segment( 3, 3 );
            const Eigen::Vector3d velocity = currentRelativeState_.segment( 3, 3 );
            const double primaryDistance = primaryPosition.norm( );
            const double coefficient = ( 1.0 + 2.0 * ppnParameterGamma ) * primaryBodyGravitationalParameterFunction_( ) *
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT / std::pow( primaryDistance, 3 );
            const Eigen::Vector3d precessionVector = primaryVelocity.cross( primaryPosition );

            currentDeSitterAcceleration_ = -coefficient * precessionVector.cross( velocity );
            currentPartialWrtVelocity_ -= coefficient * linear_algebra::getCrossProductMatrix( precessionVector );
            currentDeSitterPartialWrtCentralBodyPosition_ = 3.0 * coefficient / ( primaryDistance * primaryDistance ) *
                            precessionVector.cross( velocity ) * primaryPosition.transpose( ) +
                    coefficient * linear_algebra::getCrossProductMatrix( velocity ) *
                            linear_algebra::getCrossProductMatrix( primaryVelocity );
            currentDeSitterPartialWrtCentralBodyVelocity_ = -coefficient * linear_algebra::getCrossProductMatrix( velocity ) *
                    linear_algebra::getCrossProductMatrix( primaryPosition );
        }
        currentTime_ = currentTime;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
