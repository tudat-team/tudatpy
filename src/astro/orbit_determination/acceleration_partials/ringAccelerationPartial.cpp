/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/ringAccelerationPartial.h"

#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"

namespace tudat
{
namespace acceleration_partials
{

RingGravityPartial::RingGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::RingGravitationalAccelerationModel > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials ):
    AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::ring_gravity ),
    gravitationalParameterFunction_( accelerationModel->getGravitationalParameterFunction( ) ),
    ringRadiusFunction_( accelerationModel->getRingRadiusFunction( ) ),
    ringCache_( accelerationModel->getRingCache( ) ),
    positionFunctionOfAcceleratedBody_( std::bind( &gravitation::RingGravitationalAccelerationModel::
                                                   getCurrentPositionOfBodySubjectToAcceleration, accelerationModel ) ),
    positionFunctionOfAcceleratingBody_( std::bind( &gravitation::RingGravitationalAccelerationModel::
                                                    getCurrentPositionOfBodyExertingAcceleration, accelerationModel ) ),
    fromBodyFixedToIntegrationFrameRotation_( std::bind( &gravitation::RingGravitationalAccelerationModel::
                                                         getCurrentRotationToIntegrationFrameMatrix, accelerationModel ) ),
    accelerationFunction_( std::bind( &gravitation::RingGravitationalAccelerationModel::getAcceleration,
                                      accelerationModel ) ),
    updateFunction_( std::bind( &gravitation::RingGravitationalAccelerationModel::updateMembers,
                                accelerationModel, std::placeholders::_1 ) ),
    rotationMatrixPartials_( rotationMatrixPartials )
{

}

void RingGravityPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        // Update acceleration model
        updateFunction_( currentTime );

        // Calculate Cartesian position in frame fixed to body exerting acceleration
        Eigen::Matrix3d currentRotationToBodyFixedFrame_ = fromBodyFixedToIntegrationFrameRotation_( ).inverse( );
        bodyFixedPosition_ = currentRotationToBodyFixedFrame_ *
                ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) );

        // Calculate partial of acceleration wrt position of body undergoing acceleration.
        currentBodyFixedPartialWrtPosition_ = gravitation::computeRingHessianOfGravitationalPotential(
                bodyFixedPosition_,
                ringRadiusFunction_( ),
                gravitationalParameterFunction_( ),
                ringCache_->getEllipticIntegralB( ),
                ringCache_->getEllipticIntegralE( ),
                ringCache_->getEllipticIntegralS( ),
                ringCache_->getEllipticIntegralK( ) );

        currentPartialWrtVelocity_.setZero( );
        currentPartialWrtPosition_.setZero( );

        // Compute partial w.r.t. position in inertial frame
        currentPartialWrtPosition_ +=
                currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;

        // If rotation matrix depends on translational state, add correction partials
        if( rotationMatrixPartials_.count(
                    std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
        {
            // Compute the acceleration and body-fixed position partial, without the central term (to avoid numerical errors)
            Eigen::Vector3d nonCentralAcceleration = accelerationFunction_( );
            Eigen::Matrix3d nonCentralBodyFixedPartial = currentBodyFixedPartialWrtPosition_;

            nonCentralAcceleration -= gravitation::computeGravitationalAcceleration(
                        positionFunctionOfAcceleratedBody_( ), gravitationalParameterFunction_( ), positionFunctionOfAcceleratingBody_( ) );

            nonCentralBodyFixedPartial -=
                    currentRotationToBodyFixedFrame_ * calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody(
                        positionFunctionOfAcceleratedBody_( ), positionFunctionOfAcceleratingBody_( ), gravitationalParameterFunction_( ) ) *
                    currentRotationToBodyFixedFrame_.inverse( );

            // Compute rotation matrix partials
            std::vector< Eigen::Matrix3d > rotationPositionPartials =
                    rotationMatrixPartials_.at(
                        std::make_pair( estimatable_parameters::initial_body_state, "" ) )->
                    calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

            // Add correction terms to position and velocity partials
            for( unsigned int i = 0; i < 3; i++ )
            {
                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );

                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        currentRotationToBodyFixedFrame_.inverse( ) * nonCentralBodyFixedPartial *
                        ( rotationPositionPartials.at( i ).transpose( ) *
                          ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) ) );
                currentPartialWrtVelocity_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i + 3 ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );
            }
        }

        currentTime_ = currentTime;
    }

}

} // namespace acceleration_partials

} // namespace tudat
