/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/ellipsoidalHarmonicsAccelerationPartial.h"

#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"

namespace tudat
{
namespace acceleration_partials
{

EllipsoidalHarmonicsGravityPartial::EllipsoidalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::EllipsoidalHarmonicsGravitationalAccelerationModel > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials ):
    AccelerationPartial( acceleratedBody, acceleratingBody, accelerationModel, basic_astrodynamics::ellipsoidal_harmonic_gravity ),
    gravitationalParameterFunction_( accelerationModel->getGravitationalParameterFunction( ) ),
    gravityField_( accelerationModel->getGravityField( ) ),
    positionFunctionOfAcceleratedBody_(
            std::bind( &gravitation::EllipsoidalHarmonicsGravitationalAccelerationModel::getCurrentPositionOfBodySubjectToAcceleration,
                       accelerationModel ) ),
    positionFunctionOfAcceleratingBody_(
            std::bind( &gravitation::EllipsoidalHarmonicsGravitationalAccelerationModel::getCurrentPositionOfBodyExertingAcceleration,
                       accelerationModel ) ),
    fromBodyFixedToIntegrationFrameRotation_(
            std::bind( &gravitation::EllipsoidalHarmonicsGravitationalAccelerationModel::getCurrentRotationToIntegrationFrameMatrix,
                       accelerationModel ) ),
    accelerationFunction_(
            std::bind( &gravitation::EllipsoidalHarmonicsGravitationalAccelerationModel::getAcceleration, accelerationModel ) ),
    updateFunction_( std::bind( &gravitation::EllipsoidalHarmonicsGravitationalAccelerationModel::updateMembers,
                                accelerationModel,
                                std::placeholders::_1 ) ),
    rotationMatrixPartials_( rotationMatrixPartials )
{}

void EllipsoidalHarmonicsGravityPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        // Update acceleration model
        updateFunction_( currentTime );

        // Calculate Cartesian position in frame fixed to body exerting acceleration
        Eigen::Matrix3d currentRotationToBodyFixedFrame_ = fromBodyFixedToIntegrationFrameRotation_( ).inverse( );
        bodyFixedPosition_ =
                currentRotationToBodyFixedFrame_ * ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) );

        // The gravitational parameter is part of the field's own coefficients, so a different
        // one, such as the sum of both bodies' parameters, means scaling the whole field. The
        // acceleration model does the same.
        const double scaling = gravitationalParameterFunction_( ) / gravityField_->getGravitationalParameter( );

        currentBodyFixedPartialWrtPosition_ = scaling * gravityField_->getHessianOfPotential( bodyFixedPosition_ );

        currentPartialWrtVelocity_.setZero( );
        currentPartialWrtPosition_.setZero( );

        // Compute partial w.r.t. position in inertial frame
        currentPartialWrtPosition_ +=
                currentRotationToBodyFixedFrame_.inverse( ) * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;

        // If rotation matrix depends on translational state, add correction partials
        if( rotationMatrixPartials_.count( std::make_pair( estimatable_parameters::initial_body_state, "" ) ) > 0 )
        {
            // Work with the field minus its point-mass part, which keeps the correction terms
            // free of the large cancelling central contribution.
            Eigen::Vector3d nonCentralAcceleration = accelerationFunction_( );
            Eigen::Matrix3d nonCentralBodyFixedPartial = currentBodyFixedPartialWrtPosition_;

            nonCentralAcceleration -= gravitation::computeGravitationalAcceleration(
                    positionFunctionOfAcceleratedBody_( ), gravitationalParameterFunction_( ), positionFunctionOfAcceleratingBody_( ) );

            nonCentralBodyFixedPartial -= currentRotationToBodyFixedFrame_ *
                    calculatePartialOfPointMassGravityWrtPositionOfAcceleratedBody( positionFunctionOfAcceleratedBody_( ),
                                                                                    positionFunctionOfAcceleratingBody_( ),
                                                                                    gravitationalParameterFunction_( ) ) *
                    currentRotationToBodyFixedFrame_.inverse( );

            // Compute rotation matrix partials
            std::vector< Eigen::Matrix3d > rotationPositionPartials =
                    rotationMatrixPartials_.at( std::make_pair( estimatable_parameters::initial_body_state, "" ) )
                            ->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime );

            // Add correction terms to position and velocity partials
            for( unsigned int i = 0; i < 3; i++ )
            {
                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );

                currentPartialWrtPosition_.block( 0, i, 3, 1 ) -= currentRotationToBodyFixedFrame_.inverse( ) * nonCentralBodyFixedPartial *
                        ( rotationPositionPartials.at( i ).transpose( ) *
                          ( positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( ) ) );
                currentPartialWrtVelocity_.block( 0, i, 3, 1 ) -=
                        rotationPositionPartials.at( i + 3 ) * ( currentRotationToBodyFixedFrame_ * nonCentralAcceleration );
            }
        }

        currentTime_ = currentTime;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
