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
    ellipsoidalHarmonicsAccelerationModel_( accelerationModel ),
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

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
EllipsoidalHarmonicsGravityPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    const std::string& parameterBody = parameter->getParameterName( ).second.first;

    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
    {
        // With mutual attraction the accelerated body's parameter is part of the sum as well.
        if( parameterBody == acceleratingBody_ ||
            ( parameterBody == acceleratedBody_ && ellipsoidalHarmonicsAccelerationModel_->getIsMutualAttractionUsed( ) ) )
        {
            partialFunction = std::bind( &EllipsoidalHarmonicsGravityPartial::wrtGravitationalParameter, this, std::placeholders::_1 );
            numberOfColumns = 1;
        }
    }
    else if( parameterBody == acceleratingBody_ &&
             estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
    {
        if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) ) ) ==
            0 )
        {
            throw std::runtime_error(
                    "Error when making ellipsoidal harmonics acceleration partial, no rotation matrix partial "
                    "available for parameter " +
                    std::to_string( parameter->getParameterName( ).first ) + " of " + parameterBody + "." );
        }

        partialFunction = std::bind( &EllipsoidalHarmonicsGravityPartial::wrtRotationModelParameter,
                                     this,
                                     std::placeholders::_1,
                                     parameter->getParameterName( ).first,
                                     parameter->getSecondaryIdentifier( ) );
        numberOfColumns = 1;
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
EllipsoidalHarmonicsGravityPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfColumns = 0;

    const std::string& parameterBody = parameter->getParameterName( ).second.first;

    if( parameterBody == acceleratingBody_ &&
        estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
    {
        if( rotationMatrixPartials_.count( std::make_pair( parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) ) ) ==
            0 )
        {
            throw std::runtime_error(
                    "Error when making ellipsoidal harmonics acceleration partial, no rotation matrix partial "
                    "available for parameter " +
                    std::to_string( parameter->getParameterName( ).first ) + " of " + parameterBody + "." );
        }

        partialFunction = std::bind( &EllipsoidalHarmonicsGravityPartial::wrtRotationModelParameter,
                                     this,
                                     std::placeholders::_1,
                                     parameter->getParameterName( ).first,
                                     parameter->getSecondaryIdentifier( ) );
        numberOfColumns = parameter->getParameterSize( );
    }

    return std::make_pair( partialFunction, numberOfColumns );
}

void EllipsoidalHarmonicsGravityPartial::wrtGravitationalParameter( Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = accelerationFunction_( ) / gravitationalParameterFunction_( );
}

void EllipsoidalHarmonicsGravityPartial::wrtRotationModelParameter( Eigen::MatrixXd& partialMatrix,
                                                                    const estimatable_parameters::EstimatebleParametersEnum parameterType,
                                                                    const std::string& secondaryIdentifier )
{
    const Eigen::Vector3d relativePosition = positionFunctionOfAcceleratedBody_( ) - positionFunctionOfAcceleratingBody_( );
    const Eigen::Matrix3d rotationToIntegrationFrame = currentRotationToBodyFixedFrame_.transpose( );

    // Each column is the derivative of the rotation matrix w.r.t. one entry of the parameter.
    const std::vector< Eigen::Matrix3d > rotationMatrixPartials =
            rotationMatrixPartials_.at( std::make_pair( parameterType, secondaryIdentifier ) )
                    ->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime_ );

    partialMatrix.setZero( 3, rotationMatrixPartials.size( ) );
    for( unsigned int i = 0; i < rotationMatrixPartials.size( ); i++ )
    {
        // The rotation enters twice: it turns the body-fixed acceleration into the propagation
        // frame, and it sets where the position lies in the body-fixed frame.
        partialMatrix.col( i ) = rotationMatrixPartials[ i ] * ( currentRotationToBodyFixedFrame_ * accelerationFunction_( ) ) +
                rotationToIntegrationFrame * currentBodyFixedPartialWrtPosition_ * rotationMatrixPartials[ i ].transpose( ) *
                        relativePosition;
    }
}

void EllipsoidalHarmonicsGravityPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        // Update acceleration model
        updateFunction_( currentTime );

        // Calculate Cartesian position in frame fixed to body exerting acceleration
        currentRotationToBodyFixedFrame_ = fromBodyFixedToIntegrationFrameRotation_( ).transpose( );
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
