/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cmath>
#include "tudat/astro/observation_models/pixelCoordinatesObservationModel.h"
#include "tudat/astro/orbit_determination/observation_partials/pixelCoordinatesPartial.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{
namespace observation_partials
{

namespace
{

Eigen::Matrix3d getLeftJacobianOfSO3( const Eigen::Vector3d& rotationVector )
{
    const double angle = rotationVector.norm( );
    const Eigen::Matrix3d crossProduct = linear_algebra::getCrossProductMatrix( rotationVector );
    const Eigen::Matrix3d crossProductSquared = crossProduct * crossProduct;
    if( angle < 1.0E-4 )
    {
        // Avoid cancellation in 1-cos(angle) and angle-sin(angle) near zero.
        const double angleSquared = angle * angle;
        return Eigen::Matrix3d::Identity( ) + ( 0.5 - angleSquared / 24.0 ) * crossProduct +
                ( 1.0 / 6.0 - angleSquared / 120.0 ) * crossProductSquared;
    }
    return Eigen::Matrix3d::Identity( ) + ( 1.0 - std::cos( angle ) ) / ( angle * angle ) * crossProduct +
            ( angle - std::sin( angle ) ) / ( angle * angle * angle ) * crossProductSquared;
}

}  // namespace

Eigen::Matrix< double, 1, 3 > calculatePartialOfPixelUWrtObserverPositionCameraFrame(
        Eigen::Vector3d relativeRangeVectorCameraFrame,
        Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix )
{
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = -1.0 / relativeRangeVectorCameraFrame( 2 );
    partial( 2 ) = relativeRangeVectorCameraFrame( 0 ) / std::pow( relativeRangeVectorCameraFrame( 2 ), 2 );
    partial = focalLengthsMatrix.coeff( 0, 0 ) * partial;
    return partial;
}

Eigen::Matrix< double, 1, 3 > calculatePartialOfPixelVWrtObserverPositionCameraFrame(
        Eigen::Vector3d relativeRangeVectorCameraFrame,
        Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix )
{
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 1 ) = -1.0 / relativeRangeVectorCameraFrame( 2 );
    partial( 2 ) = relativeRangeVectorCameraFrame( 1 ) / std::pow( relativeRangeVectorCameraFrame( 2 ), 2 );
    partial = focalLengthsMatrix.coeff( 1, 1 ) * partial;
    return partial;
}

Eigen::Matrix< double, 2, 3 > calculatePartialOfPixelsWrtLinkEndPositionCameraFrame( Eigen::Vector3d relativeRangeVectorCameraFrame,
                                                                                     Eigen::DiagonalMatrix< double, 2 > focalLengthsMatrix,
                                                                                     const bool isLinkEndObserver )
{
    double partial_multiplier = ( isLinkEndObserver ) ? 1 : -1;
    Eigen::Matrix< double, 2, 3 > partial;
    partial.block( 0, 0, 1, 3 ) =
            calculatePartialOfPixelUWrtObserverPositionCameraFrame( relativeRangeVectorCameraFrame, focalLengthsMatrix );
    partial.block( 1, 0, 1, 3 ) =
            calculatePartialOfPixelVWrtObserverPositionCameraFrame( relativeRangeVectorCameraFrame, focalLengthsMatrix );
    partial *= partial_multiplier;
    return partial;
}

void PixelCoordinatesScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                      const std::vector< double >& times,
                                      const observation_models::LinkEndType fixedLinkEnd,
                                      const Eigen::VectorXd )
{
    if( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error(
                "Error, only receiver link end type is currently supported as fixed link end for pixel coordinates partial scaling." );
    }

    Eigen::Vector3d relativeRangeVector =
            linkEndStates.at( observedBodyIndex_ ).segment( 0, 3 ) - linkEndStates.at( observerIndex_ ).segment( 0, 3 );
    Eigen::Vector3d normalizedRelativeRangeVector = relativeRangeVector.normalized( );

    currentReceiverTime_ = times.at( observerIndex_ );
    currentRotationFromInertialToCameraFrame_ = rotationFromInertialToCameraFrameFunction_( currentReceiverTime_ ).toRotationMatrix( );
    currentInertialDirectionToProject_ = relativeRangeVector;
    Eigen::Matrix3d directionPartialWrtRelativePosition = Eigen::Matrix3d::Identity( );
    Eigen::Matrix3d directionPartialWrtReceiverVelocity = Eigen::Matrix3d::Zero( );
    if( correctForStellarAberration_ )
    {
        const Eigen::Vector3d beta = linkEndStates.at( observerIndex_ ).segment( 3, 3 ) / physical_constants::SPEED_OF_LIGHT;
        currentInertialDirectionToProject_ =
                observation_models::computeApparentDirectionWithStellarAberration( normalizedRelativeRangeVector, beta );

        const double directionDotVelocity = normalizedRelativeRangeVector.dot( beta );
        const double squareRoot = std::sqrt( 1.0 - beta.squaredNorm( ) + directionDotVelocity * directionDotVelocity );
        const double scaleFactor = squareRoot - directionDotVelocity;
        const Eigen::Vector3d unnormalizedApparentDirection = scaleFactor * normalizedRelativeRangeVector + beta;
        const Eigen::Matrix3d normalizationPartial =
                ( Eigen::Matrix3d::Identity( ) - currentInertialDirectionToProject_ * currentInertialDirectionToProject_.transpose( ) ) /
                unnormalizedApparentDirection.norm( );
        const Eigen::Matrix3d apparentDirectionPartialWrtUnitDirection = scaleFactor * Eigen::Matrix3d::Identity( ) +
                normalizedRelativeRangeVector * ( ( directionDotVelocity / squareRoot - 1.0 ) * beta.transpose( ) );
        directionPartialWrtRelativePosition = normalizationPartial * apparentDirectionPartialWrtUnitDirection *
                ( Eigen::Matrix3d::Identity( ) - normalizedRelativeRangeVector * normalizedRelativeRangeVector.transpose( ) ) /
                relativeRangeVector.norm( );
        directionPartialWrtReceiverVelocity = normalizationPartial *
                ( Eigen::Matrix3d::Identity( ) +
                  normalizedRelativeRangeVector *
                          ( -normalizedRelativeRangeVector + ( directionDotVelocity * normalizedRelativeRangeVector - beta ) / squareRoot )
                                  .transpose( ) ) /
                physical_constants::SPEED_OF_LIGHT;
    }

    currentRelativeRangeVectorCameraFrame_ = currentRotationFromInertialToCameraFrame_ * currentInertialDirectionToProject_;
    currentPixelLinePartialWrtCameraFramePosition_ =
            pixelLinePartialWrtCameraFramePositionFunction_( currentRelativeRangeVectorCameraFrame_ );
    hasCurrentCameraFrameQuantities_ = true;
    const Eigen::Matrix< double, 2, 3 > pixelPartialWrtInertialDirection =
            currentPixelLinePartialWrtCameraFramePosition_ * currentRotationFromInertialToCameraFrame_;
    positionScalingFactor_ = -pixelPartialWrtInertialDirection * directionPartialWrtRelativePosition;
    velocityScalingFactor_ = pixelPartialWrtInertialDirection * directionPartialWrtReceiverVelocity;

    // Only because receiver is currently the only supported fixed link end type for pixel coordinates observations.
    lightTimeCorrectionScalingFactor_ = positionScalingFactor_ * linkEndStates.at( observedBodyIndex_ ).segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT -
              linkEndStates.at( observedBodyIndex_ ).segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );

    realPositionScalingFactor_ = positionScalingFactor_ *
            ( Eigen::Matrix3d::Identity( ) +
              linkEndStates.at( observedBodyIndex_ ).segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                      ( physical_constants::SPEED_OF_LIGHT -
                        linkEndStates.at( observedBodyIndex_ ).segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );

    currentLinkEndType_ = fixedLinkEnd;
}

bool PixelCoordinatesScaling::useLinkIndependentPartials(
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier )
{
    return ( pointingCorrectionFunction_ && parameterIdentifier == pointingCorrectionParameterIdentifier_ ) ||
            ( parameterIdentifier.second.first == pointingReferenceBodyName_ &&
              rotationMatrixPartials_.count( std::make_pair( parameterIdentifier.first, parameterIdentifier.second.second ) ) != 0 );
}

std::vector< std::pair< Eigen::Matrix< double, 2, Eigen::Dynamic >, double > > PixelCoordinatesScaling::getLinkIndependentPartials(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier )
{
    if( !useLinkIndependentPartials( parameterIdentifier ) )
    {
        return {};
    }
    if( !hasCurrentCameraFrameQuantities_ )
    {
        throw std::runtime_error( "Error when computing pixel-coordinate orientation partial: scaling object has not been updated." );
    }
    if( pointingCorrectionFunction_ && parameterIdentifier == pointingCorrectionParameterIdentifier_ )
    {
        // For r_camera = Exp(correction) * r_nominal, differentiate additive changes of the rotation vector.
        const Eigen::Matrix< double, 2, Eigen::Dynamic > partial = -currentPixelLinePartialWrtCameraFramePosition_ *
                linear_algebra::getCrossProductMatrix( currentRelativeRangeVectorCameraFrame_ ) *
                getLeftJacobianOfSO3( pointingCorrectionFunction_( ) );
        return { std::make_pair( partial, currentReceiverTime_ ) };
    }

    const auto partialIterator =
            rotationMatrixPartials_.find( std::make_pair( parameterIdentifier.first, parameterIdentifier.second.second ) );
    const std::vector< Eigen::Matrix3d > rotationPartials =
            partialIterator->second->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentReceiverTime_ );
    const Eigen::Matrix< double, 2, 3 > pixelPartialWrtReferenceFrameDirection = currentPixelLinePartialWrtCameraFramePosition_ *
            currentRotationFromInertialToCameraFrame_ *
            rotationFromReferenceFrameToInertialFunction_( currentReceiverTime_ ).toRotationMatrix( );
    Eigen::Matrix< double, 2, Eigen::Dynamic > partial( 2, rotationPartials.size( ) );
    for( unsigned int i = 0; i < rotationPartials.size( ); ++i )
    {
        // R_inertial_to_camera = R_reference_to_camera * R_inertial_to_reference.
        // RotationMatrixPartial supplies dR_reference_to_inertial/dp, so transpose its derivative.
        partial.col( i ) =
                pixelPartialWrtReferenceFrameDirection * rotationPartials.at( i ).transpose( ) * currentInertialDirectionToProject_;
    }
    return { std::make_pair( partial, currentReceiverTime_ ) };
}

}  // namespace observation_partials
}  // namespace tudat
