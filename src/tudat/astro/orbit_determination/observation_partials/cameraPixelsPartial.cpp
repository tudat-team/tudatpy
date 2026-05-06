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
#include "tudat/astro/orbit_determination/observation_partials/cameraPixelsPartial.h"

namespace tudat
{
namespace observation_partials
{

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

void CameraPixelsScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                  const std::vector< double >& times,
                                  const observation_models::LinkEndType fixedLinkEnd,
                                  const Eigen::VectorXd )
{
    if( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error(
                "Error, only receiver link end type is currently supported as fixed link end for camera pixels partial scaling." );
    }

    Eigen::Vector3d relativeRangeVector =
            linkEndStates.at( observedBodyIndex_ ).segment( 0, 3 ) - linkEndStates.at( observerIndex_ ).segment( 0, 3 );
    Eigen::Vector3d normalizedRelativeRangeVector = relativeRangeVector.normalized( );

    Eigen::Quaterniond rotationFromInertialToCameraFrame = rotationFromInertialToCameraFrameFunction_( times.at( observerIndex_ ) );
    Eigen::Vector3d relativeRangeVectorCameraFrame = rotationFromInertialToCameraFrame * relativeRangeVector;
    positionScalingFactor_ =
            calculatePartialOfPixelsWrtLinkEndPositionCameraFrame( relativeRangeVectorCameraFrame, focalLengthsMatrix_, true );

    positionScalingFactor_ = positionScalingFactor_ * rotationFromInertialToCameraFrame.toRotationMatrix( );

    // Only because receiver is currently the only supported fixed link end type for camera pixels observations.
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

}  // namespace observation_partials
}  // namespace tudat
