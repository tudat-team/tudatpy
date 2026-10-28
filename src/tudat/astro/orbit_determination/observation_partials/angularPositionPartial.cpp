/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition( const Eigen::Vector3d& relativeRangeVector,
                                                                                  const bool isLinkEndReceiver )
{
    // Define multiplier of patial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = -relativeRangeVector( 1 );
    partial( 1 ) = relativeRangeVector( 0 );
    partial /= ( partialMultiplier *
                 ( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) + relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );
    return partial;
}

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition( Eigen::Vector3d relativeRangeVector,
                                                                               const bool isLinkEndReceiver )
{
    // Define multiplier of partial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    double range = relativeRangeVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    double inPlaneRange = relativeRangeVector.segment( 0, 2 ).norm( );
    partial = partialMultiplier * relativeRangeVector.transpose( ) / inPlaneRange;
    partial *= relativeRangeVector( 2 ) / ( range * range );
    partial += -partialMultiplier * ( Eigen::Vector3d::UnitZ( ) ).transpose( ) / inPlaneRange;

    return partial;
}

//! Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
//! observed object.
Eigen::Matrix< double, 2, 3 > calculatePartialOfAngularPositionWrtLinkEndPosition( Eigen::Vector3d relativeRangeVector,
                                                                                   const bool isLinkEndReceiver,
                                                                                   const bool normalizeRightAscension,
                                                                                   const Eigen::VectorXd& currentObservation )
{
    Eigen::Matrix< double, 2, 3 > angularPositionPartial;
    angularPositionPartial.block( 0, 0, 1, 3 ) =
            calculatePartialOfRightAscensionWrtLinkEndPosition( relativeRangeVector, isLinkEndReceiver );
    angularPositionPartial.block( 1, 0, 1, 3 ) = calculatePartialOfDeclinationWrtLinkEndPosition( relativeRangeVector, isLinkEndReceiver );
    if( normalizeRightAscension )
    {
        if( currentObservation.rows( ) == 0 )
        {
            throw std::runtime_error( "Error when getting angular position partial scaling for normalized observable; incompatible input" );
        }
        double rightAscension = currentObservation( 0 );
        double declination = currentObservation( 1 );
        Eigen::Matrix< double, 2, 3 > normalizedAngularPositionPartial = angularPositionPartial;

        normalizedAngularPositionPartial.block( 0, 0, 1, 3 ) = std::cos( declination ) * angularPositionPartial.block( 0, 0, 1, 3 ) -
                rightAscension * std::sin( declination ) * angularPositionPartial.block( 1, 0, 1, 3 );
        return normalizedAngularPositionPartial;
    }
    else
    {
        return angularPositionPartial;
    }
}

//! Update the scaling object to the current times and states
void AngularPositionScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                     const std::vector< double >& times,
                                     const observation_models::LinkEndType fixedLinkEnd,
                                     const Eigen::VectorXd currentObservation,
                                     const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings )
{
    Eigen::Vector3d relativeRangeVector = ( linkEndStates[ 1 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d normalizedRelativeRangeVector = relativeRangeVector.normalized( );

    // Compute common scaling factor
    velocityScalingFactor_.setZero( );
    if( !observation_models::useAngularPositionFrameSettings( ancillarySettings ) )
    {
        scalingFactor_ = calculatePartialOfAngularPositionWrtLinkEndPosition(
                relativeRangeVector, true, normalizeRightAscension_, currentObservation );
    }
    else
    {
        // Chain rule through the (aberrated) direction, rotated to the selected celestial frame
        const Eigen::Vector3d relativePosition = -relativeRangeVector;
        const Eigen::Vector3d receiverVelocity = linkEndStates[ 1 ].segment( 3, 3 );
        Eigen::Matrix3d directionPartialWrtPosition = Eigen::Matrix3d::Identity( );
        Eigen::Matrix3d directionPartialWrtVelocity = Eigen::Matrix3d::Zero( );
        if( observation_models::getPositionAngleUseStellarAberration( ancillarySettings ) )
        {
            observation_models::calculateStellarAberrationPartials(
                    relativePosition, receiverVelocity, directionPartialWrtPosition, directionPartialWrtVelocity );
        }
        const Eigen::Matrix3d rotation =
                observation_models::getPositionAngleReferenceFrameRotationFromJ2000( times[ 1 ], ancillarySettings ) *
                j2000ToGlobalFrameTransformation_.transpose( );
        const Eigen::Vector3d direction = observation_models::getAngularPositionDirectionInSelectedFrame(
                relativePosition, receiverVelocity, times[ 1 ], j2000ToGlobalFrameTransformation_, ancillarySettings );
        const Eigen::Matrix< double, 2, 3 > partialWrtReceiverToDirection = calculatePartialOfAngularPositionWrtLinkEndPosition(
                -direction, true, normalizeRightAscension_, currentObservation );
        scalingFactor_ = partialWrtReceiverToDirection * rotation * directionPartialWrtPosition;
        velocityScalingFactor_ = -partialWrtReceiverToDirection * rotation * directionPartialWrtVelocity;
    }

    // Compute scaling for receiver reference
    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 0 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ = scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) +
                  linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                          ( physical_constants::SPEED_OF_LIGHT -
                            linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }
    // Compute scaling for transmitter reference
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 1 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ = scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) +
                  linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                          ( physical_constants::SPEED_OF_LIGHT -
                            linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }

    currentLinkEndType_ = fixedLinkEnd;
}

}  // namespace observation_partials

}  // namespace tudat
