/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/relativeAngularPositionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition2(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver )
{
    // Define multiplier of patial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = - relativeRangeVector( 1 );
    partial( 1 ) = relativeRangeVector( 0 );
    partial /= ( partialMultiplier * ( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) +
                                       relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );
    return partial;
}

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition2(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver )
{
    // Define multiplier of partial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    double range = relativeRangeVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = relativeRangeVector( 0 ) * relativeRangeVector( 2 );
    partial( 1 ) = relativeRangeVector( 1 ) * relativeRangeVector( 2 );
    partial( 2 ) = - ( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) + relativeRangeVector( 1 ) * relativeRangeVector( 1 ) );
    partial *= partialMultiplier /
               ( range * range * std::sqrt( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) + relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );

    return partial;
}

//! Function to compute the derivative of (direct geometric) right ascension and declination between two sources w.r.t. position of observer.
Eigen::Matrix< double, 2, 3 > calculatePartialOfAngularPositionWrtLinkEndPosition2(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver )
{
    Eigen::Matrix< double, 2, 3 > angularPositionPartial;
    angularPositionPartial.block( 0, 0, 1, 3 ) = calculatePartialOfRightAscensionWrtLinkEndPosition2(
                relativeRangeVector, isLinkEndReceiver );
    angularPositionPartial.block( 1, 0, 1, 3 ) = calculatePartialOfDeclinationWrtLinkEndPosition2(
                relativeRangeVector, isLinkEndReceiver );
    return angularPositionPartial;
}

//! Update the scaling object to the current times and states
void RelativeAngularPositionScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                             const std::vector< double >& times,
                                             const observation_models::LinkEndType fixedLinkEnd,
                                             const Eigen::VectorXd currentObservation )
{
    if ( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating a relative angular position scaling object, fixed link end time different from receiver." );
    }

    Eigen::Vector6d firstTransmitterState = linkEndStates[ 0 ];
    Eigen::Vector6d secondTransmitterState = linkEndStates[ 1 ];
    Eigen::Vector6d receiverState = linkEndStates[ 2 ];

    Eigen::Vector3d relativeRangeVectorFirstTransmitter = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativeRangeVectorSecondTransmitter = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    Eigen::Vector3d normalizedRelativeRangeVectorFirstTransmitter = relativeRangeVectorFirstTransmitter.normalized( );
    Eigen::Vector3d normalizedRelativeRangeVectorSecondTransmitter = relativeRangeVectorSecondTransmitter.normalized( );

    // Compute common scaling factor
    scalingFactorFirstTransmitter_ = calculatePartialOfAngularPositionWrtLinkEndPosition2( relativeRangeVectorFirstTransmitter, true );
    scalingFactorSecondTransmitter_ = calculatePartialOfAngularPositionWrtLinkEndPosition2( relativeRangeVectorSecondTransmitter, true );

    // Compute scaling for receiver reference
//    if( fixedLinkEnd == observation_models::receiver )
//    {
    referenceLightTimeCorrectionScalingFirstTransmitter_ = scalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorFirstTransmitter ) );
    referenceScalingFactorFirstTransmitter_ =
            scalingFactorFirstTransmitter_ *
            ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVectorFirstTransmitter.transpose( ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorFirstTransmitter ) ) );

    referenceLightTimeCorrectionScalingSecondTransmitter_ = scalingFactorSecondTransmitter_ * linkEndStates[ 1 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorSecondTransmitter ) );
    referenceScalingFactorSecondTransmitter_ = scalingFactorSecondTransmitter_ *
            ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVectorSecondTransmitter.transpose( ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVectorSecondTransmitter ) ) );
//    }


    currentLinkEndType_ = fixedLinkEnd;

}

//! Function to calculate the observation partial(s) at required time and state
RelativeAngularPositionPartial::RelativeAngularPositionPartialReturnType RelativeAngularPositionPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector2d& currentObservation )
{
    if( linkEndOfFixedTime != relativeAngularPositionScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error relative angular position partial and scaling are inconsistent" );
    }

    RelativeAngularPositionPartialReturnType returnPartial;

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        if( positionPartialIterator_->first == observation_models::transmitter )
        {
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];
        }
        else if( positionPartialIterator_->first == observation_models::transmitter2 )
        {
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];
        }
        else if( positionPartialIterator_->first == observation_models::receiver )
        {
            currentState_  = states[ 2 ];
            currentTime_ = times[ 2 ];
        }

        // Scale position partials
        returnPartial.push_back(
                    std::make_pair(
                        relativeAngularPositionScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }

    // Add scaled light-time correction partials.
    for( unsigned int i = 0; i < lightTimeCorrectionPartialsFunctionsFirstTransmitter_.size( ); i++ )
    {
        currentLinkTimeCorrectionPartialFirstTransmitter_ = lightTimeCorrectionPartialsFunctionsFirstTransmitter_.at( i )( states, times );
        currentLinkTimeCorrectionPartialSecondTransmitter_ = lightTimeCorrectionPartialsFunctionsSecondTransmitter_.at( i )( states, times );

        if ( currentLinkTimeCorrectionPartialFirstTransmitter_.second !=
             currentLinkTimeCorrectionPartialSecondTransmitter_.second )
        {
            throw std::runtime_error( "Error when making relative angular position light time correction partials, inconsistency"
                                      " in receiver times between receiver - first transmitter and receiver - second transmitter legs." );
        }

        returnPartial.push_back(
                    std::make_pair( ( relativeAngularPositionScaler_->getLightTimePartialScalingFactorSecondTransmitter( ) *
                                    physical_constants::SPEED_OF_LIGHT * currentLinkTimeCorrectionPartialSecondTransmitter_.first )
                                    - ( relativeAngularPositionScaler_->getLightTimePartialScalingFactorFirstTransmitter( ) *
                                    physical_constants::SPEED_OF_LIGHT * currentLinkTimeCorrectionPartialFirstTransmitter_.first ) ,
                    currentLinkTimeCorrectionPartialFirstTransmitter_.second ) );
    }


    return returnPartial;
}

}

}
