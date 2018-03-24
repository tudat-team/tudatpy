/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/angularPositionPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to compute the derivative of (direct geometric) right ascension w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver )
{
    // Define multiplier of patial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = -relativeRangeVector( 1 );
    partial( 1 ) = relativeRangeVector( 0 );
    partial /= ( partialMultiplier * ( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) +
                                       relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );
    return partial;
}

//! Function to compute the derivative of (direct geometric) declination w.r.t. position of observer or observed object.
Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver )
{
    // Define multiplier of patial vector
    double partialMultiplier = ( ( isLinkEndReceiver ) ? 1.0 : -1.0 );

    // Set partial vector
    double range = relativeRangeVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = partialMultiplier * relativeRangeVector.transpose( ) / range;

    partial *= relativeRangeVector( 2 ) / ( range * range );
    partial += -partialMultiplier * ( Eigen::Vector3d::UnitZ( ) ).transpose( ) / range;

    return partial;
}

//! Function to compute the derivative of (direct geometric) right ascension and declination w.r.t. position of observer or
//! observed object.
Eigen::Matrix< double, 2, 3 > calculatePartialOfAngularPositionWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver )
{
    Eigen::Matrix< double, 2, 3 > angularPositionPartial;
    angularPositionPartial.block( 0, 0, 1, 3 ) = calculatePartialOfRightAscensionWrtLinkEndPosition(
                relativeRangeVector, isLinkEndReceiver );
    angularPositionPartial.block( 1, 0, 1, 3 ) = calculatePartialOfDeclinationWrtLinkEndPosition(
                relativeRangeVector, isLinkEndReceiver );
    return angularPositionPartial;
}

//! Update the scaling object to the current times and states
void AngularPositionScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                     const std::vector< double >& times,
                                     const observation_models::LinkEndType fixedLinkEnd,
                                     const Eigen::VectorXd currentObservation )
{
    Eigen::Vector3d relativeRangeVector = ( linkEndStates[ 1 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d normalizedRelativeRangeVector = relativeRangeVector.normalized( );

    // Compute common scaling factor
    scalingFactor_ = calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVector, true );

    // Compute scaling for receiver reference
    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 0 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ =
                scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }
    // Compute scaling for transmitter reference
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 1 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ =
                scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }

    currentLinkEndType_ = fixedLinkEnd;

}

//! Function to calculate the observation partial(s) at required time and state
AngularPositionPartial::AngularPositionPartialReturnType AngularPositionPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector2d& currentObservation )
{
    if( linkEndOfFixedTime != angularPositionScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error angular position partial and scaling are inconsistent" );
    }

    AngularPositionPartialReturnType returnPartial;

    // Iterate over all link ends
    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        if( positionPartialIterator_->first == observation_models::transmitter )
        {
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];
        }
        else if( positionPartialIterator_->first == observation_models::receiver )
        {
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];
        }

        // Scale position partials
        returnPartial.push_back(
                    std::make_pair(
                        angularPositionScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }

    // Add scaled light-time correcion partials.
    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctions_.size( ); i++ )
    {
        currentLinkTimeCorrectionPartial_ = lighTimeCorrectionPartialsFunctions_.at( i )( states, times );
        returnPartial.push_back(
                    std::make_pair( angularPositionScaler_->getLightTimePartialScalingFactor( ) *
                                    physical_constants::SPEED_OF_LIGHT * currentLinkTimeCorrectionPartial_.first,
                    currentLinkTimeCorrectionPartial_.second ) );
    }


    return returnPartial;
}

}

}
