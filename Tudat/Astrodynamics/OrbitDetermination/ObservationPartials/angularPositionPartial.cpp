#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/angularPositionPartial.h"

namespace tudat
{

namespace observation_partials
{

Eigen::Matrix< double, 1, 3 > calculatePartialOfRightAscensionWrtLinkEndPosition(
        const Eigen::Vector3d& relativeRangeVector,
        const bool isLinkEndReceiver )
{
    double partialMultiplier = TUDAT_NAN;
    if( isLinkEndReceiver )
    {
        partialMultiplier = 1.0;
    }
    else if( !isLinkEndReceiver )
    {
        partialMultiplier = -1.0;
    }

    Eigen::Matrix< double, 1, 3 > partial = Eigen::Matrix< double, 1, 3 >::Zero( );
    partial( 0 ) = -relativeRangeVector( 1 );
    partial( 1 ) = relativeRangeVector( 0 );
    partial /= ( partialMultiplier * ( relativeRangeVector( 0 ) * relativeRangeVector( 0 ) +
                                       relativeRangeVector( 1 ) * relativeRangeVector( 1 ) ) );
    return partial;
}

Eigen::Matrix< double, 1, 3 > calculatePartialOfDeclinationWrtLinkEndPosition(
        Eigen::Vector3d relativeRangeVector,
        const bool isLinkEndReceiver )
{
    double partialMultiplier = TUDAT_NAN;
    if( isLinkEndReceiver )
    {
        partialMultiplier = 1.0;
    }
    else if( !isLinkEndReceiver )
    {
        partialMultiplier = -1.0;
    }

    double range = relativeRangeVector.norm( );
    Eigen::Matrix< double, 1, 3 > partial = partialMultiplier * relativeRangeVector.transpose( ) / range;

    partial *= relativeRangeVector( 2 ) / ( range * range );
    partial += -partialMultiplier * ( Eigen::Vector3d::UnitZ( ) ).transpose( ) / range;

    return partial;
}

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

void AngularPositionScaling::update( const std::vector< basic_mathematics::Vector6d >& linkEndStates,
                                     const std::vector< double >& times,
                                     const observation_models::LinkEndType fixedLinkEnd )
{
    Eigen::Vector3d relativeRangeVector = ( linkEndStates[ 1 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d normalizedRelativeRangeVector = relativeRangeVector.normalized( );

    scalingFactor_ = calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVector, true );

    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 0 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ =
                scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 1 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ =
                scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) + linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }

}

Eigen::Matrix< double, 2, 3 > AngularPositionScaling::getScalingFactor(
        const observation_models::LinkEndType linkEndType )
{
    return referenceScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );

}

Eigen::Vector2d AngularPositionScaling::getLightTimePartialScalingFactor( )
{
    return referenceLightTimeCorrectionScaling_;
}

AngularPositionPartial::AngularPositionPartialReturnType AngularPositionPartial::calculatePartial(
        const std::vector< basic_mathematics::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime )
{
    AngularPositionPartialReturnType returnPartial;

    basic_mathematics::Vector6d currentState;
    double currentTime;

    for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
         positionPartialIterator_++ )
    {
        if( positionPartialIterator_->first == observation_models::transmitter )
        {
            currentState = states[ 0 ];
            currentTime = times[ 0 ];
        }
        else if( positionPartialIterator_->first == observation_models::receiver )
        {
            currentState = states[ 1 ];
            currentTime = times[ 1 ];
        }

        returnPartial.push_back(
                    std::make_pair(
                        angularPositionScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartial(
                              currentState, currentTime ) ), currentTime ) );
    }

    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctions_.size( ); i++ )
    {
        currentLinkTimeCorrectionPartial_ = lighTimeCorrectionPartialsFunctions_.at( i )( states, times );
        returnPartial.push_back(
                    std::make_pair( angularPositionScaler_->getLightTimePartialScalingFactor( ) * currentLinkTimeCorrectionPartial_.first,
                    currentLinkTimeCorrectionPartial_.second ) );
    }


    return returnPartial;
}

}

}
