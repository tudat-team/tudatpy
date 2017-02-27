/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayDopplerPartial.h"

namespace tudat
{

namespace observation_partials
{


//! Update the scaling object to the current times and states
void OneWayDopplerScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                 const std::vector< double >& times,
                                 const observation_models::LinkEndType fixedLinkEnd )
{
    Eigen::Vector3d lineOfSightVector = ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ).normalized( );
    Eigen::Vector3d receiverVelocity = linkEndStates.at( 1 ).segment( 3, 3 );
    Eigen::Vector3d transmitterVelocity = linkEndStates.at( 0 ).segment( 3, 3 );

    double lineOfSightVelocityReceiver = observation_models::calculateLineOfSightVelocityAsCFraction< double >(
                lineOfSightVector, receiverVelocity );
    double lineOfSightVelocityTransmitter = observation_models::calculateLineOfSightVelocityAsCFraction< double >(
                lineOfSightVector, transmitterVelocity );

    double scalingTermA = -1.0;
    double scalingTermB = 0.0;

    double currentTaylorSeriesTerm = 1.0;

    for( int i = 1; i < 3; i++ )
    {
        scalingTermB += static_cast< double >( i ) * currentTaylorSeriesTerm;
        currentTaylorSeriesTerm *= lineOfSightVelocityReceiver;
        scalingTermA -= currentTaylorSeriesTerm;
    }

    scalingTermB *= ( 1.0 - lineOfSightVelocityTransmitter );

    positionScalingFactor_ = receiverVelocity * scalingTermB + transmitterVelocity * scalingTermA;
    receiverVelocityScalingFactor_ = lineOfSightVector * scalingTermB;
    transmitterVelocityScalingFactor_ = lineOfSightVector * scalingTermA;

}

//! Function to calculate the observation partial(s) at required time and state
OneWayDopplerPartial::OneWayDopplerPartialReturnType OneWayDopplerPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime )
{
    if( linkEndOfFixedTime != oneWayDopplerScaler_->getCurrentLinkEndType( ) )
    {
        throw std::runtime_error( "Error one-way doppler partial and scaling are inconsistent" );
    }

    OneWayDopplerPartialReturnType returnPartial;

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
                        oneWayDopplerScaler_->getPositionScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartial(
                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }

    for( velocityPartialIterator_ = velocityPartialList_.begin( ); velocityPartialIterator_ != velocityPartialList_.end( );
         velocityPartialIterator_++ )
    {
        if( velocityPartialIterator_->first == observation_models::transmitter )
        {
            currentState_  = states[ 0 ];
            currentTime_ = times[ 0 ];
        }
        else if( velocityPartialIterator_->first == observation_models::receiver )
        {
            currentState_  = states[ 1 ];
            currentTime_ = times[ 1 ];
        }

        // Scale position partials
        returnPartial.push_back(
                    std::make_pair(
                        oneWayDopplerScaler_->getVelocityScalingFactor( velocityPartialIterator_->first ) *
                        ( velocityPartialIterator_->second->calculatePartial(
                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }

    return returnPartial;
}

}

}

