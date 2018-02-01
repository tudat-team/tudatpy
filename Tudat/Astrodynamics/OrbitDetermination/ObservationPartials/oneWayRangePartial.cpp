/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"

namespace tudat
{

namespace observation_partials
{


//! Update the scaling object to the current times and states
void OneWayRangeScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                 const std::vector< double >& times,
                                 const observation_models::LinkEndType fixedLinkEnd,
                                 const Eigen::VectorXd currentObservation )
{
    // Compute Euclidean distance vector
    Eigen::Vector3d rangeVector = linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 0 ].segment( 0, 3 );
    Eigen::Matrix< double, 1, 3 > rangeVectorNormalized = rangeVector.transpose( ) / rangeVector.norm( );

    // Compute scaling for receiver reference
    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = 1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 0 ].segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
    }

    // Compute scaling for transmitter reference
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ =
                1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 1 ].segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
    }

    currentLinkEndType_ = fixedLinkEnd;
}

//! Function to calculate the observation partial(s) at required time and state
OneWayRangePartial::OneWayRangePartialReturnType OneWayRangePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    if( linkEndOfFixedTime != oneWayRangeScaler_->getCurrentLinkEndType( ) )
    {
        std::cout << linkEndOfFixedTime << " " << oneWayRangeScaler_->getCurrentLinkEndType( ) << std::endl;
        throw std::runtime_error( "Error one-way range partial and scaling are inconsistent" );
    }

    OneWayRangePartialReturnType returnPartial;

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
                        oneWayRangeScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) ), currentTime_ ) );
    }

    // Add scaled light-time correcion partials.
    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctions_.size( ); i++ )
    {
        returnPartial.push_back( lighTimeCorrectionPartialsFunctions_.at( i )( states, times ) );
        returnPartial[ returnPartial.size( ) - 1 ].first *=
                physical_constants::SPEED_OF_LIGHT * oneWayRangeScaler_->getLightTimePartialScalingFactor( );
    }

    return returnPartial;
}

}

}
