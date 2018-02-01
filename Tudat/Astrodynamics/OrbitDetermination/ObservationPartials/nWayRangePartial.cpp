/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/nWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

//! Update the scaling object to the current times and states
void NWayRangeScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                               const std::vector< double >& times,
                               const observation_models::LinkEndType fixedLinkEnd,
                               const Eigen::VectorXd currentObservation )
{
    Eigen::Vector3d currentRangeVector;
    int currentIndex;

    // Define lists of link end states and tines to be used for eah one-way link.
    std::vector< Eigen::Vector6d > singleLinkEndStates;
    singleLinkEndStates.resize( 2 );
    std::vector< double > singleLinkTimes;
    singleLinkTimes.resize( 2 );

    // Find index in link ends for fixed (reference) link ends
    int fixedLinkEndIndex = observation_models::getNWayLinkIndexFromLinkEndType(
                fixedLinkEnd, numberOfLinkEnds_ );

    // Iterate over all constituent one-way range scalings
    observation_models::LinkEndType referenceLinkEnd;
    for( constituentRangeScalingIterator_ = constituentRangeScalings_.begin( );
         constituentRangeScalingIterator_ != constituentRangeScalings_.end( ); constituentRangeScalingIterator_++ )
    {
        // Set times, states and reference link end for current one-way range
        currentIndex = constituentRangeScalingIterator_->first;
        singleLinkEndStates[ 0 ] = linkEndStates.at( 2 * currentIndex );
        singleLinkEndStates[ 1 ] = linkEndStates.at( 2 * currentIndex + 1 );
        singleLinkTimes[ 0 ] = times.at( 2 * currentIndex );
        singleLinkTimes[ 1 ] = times.at( 2 * currentIndex + 1 );

        if( fixedLinkEndIndex <= constituentRangeScalingIterator_->first )
        {
            referenceLinkEnd = observation_models::transmitter;
        }
        else
        {
            referenceLinkEnd = observation_models::receiver;
        }

        // Update current one-way range scaling
        constituentRangeScalingIterator_->second->update(
                    singleLinkEndStates, singleLinkTimes, referenceLinkEnd );
    }

    // Compute scaling factors by which one-way range partias are to be multiplied before being included in the n-way range
    // partial
    for( unsigned int i = 0; i < linkEndStates.size( ) - 1; i += 2 )
    {
        currentRangeVector = ( linkEndStates[ i + 1 ] - linkEndStates[ i ] ).segment( 0, 3 ).normalized( );
        projectedRelativeVelocityRatios_[ i / 2 ] =
                ( 1.0 - currentRangeVector.dot( ( linkEndStates[ i ] ).segment( 3, 3 ) ) /
                  physical_constants::SPEED_OF_LIGHT ) /
                ( 1.0 - currentRangeVector.dot( ( linkEndStates[ i + 1 ] ).segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
    }
}

//! Function to calculate the observation partial(s) at required time and state
NWayRangePartial::NWayRangePartialReturnType NWayRangePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    NWayRangePartialReturnType completePartialSet;
    int referenceStartLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndOfFixedTime, numberOfLinkEnds_ );

    // Define states, times and reference link ends to be used for constutuent one-way range partials
    std::vector< Eigen::Vector6d > subLinkStates;
    subLinkStates.resize( 2 );
    std::vector< double > subLinkTimes;
    subLinkTimes.resize( 2 );
    observation_models::LinkEndType subLinkReference;

    double currentPartialMultiplier = TUDAT_NAN;

    for( rangePartialIterator_ = rangePartialList_.begin( ); rangePartialIterator_ != rangePartialList_.end( );
         rangePartialIterator_++ )
    {
        NWayRangePartialReturnType currentPartialSet;

        // Set link end times and states for current one-way range
        subLinkStates[ 0 ] = states[ 2 * rangePartialIterator_->first ];
        subLinkStates[ 1 ] = states[ 2 * rangePartialIterator_->first + 1 ];
        subLinkTimes[ 0 ] = times[ 2 * rangePartialIterator_->first ];
        subLinkTimes[ 1 ] = times[ 2 * rangePartialIterator_->first + 1 ];

        // Compute value by which one-way range should be scaled for inclusion into n-way range
        currentPartialMultiplier = 1.0;
        if( rangePartialIterator_->first >= referenceStartLinkEndIndex )
        {
            subLinkReference = observation_models::transmitter;

            for( int i = rangePartialIterator_->first + 1; i < numberOfLinkEnds_ - 1; i++ )
            {
                currentPartialMultiplier += nWayRangeScaler_->getProjectedRelativeVelocityRatio( i ) - 1.0;
            }
        }
        else
        {
            subLinkReference = observation_models::receiver;
            for( int i = rangePartialIterator_->first; i > 0; i-- )
            {
                currentPartialMultiplier += 1.0 / nWayRangeScaler_->getProjectedRelativeVelocityRatio( i - 1 ) - 1.0;
            }
        }

        // Compute one-way range partials
        currentPartialSet = rangePartialIterator_->second->calculatePartial( subLinkStates, subLinkTimes, subLinkReference );

        // Scale partials by required amount and add to return map.
        for( unsigned int i = 0; i < currentPartialSet.size( ); i++ )
        {
            currentPartialSet[ i ].first *= currentPartialMultiplier;
        }
        completePartialSet.insert( completePartialSet.end( ), currentPartialSet.begin( ), currentPartialSet.end( ) );
    }

    return completePartialSet;
}

}

}


