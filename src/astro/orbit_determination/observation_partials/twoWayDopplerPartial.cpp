/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/twoWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

//! Update the scaling object to the current times and states
void TwoWayDopplerScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                   const std::vector< double >& times,
                                   const observation_models::LinkEndType fixedLinkEnd,
                                   const Eigen::VectorXd currentObservation )
{
    Eigen::Vector3d currentRangeVector;

    // Define lists of link end states and tines to be used for eah one-way link.
    std::vector< Eigen::Vector6d > singleLinkEndStates;
    singleLinkEndStates.resize( 2 );
    std::vector< double > singleLinkTimes;
    singleLinkTimes.resize( 2 );

    double finiteDifferenceTimeStep = 60.0;
    double uplinkDoppler = TUDAT_NAN, downlinkDoppler = TUDAT_NAN;
    double upperturbedDoppler = TUDAT_NAN, downperturbedDoppler = TUDAT_NAN;
    double observationTime = TUDAT_NAN;
    // Find index in link ends for fixed (reference) link ends
    observation_models::LinkEndType referenceLinkEnd;
    int fixedLinkEndIndex = observation_models::getNWayLinkIndexFromLinkEndType(
                fixedLinkEnd, 3 );
    {
        singleLinkEndStates[ 0 ] = linkEndStates.at( 0 );
        singleLinkEndStates[ 1 ] = linkEndStates.at( 1 );
        singleLinkTimes[ 0 ] = times.at( 0 );
        singleLinkTimes[ 1 ] = times.at( 1 );

        if( fixedLinkEndIndex == 0 )
        {
            referenceLinkEnd = observation_models::transmitter;
            observationTime = times.at( 0 );
        }
        else
        {
            referenceLinkEnd = observation_models::receiver;
            observationTime = times.at( 1 );
        }

        upperturbedDoppler = uplinkDopplerModel_( observationTime + finiteDifferenceTimeStep, referenceLinkEnd );
        downperturbedDoppler = uplinkDopplerModel_( observationTime - finiteDifferenceTimeStep, referenceLinkEnd );
        uplinkOneWayDopplerTimeDerivative_ = ( upperturbedDoppler - downperturbedDoppler ) / ( 2.0 * finiteDifferenceTimeStep );
        uplinkDoppler = uplinkDopplerModel_( observationTime, referenceLinkEnd );

        // Update current one-way range scaling
        uplinkDopplerScaling_->update( singleLinkEndStates, singleLinkTimes, referenceLinkEnd );
        uplinkRangeScaling_->update( singleLinkEndStates, singleLinkTimes, referenceLinkEnd );
    }

    {
        singleLinkEndStates[ 0 ] = linkEndStates.at( 2 );
        singleLinkEndStates[ 1 ] = linkEndStates.at( 3 );
        singleLinkTimes[ 0 ] = times.at( 2 );
        singleLinkTimes[ 1 ] = times.at( 3 );

        if( fixedLinkEndIndex == 2 )
        {
            referenceLinkEnd = observation_models::receiver;
            observationTime = times.at( 3 );
        }
        else
        {
            referenceLinkEnd = observation_models::transmitter;
            observationTime = times.at( 2 );
        }

        upperturbedDoppler = downlinkDopplerModel_( observationTime + finiteDifferenceTimeStep, referenceLinkEnd );
        downperturbedDoppler = downlinkDopplerModel_( observationTime - finiteDifferenceTimeStep, referenceLinkEnd );
        downlinkOneWayDopplerTimeDerivative_ = ( upperturbedDoppler - downperturbedDoppler ) / ( 2.0 * finiteDifferenceTimeStep );
        downlinkDoppler = downlinkDopplerModel_( observationTime, referenceLinkEnd );

        // Update current one-way range scaling
        downlinkDopplerScaling_->update( singleLinkEndStates, singleLinkTimes, referenceLinkEnd );
        downlinkRangeScaling_->update( singleLinkEndStates, singleLinkTimes, referenceLinkEnd );
    }

    projectedRelativeVelocityRatios_[ 0 ] = downlinkDoppler + 1.0;
    projectedRelativeVelocityRatios_[ 1 ] = uplinkDoppler + 1.0;
}


//! Function to calculate the observation partial(s) at required time and state
TwoWayDopplerPartial::TwoWayDopplerPartialReturnType TwoWayDopplerPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings,
        const Eigen::Vector1d& currentObservation )
{
    TwoWayDopplerPartialReturnType completePartialSet;
    int referenceStartLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndOfFixedTime, numberOfLinkEnds_ );

    // Define states, times and reference link ends to be used for constutuent one-way range partials
    std::vector< Eigen::Vector6d > subLinkStates;
    subLinkStates.resize( 2 );
    std::vector< double > subLinkTimes;
    subLinkTimes.resize( 2 );
    observation_models::LinkEndType subLinkReference;

    double currentPartialMultiplier = TUDAT_NAN;
    double frequencyScaling = 1.0;
    if( scalingFactorFunction_ != nullptr )
    {
        frequencyScaling = scalingFactorFunction_( linkEndOfFixedTime, states, times, ancillarySettings );
    }

    for( dopplerPartialIterator_ = dopplerPartialList_.begin( ); dopplerPartialIterator_ != dopplerPartialList_.end( );
         dopplerPartialIterator_++ )
    {
        TwoWayDopplerPartialReturnType currentPartialSet;

        // Set link end times and states for current one-way range
        subLinkStates[ 0 ] = states[ 2 * dopplerPartialIterator_->first ];
        subLinkStates[ 1 ] = states[ 2 * dopplerPartialIterator_->first + 1 ];
        subLinkTimes[ 0 ] = times[ 2 * dopplerPartialIterator_->first ];
        subLinkTimes[ 1 ] = times[ 2 * dopplerPartialIterator_->first + 1 ];

        // Compute value by which one-way range should be scaled for inclusion into n-way range
        currentPartialMultiplier = frequencyScaling * twoWayDopplerScaler_->getProjectedRelativeVelocityRatio( dopplerPartialIterator_->first );

        if( dopplerPartialIterator_->first >= referenceStartLinkEndIndex )
        {
            subLinkReference = observation_models::transmitter;
        }
        else
        {
            subLinkReference = observation_models::receiver;
        }


        // Compute one-way range partials
        currentPartialSet = dopplerPartialIterator_->second->calculatePartial(
                subLinkStates, subLinkTimes, subLinkReference, ancillarySettings );

        // Scale partials by required amount and add to return map.
        for( unsigned int i = 0; i < currentPartialSet.size( ); i++ )
        {
            currentPartialSet[ i ].first *= currentPartialMultiplier;
        }
        completePartialSet.insert( completePartialSet.end( ), currentPartialSet.begin( ), currentPartialSet.end( ) );


        if( rangePartialList_.count( dopplerPartialIterator_->first ) > 0 )
        {
            if( ( linkEndOfFixedTime == observation_models::transmitter && dopplerPartialIterator_->first == 1 )||
                    ( linkEndOfFixedTime == observation_models::receiver && dopplerPartialIterator_->first == 0 ) )
            {
                currentPartialSet = rangePartialList_.at( dopplerPartialIterator_->first )->calculatePartial(
                            subLinkStates, subLinkTimes, subLinkReference, ancillarySettings );

                for( unsigned int i = 0; i < currentPartialSet.size( ); i++ )
                {
                    currentPartialSet[ i ].first *= ( currentPartialMultiplier ) *
                            twoWayDopplerScaler_->getRelevantOneWayDopplerTimePartial( linkEndOfFixedTime ) / physical_constants::SPEED_OF_LIGHT;
                }
                completePartialSet.insert( completePartialSet.end( ), currentPartialSet.begin( ), currentPartialSet.end( ) );
            }
        }
    }

    return completePartialSet;
}

}

}


