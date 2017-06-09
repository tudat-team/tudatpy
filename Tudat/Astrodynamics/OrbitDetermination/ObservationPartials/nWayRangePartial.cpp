#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/nWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace observation_partials
{

void NWayRangeScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                               const std::vector< double >& times,
                               const observation_models::LinkEndType fixedLinkEnd,
                               const Eigen::VectorXd currentObservation )
{
    Eigen::Vector3d currentRangeVector;
    int currentIndex;

    std::vector< Eigen::Vector6d > singleLinkEndStates;
    singleLinkEndStates.resize( 2 );
    const std::vector< double > singleLinkTimes;

    int fixedLinkEndIndex = observation_models::getNWayLinkIndexFromLinkEndType(
                fixedLinkEnd, numberOfLinkEnds_ );
    observation_models::LinkEndType referenceLinkEnd;

    for( constituentRangeScalingIterator_ = constituentRangeScalings_.begin( );
         constituentRangeScalingIterator_ != constituentRangeScalings_.end( ); constituentRangeScalingIterator_++ )
    {
        currentIndex = constituentRangeScalingIterator_->first;
        singleLinkEndStates[ 0 ] = linkEndStates.at( 2 * currentIndex );
        singleLinkEndStates[ 1 ] = linkEndStates.at( 2 * currentIndex + 1 );

        if( fixedLinkEndIndex <= constituentRangeScalingIterator_->first )
        {
            referenceLinkEnd = observation_models::transmitter;
        }
        else
        {
            referenceLinkEnd = observation_models::receiver;
        }
        constituentRangeScalingIterator_->second->update(
                    singleLinkEndStates, singleLinkTimes, referenceLinkEnd );
    }

    for( unsigned int i = 0; i < linkEndStates.size( ) - 1; i += 2 )
    {
        currentRangeVector = ( linkEndStates[ i + 1 ] - linkEndStates[ i ] ).segment( 0, 3 ).normalized( );
        projectedRelativeVelocityRatios_[ i / 2 ] = ( 1.0 - currentRangeVector.dot( ( linkEndStates[ i ] ).segment( 3, 3 ) ) /
                                                      physical_constants::SPEED_OF_LIGHT ) / ( 1.0 - currentRangeVector.dot( ( linkEndStates[ i + 1 ] ).segment( 3, 3 ) ) /
                                                                                               physical_constants::SPEED_OF_LIGHT );
    }
}



NWayRangePartial::NWayRangePartialReturnType NWayRangePartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
{
    NWayRangePartialReturnType completePartialSet;
    int referenceStartLinkEndIndex = getNWayLinkIndexFromLinkEndType( linkEndOfFixedTime, numberOfLinkEnds_ );

    std::vector< Eigen::Vector6d > subLinkStates;
    subLinkStates.resize( 2 );
    std::vector< double > subLinkTimes;
    subLinkTimes.resize( 2 );
    observation_models::LinkEndType subLinkReference;

    double additionalLinksMultiplier;

    //int counter = 0;
    for( rangePartialIterator_ = rangePartialList_.begin( ); rangePartialIterator_ != rangePartialList_.end( );
         rangePartialIterator_++ )
    {
        additionalLinksMultiplier = 1.0;
        NWayRangePartialReturnType currentPartialSet;

        subLinkStates[ 0 ] = states[ 2 * rangePartialIterator_->first ];
        subLinkStates[ 1 ] = states[ 2 * rangePartialIterator_->first + 1 ];
        subLinkTimes[ 0 ] = times[ 2 * rangePartialIterator_->first ];
        subLinkTimes[ 1 ] = times[ 2 * rangePartialIterator_->first + 1 ];

        //if( ! ( rangePartialIterator_->second->getParameterIdentifier( ).first == estimatable_parameters::observation_bias ) )
        {
            if( rangePartialIterator_->first >= referenceStartLinkEndIndex )
            {
                subLinkReference = observation_models::transmitter;

                for( int i = rangePartialIterator_->first + 1; i < numberOfLinkEnds_ - 1; i++ )
                {
                    additionalLinksMultiplier += nWayRangeScaler_->getProjectedRelativeVelocityRatio( i ) - 1.0;
                }
            }
            else
            {
                subLinkReference = observation_models::receiver;
                for( int i = rangePartialIterator_->first; i > 0; i-- )
                {
                    additionalLinksMultiplier += 1.0 / nWayRangeScaler_->getProjectedRelativeVelocityRatio( i - 1 ) - 1.0;
                }
            }
        }

        currentPartialSet = rangePartialIterator_->second->calculatePartial( subLinkStates, subLinkTimes, subLinkReference );

        for( unsigned int i = 0; i < currentPartialSet.size( ); i++ )
        {
            currentPartialSet[ i ].first *= additionalLinksMultiplier;
        }

        completePartialSet.insert( completePartialSet.end( ), currentPartialSet.begin( ), currentPartialSet.end( ) );
    }

    return completePartialSet;
}

}

}


