#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"

namespace tudat
{

namespace observation_partials
{



void OneWayRangeScaling::update( const std::vector< basic_mathematics::Vector6d >& linkEndStates,
                                 const std::vector< double >& times,
                                 const observation_models::LinkEndType fixedLinkEnd )
{
    Eigen::Vector3d rangeVector = linkEndStates[ 1 ].segment( 0, 3 ) - linkEndStates[ 0 ].segment( 0, 3 );
    Eigen::Matrix< double, 1, 3 > rangeVectorNormalized = rangeVector.transpose( ) / rangeVector.norm( );
    baseScalingFactor_ = rangeVectorNormalized;
    receiverReferenceScalingFactor_ =  rangeVectorNormalized /
            ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 0 ].segment( 3, 3 ) ) /
              physical_constants::SPEED_OF_LIGHT );
    transmitterReferenceScalingFactor_ =  rangeVectorNormalized /
            ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 1 ].segment( 3, 3 ) ) /
              physical_constants::SPEED_OF_LIGHT );
}

Eigen::Matrix< double, 1, 3 > OneWayRangeScaling::getScalingFactor(
        const observation_models::LinkEndType linkEndType, const observation_models::LinkEndType referenceTimeLinkEnd  )
{
    Eigen::Matrix< double, 1, 3 > scaling;

    if( referenceTimeLinkEnd == observation_models::transmitter )
    {
        scaling = transmitterReferenceScalingFactor_;
    }
    else if( referenceTimeLinkEnd == observation_models::receiver )
    {
        scaling = receiverReferenceScalingFactor_;
    }
    else
    {
        std::cerr<<"Error when getting one-way range scaling "<<linkEndType<<std::endl;
    }

    if( linkEndType == observation_models::receiver )
    {
        scaling *= -1.0;
    }

   // std::cout<<"1-way scal. "<<scaling<<std::endl;

    return scaling;
}

OneWayRangePartial::OneWayRangePartialReturnType OneWayRangePartial::calculatePartial(
        const std::vector< basic_mathematics::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime )
{
    OneWayRangePartialReturnType returnPartial;

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
                        oneWayRangeScaler_->getScalingFactor( positionPartialIterator_->first, linkEndOfFixedTime ) *
                        ( positionPartialIterator_->second->calculatePartial(
                              currentState, currentTime ) ), currentTime ) );
    }

    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctions_.size( ); i++ )
    {
        returnPartial.push_back( lighTimeCorrectionPartialsFunctions_.at( i )( states, times ) );
        returnPartial[ returnPartial.size( ) - 1 ].first *= -1.0 * physical_constants::SPEED_OF_LIGHT;
    }

    return returnPartial;
}

}

}
