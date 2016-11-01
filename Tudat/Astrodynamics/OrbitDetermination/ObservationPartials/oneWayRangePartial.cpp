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

    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = 1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 0 ].segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
    }
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ =
                1.0 / ( 1.0 - rangeVectorNormalized.transpose( ).dot( linkEndStates[ 1 ].segment( 3, 3 ) ) /
                physical_constants::SPEED_OF_LIGHT );
        referenceScalingFactor_ =  rangeVectorNormalized * referenceLightTimeCorrectionScaling_;
    }
}

Eigen::Matrix< double, 1, 3 > OneWayRangeScaling::getScalingFactor(
        const observation_models::LinkEndType linkEndType )
{
    return referenceScalingFactor_ * ( ( linkEndType == observation_models::transmitter ) ? ( -1.0 ) : ( 1.0 ) );
}

double OneWayRangeScaling::getLightTimePartialScalingFactor( )
{
   return referenceLightTimeCorrectionScaling_;
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
                        oneWayRangeScaler_->getScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartial(
                              currentState, currentTime ) ), currentTime ) );
    }

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
