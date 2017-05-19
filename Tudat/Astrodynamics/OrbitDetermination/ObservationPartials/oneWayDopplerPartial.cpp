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

Eigen::Vector3d computePartialOfUnitVectorWrtLinkEndTime(
        const Eigen::Vector3d& vectorToReceiver,
        const Eigen::Vector3d& unitVectorToReceiver,
        const double linkEndDistance,
        const Eigen::Vector3d linkEndVelocity )
{
    return linkEndVelocity / linkEndDistance - vectorToReceiver * unitVectorToReceiver.dot( linkEndVelocity ) /
            ( linkEndDistance * linkEndDistance );
}

double computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
        const Eigen::Vector3d& vectorToReceiver,
        const Eigen::Vector3d& projectedLinkEndVelocity,
        const Eigen::Vector3d& variableLinkEndVelocity,
        const Eigen::Vector3d& linkEndAcceleration,
        const bool linkEndIsReceiver )
{
     Eigen::Vector3d normalizedVector = vectorToReceiver.normalized( );
     double distance = vectorToReceiver.norm( );

     return static_cast< double >( linkEndIsReceiver ? 1.0 : -1.0 ) * computePartialOfUnitVectorWrtLinkEndTime(
                 vectorToReceiver, normalizedVector, distance, variableLinkEndVelocity ).dot( projectedLinkEndVelocity ) +
             normalizedVector.dot( linkEndAcceleration );
}

//! Update the scaling object to the current times and states
void OneWayDopplerScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                   const std::vector< double >& times,
                                   const observation_models::LinkEndType fixedLinkEnd,
                                   const Eigen::VectorXd currentObservation )
{
    double distance = ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ).norm( );
    Eigen::Vector3d lineOfSightVector = ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ).normalized( );
    Eigen::Vector3d receiverVelocity = linkEndStates.at( 1 ).segment( 3, 3 );
    Eigen::Vector3d transmitterVelocity = linkEndStates.at( 0 ).segment( 3, 3 );

    double lineOfSightVelocityReceiver = observation_models::calculateLineOfSightVelocityAsCFraction< double >(
                lineOfSightVector, receiverVelocity );
    double lineOfSightVelocityTransmitter = observation_models::calculateLineOfSightVelocityAsCFraction< double >(
                lineOfSightVector, transmitterVelocity );


    double transmitterPartialScalingTerm = 0.0;
    double receiverPartialScalingTerm = -1.0;

    double currentTaylorSeriesTerm = 1.0;

    for( int i = 1; i <= 3; i++ )
    {
        transmitterPartialScalingTerm += static_cast< double >( i ) * currentTaylorSeriesTerm;
        currentTaylorSeriesTerm *= lineOfSightVelocityTransmitter;
        receiverPartialScalingTerm -= currentTaylorSeriesTerm;
    }

    transmitterPartialScalingTerm *= ( 1.0 - lineOfSightVelocityReceiver );

    positionScalingFactor_ =
            ( receiverVelocity.transpose( ) * receiverPartialScalingTerm + transmitterVelocity.transpose( ) * transmitterPartialScalingTerm ) /
            physical_constants::SPEED_OF_LIGHT *
            ( Eigen::Matrix3d::Identity( ) - lineOfSightVector * lineOfSightVector.transpose( ) ) / distance;

    if( fixedLinkEnd == observation_models::receiver )
    {
        lightTimeEffectPositionScalingFactor_ =
                - transmitterPartialScalingTerm * lineOfSightVector.transpose( ) * computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                    ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
                    transmitterVelocity + receiverVelocity, transmitterVelocity,  transmitterAccelerationFunction_( times.at( 0 ) ), false )
                / ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );
//        lightTimeEffectPositionScalingFactor_ -=
//                receiverPartialScalingTerm * lineOfSightVector.transpose( ) * computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
//                    ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
//                    receiverVelocity, transmitterVelocity,  Eigen::Vector3d::Zero( ), false )
//                / ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );
    }
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        lightTimeEffectPositionScalingFactor_ =
                receiverPartialScalingTerm * lineOfSightVector.transpose( ) * computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                    ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
                    receiverVelocity + transmitterVelocity, receiverVelocity, receiverAccelerationFunction_( times.at( 1 ) ), true )
                / ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );
//        lightTimeEffectPositionScalingFactor_ +=
//                transmitterPartialScalingTerm * lineOfSightVector.transpose( ) * computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
//                    ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
//                    transmitterVelocity, receiverVelocity, Eigen::Vector3d::Zero( ), true )
//                / ( physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT );
    }


    positionScalingFactor_ += lightTimeEffectPositionScalingFactor_;
    receiverVelocityScalingFactor_ = -lineOfSightVector.transpose( ) * receiverPartialScalingTerm / physical_constants::SPEED_OF_LIGHT;
    transmitterVelocityScalingFactor_ = -lineOfSightVector.transpose( ) * transmitterPartialScalingTerm / physical_constants::SPEED_OF_LIGHT;

    currentLinkEndType_ = fixedLinkEnd;

}

//! Function to calculate the observation partial(s) at required time and state
OneWayDopplerPartial::OneWayDopplerPartialReturnType OneWayDopplerPartial::calculatePartial(
        const std::vector< Eigen::Vector6d >& states,
        const std::vector< double >& times,
        const observation_models::LinkEndType linkEndOfFixedTime,
        const Eigen::Vector1d& currentObservation )
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
                        ( positionPartialIterator_->second->calculatePartialOfPosition(
                              currentState_ , currentTime_ ) ) +
                        oneWayDopplerScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
                        ( positionPartialIterator_->second->calculatePartialOfVelocity(
                              currentState_ , currentTime_ ) ), currentTime_ ) );

//        std::cout<<"Position: "<<oneWayDopplerScaler_->getPositionScalingFactor( positionPartialIterator_->first ) *
//                   ( positionPartialIterator_->second->calculatePartialOfPosition(
//                         currentState_ , currentTime_ ) )<<std::endl;
//        std::cout<<"Velocity: "<<oneWayDopplerScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
//                   ( positionPartialIterator_->second->calculatePartialOfVelocity(
//                         currentState_ , currentTime_ ) )<<std::endl;
//        std::cout<<"Ratio: "<<
//                   ( oneWayDopplerScaler_->getPositionScalingFactor( positionPartialIterator_->first ) *
//                                      ( positionPartialIterator_->second->calculatePartialOfPosition(
//                                            currentState_ , currentTime_ ) ) ).cwiseQuotient
//                   ( oneWayDopplerScaler_->getVelocityScalingFactor( positionPartialIterator_->first ) *
//                   ( positionPartialIterator_->second->calculatePartialOfVelocity(
//                         currentState_ , currentTime_ ) ) )<<std::endl<<std::endl;
//        std::cout<<"Pos. part: "<<positionPartialIterator_->second->calculatePartialOfPosition(
//                       currentState_ , currentTime_ )<<std::endl;
//        std::cout<<"Vel. part: "<<positionPartialIterator_->second->calculatePartialOfVelocity(
//                       currentState_ , currentTime_ )<<std::endl;
    }

    return returnPartial;
}

}

}

