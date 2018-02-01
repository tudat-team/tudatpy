/*    Copyright (c) 2010-2018, Delft University of Technology
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

//! Function to computed the derivative of the unit vector from transmitter to receiver w.r.t. the observation time
Eigen::Vector3d computePartialOfUnitVectorWrtLinkEndTime(
        const Eigen::Vector3d& vectorToReceiver,
        const Eigen::Vector3d& unitVectorToReceiver,
        const double linkEndDistance,
        const Eigen::Vector3d linkEndVelocity )
{
    return linkEndVelocity / linkEndDistance - vectorToReceiver * unitVectorToReceiver.dot( linkEndVelocity ) /
            ( linkEndDistance * linkEndDistance );
}

//! Function to computed the derivative of velocity component along line-of-sight vector w.r.t. the observation time
double computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
        const Eigen::Vector3d& vectorToReceiver,
        const Eigen::Vector3d& projectedLinkEndVelocity,
        const Eigen::Vector3d& variableLinkEndVelocity,
        const Eigen::Vector3d& projectedLinkEndAcceleration,
        const bool linkEndIsReceiver,
        const bool projectedLinkEndIsVariableLinkEnd )
{
     Eigen::Vector3d normalizedVector = vectorToReceiver.normalized( );
     double distance = vectorToReceiver.norm( );

     return static_cast< double >( linkEndIsReceiver ? 1.0 : -1.0 ) * computePartialOfUnitVectorWrtLinkEndTime(
                 vectorToReceiver, normalizedVector, distance, variableLinkEndVelocity ).dot( projectedLinkEndVelocity ) +
             static_cast< double >( projectedLinkEndIsVariableLinkEnd ) * normalizedVector.dot( projectedLinkEndAcceleration );
}

//! Update the scaling object to the current times and states
void OneWayDopplerDirectFirstOrderProperTimeComponentScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
             const std::vector< double >& times,
             const observation_models::LinkEndType fixedLinkEnd,
             const Eigen::VectorXd currentObservation )
{
    // Get relative state
    Eigen::Vector6d relativeState = properTimeRateModel_->getComputationPointRelativeState(
                times, linkEndStates );
    currentDistance_ = relativeState.segment( 0, 3 ).norm( );
    currentGravitationalParameter_ = properTimeRateModel_->getGravitationalParameter( );

    currentLinkEndTime_ == ( linkEndWithPartial_ == observation_models::transmitter ) ? ( times.at( 0 ) ) : ( times.at( 1 ) );

    // Compute partials w.r.t. position and velocity
    if( computeStatePartials_ )
    {
        partialWrPosition_ = -physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                ( 1.0 + relativity::equivalencePrincipleLpiViolationParameter ) *
                currentGravitationalParameter_ / ( currentDistance_ * currentDistance_ ) *
                ( relativeState.segment( 0, 3 ).normalized( ) ).transpose( );
        partialWrtVelocity_ = physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                ( relativeState.segment( 3, 3 ) ).transpose( );
    }
}

//! Function to retrieve the scaling factor for the derivative w.r.t. the position of a given link end
Eigen::Matrix< double, 1, 3 > OneWayDopplerDirectFirstOrderProperTimeComponentScaling::getPositionScalingFactor(
        const observation_models::LinkEndType linkEndType )
{
    if( computeStatePartials_ )
    {
        return partialWrPosition_ * ( ( linkEndType == properTimeRateModel_->getComputationPointLinkEndType( ) ) ?
                                          ( 1.0 ) : ( 0.0 ) );
    }
    else
    {
        return Eigen::Matrix< double, 1, 3 >::Zero( );
    }
}


//! Function to retrieve the scaling factor for the derivative w.r.t. the velocity of a given link end
Eigen::Matrix< double, 1, 3 >
OneWayDopplerDirectFirstOrderProperTimeComponentScaling::getVelocityScalingFactor(
        const observation_models::LinkEndType linkEndType )
{
    if( computeStatePartials_ )
    {
        return partialWrtVelocity_ * ( ( linkEndType == properTimeRateModel_->getComputationPointLinkEndType( ) ) ?
                                           ( 1.0 ) : ( 0.0 ) );
    }
    else
    {
        return Eigen::Matrix< double, 1, 3 >::Zero( );
    }
}

//! Function to get the direct partial derivative, and associated time, of proper time
std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double >
OneWayDopplerDirectFirstOrderProperTimeComponentScaling::getProperTimeParameterPartial(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterType )
{
    if( parameterType.first == estimatable_parameters::equivalence_principle_lpi_violation_parameter )
    {
        return std::make_pair( ( Eigen::Matrix< double, 1, 1 >( ) << getEquivalencePrincipleViolationParameterPartial( ) ).finished( ),
                               currentLinkEndTime_ );
    }
    else
    {
        return std::make_pair( Eigen::Matrix< double, 1, 0 >::Zero( ), TUDAT_NAN );
    }
}

//! Function to get the size of the direct dependency of proper time rate on parameter
int OneWayDopplerDirectFirstOrderProperTimeComponentScaling::getParameterDependencySize(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterType )
{
    if( parameterType.first == estimatable_parameters::equivalence_principle_lpi_violation_parameter )
    {
        return 1;
    }
    else
    {
        return 0;
    }
}

//! Update the scaling object to the current times and states
void OneWayDopplerScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                   const std::vector< double >& times,
                                   const observation_models::LinkEndType fixedLinkEnd,
                                   const Eigen::VectorXd currentObservation )
{
    // Compute geometry.
    double distance = ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ).norm( );
    Eigen::Vector3d lineOfSightVector = ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ).normalized( );
    Eigen::Vector3d receiverVelocity = linkEndStates.at( 1 ).segment( 3, 3 );
    Eigen::Vector3d transmitterVelocity = linkEndStates.at( 0 ).segment( 3, 3 );

    double lineOfSightVelocityReceiver = observation_models::calculateLineOfSightVelocityAsCFraction< double >(
                lineOfSightVector, receiverVelocity );
    double lineOfSightVelocityTransmitter = observation_models::calculateLineOfSightVelocityAsCFraction< double >(
                lineOfSightVector, transmitterVelocity );

    // Compute scaling terms of transmitter/receiver
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

    // Compute position partial scaling term,
    positionScalingFactor_ =
            ( receiverVelocity.transpose( ) * receiverPartialScalingTerm +
              transmitterVelocity.transpose( ) * transmitterPartialScalingTerm ) /
            physical_constants::SPEED_OF_LIGHT *
            ( Eigen::Matrix3d::Identity( ) - lineOfSightVector * lineOfSightVector.transpose( ) ) / distance;

    if( fixedLinkEnd == observation_models::receiver )
    {
        lightTimeEffectPositionScalingFactor_ =
                -1.0 / ( ( physical_constants::SPEED_OF_LIGHT - lineOfSightVector.dot( transmitterVelocity ) ) *
                         physical_constants::SPEED_OF_LIGHT ) *
                ( transmitterPartialScalingTerm *
                  computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                      ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
                      transmitterVelocity, transmitterVelocity, transmitterAccelerationFunction_( times.at( 0 ) ), false, true ) +
                  receiverPartialScalingTerm *
                  computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                      ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
                      receiverVelocity, transmitterVelocity,  Eigen::Vector3d::Zero( ), false, false ) );
    }
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        lightTimeEffectPositionScalingFactor_ =
                1.0 / ( ( physical_constants::SPEED_OF_LIGHT - lineOfSightVector.dot( receiverVelocity ) ) *
                   physical_constants::SPEED_OF_LIGHT ) *
                ( receiverPartialScalingTerm *
                computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                    ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
                    receiverVelocity, receiverVelocity, receiverAccelerationFunction_( times.at( 1 ) ), true, true ) +
                transmitterPartialScalingTerm *
                computePartialOfProjectedLinkEndVelocityWrtAssociatedTime(
                    ( linkEndStates.at( 1 ) - linkEndStates.at( 0 ) ).segment( 0, 3 ),
                    transmitterVelocity, receiverVelocity, Eigen::Vector3d::Zero( ), true, true ) );
    }


    positionScalingFactor_ += lineOfSightVector.transpose( ) * lightTimeEffectPositionScalingFactor_;

    // Compute velocity scaling terms
    receiverVelocityScalingFactor_ = -lineOfSightVector.transpose( ) * receiverPartialScalingTerm /
            physical_constants::SPEED_OF_LIGHT;
    transmitterVelocityScalingFactor_ = -lineOfSightVector.transpose( ) * transmitterPartialScalingTerm /
            physical_constants::SPEED_OF_LIGHT;

    // Update proper time scaling objects.
    currentLinkEndType_ = fixedLinkEnd;

    if( transmitterProperTimePartials_ != NULL )
    {
        transmitterProperTimePartials_->update(
                    linkEndStates, times, fixedLinkEnd, currentObservation );
    }

    if( receiverProperTimePartials_ != NULL )
    {
        receiverProperTimePartials_->update(
                    linkEndStates, times, fixedLinkEnd, currentObservation );
    }

}

//! Function to retrieve the position scaling factor for specific link end
Eigen::Matrix< double, 1, 3 > OneWayDopplerScaling::getPositionScalingFactor( const observation_models::LinkEndType linkEndType )
{
    Eigen::Matrix< double, 1, 3 > scalingFactor =
            positionScalingFactor_ * ( ( linkEndType == observation_models::receiver ) ? ( 1.0 ) : ( -1.0 ) );

    if( transmitterProperTimePartials_ != NULL )
    {
        scalingFactor += transmitterProperTimePartials_->getPositionScalingFactor( linkEndType );
    }

    if( receiverProperTimePartials_ != NULL )
    {
        scalingFactor -= receiverProperTimePartials_->getPositionScalingFactor( linkEndType );
    }

    return scalingFactor;
}

//! Function to retrieve the velocity scaling factor for specific link end
Eigen::Matrix< double, 1, 3 > OneWayDopplerScaling::getVelocityScalingFactor( const observation_models::LinkEndType linkEndType )
{
    Eigen::Matrix< double, 1, 3 > scalingFactor =
            ( ( linkEndType == observation_models::receiver ) ? ( transmitterVelocityScalingFactor_ ) :
                                                                ( receiverVelocityScalingFactor_ ) );

    if( transmitterProperTimePartials_ != NULL )
    {
        scalingFactor += transmitterProperTimePartials_->getVelocityScalingFactor( linkEndType );
    }

    if( receiverProperTimePartials_ != NULL )
    {
        scalingFactor -= receiverProperTimePartials_->getVelocityScalingFactor( linkEndType );
    }

    return scalingFactor;
}

//! Function to get the size of the direct dependency of proper time rate on parameter
int OneWayDopplerScaling::getProperTimeParameterDependencySize(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterType )
{
    int transmitterDependencySize = 0;
    if( transmitterProperTimePartials_ != NULL )
    {
        transmitterDependencySize = transmitterProperTimePartials_->getParameterDependencySize( parameterType );
    }

    int receiverDependencySize = 0;
    if( receiverProperTimePartials_ != NULL )
    {
        receiverDependencySize = receiverProperTimePartials_->getParameterDependencySize( parameterType );
    }


    int totalDependencySize = 0;
    if( transmitterDependencySize == receiverDependencySize )
    {
        totalDependencySize = transmitterDependencySize;
    }
    else if( transmitterDependencySize == 0 )
    {
        totalDependencySize = receiverDependencySize;
    }
    else if( receiverDependencySize == 0 )
    {
        totalDependencySize = transmitterDependencySize;
    }
    else
    {
        throw std::runtime_error( "Error, proper time parameter dependency has inconsistent size" );
    }
    return totalDependencySize;
}

//! Function to get the direct partial derivatives, and associated times, of proper time components of Doppler partials
std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > OneWayDopplerScaling::getProperTimeParameterPartial(
        const estimatable_parameters::EstimatebleParameterIdentifier parameterType  )
{
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > transmitterPartial =
            transmitterProperTimePartials_->getProperTimeParameterPartial( parameterType );
    std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > receiverPartial =
            receiverProperTimePartials_->getProperTimeParameterPartial( parameterType );

    std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > totalPartial;
    if( transmitterPartial.first.rows( ) == receiverPartial.first.rows( ) )
    {
        totalPartial.push_back( transmitterPartial );
        totalPartial.push_back( receiverPartial );

        totalPartial[ totalPartial.size( ) - 1 ].first *= -1.0;
    }
    else if( transmitterPartial.first.rows( ) == 0 )
    {
        totalPartial.push_back( receiverPartial );
        totalPartial[ totalPartial.size( ) - 1 ].first *= -1.0;
    }
    else if( receiverPartial.first.rows( ) == 0 )
    {
        totalPartial.push_back( transmitterPartial );
    }
    else
    {
        throw std::runtime_error( "Error, proper time parameter partials have inconsistent size" );
    }
    return totalPartial;
}

//! Constructor
OneWayDopplerPartial::OneWayDopplerPartial(
        const boost::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler,
        const std::map< observation_models::LinkEndType, boost::shared_ptr< CartesianStatePartial > >& positionPartialList,
        const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lighTimeCorrectionPartials ):
    ObservationPartial< 1 >( parameterIdentifier ), oneWayDopplerScaler_( oneWayDopplerScaler ),
    positionPartialList_( positionPartialList )
{
    std::pair< boost::function< SingleOneWayDopplerPartialReturnType(
                const std::vector< Eigen::Vector6d >&, const std::vector< double >& ) >, bool > lightTimeCorrectionPartial;

    // Create light time correction partial functions
    for( unsigned int i = 0; i < lighTimeCorrectionPartials.size( ); i++ )
    {
        lightTimeCorrectionPartial = getLightTimeParameterPartialFunction(
                    parameterIdentifier_, lighTimeCorrectionPartials.at( i ) );
        if( lightTimeCorrectionPartial.second != 0 )
        {
            lighTimeCorrectionPartialsFunctions_.push_back( lightTimeCorrectionPartial.first );
        }
    }

    if( oneWayDopplerScaler->getProperTimeParameterDependencySize( parameterIdentifier ) > 0 )
    {
        addProperTimeParameterPartials_ = true;
    }
    else
    {
        addProperTimeParameterPartials_ = false;
    }
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
    }

    // Add scaled light-time correcion partials.
    for( unsigned int i = 0; i < lighTimeCorrectionPartialsFunctions_.size( ); i++ )
    {

        returnPartial.push_back( lighTimeCorrectionPartialsFunctions_.at( i )( states, times ) );
        returnPartial[ returnPartial.size( ) - 1 ].first *=
                physical_constants::SPEED_OF_LIGHT * oneWayDopplerScaler_->getLightTimeCorrectionPartialScaling( );
    }

    if( addProperTimeParameterPartials_ )
    {
        OneWayDopplerPartialReturnType properTimeReturnPartials = oneWayDopplerScaler_->getProperTimeParameterPartial(
                    parameterIdentifier_ );
        for( unsigned int i = 0; i < properTimeReturnPartials.size( ); i++ )
        {
            returnPartial.push_back( properTimeReturnPartials.at( i ) );
        }
    }

    return returnPartial;
}

}

}

