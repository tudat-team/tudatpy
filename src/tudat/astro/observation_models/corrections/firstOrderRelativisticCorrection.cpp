/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/relativity/relativisticLightTimeCorrection.h"
#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"

namespace tudat
{

namespace observation_models
{

//! Function to calculate first order relativistic light time correction due to set of gravitating point masses.
double FirstOrderLightTimeCorrectionCalculator::calculateLightTimeCorrectionWithMultiLegLinkEndStates(
        const std::vector< Eigen::Vector6d >& linkEndsStates,
        const std::vector< double >& linkEndsTimes,
        const unsigned int currentMultiLegTransmitterIndex,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings )
{
    // Retrieve state and time of receiver and transmitter
    Eigen::Vector6d transmitterState, receiverState;
    double transmissionTime, receptionTime;
    getTransmissionReceptionTimesAndStates( linkEndsStates,
                                            linkEndsTimes,
                                            currentMultiLegTransmitterIndex,
                                            transmitterState,
                                            receiverState,
                                            transmissionTime,
                                            receptionTime );

    // Retrieve ppn parameter gamma.
    double ppnParameterGamma = ppnParameterGammaFunction_( );

    // Initialize correction to zero.
    currentTotalLightTimeCorrection_ = 0.0;

    double evaluationTime = TUDAT_NAN;
    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        evaluationTime = transmissionTime + lightTimeEvaluationContribution_.at( i ) * ( receptionTime - transmissionTime );
        // Calculate correction due to current body and add to total.
        currentLighTimeCorrectionComponents_.at( i ) = relativity::calculateFirstOrderLightTimeCorrectionFromCentralBody(
                perturbingBodyGravitationalParameterFunctions_.at( i )( ),
                transmitterState.segment( 0, 3 ),
                receiverState.segment( 0, 3 ),
                perturbingBodyStateFunctions_.at( i )( evaluationTime ).segment( 0, 3 ),
                ppnParameterGamma,
                useBending_ );
        currentTotalLightTimeCorrection_ += currentLighTimeCorrectionComponents_.at( i );
    }

    return currentTotalLightTimeCorrection_;
}

//! Partial derivative with respect to link-end time at fixed link-end states.
double FirstOrderLightTimeCorrectionCalculator::calculateLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime,
        const LinkEndType linkEndAtWhichPartialIsEvaluated,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings )
{
    double timePartial = 0.0;
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        timePartial += calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
                transmitterState, receiverState, transmissionTime, receptionTime, linkEndAtWhichPartialIsEvaluated, i );
    }
    return timePartial;
}

//! Partial derivative of a single body's contribution with respect to link-end time at fixed link-end states.
double FirstOrderLightTimeCorrectionCalculator::calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndTime(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime,
        const LinkEndType linkEndAtWhichPartialIsEvaluated,
        const unsigned int bodyIndex )
{
    if( linkEndAtWhichPartialIsEvaluated != transmitter && linkEndAtWhichPartialIsEvaluated != receiver )
    {
        throw std::runtime_error( "Error when computing relativistic light-time time partial: invalid link end." );
    }

    const double receptionTimeWeight = lightTimeEvaluationContribution_.at( bodyIndex );
    const double timeWeight = ( linkEndAtWhichPartialIsEvaluated == receiver ) ? receptionTimeWeight : 1.0 - receptionTimeWeight;
    if( timeWeight == 0.0 )
    {
        return 0.0;
    }

    const double ppnParameterGamma = ppnParameterGammaFunction_( );
    const double evaluationTime = transmissionTime + receptionTimeWeight * ( receptionTime - transmissionTime );
    const Eigen::Vector6d perturbingBodyState = perturbingBodyStateFunctions_.at( bodyIndex )( evaluationTime );
    const double gravitationalParameter = perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( );
    const Eigen::RowVector3d transmitterGradient =
            relativity::calculateFirstOrderCentralBodyLightTimeCorrectionGradient( gravitationalParameter,
                                                                                   transmitterState.segment( 0, 3 ),
                                                                                   receiverState.segment( 0, 3 ),
                                                                                   perturbingBodyState.segment( 0, 3 ),
                                                                                   false,
                                                                                   ppnParameterGamma,
                                                                                   useBending_ );
    const Eigen::RowVector3d receiverGradient =
            relativity::calculateFirstOrderCentralBodyLightTimeCorrectionGradient( gravitationalParameter,
                                                                                   transmitterState.segment( 0, 3 ),
                                                                                   receiverState.segment( 0, 3 ),
                                                                                   perturbingBodyState.segment( 0, 3 ),
                                                                                   true,
                                                                                   ppnParameterGamma,
                                                                                   useBending_ );

    // Translational invariance gives the body-position gradient as minus the sum of the link-end gradients.
    return -timeWeight * ( transmitterGradient + receiverGradient ).dot( perturbingBodyState.segment( 3, 3 ) );
}

//! Function to compute the partial derivative of the light-time correction w.r.t. link end position
Eigen::Matrix< double, 3, 1 > FirstOrderLightTimeCorrectionCalculator::calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime,
        const LinkEndType linkEndAtWhichPartialIsEvaluated,
        const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings )
{
    // Initialize correction to zero.
    Eigen::Matrix< double, 3, 1 > currentTotalLightTimeCorrectionPartial_ = Eigen::Matrix< double, 3, 1 >::Zero( );

    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        currentTotalLightTimeCorrectionPartial_ += calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                transmitterState, receiverState, transmissionTime, receptionTime, linkEndAtWhichPartialIsEvaluated, i );
    }

    return currentTotalLightTimeCorrectionPartial_;
}

//! Function to compute the partial derivative of a single body's light-time correction w.r.t. link end position
Eigen::Matrix< double, 3, 1 >
FirstOrderLightTimeCorrectionCalculator::calculateSingleBodyLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
        const Eigen::Vector6d& transmitterState,
        const Eigen::Vector6d& receiverState,
        const double transmissionTime,
        const double receptionTime,
        const LinkEndType linkEndAtWhichPartialIsEvaluated,
        const unsigned int bodyIndex )
{
    const double evaluationTime =
            transmissionTime + lightTimeEvaluationContribution_.at( bodyIndex ) * ( receptionTime - transmissionTime );
    return relativity::calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
                   perturbingBodyGravitationalParameterFunctions_.at( bodyIndex )( ),
                   transmitterState.segment( 0, 3 ),
                   receiverState.segment( 0, 3 ),
                   perturbingBodyStateFunctions_.at( bodyIndex )( evaluationTime ).segment( 0, 3 ),
                   ( linkEndAtWhichPartialIsEvaluated == receiver ),
                   ppnParameterGammaFunction_( ),
                   useBending_ )
            .transpose( );
}

}  // namespace observation_models

}  // namespace tudat
