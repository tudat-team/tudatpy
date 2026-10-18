/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include <limits>
#include <stdexcept>

#include "tudat/astro/basic_astro/physicalConstants.h"

#include "tudat/astro/relativity/relativisticLightTimeCorrection.h"
#include "tudat/astro/observation_models/corrections/firstOrderRelativisticCorrection.h"

namespace tudat
{

namespace observation_models
{

namespace
{

// Evaluate each moving gravitating body at its closest approach to the straight
// photon segment. The bounded iteration uses local relative position/velocity;
// normal planetary motion gives a convex distance function over a light-time leg.
Eigen::Vector3d bodyPositionAtClosestPhotonApproach(
        const std::function< Eigen::Vector6d( const double ) >& bodyStateFunction,
        const Eigen::Vector3d& transmitterPosition,
        const Eigen::Vector3d& receiverPosition,
        const double transmissionTime,
        const double receptionTime )
{
    const double duration = receptionTime - transmissionTime;
    if( !std::isfinite( transmissionTime ) || !std::isfinite( receptionTime ) ||
        !std::isfinite( duration ) || duration < 0.0 ||
        !transmitterPosition.allFinite( ) || !receiverPosition.allFinite( ) )
    {
        throw std::invalid_argument( "Invalid light-time leg for closest-approach body evaluation." );
    }
    const auto checkedBodyState = [ & ]( const double epoch ) -> Eigen::Vector6d {
        const Eigen::Vector6d state = bodyStateFunction( epoch );
        if( !std::isfinite( epoch ) || !state.allFinite( ) )
        {
            throw std::runtime_error( "Nonfinite body state in closest-approach light-time correction." );
        }
        return state;
    };
    if( duration == 0.0 )
    {
        return checkedBodyState( transmissionTime ).head< 3 >( );
    }
    const Eigen::Vector3d displacement = receiverPosition - transmitterPosition;
    if( !displacement.allFinite( ) )
    {
        throw std::invalid_argument( "Nonfinite photon displacement in closest-approach light-time correction." );
    }
    // A zero spatial segment has no photon direction. Retain a defined,
    // symmetric midpoint evaluation rather than dividing by a zero path length.
    if( displacement.isZero( 0.0 ) )
    {
        return checkedBodyState( transmissionTime + 0.5 * duration ).head< 3 >( );
    }
    const Eigen::Vector3d photonVelocity = displacement / duration;
    const double tolerance = std::max(
            1.0e-9, 2.0 * std::numeric_limits< double >::epsilon( ) *
                    std::max( std::abs( transmissionTime ), std::abs( receptionTime ) ) );
    double offset = 0.5 * duration;
    for( unsigned int iteration = 0; iteration < 20; ++iteration )
    {
        const double epoch = transmissionTime + offset;
        const Eigen::Vector6d bodyState = checkedBodyState( epoch );
        const Eigen::Vector3d relativePosition = transmitterPosition + offset * photonVelocity - bodyState.head< 3 >( );
        const Eigen::Vector3d relativeVelocity = photonVelocity - bodyState.tail< 3 >( );
        const double speedSquared = relativeVelocity.squaredNorm( );
        const double numerator = relativePosition.dot( relativeVelocity );
        if( !relativePosition.allFinite( ) || !relativeVelocity.allFinite( ) ||
            !std::isfinite( speedSquared ) || speedSquared <= 0.0 || !std::isfinite( numerator ) )
        {
            throw std::runtime_error( "Degenerate relative motion in closest-approach light-time correction." );
        }
        const double step = numerator / speedSquared;
        if( !std::isfinite( step ) )
        {
            throw std::runtime_error( "Nonfinite closest-approach light-time iteration." );
        }
        const double candidateOffset = offset - step;
        if( !std::isfinite( candidateOffset ) )
        {
            throw std::runtime_error( "Nonfinite closest-approach light-time epoch update." );
        }
        const double updatedOffset = std::clamp( candidateOffset, 0.0, duration );
        const double updatedEpoch = transmissionTime + updatedOffset;
        if( std::abs( updatedOffset - offset ) <= tolerance || updatedEpoch == epoch )
        {
            return checkedBodyState( updatedEpoch ).head< 3 >( );
        }
        offset = updatedOffset;
    }
    throw std::runtime_error( "Closest-approach light-time body epoch did not converge within 20 updates." );
}

}  // namespace

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

    // Iterate over all gravitating bodies independently of their link-end names.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        const Eigen::Vector3d perturbingBodyPosition = bodyPositionAtClosestPhotonApproach(
                perturbingBodyStateFunctions_.at( i ), transmitterState.head< 3 >( ), receiverState.head< 3 >( ),
                transmissionTime, receptionTime );
        // Calculate correction due to current body and add to total.
        currentLighTimeCorrectionComponents_.at( i ) = relativity::calculateFirstOrderLightTimeCorrectionFromCentralBody(
                perturbingBodyGravitationalParameterFunctions_.at( i )( ),
                transmitterState.segment( 0, 3 ),
                receiverState.segment( 0, 3 ),
                perturbingBodyPosition,
                ppnParameterGamma,
                useBending_ );
        currentTotalLightTimeCorrection_ += currentLighTimeCorrectionComponents_.at( i );
    }

    return currentTotalLightTimeCorrection_;
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
    // Retrieve ppn parameter gamma.
    double ppnParameterGamma = ppnParameterGammaFunction_( );

    // Initialize correction to zero.
    Eigen::Matrix< double, 3, 1 > currentTotalLightTimeCorrectionPartial_ = Eigen::Matrix< double, 3, 1 >::Zero( );

    // Iterate over all gravitating bodies.
    for( unsigned int i = 0; i < perturbingBodyStateFunctions_.size( ); i++ )
    {
        const Eigen::Vector3d perturbingBodyPosition = bodyPositionAtClosestPhotonApproach(
                perturbingBodyStateFunctions_.at( i ), transmitterState.head< 3 >( ), receiverState.head< 3 >( ),
                transmissionTime, receptionTime );
        // Preserve the existing fixed-body-epoch gradient approximation. This
        // does not differentiate the selected closest-approach epoch itself.

        // Calculate correction due to current body and add to total.
        currentTotalLightTimeCorrectionPartial_ += relativity::calculateFirstOrderCentralBodyLightTimeCorrectionGradient(
                perturbingBodyGravitationalParameterFunctions_.at( i )( ),
                transmitterState.segment( 0, 3 ),
                receiverState.segment( 0, 3 ),
                perturbingBodyPosition,
                ( linkEndAtWhichPartialIsEvaluated == receiver ),
                ppnParameterGamma,
                useBending_ );
    }

    return currentTotalLightTimeCorrectionPartial_;
}

}  // namespace observation_models

}  // namespace tudat
