/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
#define TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H

#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/basic_astro/physicalConstants.h"
#include "tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h"
#include "tudat/astro/observation_models/stellarAberrationCorrection.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Calculate the partials of position angle and separation distance w.r.t. the two (observed) unit lines of sight.
/*!
 * \param firstLineOfSight Unit direction from receiver to first transmitter, as used in the observable.
 * \param secondLineOfSight Unit direction from receiver to second transmitter, as used in the observable.
 * \param northPoleDirection Unit celestial north pole defining the position angle, in the same frame as the lines of sight.
 * \param calculatePositionAngle Whether the position-angle row is computed (it is left zero otherwise).
 * \return Partials of [position angle; separation] w.r.t. the first and second line of sight, respectively.
 */
std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > calculatePositionAngleAndSeparationPartialWrtLinesOfSight(
        const Eigen::Vector3d& firstLineOfSight,
        const Eigen::Vector3d& secondLineOfSight,
        const Eigen::Vector3d& northPoleDirection,
        const bool calculatePositionAngle = true );

//! Derived class for scaling three-dimensional position partials to position angle and/or separation distance partials
/*!
 * Derived class for scaling three-dimensional position partials to position angle and/or separation distance partials. The link
 * end states are those of the first transmitter, second transmitter and receiver (in this order). The north pole defining the
 * position angle, and the use of stellar aberration, are taken from the ancillary settings of the observation. If stellar
 * aberration is used, the observable also depends on the receiver velocity.
 */
template< int ObservationSize >
class PositionAngleAndSeparationScaling : public DirectPositionPartialScaling< ObservationSize >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param observableType Observable type: position_angle or separation_distance (ObservationSize 1), or
     * position_angle_and_separation (ObservationSize 2).
     * \param j2000ToGlobalFrameTransformation Rotation from ICRF/J2000 to the frame of the link-end states.
     */
    PositionAngleAndSeparationScaling( const observation_models::ObservableType observableType,
                                       const Eigen::Matrix3d& j2000ToGlobalFrameTransformation ):
        DirectPositionPartialScaling< ObservationSize >( observableType ),
        j2000ToGlobalFrameTransformation_( j2000ToGlobalFrameTransformation ),
        calculatePositionAngle_( observableType != observation_models::separation_distance ),
        firstComponent_( observableType == observation_models::separation_distance ? 1 : 0 )
    {
        if( observation_models::getObservableSize( observableType ) != ObservationSize )
        {
            throw std::runtime_error( "Error when creating position angle/separation distance partial scaling, observable " +
                                      observation_models::getObservableName( observableType ) + " is incompatible with size " +
                                      std::to_string( ObservationSize ) );
        }
    }

    //! Update the scaling object to the current times and states
    /*!
     *  Update the scaling object to the current times and states
     *  \param linkEndStates List of states at each link end during observation (first transmitter, second transmitter, receiver)
     *  \param times List of times at each link end during observation.
     *  \param fixedLinkEnd Link end at which observation time is defined, i.e. link end for which associated time
     *  is kept constant when computing observable.
     *  \param currentObservation Value of observation for which partial scaling is to be computed
     *  \param ancillarySettings Ancillary settings of the observation, defining north pole and use of stellar aberration
     */
    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation,
                 const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr ) override
    {
        if( fixedLinkEnd != observation_models::receiver )
        {
            throw std::runtime_error(
                    "Error when updating position angle/separation distance scaling object, fixed link end time different from receiver." );
        }

        Eigen::Vector3d northPoleDirection = Eigen::Vector3d::UnitZ( );
        if( calculatePositionAngle_ )
        {
            northPoleDirection = j2000ToGlobalFrameTransformation_ *
                    observation_models::getPositionAngleReferencePoleInJ2000( times.at( 2 ), ancillarySettings );
        }
        useStellarAberration_ = observation_models::getPositionAngleUseStellarAberration( ancillarySettings );
        const Eigen::Vector3d receiverVelocity = linkEndStates.at( 2 ).segment< 3 >( 3 );

        // Compute observed lines of sight, and their partials w.r.t. the relative positions and receiver velocity
        std::vector< Eigen::Vector3d > linesOfSight( 2 );
        std::vector< Eigen::Matrix3d > lineOfSightPartialsWrtRelativePosition( 2 );
        std::vector< Eigen::Matrix3d > lineOfSightPartialsWrtReceiverVelocity( 2, Eigen::Matrix3d::Zero( ) );
        for( unsigned int i = 0; i < 2; i++ )
        {
            const Eigen::Vector3d relativePosition = linkEndStates.at( i ).segment< 3 >( 0 ) - linkEndStates.at( 2 ).segment< 3 >( 0 );
            if( useStellarAberration_ )
            {
                linesOfSight[ i ] =
                        observation_models::calculateApparentDirectionWithStellarAberration( relativePosition, receiverVelocity );
                observation_models::calculateStellarAberrationPartials( relativePosition,
                                                                        receiverVelocity,
                                                                        lineOfSightPartialsWrtRelativePosition[ i ],
                                                                        lineOfSightPartialsWrtReceiverVelocity[ i ] );
            }
            else
            {
                linesOfSight[ i ] = relativePosition.normalized( );
                lineOfSightPartialsWrtRelativePosition[ i ] =
                        ( Eigen::Matrix3d::Identity( ) - linesOfSight[ i ] * linesOfSight[ i ].transpose( ) ) / relativePosition.norm( );
            }
        }

        observation_models::checkPositionAngleGeometry< double >(
                linesOfSight[ 0 ], linesOfSight[ 1 ], northPoleDirection, calculatePositionAngle_ );
        const std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > lineOfSightPartials =
                calculatePositionAngleAndSeparationPartialWrtLinesOfSight(
                        linesOfSight[ 0 ], linesOfSight[ 1 ], northPoleDirection.normalized( ), calculatePositionAngle_ );

        // Compute scaling factors, with and without light-time correction, for each transmitter
        velocityScalingFactor_ = Eigen::Matrix< double, ObservationSize, 3 >::Zero( );
        for( unsigned int i = 0; i < 2; i++ )
        {
            const Eigen::Matrix< double, 2, 3 >& currentLineOfSightPartial =
                    ( i == 0 ) ? lineOfSightPartials.first : lineOfSightPartials.second;
            const Eigen::Vector3d transmitterVelocity = linkEndStates.at( i ).segment< 3 >( 3 );
            const Eigen::Vector3d normalizedRelativePosition =
                    ( linkEndStates.at( i ).segment< 3 >( 0 ) - linkEndStates.at( 2 ).segment< 3 >( 0 ) ).normalized( );
            const double lightTimeDenominator = physical_constants::SPEED_OF_LIGHT + normalizedRelativePosition.dot( transmitterVelocity );

            fixedTimeScalingFactors_[ i ] = ( currentLineOfSightPartial * lineOfSightPartialsWrtRelativePosition[ i ] )
                                                    .block( firstComponent_, 0, ObservationSize, 3 );
            scalingFactors_[ i ] = fixedTimeScalingFactors_[ i ] *
                    ( Eigen::Matrix3d::Identity( ) - transmitterVelocity * normalizedRelativePosition.transpose( ) / lightTimeDenominator );
            lightTimeCorrectionScalingFactors_[ i ] = -fixedTimeScalingFactors_[ i ] * transmitterVelocity / lightTimeDenominator;
            velocityScalingFactor_ += ( currentLineOfSightPartial * lineOfSightPartialsWrtReceiverVelocity[ i ] )
                                              .block( firstComponent_, 0, ObservationSize, 3 );
        }

        currentLinkEndType_ = fixedLinkEnd;
    }

    //! Function to retrieve the position scaling factor (with light-time correction) for specific link end
    Eigen::Matrix< double, ObservationSize, 3 > getPositionScalingFactor( const observation_models::LinkEndType linkEndType ) override
    {
        return getLinkEndScalingFactor( scalingFactors_, linkEndType );
    }

    //! Function to retrieve the position scaling factor (without light-time correction) for specific link end
    Eigen::Matrix< double, ObservationSize, 3 > getFixedTimePositionScalingFactor(
            const observation_models::LinkEndType linkEndType ) override
    {
        return getLinkEndScalingFactor( fixedTimeScalingFactors_, linkEndType );
    }

    //! Function to retrieve the velocity scaling factor for specific link end (non-zero only for receiver, with stellar aberration)
    Eigen::Matrix< double, ObservationSize, 3 > getVelocityScalingFactor( const observation_models::LinkEndType linkEndType ) override
    {
        if( linkEndType == observation_models::receiver )
        {
            return velocityScalingFactor_;
        }
        return Eigen::Matrix< double, ObservationSize, 3 >::Zero( );
    }

    Eigen::Matrix< double, ObservationSize, 3 > getFixedTimeVelocityScalingFactor(
            const observation_models::LinkEndType linkEndType ) override
    {
        return getVelocityScalingFactor( linkEndType );
    }

    bool isVelocityScalingNonZero( const observation_models::LinkEndType linkEndType ) override
    {
        return useStellarAberration_ && linkEndType == observation_models::receiver;
    }

    Eigen::Matrix< double, ObservationSize, 1 > getLightTimePartialScalingFactor( ) override
    {
        throw std::runtime_error(
                "Error when getting position angle/separation distance light-time partial scaling, leg index must be provided." );
    }

    //! Function to retrieve the factor by which the light-time partials of the given leg (0: first transmitter, 1: second
    //! transmitter) are to be scaled
    Eigen::Matrix< double, ObservationSize, 1 > getLegLightTimePartialScalingFactor( const int legIndex ) override
    {
        return lightTimeCorrectionScalingFactors_.at( legIndex );
    }

    observation_models::LinkEndType getCurrentLinkEndType( ) override
    {
        return currentLinkEndType_;
    }

private:
    Eigen::Matrix< double, ObservationSize, 3 > getLinkEndScalingFactor(
            const std::vector< Eigen::Matrix< double, ObservationSize, 3 > >& transmitterScalingFactors,
            const observation_models::LinkEndType linkEndType )
    {
        switch( linkEndType )
        {
            case observation_models::transmitter:
                return transmitterScalingFactors.at( 0 );
            case observation_models::transmitter2:
                return transmitterScalingFactors.at( 1 );
            case observation_models::receiver:
                return -transmitterScalingFactors.at( 0 ) - transmitterScalingFactors.at( 1 );
            default:
                throw std::runtime_error(
                        "Error when getting position angle/separation distance scaling factor, incorrect link end type: " +
                        observation_models::getLinkEndTypeString( linkEndType ) );
        }
    }

    //! Rotation from ICRF/J2000 to the frame of the link-end states.
    Eigen::Matrix3d j2000ToGlobalFrameTransformation_;

    //! Whether the position angle is part of the observable (false for separation-distance-only observables).
    bool calculatePositionAngle_;

    //! Index in [position angle; separation] of the first component of the observable.
    int firstComponent_;

    //! Whether stellar aberration was used in the last call to update( )
    bool useStellarAberration_ = false;

    //! Scaling factors, per transmitter, with light-time correction
    std::vector< Eigen::Matrix< double, ObservationSize, 3 > > scalingFactors_ =
            std::vector< Eigen::Matrix< double, ObservationSize, 3 > >( 2 );

    //! Scaling factors, per transmitter, without light-time correction
    std::vector< Eigen::Matrix< double, ObservationSize, 3 > > fixedTimeScalingFactors_ =
            std::vector< Eigen::Matrix< double, ObservationSize, 3 > >( 2 );

    //! Light-time correction partial scaling factors, per transmitter
    std::vector< Eigen::Matrix< double, ObservationSize, 1 > > lightTimeCorrectionScalingFactors_ =
            std::vector< Eigen::Matrix< double, ObservationSize, 1 > >( 2 );

    //! Receiver velocity scaling factor (due to stellar aberration)
    Eigen::Matrix< double, ObservationSize, 3 > velocityScalingFactor_ = Eigen::Matrix< double, ObservationSize, 3 >::Zero( );

    //! Fixed link end for last computation of update() function.
    observation_models::LinkEndType currentLinkEndType_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
