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

#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

#include <Eigen/Core>

#include "tudat/astro/observation_models/positionAngleAndSeparationObservationModel.h"
#include "tudat/astro/observation_models/stellarAberrationCorrection.h"
#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/positionPartials.h"
#include "tudat/math/basic/linearAlgebra.h"

namespace tudat
{

namespace observation_partials
{

//! Partials of position angle and separation distance w.r.t. the angular positions of both transmitters and the receiver velocity.
struct PositionAngleAndSeparationGeometryPartials {
    //! Partial of [position angle; separation] w.r.t. [right ascension; declination] of the first transmitter, as seen by receiver.
    Eigen::Matrix2d partialWrtFirstAngularPosition;

    //! Partial of [position angle; separation] w.r.t. [right ascension; declination] of the second transmitter, as seen by receiver.
    Eigen::Matrix2d partialWrtSecondAngularPosition;

    //! Partial of [position angle; separation] w.r.t. the inertial receiver velocity (zero for astrometric directions).
    Eigen::Matrix< double, 2, 3 > partialWrtReceiverVelocity;
};

//! Calculate the partials of position angle and separation distance w.r.t. the two (observed) unit lines of sight.
/*!
 * \param firstLineOfSight Unit direction from receiver to first transmitter, as used in the observable.
 * \param secondLineOfSight Unit direction from receiver to second transmitter, as used in the observable.
 * \param northPoleDirection Celestial north pole defining the position angle, in the same frame as the lines of sight.
 * \param calculatePositionAngle Whether the position-angle row is computed (it is left zero otherwise).
 * \return Partials of [position angle; separation] w.r.t. the first and second line of sight, respectively.
 */
inline std::pair< Eigen::Matrix< double, 2, 3 >, Eigen::Matrix< double, 2, 3 > > calculatePositionAngleAndSeparationPartialWrtLinesOfSight(
        const Eigen::Vector3d& firstLineOfSight,
        const Eigen::Vector3d& secondLineOfSight,
        const Eigen::Vector3d& northPoleDirection,
        const bool calculatePositionAngle = true )
{
    Eigen::Matrix< double, 2, 3 > firstUnitVectorPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
    Eigen::Matrix< double, 2, 3 > secondUnitVectorPartial = Eigen::Matrix< double, 2, 3 >::Zero( );
    if( calculatePositionAngle )
    {
        const Eigen::Vector3d unnormalisedEastDirection = northPoleDirection.cross( firstLineOfSight );
        const double eastDirectionNorm = unnormalisedEastDirection.norm( );
        const Eigen::Vector3d eastDirection = unnormalisedEastDirection / eastDirectionNorm;
        const Eigen::Vector3d northDirection = firstLineOfSight.cross( eastDirection );

        const double eastProjection = secondLineOfSight.dot( eastDirection );
        const double northProjection = secondLineOfSight.dot( northDirection );
        const double positionAngleDenominator = eastProjection * eastProjection + northProjection * northProjection;

        const Eigen::Matrix3d eastDirectionPartial = ( Eigen::Matrix3d::Identity( ) - eastDirection * eastDirection.transpose( ) ) *
                linear_algebra::getCrossProductMatrix( northPoleDirection ) / eastDirectionNorm;
        const Eigen::Matrix3d northDirectionPartial = -linear_algebra::getCrossProductMatrix( eastDirection ) +
                linear_algebra::getCrossProductMatrix( firstLineOfSight ) * eastDirectionPartial;

        firstUnitVectorPartial.row( 0 ) = ( northProjection * secondLineOfSight.transpose( ) * eastDirectionPartial -
                                            eastProjection * secondLineOfSight.transpose( ) * northDirectionPartial ) /
                positionAngleDenominator;
        secondUnitVectorPartial.row( 0 ) =
                ( northProjection * eastDirection.transpose( ) - eastProjection * northDirection.transpose( ) ) / positionAngleDenominator;
    }

    const double cosineOfSeparation = firstLineOfSight.dot( secondLineOfSight );
    const double sineOfSeparation = firstLineOfSight.cross( secondLineOfSight ).norm( );
    firstUnitVectorPartial.row( 1 ) = ( cosineOfSeparation * firstLineOfSight - secondLineOfSight ).transpose( ) / sineOfSeparation;
    secondUnitVectorPartial.row( 1 ) = ( cosineOfSeparation * secondLineOfSight - firstLineOfSight ).transpose( ) / sineOfSeparation;

    return std::make_pair( firstUnitVectorPartial, secondUnitVectorPartial );
}

//! Calculate the partials of position angle and separation distance w.r.t. both angular positions and the receiver velocity.
/*!
 * The angular positions are the right ascension and declination (in the frame of the link-end states) of the astrometric
 * (light-time-corrected geometric) directions from the receiver to each transmitter. If stellar aberration is applied, the partials
 * include the derivative of the aberrated directions w.r.t. the astrometric directions, and w.r.t. the receiver velocity.
 * \param linkEndStates States of first transmitter, second transmitter and receiver (in this order), in the global frame.
 * \param northPoleDirection Celestial north pole defining the position angle, in the global frame.
 * \param calculatePositionAngle Whether the position-angle partials are computed (they are left zero otherwise).
 * \param applyStellarAberration Whether the observable uses directions aberrated with the receiver velocity.
 * \return Partials of [position angle; separation] w.r.t. both angular positions and the receiver velocity.
 */
inline PositionAngleAndSeparationGeometryPartials calculatePositionAngleAndSeparationPartialWrtAngularPositions(
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const Eigen::Vector3d& northPoleDirection,
        const bool calculatePositionAngle = true,
        const bool applyStellarAberration = false )
{
    if( linkEndStates.size( ) != 3 )
    {
        throw std::runtime_error( "Position-angle partials require states for two transmitters and one receiver." );
    }

    const Eigen::Vector3d firstRelativePosition = linkEndStates.at( 0 ).segment< 3 >( 0 ) - linkEndStates.at( 2 ).segment< 3 >( 0 );
    const Eigen::Vector3d secondRelativePosition = linkEndStates.at( 1 ).segment< 3 >( 0 ) - linkEndStates.at( 2 ).segment< 3 >( 0 );
    const Eigen::Vector3d receiverVelocity = linkEndStates.at( 2 ).segment< 3 >( 3 );

    const Eigen::Vector3d firstAstrometricLineOfSight = firstRelativePosition.normalized( );
    const Eigen::Vector3d secondAstrometricLineOfSight = secondRelativePosition.normalized( );

    // Retrieve directions as used in the observable, and their partials w.r.t. the astrometric directions and receiver velocity
    Eigen::Vector3d firstLineOfSight = firstAstrometricLineOfSight;
    Eigen::Vector3d secondLineOfSight = secondAstrometricLineOfSight;
    Eigen::Matrix3d firstLineOfSightPartialWrtAstrometricDirection = Eigen::Matrix3d::Identity( );
    Eigen::Matrix3d secondLineOfSightPartialWrtAstrometricDirection = Eigen::Matrix3d::Identity( );
    Eigen::Matrix3d firstLineOfSightPartialWrtReceiverVelocity = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d secondLineOfSightPartialWrtReceiverVelocity = Eigen::Matrix3d::Zero( );
    if( applyStellarAberration )
    {
        firstLineOfSight =
                observation_models::calculateApparentDirectionWithStellarAberration( firstAstrometricLineOfSight, receiverVelocity );
        secondLineOfSight =
                observation_models::calculateApparentDirectionWithStellarAberration( secondAstrometricLineOfSight, receiverVelocity );
        observation_models::calculateStellarAberrationPartials( firstAstrometricLineOfSight,
                                                                receiverVelocity,
                                                                firstLineOfSightPartialWrtAstrometricDirection,
                                                                firstLineOfSightPartialWrtReceiverVelocity );
        observation_models::calculateStellarAberrationPartials( secondAstrometricLineOfSight,
                                                                receiverVelocity,
                                                                secondLineOfSightPartialWrtAstrometricDirection,
                                                                secondLineOfSightPartialWrtReceiverVelocity );
    }

    observation_models::checkPositionAngleGeometry< double >(
            firstLineOfSight, secondLineOfSight, northPoleDirection, calculatePositionAngle );

    const auto lineOfSightPartials = calculatePositionAngleAndSeparationPartialWrtLinesOfSight(
            firstLineOfSight, secondLineOfSight, northPoleDirection.normalized( ), calculatePositionAngle );

    const double singularityTolerance = 100.0 * std::numeric_limits< double >::epsilon( );
    const auto getUnitVectorPartialWrtAngularPosition = [ singularityTolerance ]( const Eigen::Vector3d& lineOfSight ) {
        const double horizontalNorm = lineOfSight.head< 2 >( ).norm( );
        if( horizontalNorm <= singularityTolerance )
        {
            throw std::runtime_error( "Position-angle partials cannot be formed from right ascension at a pole of the state frame." );
        }

        Eigen::Matrix< double, 3, 2 > partial;
        partial.col( 0 ) << -lineOfSight.y( ), lineOfSight.x( ), 0.0;
        partial.col( 1 ) << -lineOfSight.x( ) * lineOfSight.z( ) / horizontalNorm, -lineOfSight.y( ) * lineOfSight.z( ) / horizontalNorm,
                horizontalNorm;
        return partial;
    };

    PositionAngleAndSeparationGeometryPartials geometryPartials;
    geometryPartials.partialWrtFirstAngularPosition = lineOfSightPartials.first * firstLineOfSightPartialWrtAstrometricDirection *
            getUnitVectorPartialWrtAngularPosition( firstAstrometricLineOfSight );
    geometryPartials.partialWrtSecondAngularPosition = lineOfSightPartials.second * secondLineOfSightPartialWrtAstrometricDirection *
            getUnitVectorPartialWrtAngularPosition( secondAstrometricLineOfSight );
    geometryPartials.partialWrtReceiverVelocity = lineOfSightPartials.first * firstLineOfSightPartialWrtReceiverVelocity +
            lineOfSightPartials.second * secondLineOfSightPartialWrtReceiverVelocity;
    return geometryPartials;
}

//! Scaling object that updates the independent angular-position scalers for both light-time legs.
class PositionAngleAndSeparationPartialScaling : public PositionPartialScaling
{
public:
    PositionAngleAndSeparationPartialScaling( const std::shared_ptr< PositionPartialScaling > firstPartialScaling,
                                              const std::shared_ptr< PositionPartialScaling > secondPartialScaling ):
        firstPartialScaling_( firstPartialScaling ), secondPartialScaling_( secondPartialScaling )
    {}

    void update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd ) override
    {
        if( fixedLinkEnd != observation_models::receiver )
        {
            throw std::runtime_error( "Position-angle and separation observations must be referenced at the receiver." );
        }
        if( linkEndStates.size( ) != 3 || times.size( ) != 3 )
        {
            throw std::runtime_error( "Position-angle partial scaling requires three link-end states and times." );
        }

        firstPartialScaling_->update( { linkEndStates.at( 0 ), linkEndStates.at( 2 ) },
                                      { times.at( 0 ), times.at( 2 ) },
                                      observation_models::receiver,
                                      Eigen::Vector2d::Constant( TUDAT_NAN ) );
        secondPartialScaling_->update( { linkEndStates.at( 1 ), linkEndStates.at( 2 ) },
                                       { times.at( 1 ), times.at( 2 ) },
                                       observation_models::receiver,
                                       Eigen::Vector2d::Constant( TUDAT_NAN ) );
    }

private:
    std::shared_ptr< PositionPartialScaling > firstPartialScaling_;
    std::shared_ptr< PositionPartialScaling > secondPartialScaling_;
};

//! Class to compute partials of position angle and/or separation distance observables.
/*!
 * The partials are computed from the partials of the angular positions (right ascension and declination) of both transmitters,
 * as seen from the receiver. If the observation ancillary settings select aberrated directions, the contribution of the stellar
 * aberration to the partials is included, both through the change in the astrometric directions and through the dependency on the
 * receiver velocity.
 */
template< int ObservationSize >
class PositionAngleAndSeparationPartial : public ObservationPartial< ObservationSize >
{
public:
    //! Constructor
    /*!
     * \param firstAngularPositionPartial Angular-position partial for the link from first transmitter to receiver (may be nullptr).
     * \param secondAngularPositionPartial Angular-position partial for the link from second transmitter to receiver (may be nullptr).
     * \param j2000ToGlobalFrameTransformation Rotation from ICRF/J2000 to the frame of the link-end states.
     * \param observableType Observable type: position_angle or separation_distance (ObservationSize 1), or
     * position_angle_and_separation (ObservationSize 2).
     * \param receiverStatePartial Partial of the receiver state w.r.t. the parameter, used for the receiver-velocity dependency of
     * stellar aberration (nullptr if the receiver state does not depend on the parameter).
     */
    PositionAngleAndSeparationPartial( const std::shared_ptr< ObservationPartial< 2 > > firstAngularPositionPartial,
                                       const std::shared_ptr< ObservationPartial< 2 > > secondAngularPositionPartial,
                                       const Eigen::Matrix3d& j2000ToGlobalFrameTransformation,
                                       const observation_models::ObservableType observableType,
                                       const std::shared_ptr< CartesianStatePartial > receiverStatePartial = nullptr ):
        ObservationPartial< ObservationSize >( getParameterIdentifier( firstAngularPositionPartial, secondAngularPositionPartial ) ),
        firstAngularPositionPartial_( firstAngularPositionPartial ), secondAngularPositionPartial_( secondAngularPositionPartial ),
        j2000ToGlobalFrameTransformation_( j2000ToGlobalFrameTransformation ),
        calculatePositionAngle_( observableType != observation_models::separation_distance ),
        firstComponent_( observableType == observation_models::separation_distance ? 1 : 0 ), receiverStatePartial_( receiverStatePartial )
    {
        if( ( ObservationSize == 2 && observableType != observation_models::position_angle_and_separation ) ||
            ( ObservationSize == 1 && observableType != observation_models::position_angle &&
              observableType != observation_models::separation_distance ) )
        {
            throw std::runtime_error( "Invalid observable type " + observation_models::getObservableName( observableType ) +
                                      " requested for a position-angle and separation partial of size " +
                                      std::to_string( ObservationSize ) + "." );
        }
    }

    std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, ObservationSize, 1 >& =
                    Eigen::Matrix< double, ObservationSize, 1 >::Constant( TUDAT_NAN ) ) override
    {
        if( linkEndOfFixedTime != observation_models::receiver )
        {
            throw std::runtime_error( "Position-angle and separation partials must be referenced at the receiver." );
        }
        if( states.size( ) != 3 || times.size( ) != 3 )
        {
            throw std::runtime_error( "Position-angle and separation partials require three link-end states and times." );
        }

        // Match the observation model's frame convention at the receiver epoch.
        // The epoch is held fixed with respect to the estimated parameters here.
        Eigen::Vector3d referencePoleDirection = Eigen::Vector3d::UnitZ( );
        if( calculatePositionAngle_ )
        {
            referencePoleDirection = j2000ToGlobalFrameTransformation_ *
                    observation_models::getPositionAngleReferencePoleInJ2000( times.at( 2 ), ancillarySettings );
        }
        const bool applyStellarAberration = observation_models::getPositionAngleDirectionType( ancillarySettings ) ==
                observation_models::aberrated_position_angle_direction;
        const PositionAngleAndSeparationGeometryPartials geometryPartials = calculatePositionAngleAndSeparationPartialWrtAngularPositions(
                states, referencePoleDirection, calculatePositionAngle_, applyStellarAberration );
        const Eigen::Matrix< double, ObservationSize, 2 > firstTransformation =
                geometryPartials.partialWrtFirstAngularPosition.block( firstComponent_, 0, ObservationSize, 2 );
        const Eigen::Matrix< double, ObservationSize, 2 > secondTransformation =
                geometryPartials.partialWrtSecondAngularPosition.block( firstComponent_, 0, ObservationSize, 2 );

        std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > > partials;
        appendTransformedPartials( partials,
                                   firstAngularPositionPartial_,
                                   firstTransformation,
                                   { states.at( 0 ), states.at( 2 ) },
                                   { times.at( 0 ), times.at( 2 ) },
                                   ancillarySettings );
        appendTransformedPartials( partials,
                                   secondAngularPositionPartial_,
                                   secondTransformation,
                                   { states.at( 1 ), states.at( 2 ) },
                                   { times.at( 1 ), times.at( 2 ) },
                                   ancillarySettings );

        // Add dependency of stellar aberration on receiver velocity (receiver time is fixed).
        if( applyStellarAberration && receiverStatePartial_ != nullptr )
        {
            const Eigen::Matrix< double, ObservationSize, 3 > receiverVelocityTransformation =
                    geometryPartials.partialWrtReceiverVelocity.block( firstComponent_, 0, ObservationSize, 3 );
            partials.push_back( std::make_pair(
                    receiverVelocityTransformation * receiverStatePartial_->calculatePartialOfVelocity( states.at( 2 ), times.at( 2 ) ),
                    times.at( 2 ) ) );
        }
        return partials;
    }

private:
    static estimatable_parameters::EstimatebleParameterIdentifier getParameterIdentifier(
            const std::shared_ptr< ObservationPartial< 2 > >& firstPartial,
            const std::shared_ptr< ObservationPartial< 2 > >& secondPartial )
    {
        if( firstPartial != nullptr )
        {
            return firstPartial->getParameterIdentifier( );
        }
        if( secondPartial != nullptr )
        {
            return secondPartial->getParameterIdentifier( );
        }
        throw std::runtime_error( "Cannot create a position-angle partial from two null angular-position partials." );
    }

    static void appendTransformedPartials(
            std::vector< std::pair< Eigen::Matrix< double, ObservationSize, Eigen::Dynamic >, double > >& outputPartials,
            const std::shared_ptr< ObservationPartial< 2 > >& angularPositionPartial,
            const Eigen::Matrix< double, ObservationSize, 2 >& transformation,
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const std::shared_ptr< observation_models::ObservationAncillarySimulationSettings >& ancillarySettings )
    {
        if( angularPositionPartial == nullptr )
        {
            return;
        }
        const auto angularPartials = angularPositionPartial->calculatePartial(
                states, times, observation_models::receiver, ancillarySettings, Eigen::Vector2d::Constant( TUDAT_NAN ) );
        for( const auto& angularPartial : angularPartials )
        {
            outputPartials.push_back( std::make_pair( transformation * angularPartial.first, angularPartial.second ) );
        }
    }

    std::shared_ptr< ObservationPartial< 2 > > firstAngularPositionPartial_;
    std::shared_ptr< ObservationPartial< 2 > > secondAngularPositionPartial_;
    Eigen::Matrix3d j2000ToGlobalFrameTransformation_;

    //! Whether the position angle is part of the observable (false for separation-distance-only observables).
    bool calculatePositionAngle_;

    //! Index in [position angle; separation] of the first component of the observable.
    int firstComponent_;

    std::shared_ptr< CartesianStatePartial > receiverStatePartial_;
};

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_POSITIONANGLEANDSEPARATIONPARTIAL_H
