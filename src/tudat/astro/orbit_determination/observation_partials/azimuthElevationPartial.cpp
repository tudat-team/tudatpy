/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/azimuthElevationPartial.h"

#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"
#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{

namespace observation_partials
{

Eigen::Matrix2d calculatePartialOfAzimuthElevationWrtAngularPosition( const Eigen::Vector2d& angularPosition,
                                                                      const Eigen::Matrix3d& rotationFromInertialToTopocentricFrame,
                                                                      const bool invertLineOfSight )
{
    const double rightAscension = angularPosition.x( );
    const double declination = angularPosition.y( );
    const double cosineDeclination = std::cos( declination );
    const double lineOfSightSign = invertLineOfSight ? -1.0 : 1.0;

    Eigen::Vector3d inertialLineOfSight = ( Eigen::Vector3d( ) << cosineDeclination * std::cos( rightAscension ),
                                            cosineDeclination * std::sin( rightAscension ),
                                            std::sin( declination ) )
                                                  .finished( );
    Eigen::Vector3d topocentricRelativePosition = lineOfSightSign * rotationFromInertialToTopocentricFrame * inertialLineOfSight;

    const double east = topocentricRelativePosition.x( );
    const double north = topocentricRelativePosition.y( );
    const double up = topocentricRelativePosition.z( );
    const double horizontalRangeSquared = east * east + north * north;
    const double horizontalRange = std::sqrt( horizontalRangeSquared );
    const double rangeSquared = horizontalRangeSquared + up * up;

    // First map the topocentric line of sight to azimuth/elevation.
    Eigen::Matrix< double, 2, 3 > partialWrtTopocentricPosition = Eigen::Matrix< double, 2, 3 >::Zero( );
    partialWrtTopocentricPosition( 0, 0 ) = north / horizontalRangeSquared;
    partialWrtTopocentricPosition( 0, 1 ) = -east / horizontalRangeSquared;

    partialWrtTopocentricPosition( 1, 0 ) = -east * up / ( rangeSquared * horizontalRange );
    partialWrtTopocentricPosition( 1, 1 ) = -north * up / ( rangeSquared * horizontalRange );
    partialWrtTopocentricPosition( 1, 2 ) = horizontalRange / rangeSquared;

    // Then map right ascension/declination perturbations to the inertial unit line of sight.
    Eigen::Matrix< double, 3, 2 > partialOfUnitVectorWrtAngularPosition;
    partialOfUnitVectorWrtAngularPosition.col( 0 ) =
            ( Eigen::Vector3d( ) << -cosineDeclination * std::sin( rightAscension ), cosineDeclination * std::cos( rightAscension ), 0.0 )
                    .finished( );
    partialOfUnitVectorWrtAngularPosition.col( 1 ) = ( Eigen::Vector3d( ) << -std::sin( declination ) * std::cos( rightAscension ),
                                                       -std::sin( declination ) * std::sin( rightAscension ),
                                                       cosineDeclination )
                                                             .finished( );

    return partialWrtTopocentricPosition * lineOfSightSign * rotationFromInertialToTopocentricFrame * partialOfUnitVectorWrtAngularPosition;
}

Eigen::Vector2d calculatePartialOfAzimuthElevationWrtStationTime(
        const Eigen::Vector3d& relativeRangeVectorFromStationToTarget,
        const std::shared_ptr< ground_stations::PointingAnglesCalculator >& pointingAnglesCalculator,
        const double stationTime )
{
    const double timeStep = 1.0E-3;
    // The station frame is time-dependent through body rotation; evaluate this term numerically.
    std::pair< double, double > backwardPointingAngles =
            pointingAnglesCalculator->calculatePointingAngles( relativeRangeVectorFromStationToTarget, stationTime - timeStep );
    std::pair< double, double > forwardPointingAngles =
            pointingAnglesCalculator->calculatePointingAngles( relativeRangeVectorFromStationToTarget, stationTime + timeStep );

    double azimuthDifference = forwardPointingAngles.second - backwardPointingAngles.second;
    if( azimuthDifference > mathematical_constants::PI )
    {
        azimuthDifference -= 2.0 * mathematical_constants::PI;
    }
    else if( azimuthDifference < -mathematical_constants::PI )
    {
        azimuthDifference += 2.0 * mathematical_constants::PI;
    }

    return ( Eigen::Vector2d( ) << azimuthDifference, forwardPointingAngles.first - backwardPointingAngles.first ).finished( ) /
            ( 2.0 * timeStep );
}

void AzimuthElevationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                      const std::vector< double >& times,
                                      const observation_models::LinkEndType fixedLinkEnd,
                                      const Eigen::VectorXd currentObservation )
{
    const int stationIndex = ( stationLinkEndType_ == observation_models::transmitter ) ? 0 : 1;
    const int targetIndex = ( stationLinkEndType_ == observation_models::transmitter ) ? 1 : 0;

    // Reconstruct the inertial-to-topocentric rotation used by the pointing calculator at station time.
    Eigen::Matrix3d rotationFromInertialToTopocentricFrame;
    rotationFromInertialToTopocentricFrame.col( 0 ) =
            pointingAnglesCalculator_->convertVectorFromInertialToTopocentricFrame( Eigen::Vector3d::UnitX( ), times.at( stationIndex ) );
    rotationFromInertialToTopocentricFrame.col( 1 ) =
            pointingAnglesCalculator_->convertVectorFromInertialToTopocentricFrame( Eigen::Vector3d::UnitY( ), times.at( stationIndex ) );
    rotationFromInertialToTopocentricFrame.col( 2 ) =
            pointingAnglesCalculator_->convertVectorFromInertialToTopocentricFrame( Eigen::Vector3d::UnitZ( ), times.at( stationIndex ) );

    Eigen::Vector3d relativeRangeVector = ( linkEndStates[ 1 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d normalizedRelativeRangeVector = relativeRangeVector.normalized( );
    Eigen::Vector3d angularPositionUnitVector = -normalizedRelativeRangeVector;
    Eigen::Vector2d angularPosition = ( Eigen::Vector2d( ) << std::atan2( angularPositionUnitVector.y( ), angularPositionUnitVector.x( ) ),
                                        std::asin( angularPositionUnitVector.z( ) ) )
                                              .finished( );

    // Use the angular-position scaling architecture and only replace d(RA/DEC) by d(Az/El)/d(RA/DEC).
    Eigen::Matrix< double, 2, 3 > angularPositionScalingFactor =
            calculatePartialOfAngularPositionWrtLinkEndPosition( relativeRangeVector, true );
    scalingFactor_ =
            calculatePartialOfAzimuthElevationWrtAngularPosition(
                    angularPosition, rotationFromInertialToTopocentricFrame, stationLinkEndType_ == observation_models::transmitter ) *
            angularPositionScalingFactor;

    if( fixedLinkEnd == observation_models::receiver )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 0 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ = scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) +
                  linkEndStates[ 0 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                          ( physical_constants::SPEED_OF_LIGHT -
                            linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }
    else if( fixedLinkEnd == observation_models::transmitter )
    {
        referenceLightTimeCorrectionScaling_ = scalingFactor_ * linkEndStates[ 1 ].segment( 3, 3 ) /
                ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) );
        referenceScalingFactor_ = scalingFactor_ *
                ( Eigen::Matrix3d::Identity( ) +
                  linkEndStates[ 1 ].segment( 3, 3 ) * normalizedRelativeRangeVector.transpose( ) /
                          ( physical_constants::SPEED_OF_LIGHT -
                            linkEndStates[ 1 ].segment( 3, 3 ).dot( normalizedRelativeRangeVector ) ) );
    }

    if( fixedLinkEnd != stationLinkEndType_ )
    {
        // If the fixed-time link end is not the station, the light-time correction changes station time and frame orientation.
        const Eigen::Vector3d relativeRangeVectorFromStationToTarget =
                ( linkEndStates.at( targetIndex ) - linkEndStates.at( stationIndex ) ).segment( 0, 3 );
        const Eigen::Vector3d referenceVelocity = ( fixedLinkEnd == observation_models::receiver ) ? linkEndStates.at( 0 ).segment( 3, 3 )
                                                                                                   : linkEndStates.at( 1 ).segment( 3, 3 );
        const double lightTimeDenominator = physical_constants::SPEED_OF_LIGHT - referenceVelocity.dot( normalizedRelativeRangeVector );
        const double stationTimeSign = ( stationLinkEndType_ == observation_models::receiver ) ? 1.0 : -1.0;
        const Eigen::Vector2d stationTimeScaling = calculatePartialOfAzimuthElevationWrtStationTime(
                relativeRangeVectorFromStationToTarget, pointingAnglesCalculator_, times.at( stationIndex ) );

        referenceLightTimeCorrectionScaling_ += stationTimeSign * stationTimeScaling / lightTimeDenominator;
        referenceScalingFactor_ += stationTimeSign * stationTimeScaling * normalizedRelativeRangeVector.transpose( ) / lightTimeDenominator;
    }

    currentLinkEndType_ = fixedLinkEnd;
}

}  // namespace observation_partials

}  // namespace tudat
