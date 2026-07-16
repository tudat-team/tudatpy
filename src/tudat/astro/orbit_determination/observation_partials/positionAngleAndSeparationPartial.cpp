/*    Copyright (c) 2010-2024, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/observation_partials/positionAngleAndSeparationPartial.h"

namespace tudat
{

namespace observation_partials
{

PositionAngleScaling::PositionAngleScaling( ): DirectPositionPartialScaling< 1 >( observation_models::position_angle )
{
    psScaling_ = std::make_shared< PositionAngleAndSeparationScaling >( );
}

PositionAngleScaling::~PositionAngleScaling( ) {}

void PositionAngleScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                   const std::vector< double >& times,
                                   const observation_models::LinkEndType fixedLinkEnd,
                                   const Eigen::VectorXd currentObservation )
{
    // Delegate to the PS scaling
    psScaling_->update( linkEndStates, times, fixedLinkEnd, currentObservation );

    // Extract position angle row (first row) from the PS scaling
    referenceScalingFactorFirstTransmitter_ = psScaling_->getPositionScalingFactor( observation_models::transmitter ).row( 0 );
    referenceScalingFactorSecondTransmitter_ = psScaling_->getPositionScalingFactor( observation_models::transmitter2 ).row( 0 );
    referenceLightTimeCorrectionScaling_ = psScaling_->getLightTimePartialScalingFactor( ).segment( 0, 1 );
    currentLinkEndType_ = fixedLinkEnd;
}

SeparationScaling::SeparationScaling( ): DirectPositionPartialScaling< 1 >( observation_models::separation_distance )
{
    psScaling_ = std::make_shared< PositionAngleAndSeparationScaling >( );
}

SeparationScaling::~SeparationScaling( ) {}

void SeparationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                const std::vector< double >& times,
                                const observation_models::LinkEndType fixedLinkEnd,
                                const Eigen::VectorXd currentObservation )
{
    // Delegate to the PS scaling
    psScaling_->update( linkEndStates, times, fixedLinkEnd, currentObservation );

    // Extract separation row (second row) from the PS scaling
    referenceScalingFactorFirstTransmitter_ = psScaling_->getPositionScalingFactor( observation_models::transmitter ).row( 1 );
    referenceScalingFactorSecondTransmitter_ = psScaling_->getPositionScalingFactor( observation_models::transmitter2 ).row( 1 );
    referenceLightTimeCorrectionScaling_ = psScaling_->getLightTimePartialScalingFactor( ).segment( 1, 1 );
    currentLinkEndType_ = fixedLinkEnd;
}

void PositionAngleAndSeparationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                                const std::vector< double >& times,
                                                const observation_models::LinkEndType fixedLinkEnd,
                                                const Eigen::VectorXd currentObservation )
{
    if( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating position angle and separation scaling object, fixed link end must be receiver." );
    }

    Eigen::Vector3d relativePositionVector1 = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativePositionVector2 = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    // Unit vectors from receiver to each transmitter
    Eigen::Vector3d unitVectorToTransmitter1 = relativePositionVector1 / relativePositionVector1.norm( );
    Eigen::Vector3d unitVectorToTransmitter2 = relativePositionVector2 / relativePositionVector2.norm( );
    double normRelativePosition1 = relativePositionVector1.norm( );
    double normRelativePosition2 = relativePositionVector2.norm( );

    // ===== Position angle: direct vector formulation =====
    // Tangent plane basis at u1: eastDirection = northPoleDirection × u1 / ||northPoleDirection × u1||, northDirection = u1 × eastDirection
    Eigen::Vector3d northPoleDirection;
    northPoleDirection << 0.0, 0.0, 1.0;
    Eigen::Vector3d eastDirection = northPoleDirection.cross( unitVectorToTransmitter1 );
    double normEastDirection = eastDirection.norm( );
    eastDirection = eastDirection / normEastDirection;
    Eigen::Vector3d northDirection = unitVectorToTransmitter1.cross( eastDirection );

    double eastComponent = unitVectorToTransmitter2.dot( eastDirection );
    double northComponent = unitVectorToTransmitter2.dot( northDirection );
    double denominatorPositionAngle = northComponent * northComponent + eastComponent * eastComponent;

    // ∂θ/∂u2 = (northComponent·eastDirection - eastComponent·northDirection) / denominator
    Eigen::Vector3d dPositionAngle_dUnitVector2 =
            ( northComponent * eastDirection - eastComponent * northDirection ) / denominatorPositionAngle;

    // ∂θ/∂u1 via chain rule through eastDirection and northDirection
    Eigen::Matrix3d skewNorthPoleDirection;
    skewNorthPoleDirection << 0.0, -northPoleDirection( 2 ), northPoleDirection( 1 ), northPoleDirection( 2 ), 0.0,
            -northPoleDirection( 0 ), -northPoleDirection( 1 ), northPoleDirection( 0 ), 0.0;
    Eigen::Matrix3d projectionEastDirection = Eigen::Matrix3d::Identity( ) - eastDirection * eastDirection.transpose( );
    Eigen::Matrix3d dEastDirection_dUnitVector1 = projectionEastDirection * skewNorthPoleDirection / normEastDirection;

    Eigen::Matrix3d skewUnitVector1;
    skewUnitVector1 << 0.0, -unitVectorToTransmitter1( 2 ), unitVectorToTransmitter1( 1 ), unitVectorToTransmitter1( 2 ), 0.0,
            -unitVectorToTransmitter1( 0 ), -unitVectorToTransmitter1( 1 ), unitVectorToTransmitter1( 0 ), 0.0;
    Eigen::Matrix3d skewEastDirection;
    skewEastDirection << 0.0, -eastDirection( 2 ), eastDirection( 1 ), eastDirection( 2 ), 0.0, -eastDirection( 0 ), -eastDirection( 1 ),
            eastDirection( 0 ), 0.0;
    Eigen::Matrix3d dNorthDirection_dUnitVector1 = -skewEastDirection + skewUnitVector1 * dEastDirection_dUnitVector1;

    // ∂θ/∂u1 = (northComponent·∂eastComponent/∂u1 - eastComponent·∂northComponent/∂u1) / denominator  where ∂eastComponent/∂u1 and
    // ∂northComponent/∂u1 are vectors
    Eigen::Vector3d dEastComponent_dUnitVector1 = dEastDirection_dUnitVector1.transpose( ) * unitVectorToTransmitter2;
    Eigen::Vector3d dNorthComponent_dUnitVector1 = dNorthDirection_dUnitVector1.transpose( ) * unitVectorToTransmitter2;
    Eigen::Vector3d dPositionAngle_dUnitVector1 =
            ( northComponent * dEastComponent_dUnitVector1 - eastComponent * dNorthComponent_dUnitVector1 ) / denominatorPositionAngle;

    // ===== Separation: vector formulation =====
    double dotProductUnitVectors = unitVectorToTransmitter1.dot( unitVectorToTransmitter2 );
    Eigen::Vector3d crossProductUnitVectors = unitVectorToTransmitter1.cross( unitVectorToTransmitter2 );
    double normCrossProduct = crossProductUnitVectors.norm( );
    double cosSeparationAngle = dotProductUnitVectors;
    double sinSeparationAngle = normCrossProduct;

    Eigen::Vector3d unitCrossProduct = crossProductUnitVectors / normCrossProduct;
    Eigen::Vector3d dSinSeparation_dUnitVector1 = -unitCrossProduct.cross( unitVectorToTransmitter2 );
    Eigen::Vector3d dSinSeparation_dUnitVector2 = unitCrossProduct.cross( unitVectorToTransmitter1 );
    Eigen::Vector3d dSeparation_dUnitVector1 =
            cosSeparationAngle * dSinSeparation_dUnitVector1 - sinSeparationAngle * unitVectorToTransmitter2;
    Eigen::Vector3d dSeparation_dUnitVector2 =
            cosSeparationAngle * dSinSeparation_dUnitVector2 - sinSeparationAngle * unitVectorToTransmitter1;

    // Chain rule: ∂/∂r = ∂/∂u · (I - u·uᵀ)/r
    Eigen::Matrix3d projectionMatrix1 = Eigen::Matrix3d::Identity( ) - unitVectorToTransmitter1 * unitVectorToTransmitter1.transpose( );
    Eigen::Matrix3d projectionMatrix2 = Eigen::Matrix3d::Identity( ) - unitVectorToTransmitter2 * unitVectorToTransmitter2.transpose( );

    Eigen::Matrix< double, 1, 3 > dPositionAngle_dRelativePosition1 =
            dPositionAngle_dUnitVector1.transpose( ) * projectionMatrix1 / normRelativePosition1;
    Eigen::Matrix< double, 1, 3 > dPositionAngle_dRelativePosition2 =
            dPositionAngle_dUnitVector2.transpose( ) * projectionMatrix2 / normRelativePosition2;
    Eigen::Matrix< double, 1, 3 > dSeparation_dRelativePosition1 =
            dSeparation_dUnitVector1.transpose( ) * projectionMatrix1 / normRelativePosition1;
    Eigen::Matrix< double, 1, 3 > dSeparation_dRelativePosition2 =
            dSeparation_dUnitVector2.transpose( ) * projectionMatrix2 / normRelativePosition2;

    // Assemble 2x3 matrices
    referenceScalingFactorFirstTransmitter_.row( 0 ) = dPositionAngle_dRelativePosition1;
    referenceScalingFactorFirstTransmitter_.row( 1 ) = dSeparation_dRelativePosition1;
    referenceScalingFactorSecondTransmitter_.row( 0 ) = dPositionAngle_dRelativePosition2;
    referenceScalingFactorSecondTransmitter_.row( 1 ) = dSeparation_dRelativePosition2;

    // Light-time correction scaling
    Eigen::Vector3d normalizedRelativePosition1 = relativePositionVector1.normalized( );
    referenceLightTimeCorrectionScaling_ = referenceScalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normalizedRelativePosition1 ) );

    currentLinkEndType_ = fixedLinkEnd;
}

}  // namespace observation_partials

}  // namespace tudat
