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

    Eigen::Vector3d r1 = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d r2 = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    // Unit vectors from receiver to each transmitter
    Eigen::Vector3d u1 = r1 / r1.norm( );
    Eigen::Vector3d u2 = r2 / r2.norm( );
    double r1norm = r1.norm( );
    double r2norm = r2.norm( );

    // ===== Position angle: direct vector formulation =====
    // Tangent plane basis at u1: e_perp = z_hat × u1 / ||z_hat × u1||, e_para = u1 × e_perp
    Eigen::Vector3d zHat;
    zHat << 0.0, 0.0, 1.0;
    Eigen::Vector3d ePerp = zHat.cross( u1 );
    double ePerpNorm = ePerp.norm( );
    ePerp = ePerp / ePerpNorm;
    Eigen::Vector3d ePara = u1.cross( ePerp );

    double y = u2.dot( ePerp );
    double x = u2.dot( ePara );
    double denomPA = x * x + y * y;

    // ∂θ/∂u2 = (x·e_perp - y·e_para) / denom
    Eigen::Vector3d dTheta_dU2 = ( x * ePerp - y * ePara ) / denomPA;

    // ∂θ/∂u1 via chain rule through e_perp and e_para
    Eigen::Matrix3d skewZ;
    skewZ << 0.0, -zHat( 2 ), zHat( 1 ), zHat( 2 ), 0.0, -zHat( 0 ), -zHat( 1 ), zHat( 0 ), 0.0;
    Eigen::Matrix3d projPerp = Eigen::Matrix3d::Identity( ) - ePerp * ePerp.transpose( );
    Eigen::Matrix3d dePerp_dU1 = -projPerp * skewZ / ePerpNorm;

    Eigen::Matrix3d skewU1;
    skewU1 << 0.0, -u1( 2 ), u1( 1 ), u1( 2 ), 0.0, -u1( 0 ), -u1( 1 ), u1( 0 ), 0.0;
    Eigen::Matrix3d skewEPerp;
    skewEPerp << 0.0, -ePerp( 2 ), ePerp( 1 ), ePerp( 2 ), 0.0, -ePerp( 0 ), -ePerp( 1 ), ePerp( 0 ), 0.0;
    Eigen::Matrix3d dePara_dU1 = skewEPerp + skewU1 * dePerp_dU1;

    // ∂θ/∂u1 = (x·∂y/∂u1 - y·∂x/∂u1) / denom  where ∂y/∂u1 and ∂x/∂u1 are vectors
    Eigen::Vector3d dY_dU1_vec = dePerp_dU1.transpose( ) * u2;
    Eigen::Vector3d dX_dU1_vec = dePara_dU1.transpose( ) * u2;
    Eigen::Vector3d dTheta_dU1 = ( x * dY_dU1_vec - y * dX_dU1_vec ) / denomPA;

    // ===== Separation: vector formulation =====
    double dotProduct = u1.dot( u2 );
    Eigen::Vector3d crossProduct = u1.cross( u2 );
    double crossNorm = crossProduct.norm( );
    double cosTheta = dotProduct;
    double sinTheta = crossNorm;

    Eigen::Vector3d crossHat = crossProduct / crossNorm;
    Eigen::Vector3d dSin_dU1 = -crossHat.cross( u2 );
    Eigen::Vector3d dSin_dU2 = crossHat.cross( u1 );
    Eigen::Vector3d dSep_dU1 = cosTheta * dSin_dU1 - sinTheta * u2;
    Eigen::Vector3d dSep_dU2 = cosTheta * dSin_dU2 - sinTheta * u1;

    // Chain rule: ∂/∂r = ∂/∂u · (I - u·uᵀ)/r
    Eigen::Matrix3d proj1 = Eigen::Matrix3d::Identity( ) - u1 * u1.transpose( );
    Eigen::Matrix3d proj2 = Eigen::Matrix3d::Identity( ) - u2 * u2.transpose( );

    Eigen::Matrix< double, 1, 3 > dPA_dr1 = dTheta_dU1.transpose( ) * proj1 / r1norm;
    Eigen::Matrix< double, 1, 3 > dPA_dr2 = dTheta_dU2.transpose( ) * proj2 / r2norm;
    Eigen::Matrix< double, 1, 3 > dSep_dr1 = dSep_dU1.transpose( ) * proj1 / r1norm;
    Eigen::Matrix< double, 1, 3 > dSep_dr2 = dSep_dU2.transpose( ) * proj2 / r2norm;

    // Assemble 2x3 matrices
    referenceScalingFactorFirstTransmitter_.row( 0 ) = dPA_dr1;
    referenceScalingFactorFirstTransmitter_.row( 1 ) = dSep_dr1;
    referenceScalingFactorSecondTransmitter_.row( 0 ) = dPA_dr2;
    referenceScalingFactorSecondTransmitter_.row( 1 ) = dSep_dr2;

    // Light-time correction scaling
    Eigen::Vector3d normRelRange1 = r1.normalized( );
    referenceLightTimeCorrectionScaling_ = referenceScalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normRelRange1 ) );

    currentLinkEndType_ = fixedLinkEnd;
}

}  // namespace observation_partials

}  // namespace tudat
