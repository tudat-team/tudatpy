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

void PositionAngleScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                   const std::vector< double >& times,
                                   const observation_models::LinkEndType fixedLinkEnd,
                                   const Eigen::VectorXd currentObservation )
{
    if( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating position angle scaling object, fixed link end must be receiver." );
    }

    // States: [0]=transmitter, [1]=transmitter2, [2]=receiver
    Eigen::Vector3d r1 = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d r2 = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    // Unit vectors from receiver to each transmitter
    Eigen::Vector3d u1 = r1 / r1.norm( );
    Eigen::Vector3d u2 = r2 / r2.norm( );
    double r1norm = r1.norm( );
    double r2norm = r2.norm( );

    // Tangent plane basis at u1: e_perp = z_hat × u1 / ||z_hat × u1||, e_para = u1 × e_perp
    Eigen::Vector3d zHat;
    zHat << 0.0, 0.0, 1.0;
    Eigen::Vector3d ePerp = zHat.cross( u1 );
    double ePerpNorm = ePerp.norm( );
    ePerp = ePerp / ePerpNorm;
    Eigen::Vector3d ePara = u1.cross( ePerp );

    // Position angle: θ = atan2(u2 · e_perp, u2 · e_para)
    double y = u2.dot( ePerp );
    double x = u2.dot( ePara );
    double denom = x * x + y * y;

    // ∂θ/∂u2 = (x·e_perp - y·e_para) / denom
    Eigen::Vector3d dTheta_dU2 = ( x * ePerp - y * ePara ) / denom;

    // ∂θ/∂u1: e_perp depends on u1, e_para depends on u1
    // ∂e_perp/∂u1 = -(1/||z×u1||) * [u1]× * [z]× + (z×u1)/||z×u1||³ * ((z×u1)·[u1]×) ... simplified:
    // Let's compute numerically: ∂θ/∂u1 = -∂θ/∂u2 (from rotation invariance) plus correction from e_perp/e_para dependence
    // More precisely, using the chain rule through e_perp and e_para:
    // ∂e_perp/∂u1 = -1/||z×u1|| * (I - e_perp·e_perpᵀ) · [z]×   (projection of z-cross onto tangent plane)
    // ∂e_para/∂u1 = [u1]× · ∂e_perp/∂u1 + [e_perp]×
    // This gets complex. Let's use the direct approach:
    // θ = atan2(u2·(z×u1)/||z×u1||, u2·(u1×(z×u1)/||z×u1||))
    // A simpler approach: use the full chain rule with numerical derivatives as fallback

    // Direct computation using the formula:
    // ∂θ/∂r1 = ∂θ/∂u1 · ∂u1/∂r1  where ∂u1/∂r1 = (I - u1·u1ᵀ)/r1
    // ∂θ/∂r2 = ∂θ/∂u2 · ∂u2/∂r2  where ∂u2/∂r2 = (I - u2·u2ᵀ)/r2

    // For ∂θ/∂u1, note that θ depends on u1 through e_perp and e_para:
    // ∂θ/∂u1 = u2·∂e_perp/∂u1 · (x/denom) ... actually let's just compute it directly:
    // ∂(u2·e_perp)/∂u1 = u2ᵀ · ∂e_perp/∂u1
    // ∂(u2·e_para)/∂u1 = u2ᵀ · ∂e_para/∂u1

    // ∂e_perp/∂u1 = -1/||z×u1|| * (I - e_perp·e_perpᵀ) · [z]×
    Eigen::Matrix3d skewZ;
    skewZ << 0.0, -zHat( 2 ), zHat( 1 ), zHat( 2 ), 0.0, -zHat( 0 ), -zHat( 1 ), zHat( 0 ), 0.0;
    Eigen::Matrix3d projPerp = Eigen::Matrix3d::Identity( ) - ePerp * ePerp.transpose( );
    Eigen::Matrix3d dePerp_dU1 = -projPerp * skewZ / ePerpNorm;

    // ∂e_para/∂u1 = ∂(u1 × e_perp)/∂u1 = [e_perp]× + [u1]× · ∂e_perp/∂u1
    Eigen::Matrix3d skewU1;
    skewU1 << 0.0, -u1( 2 ), u1( 1 ), u1( 2 ), 0.0, -u1( 0 ), -u1( 1 ), u1( 0 ), 0.0;
    Eigen::Matrix3d skewEPerp;
    skewEPerp << 0.0, -ePerp( 2 ), ePerp( 1 ), ePerp( 2 ), 0.0, -ePerp( 0 ), -ePerp( 1 ), ePerp( 0 ), 0.0;
    Eigen::Matrix3d dePara_dU1 = skewEPerp + skewU1 * dePerp_dU1;

    // ∂θ/∂u1 = ∂(atan2(y,x))/∂u1 = (x·∂y/∂u1 - y·∂x/∂u1) / denom
    Eigen::Vector3d dY_dU1_vec = dePerp_dU1.transpose( ) * u2;
    Eigen::Vector3d dX_dU1_vec = dePara_dU1.transpose( ) * u2;
    Eigen::Vector3d dTheta_dU1 = ( x * dY_dU1_vec - y * dX_dU1_vec ) / denom;

    // Chain rule: ∂θ/∂r = ∂θ/∂u · ∂u/∂r = ∂θ/∂u · (I - u·uᵀ)/r
    Eigen::Matrix3d proj1 = Eigen::Matrix3d::Identity( ) - u1 * u1.transpose( );
    Eigen::Matrix3d proj2 = Eigen::Matrix3d::Identity( ) - u2 * u2.transpose( );

    referenceScalingFactorFirstTransmitter_ = dTheta_dU1.transpose( ) * proj1 / r1norm;
    referenceScalingFactorSecondTransmitter_ = dTheta_dU2.transpose( ) * proj2 / r2norm;

    // Light-time correction scaling
    Eigen::Vector3d normRelRange1 = r1.normalized( );
    referenceLightTimeCorrectionScaling_( 0, 0 ) =
            ( referenceScalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) ).value( ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normRelRange1 ) );

    currentLinkEndType_ = fixedLinkEnd;
}

void SeparationScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                                const std::vector< double >& times,
                                const observation_models::LinkEndType fixedLinkEnd,
                                const Eigen::VectorXd currentObservation )
{
    if( fixedLinkEnd != observation_models::receiver )
    {
        throw std::runtime_error( "Error when updating separation scaling object, fixed link end must be receiver." );
    }

    Eigen::Vector3d relativeRange1 = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativeRange2 = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    // Unit vectors from receiver to each transmitter
    Eigen::Vector3d u1 = relativeRange1 / relativeRange1.norm( );
    Eigen::Vector3d u2 = relativeRange2 / relativeRange2.norm( );

    double dotProduct = u1.dot( u2 );
    Eigen::Vector3d crossProduct = u1.cross( u2 );
    double crossNorm = crossProduct.norm( );

    double cosTheta = dotProduct;
    double sinTheta = crossNorm;

    // θ = atan2(sinθ, cosθ) → dθ = cosθ·dsinθ - sinθ·dcosθ  (since denom = sin²+cos²=1)
    // d(cosθ)/d(u1) = u2,  d(cosθ)/d(u2) = u1
    // d(sinθ)/d(u₁) = (u₁×u₂)̂ × (-u₂) = -(u₁×u₂)̂ × u₂
    // d(sinθ)/d(u₂) = (u₁×u₂)̂ × u₁
    Eigen::Vector3d crossHat = crossProduct / crossNorm;

    Eigen::Vector3d dSin_dU1 = -crossHat.cross( u2 );
    Eigen::Vector3d dSin_dU2 = crossHat.cross( u1 );

    // dθ/du₁ = cosθ·d(sinθ)/du₁ - sinθ·d(cosθ)/du₁
    Eigen::Vector3d dSep_dU1 = cosTheta * dSin_dU1 - sinTheta * u2;
    // dθ/du₂ = cosθ·d(sinθ)/du₂ - sinθ·d(cosθ)/du₂
    Eigen::Vector3d dSep_dU2 = cosTheta * dSin_dU2 - sinTheta * u1;

    // Chain rule: du/dr = (I - u·uᵀ) / r
    double r1 = relativeRange1.norm( );
    double r2 = relativeRange2.norm( );

    Eigen::Matrix3d proj1 = Eigen::Matrix3d::Identity( ) - u1 * u1.transpose( );
    Eigen::Matrix3d proj2 = Eigen::Matrix3d::Identity( ) - u2 * u2.transpose( );

    referenceScalingFactorFirstTransmitter_ = dSep_dU1.transpose( ) * proj1 / r1;
    referenceScalingFactorSecondTransmitter_ = dSep_dU2.transpose( ) * proj2 / r2;

    Eigen::Vector3d normRelRange1 = relativeRange1.normalized( );
    referenceLightTimeCorrectionScaling_( 0, 0 ) =
            ( referenceScalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) ).value( ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normRelRange1 ) );

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
