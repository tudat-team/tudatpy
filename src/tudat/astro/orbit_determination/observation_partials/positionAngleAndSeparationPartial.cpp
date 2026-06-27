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
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"

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
    Eigen::Vector3d relativeRange1 = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativeRange2 = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    // Compute RA/Dec for both transmitters
    double ra1 = 2.0 *
            std::atan( relativeRange1[ 1 ] /
                       ( std::sqrt( relativeRange1[ 0 ] * relativeRange1[ 0 ] + relativeRange1[ 1 ] * relativeRange1[ 1 ] ) +
                         relativeRange1[ 0 ] ) );
    double dec1 = mathematical_constants::PI / 2.0 - std::acos( relativeRange1[ 2 ] / relativeRange1.norm( ) );

    double ra2 = 2.0 *
            std::atan( relativeRange2[ 1 ] /
                       ( std::sqrt( relativeRange2[ 0 ] * relativeRange2[ 0 ] + relativeRange2[ 1 ] * relativeRange2[ 1 ] ) +
                         relativeRange2[ 0 ] ) );
    double dec2 = mathematical_constants::PI / 2.0 - std::acos( relativeRange2[ 2 ] / relativeRange2.norm( ) );

    double dRA = ra2 - ra1;
    double sinDRA = std::sin( dRA );
    double cosDRA = std::cos( dRA );
    double sinD1 = std::sin( dec1 );
    double cosD1 = std::cos( dec1 );
    double sinD2 = std::sin( dec2 );
    double cosD2 = std::cos( dec2 );

    double N = sinDRA * cosD2;
    double D = cosD1 * sinD2 - sinD1 * cosD2 * cosDRA;
    double denom = N * N + D * D;

    // dPA/d(alpha1, delta1, alpha2, delta2)
    double dPA_dRA1 = -( D * cosDRA * cosD2 + N * sinD1 * cosD2 * sinDRA ) / denom;
    double dPA_dDec1 = -( N * ( -sinD1 * sinD2 - cosD1 * cosD2 * cosDRA ) ) / denom;
    double dPA_dRA2 = ( D * cosDRA * cosD2 + N * sinD1 * cosD2 * sinDRA ) / denom;
    double dPA_dDec2 = ( D * ( -sinDRA * sinD2 ) + N * ( cosD1 * cosD2 + sinD1 * sinD2 * cosDRA ) ) / denom;

    // Get d(RA,Dec)/d(position) for each link end
    Eigen::Matrix< double, 1, 3 > dRA1_dr1 = calculatePartialOfRightAscensionWrtLinkEndPosition( relativeRange1, true );
    Eigen::Matrix< double, 1, 3 > dDec1_dr1 = calculatePartialOfDeclinationWrtLinkEndPosition( relativeRange1, true );
    Eigen::Matrix< double, 1, 3 > dRA2_dr2 = calculatePartialOfRightAscensionWrtLinkEndPosition( relativeRange2, true );
    Eigen::Matrix< double, 1, 3 > dDec2_dr2 = calculatePartialOfDeclinationWrtLinkEndPosition( relativeRange2, true );

    // Chain rule: dPA/dr = dPA/dRA * dRA/dr + dPA/dDec * dDec/dr
    referenceScalingFactorFirstTransmitter_ = dPA_dRA1 * dRA1_dr1 + dPA_dDec1 * dDec1_dr1;
    referenceScalingFactorSecondTransmitter_ = dPA_dRA2 * dRA2_dr2 + dPA_dDec2 * dDec2_dr2;

    // Light-time correction scaling
    Eigen::Vector3d normRelRange1 = relativeRange1.normalized( );
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

    Eigen::Vector3d relativeRange1 = ( linkEndStates[ 2 ] - linkEndStates[ 0 ] ).segment( 0, 3 );
    Eigen::Vector3d relativeRange2 = ( linkEndStates[ 2 ] - linkEndStates[ 1 ] ).segment( 0, 3 );

    double ra1 = 2.0 *
            std::atan( relativeRange1[ 1 ] /
                       ( std::sqrt( relativeRange1[ 0 ] * relativeRange1[ 0 ] + relativeRange1[ 1 ] * relativeRange1[ 1 ] ) +
                         relativeRange1[ 0 ] ) );
    double dec1 = mathematical_constants::PI / 2.0 - std::acos( relativeRange1[ 2 ] / relativeRange1.norm( ) );

    double ra2 = 2.0 *
            std::atan( relativeRange2[ 1 ] /
                       ( std::sqrt( relativeRange2[ 0 ] * relativeRange2[ 0 ] + relativeRange2[ 1 ] * relativeRange2[ 1 ] ) +
                         relativeRange2[ 0 ] ) );
    double dec2 = mathematical_constants::PI / 2.0 - std::acos( relativeRange2[ 2 ] / relativeRange2.norm( ) );

    double dRA = ra2 - ra1;
    double sinDRA = std::sin( dRA );
    double cosDRA = std::cos( dRA );
    double sinD1 = std::sin( dec1 );
    double cosD1 = std::cos( dec1 );
    double sinD2 = std::sin( dec2 );
    double cosD2 = std::cos( dec2 );

    // --- Position Angle partials ---
    double N = sinDRA * cosD2;
    double D = cosD1 * sinD2 - sinD1 * cosD2 * cosDRA;
    double denomPA = N * N + D * D;

    double dPA_dRA1 = -( D * cosDRA * cosD2 + N * sinD1 * cosD2 * sinDRA ) / denomPA;
    double dPA_dDec1 = -( N * ( -sinD1 * sinD2 - cosD1 * cosD2 * cosDRA ) ) / denomPA;
    double dPA_dRA2 = ( D * cosDRA * cosD2 + N * sinD1 * cosD2 * sinDRA ) / denomPA;
    double dPA_dDec2 = ( D * ( -sinDRA * sinD2 ) + N * ( cosD1 * cosD2 + sinD1 * sinD2 * cosDRA ) ) / denomPA;

    // --- Separation partials ---

    // Unit vectors from receiver to each transmitter
    Eigen::Vector3d u1 = relativeRange1 / relativeRange1.norm( );
    Eigen::Vector3d u2 = relativeRange2 / relativeRange2.norm( );

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

    double r1 = relativeRange1.norm( );
    double r2 = relativeRange2.norm( );

    Eigen::Matrix3d proj1 = Eigen::Matrix3d::Identity( ) - u1 * u1.transpose( );
    Eigen::Matrix3d proj2 = Eigen::Matrix3d::Identity( ) - u2 * u2.transpose( );

    Eigen::Matrix< double, 1, 3 > dSep_dr1 = dSep_dU1.transpose( ) * proj1 / r1;
    Eigen::Matrix< double, 1, 3 > dSep_dr2 = dSep_dU2.transpose( ) * proj2 / r2;

    // Get d(RA,Dec)/d(position) for PA chain rule
    Eigen::Matrix< double, 1, 3 > dRA1_dr1 = calculatePartialOfRightAscensionWrtLinkEndPosition( relativeRange1, true );
    Eigen::Matrix< double, 1, 3 > dDec1_dr1 = calculatePartialOfDeclinationWrtLinkEndPosition( relativeRange1, true );
    Eigen::Matrix< double, 1, 3 > dRA2_dr2 = calculatePartialOfRightAscensionWrtLinkEndPosition( relativeRange2, true );
    Eigen::Matrix< double, 1, 3 > dDec2_dr2 = calculatePartialOfDeclinationWrtLinkEndPosition( relativeRange2, true );

    // Chain rule for transmitter 1 (PA uses RA/Dec; Sep uses vector)
    Eigen::Matrix< double, 1, 3 > dPA_dr1 = dPA_dRA1 * dRA1_dr1 + dPA_dDec1 * dDec1_dr1;

    // Chain rule for transmitter 2
    Eigen::Matrix< double, 1, 3 > dPA_dr2 = dPA_dRA2 * dRA2_dr2 + dPA_dDec2 * dDec2_dr2;

    // Assemble 2x3 matrices
    referenceScalingFactorFirstTransmitter_.row( 0 ) = dPA_dr1;
    referenceScalingFactorFirstTransmitter_.row( 1 ) = dSep_dr1;
    referenceScalingFactorSecondTransmitter_.row( 0 ) = dPA_dr2;
    referenceScalingFactorSecondTransmitter_.row( 1 ) = dSep_dr2;

    // Light-time correction scaling
    Eigen::Vector3d normRelRange1 = relativeRange1.normalized( );
    referenceLightTimeCorrectionScaling_ = referenceScalingFactorFirstTransmitter_ * linkEndStates[ 0 ].segment( 3, 3 ) /
            ( physical_constants::SPEED_OF_LIGHT - linkEndStates[ 0 ].segment( 3, 3 ).dot( normRelRange1 ) );

    currentLinkEndType_ = fixedLinkEnd;
}

}  // namespace observation_partials

}  // namespace tudat
