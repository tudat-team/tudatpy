/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/fullTwoBodySphericalHarmonicGravityPartial.h"

#include "tudat/math/basic/coordinateConversions.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicPartialFunctions.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/sphericalHarmonicTransformations.h"
#include "tudat/math/basic/sphericalHarmonics.h"

namespace tudat
{

namespace acceleration_partials
{

namespace detail
{

std::array< Eigen::Matrix3d, 4 > getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternionForFullTwoBodyTorque(
        const Eigen::Quaterniond& rotationFromInertialToBodyFixedFrame );

Eigen::Matrix4d getLeftQuaternionMultiplicationMatrix( const Eigen::Vector4d& quaternion );

Eigen::Matrix4d getRightQuaternionMultiplicationMatrix( const Eigen::Vector4d& quaternion );

void computeDerivativeOfWignerDMatricesWrtRelativeQuaternion( const Eigen::Quaterniond& relativeRotationFromBody2ToBody1,
                                                              const std::shared_ptr< basic_mathematics::WignerDMatricesCache >& wignerCache,
                                                              std::array< std::vector< Eigen::MatrixXcd >, 4 >& derivatives );

Eigen::Matrix< double, 3, 2 > computeCurrentBodyFixedBasisFunctionGradients(
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache >& sphericalHarmonicsCache,
        const double cosineOfLatitude,
        const double preMultiplier,
        const double equatorialRadiusRatioPower,
        const int totalDegree,
        const int totalOrder );

}  // namespace detail

//! Constructor.
/*!
 * Initializes caches required for analytical derivatives of the full two-body acceleration model
 * (Dirkx et al. (2019), Eqs. (47)-(49), used in Eq. (55)).
 */
FullTwoBodySphericalHarmonicsGravityPartial::FullTwoBodySphericalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > accelerationModel ):
    AccelerationPartial( acceleratedBody,
                         acceleratingBody,
                         accelerationModel,
                         basic_astrodynamics::full_two_body_spherical_harmonic_gravity ),
    accelerationModel_( accelerationModel ), effectiveMutualPotentialField_( accelerationModel_->getEffectiveMutualPotentialField( ) ),
    sphericalHarmonicsCache_(
            std::make_shared< basic_mathematics::SphericalHarmonicsCache >( *accelerationModel_->getSphericalHarmonicsCache( ) ) ),
    coefficientCombinationsToUse_( effectiveMutualPotentialField_->getCoefficientCombinationsToUse( ) )
{
    // Cache setup supports derivatives of the Eq. (49) expansion used in translational Eq. (55),
    // with effective coefficients defined through Eqs. (47)-(48).
    sphericalHarmonicsCache_->getLegendreCache( ).setComputeSecondDerivatives( 1 );

    const int numberOfEffectiveCoefficients = effectiveMutualPotentialField_->getTotalVectorSize( );
    currentPartialsWrtEffectiveCoefficients_.resize( numberOfEffectiveCoefficients, Eigen::Matrix< double, 3, 2 >::Zero( ) );
    currentBodyFixedPartialsWrtEffectiveCoefficients_.resize( numberOfEffectiveCoefficients, Eigen::Matrix< double, 3, 2 >::Zero( ) );
    currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.resize( numberOfEffectiveCoefficients, Eigen::Matrix2d::Zero( ) );

    const Eigen::MatrixXd& transformedCosineCoefficients = effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& transformedSineCoefficients = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );
    body2BasisCosineCoefficients_.setZero( transformedCosineCoefficients.rows( ), transformedCosineCoefficients.cols( ) );
    body2BasisSineCoefficients_.setZero( transformedSineCoefficients.rows( ), transformedSineCoefficients.cols( ) );
    transformedCosineBody2CoefficientPartialsScratch_.setZero( transformedCosineCoefficients.rows( ),
                                                               transformedCosineCoefficients.cols( ) );
    transformedSineBody2CoefficientPartialsScratch_.setZero( transformedSineCoefficients.rows( ), transformedSineCoefficients.cols( ) );
    partialOfTransformedCosineCoefficientsBody2Scratch_.setZero( transformedCosineCoefficients.rows( ),
                                                                 transformedCosineCoefficients.cols( ) );
    partialOfTransformedSineCoefficientsBody2Scratch_.setZero( transformedSineCoefficients.rows( ), transformedSineCoefficients.cols( ) );
    currentPartialWrtQuaternionOfBody1_.setZero( );
    currentPartialWrtQuaternionOfBody2_.setZero( );
}

void FullTwoBodySphericalHarmonicsGravityPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        // Step 1: synchronize underlying acceleration model state (effective coefficients and geometry).
        accelerationModel_->updateMembers( currentTime );

        // Step 2: cache current rotations/positions and scalar factors used in Eq. (55) derivatives.
        currentRotationToBodyFixedFrame_ = accelerationModel_->getCurrentRotationFromInertialToBody1( ).toRotationMatrix( );
        currentRotationToInertialFrame_ = currentRotationToBodyFixedFrame_.transpose( );
        currentBodyFixedRelativePosition_ = accelerationModel_->getCurrentBodyFixedRelativePosition( );

        currentSphericalPosition_ = coordinate_conversions::convertCartesianToSpherical( currentBodyFixedRelativePosition_ );
        currentSphericalPosition_( 1 ) = mathematical_constants::PI / 2.0 - currentSphericalPosition_( 1 );
        currentDistance_ = currentBodyFixedRelativePosition_.norm( );
        currentGravitationalParameter_ = accelerationModel_->getCurrentGravitationalParameter( );

        currentRadius1Powers_ = accelerationModel_->getRadius1Powers( );
        currentRadius2Powers_ = accelerationModel_->getRadius2Powers( );

        currentBodyFixedSphericalGradient_ =
                coordinate_conversions::getSphericalToCartesianGradientMatrix( currentBodyFixedRelativePosition_ ).inverse( ) *
                accelerationModel_->getMutualPotentialGradient( );

        // Step 3: update analytical partial blocks w.r.t. position and effective coefficients.
        updateCurrentPositionPartial( );
        updateCurrentPartialsWrtEffectiveCoefficients( );

        // Step 4: update Jacobian of effective coefficients (Eqs. (47)-(48)) w.r.t transformed body-2 coefficients.
        effectiveMutualPotentialField_->computePartialsOfFullCoefficientsWrtTransformedCoefficients(
                currentEffectiveCoefficientsWrtTransformedBody2Coefficients_ );
        updateCurrentOrientationPartials( );

        currentPartialWrtVelocity_.setZero( );
        currentTime_ = currentTime;
    }
}

//! Update \partial a / \partial r in body-fixed and inertial frames.
/*!
 * Computes the Hessian of the effective potential summation (Eq. (49)) and maps it to Cartesian coordinates,
 * yielding the Jacobian entering derivatives of Eq. (55).
 */
void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentPositionPartial( )
{
    const double sineOfLatitude = std::sin( currentSphericalPosition_( 1 ) );
    const double cosineOfLatitude = std::cos( currentSphericalPosition_( 1 ) );
    const double preMultiplier = currentGravitationalParameter_ / currentDistance_;
    sphericalHarmonicsCache_->update( TUDAT_NAN, sineOfLatitude, currentSphericalPosition_( 2 ), TUDAT_NAN );

    Eigen::Matrix3d sphericalHessian = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d currentSphericalHessianContribution = Eigen::Matrix3d::Zero( );

    // Accumulate spherical Hessian contribution for each selected (l1,m1,l2,m2) term from Eq. (49).
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

        const int totalDegree = degreeOfBody1 + degreeOfBody2;
        const double equatorialRadiusRatioPower = currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );

        // Expand to signed-order variants consistent with real effective coefficients (Eqs. (47)-(49)).
        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            const bool computeTerm =
                    gravitation::getSignedOrdersForCombinationCase( j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 );
            if( computeTerm )
            {
                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const double effectiveCosineCoefficient = effectiveMutualPotentialField_->getEffectiveCosineCoefficient(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double effectiveSineCoefficient = effectiveMutualPotentialField_->getEffectiveSineCoefficient(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                computePotentialSphericalHessian(
                        currentDistance_,
                        equatorialRadiusRatioPower,
                        sphericalHarmonicsCache_->getCosineOfMultipleLongitude( totalOrder ),
                        sphericalHarmonicsCache_->getSineOfMultipleLongitude( totalOrder ),
                        cosineOfLatitude,
                        sineOfLatitude,
                        preMultiplier,
                        totalDegree,
                        totalOrder,
                        effectiveCosineCoefficient,
                        effectiveSineCoefficient,
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder ),
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomialDerivative( totalDegree, totalOrder ),
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomialSecondDerivative( totalDegree, totalOrder ),
                        currentSphericalHessianContribution );
                // Eq. (49): per-term contribution to the spherical Hessian of the mutual potential.
                sphericalHessian += currentSphericalHessianContribution;
            }
        }
    }

    const Eigen::Matrix3d sphericalToCartesianGradientMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( currentBodyFixedRelativePosition_ );
    currentBodyFixedPartialWrtPosition_ =
            sphericalToCartesianGradientMatrix * sphericalHessian * sphericalToCartesianGradientMatrix.transpose( );
    currentBodyFixedPartialWrtPosition_ += coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
            currentBodyFixedSphericalGradient_, currentBodyFixedRelativePosition_ );
    // Eq. (55): Cartesian Jacobian of acceleration from the gradient/Hessian of the mutual potential.
    currentPartialWrtPosition_ = currentRotationToInertialFrame_ * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;
}

//! Update \partial a / \partial C_eff and \partial a / \partial S_eff for all effective indices.
/*!
 * Evaluates basis-function gradients of Eq. (49) to obtain coefficient derivatives used in
 * chain-rule mappings to body-1/body-2 coefficient blocks.
 */
void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentPartialsWrtEffectiveCoefficients( )
{
    for( unsigned int i = 0; i < currentPartialsWrtEffectiveCoefficients_.size( ); i++ )
    {
        currentPartialsWrtEffectiveCoefficients_.at( i ).setZero( );
    }

    const double sineOfLatitude = std::sin( currentSphericalPosition_( 1 ) );
    const double cosineOfLatitude = std::cos( currentSphericalPosition_( 1 ) );
    const double preMultiplier = currentGravitationalParameter_ / currentDistance_;
    sphericalHarmonicsCache_->update( TUDAT_NAN, sineOfLatitude, currentSphericalPosition_( 2 ), TUDAT_NAN );

    // Loop over all selected interaction tuples and accumulate signed-order contributions.
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

        const int totalDegree = degreeOfBody1 + degreeOfBody2;
        const double equatorialRadiusRatioPower = currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );

        for( int j = 0; j < 4; j++ )
        {
            int signedOrderOfBody1 = 0;
            int signedOrderOfBody2 = 0;
            const bool computeTerm =
                    gravitation::getSignedOrdersForCombinationCase( j, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 );
            if( computeTerm )
            {
                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const int effectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double legendrePolynomial =
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
                const double legendrePolynomialDerivative =
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomialDerivative( totalDegree, totalOrder );

                // Eq. (49): basis-gradient contribution for one signed (l,m) term.
                currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ).block( 0, 0, 3, 1 ) +=
                        basic_mathematics::computePotentialGradient( currentDistance_,
                                                                     equatorialRadiusRatioPower,
                                                                     sphericalHarmonicsCache_->getCosineOfMultipleLongitude( totalOrder ),
                                                                     sphericalHarmonicsCache_->getSineOfMultipleLongitude( totalOrder ),
                                                                     cosineOfLatitude,
                                                                     preMultiplier,
                                                                     totalDegree,
                                                                     totalOrder,
                                                                     1.0,
                                                                     0.0,
                                                                     legendrePolynomial,
                                                                     legendrePolynomialDerivative );

                currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ).block( 0, 1, 3, 1 ) +=
                        basic_mathematics::computePotentialGradient( currentDistance_,
                                                                     equatorialRadiusRatioPower,
                                                                     sphericalHarmonicsCache_->getCosineOfMultipleLongitude( totalOrder ),
                                                                     sphericalHarmonicsCache_->getSineOfMultipleLongitude( totalOrder ),
                                                                     cosineOfLatitude,
                                                                     preMultiplier,
                                                                     totalDegree,
                                                                     totalOrder,
                                                                     0.0,
                                                                     1.0,
                                                                     legendrePolynomial,
                                                                     legendrePolynomialDerivative );
            }
        }
    }

    Eigen::Matrix< double, 3, 2 > currentBodyFixedPartialsWrtEffectiveCoefficients = Eigen::Matrix< double, 3, 2 >::Zero( );
    for( unsigned int i = 0; i < currentPartialsWrtEffectiveCoefficients_.size( ); i++ )
    {
        if( currentPartialsWrtEffectiveCoefficients_.at( i ).norm( ) > 0.0 )
        {
            currentBodyFixedPartialsWrtEffectiveCoefficients.block( 0, 0, 3, 1 ) =
                    coordinate_conversions::convertSphericalToCartesianGradient(
                            currentPartialsWrtEffectiveCoefficients_.at( i ).block( 0, 0, 3, 1 ), currentBodyFixedRelativePosition_ );
            currentBodyFixedPartialsWrtEffectiveCoefficients.block( 0, 1, 3, 1 ) =
                    coordinate_conversions::convertSphericalToCartesianGradient(
                            currentPartialsWrtEffectiveCoefficients_.at( i ).block( 0, 1, 3, 1 ), currentBodyFixedRelativePosition_ );
            currentBodyFixedPartialsWrtEffectiveCoefficients_.at( i ) = currentBodyFixedPartialsWrtEffectiveCoefficients;
            currentPartialsWrtEffectiveCoefficients_.at( i ) =
                    currentRotationToInertialFrame_ * currentBodyFixedPartialsWrtEffectiveCoefficients;
        }
        else
        {
            currentBodyFixedPartialsWrtEffectiveCoefficients_.at( i ).setZero( );
        }
    }
}

//! Update \partial a / \partial q for both body orientation quaternions.
/*!
 * Implements the quaternion chain rule in Dirkx et al. (2019), Eqs. (70)-(79): body-1 orientation affects
 * the output-frame rotation, the body-fixed relative position, and the relative coefficient rotation, while
 * body-2 orientation affects only the relative coefficient rotation.
 */
void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentOrientationPartials( )
{
    const std::array< Eigen::Matrix3d, 4 > derivativeOfBody1RotationFromBodyFixedToInertialWrtQuaternion =
            detail::getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternionForFullTwoBodyTorque(
                    accelerationModel_->getCurrentRotationFromInertialToBody1( ) );

    const Eigen::MatrixXd& cosineCoefficientsOfBody2 = effectiveMutualPotentialField_->getCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& sineCoefficientsOfBody2 = effectiveMutualPotentialField_->getSineCoefficientsOfBody2( );
    const Eigen::Vector4d quaternionVectorOfBody1 =
            linear_algebra::convertQuaternionToVectorFormat( accelerationModel_->getCurrentRotationFromInertialToBody1( ) );
    const Eigen::Matrix3d currentRotationFromInertialToBody2 =
            accelerationModel_->getCurrentRotationFromBody2ToBody1( ).toRotationMatrix( ).transpose( ) * currentRotationToBodyFixedFrame_;
    const Eigen::Vector4d quaternionVectorOfBody2 =
            linear_algebra::convertQuaternionToVectorFormat( Eigen::Quaterniond( currentRotationFromInertialToBody2 ) );
    const Eigen::Vector4d conjugatedQuaternionVectorOfBody2 = ( Eigen::Vector4d( ) << quaternionVectorOfBody2( 0 ),
                                                                -quaternionVectorOfBody2( 1 ),
                                                                -quaternionVectorOfBody2( 2 ),
                                                                -quaternionVectorOfBody2( 3 ) )
                                                                      .finished( );

    const Eigen::Matrix4d partialOfInertialToBodyQuaternionWrtBodyToInertialQuaternion =
            Eigen::Vector4d( 1.0, -1.0, -1.0, -1.0 ).asDiagonal( );
    const Eigen::Matrix4d partialOfRelativeQuaternionWrtQuaternionOfBody1 =
            detail::getRightQuaternionMultiplicationMatrix( conjugatedQuaternionVectorOfBody2 ) *
            partialOfInertialToBodyQuaternionWrtBodyToInertialQuaternion;
    const Eigen::Matrix4d partialOfRelativeQuaternionWrtQuaternionOfBody2 =
            detail::getLeftQuaternionMultiplicationMatrix( quaternionVectorOfBody1 );

    detail::computeDerivativeOfWignerDMatricesWrtRelativeQuaternion(
            accelerationModel_->getCurrentRotationFromBody2ToBody1( ),
            effectiveMutualPotentialField_->getTransformationCache( )->getWignerDMatricesCache( ),
            derivativeOfWignerDMatricesWrtRelativeQuaternionScratch_ );

    const double cosineOfLatitude = std::cos( currentSphericalPosition_( 1 ) );
    const double preMultiplier = currentGravitationalParameter_ / currentDistance_;
    std::array< Eigen::Vector3d, 4 > partialOfMutualPotentialGradientWrtRelativeQuaternion;

    for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
    {
        partialOfMutualPotentialGradientWrtRelativeQuaternion.at( quaternionIndex ).setZero( );

        basic_mathematics::transformSphericalHarmonicCoefficientsWithWignerD(
                cosineCoefficientsOfBody2,
                sineCoefficientsOfBody2,
                derivativeOfWignerDMatricesWrtRelativeQuaternionScratch_.at( quaternionIndex ),
                partialOfTransformedCosineCoefficientsBody2Scratch_,
                partialOfTransformedSineCoefficientsBody2Scratch_,
                accelerationModel_->getAreCoefficientsNormalized( ) );

        for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );
            const int totalDegree = degreeOfBody1 + degreeOfBody2;
            const double equatorialRadiusRatioPower = currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );

            for( int variant = 0; variant < 4; variant++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                if( !gravitation::getSignedOrdersForCombinationCase(
                            variant, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                {
                    continue;
                }

                const int transformedOrderOfBody2 = std::abs( signedOrderOfBody2 );
                const int totalOrder = std::abs( signedOrderOfBody1 + signedOrderOfBody2 );
                const int effectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

                const Eigen::Vector2d partialOfTransformedCoefficients =
                        ( Eigen::Vector2d( ) << partialOfTransformedCosineCoefficientsBody2Scratch_( degreeOfBody2,
                                                                                                     transformedOrderOfBody2 ),
                          partialOfTransformedSineCoefficientsBody2Scratch_( degreeOfBody2, transformedOrderOfBody2 ) )
                                .finished( );
                partialOfMutualPotentialGradientWrtRelativeQuaternion.at( quaternionIndex ) +=
                        detail::computeCurrentBodyFixedBasisFunctionGradients( currentBodyFixedRelativePosition_,
                                                                               sphericalHarmonicsCache_,
                                                                               cosineOfLatitude,
                                                                               preMultiplier,
                                                                               equatorialRadiusRatioPower,
                                                                               totalDegree,
                                                                               totalOrder ) *
                        ( currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.at( effectiveIndex ) *
                          partialOfTransformedCoefficients );
            }
        }
    }

    currentPartialWrtQuaternionOfBody1_.setZero( );
    currentPartialWrtQuaternionOfBody2_.setZero( );
    for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
    {
        const Eigen::Vector3d partialOfBodyFixedPositionWrtQuaternionOfBody1 =
                derivativeOfBody1RotationFromBodyFixedToInertialWrtQuaternion.at( quaternionIndex ).transpose( ) *
                accelerationModel_->getCurrentRelativePosition( );

        Eigen::Vector3d coefficientContributionWrtQuaternionOfBody1 = Eigen::Vector3d::Zero( );
        Eigen::Vector3d coefficientContributionWrtQuaternionOfBody2 = Eigen::Vector3d::Zero( );
        for( int relativeQuaternionIndex = 0; relativeQuaternionIndex < 4; relativeQuaternionIndex++ )
        {
            coefficientContributionWrtQuaternionOfBody1 +=
                    partialOfMutualPotentialGradientWrtRelativeQuaternion.at( relativeQuaternionIndex ) *
                    partialOfRelativeQuaternionWrtQuaternionOfBody1( relativeQuaternionIndex, quaternionIndex );
            coefficientContributionWrtQuaternionOfBody2 +=
                    partialOfMutualPotentialGradientWrtRelativeQuaternion.at( relativeQuaternionIndex ) *
                    partialOfRelativeQuaternionWrtQuaternionOfBody2( relativeQuaternionIndex, quaternionIndex );
        }

        currentPartialWrtQuaternionOfBody1_.col( quaternionIndex ) =
                derivativeOfBody1RotationFromBodyFixedToInertialWrtQuaternion.at( quaternionIndex ) *
                        accelerationModel_->getMutualPotentialGradient( ) +
                currentRotationToInertialFrame_ *
                        ( currentBodyFixedPartialWrtPosition_ * partialOfBodyFixedPositionWrtQuaternionOfBody1 +
                          coefficientContributionWrtQuaternionOfBody1 );
        currentPartialWrtQuaternionOfBody2_.col( quaternionIndex ) =
                currentRotationToInertialFrame_ * coefficientContributionWrtQuaternionOfBody2;
    }
}

//! Compute derivatives of transformed body-2 coefficients w.r.t. one original body-2 coefficient.
/*!
 * Uses the same spherical-harmonic rotation machinery as the model to obtain linearized transformed
 * coefficient mappings required for chain-rule terms in Eq. (47)-(48).
 */
void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentTransformedBody2CoefficientPartials(
        const int degree,
        const int order,
        const bool wrtCosineCoefficient,
        Eigen::MatrixXd& transformedCosinePartials,
        Eigen::MatrixXd& transformedSinePartials )
{
    const Eigen::MatrixXd& transformedCosineCoefficients = effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& transformedSineCoefficients = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );

    body2BasisCosineCoefficients_.setZero( transformedCosineCoefficients.rows( ), transformedCosineCoefficients.cols( ) );
    body2BasisSineCoefficients_.setZero( transformedSineCoefficients.rows( ), transformedSineCoefficients.cols( ) );

    if( degree < body2BasisCosineCoefficients_.rows( ) && order < body2BasisCosineCoefficients_.cols( ) )
    {
        if( wrtCosineCoefficient )
        {
            body2BasisCosineCoefficients_( degree, order ) = 1.0;
        }
        else
        {
            body2BasisSineCoefficients_( degree, order ) = 1.0;
        }
    }

    effectiveMutualPotentialField_->getTransformationCache( )->transformCoefficientsAtDegree(
            body2BasisCosineCoefficients_,
            body2BasisSineCoefficients_,
            transformedCosinePartials,
            transformedSinePartials,
            accelerationModel_->getAreCoefficientsNormalized( ) );
    // Eqs. (47)-(48): these transformed-coefficient partials are the chain-rule inputs for effective coefficients.
}

//! Compute partial of acceleration w.r.t. a cosine coefficient block of body 1.
void FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody1(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::MatrixXd& transformedCosineCoefficientsOfBody2 =
            effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& transformedSineCoefficientsOfBody2 = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( j ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( blockIndices.at( i ).first != degreeOfBody1 )
            {
                continue;
            }

            for( int k = 0; k < 4; k++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                if( !gravitation::getSignedOrdersForCombinationCase(
                            k, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                {
                    continue;
                }

                if( blockIndices.at( i ).second != std::abs( signedOrderOfBody1 ) )
                {
                    continue;
                }

                const int transformedOrderOfBody2 = std::abs( signedOrderOfBody2 );
                if( degreeOfBody2 >= transformedCosineCoefficientsOfBody2.rows( ) ||
                    transformedOrderOfBody2 >= transformedCosineCoefficientsOfBody2.cols( ) )
                {
                    continue;
                }

                const int effectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double coefficientMultiplier = effectiveMutualPotentialField_->getMultiplier(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double signOrderOfBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;

                Eigen::Vector2d effectiveCoefficientsPartial;
                effectiveCoefficientsPartial( 0 ) =
                        coefficientMultiplier * transformedCosineCoefficientsOfBody2( degreeOfBody2, transformedOrderOfBody2 );
                effectiveCoefficientsPartial( 1 ) = coefficientMultiplier * signOrderOfBody2 * signTotalOrder *
                        transformedSineCoefficientsOfBody2( degreeOfBody2, transformedOrderOfBody2 );
                partialMatrix.block( 0, i, 3, 1 ) +=
                        currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) * effectiveCoefficientsPartial;
            }
        }
    }
}

//! Compute partial of acceleration w.r.t. a sine coefficient block of body 1.
void FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody1( const std::vector< std::pair< int, int > >& blockIndices,
                                                                                  Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::MatrixXd& transformedCosineCoefficientsOfBody2 =
            effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& transformedSineCoefficientsOfBody2 = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( j ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( blockIndices.at( i ).first != degreeOfBody1 )
            {
                continue;
            }

            for( int k = 0; k < 4; k++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                if( !gravitation::getSignedOrdersForCombinationCase(
                            k, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                {
                    continue;
                }

                if( blockIndices.at( i ).second != std::abs( signedOrderOfBody1 ) )
                {
                    continue;
                }

                const int transformedOrderOfBody2 = std::abs( signedOrderOfBody2 );
                if( degreeOfBody2 >= transformedCosineCoefficientsOfBody2.rows( ) ||
                    transformedOrderOfBody2 >= transformedCosineCoefficientsOfBody2.cols( ) )
                {
                    continue;
                }

                const int effectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double coefficientMultiplier = effectiveMutualPotentialField_->getMultiplier(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const double signOrderOfBody1 = ( signedOrderOfBody1 < 0 ) ? -1.0 : 1.0;
                const double signOrderOfBody2 = ( signedOrderOfBody2 < 0 ) ? -1.0 : 1.0;
                const double signTotalOrder = ( ( signedOrderOfBody1 + signedOrderOfBody2 ) < 0 ) ? -1.0 : 1.0;

                Eigen::Vector2d effectiveCoefficientsPartial;
                effectiveCoefficientsPartial( 0 ) = -coefficientMultiplier * signOrderOfBody1 * signOrderOfBody2 *
                        transformedSineCoefficientsOfBody2( degreeOfBody2, transformedOrderOfBody2 );
                effectiveCoefficientsPartial( 1 ) = coefficientMultiplier * signOrderOfBody1 * signTotalOrder *
                        transformedCosineCoefficientsOfBody2( degreeOfBody2, transformedOrderOfBody2 );
                partialMatrix.block( 0, i, 3, 1 ) +=
                        currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) * effectiveCoefficientsPartial;
            }
        }
    }
}

//! Compute partial of acceleration w.r.t. a cosine coefficient block of body 2.
void FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody2(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        const int degree = blockIndices.at( i ).first;
        const int order = blockIndices.at( i ).second;

        updateCurrentTransformedBody2CoefficientPartials(
                degree, order, true, transformedCosineBody2CoefficientPartialsScratch_, transformedSineBody2CoefficientPartialsScratch_ );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( j ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( degreeOfBody2 != degree )
            {
                continue;
            }

            for( int k = 0; k < 4; k++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                if( !gravitation::getSignedOrdersForCombinationCase(
                            k, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                {
                    continue;
                }

                const int transformedOrderOfBody2 = std::abs( signedOrderOfBody2 );
                if( degreeOfBody2 >= transformedCosineBody2CoefficientPartialsScratch_.rows( ) ||
                    transformedOrderOfBody2 >= transformedCosineBody2CoefficientPartialsScratch_.cols( ) )
                {
                    continue;
                }

                const int effectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const Eigen::Vector2d transformedCoefficientPartials =
                        ( Eigen::Vector2d( ) << transformedCosineBody2CoefficientPartialsScratch_( degreeOfBody2, transformedOrderOfBody2 ),
                          transformedSineBody2CoefficientPartialsScratch_( degreeOfBody2, transformedOrderOfBody2 ) )
                                .finished( );
                partialMatrix.block( 0, i, 3, 1 ) += currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) *
                        ( currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.at( effectiveIndex ) *
                          transformedCoefficientPartials );
            }
        }
    }
}

//! Compute partial of acceleration w.r.t. a sine coefficient block of body 2.
void FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody2( const std::vector< std::pair< int, int > >& blockIndices,
                                                                                  Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        const int degree = blockIndices.at( i ).first;
        const int order = blockIndices.at( i ).second;

        updateCurrentTransformedBody2CoefficientPartials(
                degree, order, false, transformedCosineBody2CoefficientPartialsScratch_, transformedSineBody2CoefficientPartialsScratch_ );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( j ) );
            const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( j ) );
            const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( j ) );

            if( degreeOfBody2 != degree )
            {
                continue;
            }

            for( int k = 0; k < 4; k++ )
            {
                int signedOrderOfBody1 = 0;
                int signedOrderOfBody2 = 0;
                if( !gravitation::getSignedOrdersForCombinationCase(
                            k, orderOfBody1, orderOfBody2, signedOrderOfBody1, signedOrderOfBody2 ) )
                {
                    continue;
                }

                const int transformedOrderOfBody2 = std::abs( signedOrderOfBody2 );
                if( degreeOfBody2 >= transformedCosineBody2CoefficientPartialsScratch_.rows( ) ||
                    transformedOrderOfBody2 >= transformedCosineBody2CoefficientPartialsScratch_.cols( ) )
                {
                    continue;
                }

                const int effectiveIndex = effectiveMutualPotentialField_->getEffectiveIndex(
                        degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
                const Eigen::Vector2d transformedCoefficientPartials =
                        ( Eigen::Vector2d( ) << transformedCosineBody2CoefficientPartialsScratch_( degreeOfBody2, transformedOrderOfBody2 ),
                          transformedSineBody2CoefficientPartialsScratch_( degreeOfBody2, transformedOrderOfBody2 ) )
                                .finished( );
                partialMatrix.block( 0, i, 3, 1 ) += currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) *
                        ( currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.at( effectiveIndex ) *
                          transformedCoefficientPartials );
            }
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtGravitationalParameter( Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = accelerationModel_->getAcceleration( ) / currentGravitationalParameter_;
}

//! Return scalar-parameter partials for gravitational parameters used by the acceleration model.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicsGravityPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter &&
        ( parameter->getParameterName( ).second.first == acceleratingBody_ ||
          ( parameter->getParameterName( ).second.first == acceleratedBody_ && accelerationModel_->getIsMutualAttractionUsed( ) ) ) )
    {
        partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtGravitationalParameter, this, std::placeholders::_1 );
        numberOfRows = 1;
    }

    return std::make_pair( partialFunction, numberOfRows );
}

//! Return vector-parameter partials for spherical-harmonic coefficient blocks.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicsGravityPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    int numberOfRows = 0;

    if( parameter->getParameterName( ).first == estimatable_parameters::spherical_harmonics_cosine_coefficient_block )
    {
        std::shared_ptr< estimatable_parameters::SphericalHarmonicsCosineCoefficients > coefficientsParameter =
                std::dynamic_pointer_cast< estimatable_parameters::SphericalHarmonicsCosineCoefficients >( parameter );
        if( coefficientsParameter == nullptr )
        {
            throw std::runtime_error( "Error when creating mutual extended-body SH partial w.r.t. cosine coefficients, cast failed." );
        }

        if( parameter->getParameterName( ).second.first == acceleratedBody_ )
        {
            partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody1,
                                         this,
                                         coefficientsParameter->getBlockIndices( ),
                                         std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
        else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
        {
            partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody2,
                                         this,
                                         coefficientsParameter->getBlockIndices( ),
                                         std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::spherical_harmonics_sine_coefficient_block )
    {
        std::shared_ptr< estimatable_parameters::SphericalHarmonicsSineCoefficients > coefficientsParameter =
                std::dynamic_pointer_cast< estimatable_parameters::SphericalHarmonicsSineCoefficients >( parameter );
        if( coefficientsParameter == nullptr )
        {
            throw std::runtime_error( "Error when creating mutual extended-body SH partial w.r.t. sine coefficients, cast failed." );
        }

        if( parameter->getParameterName( ).second.first == acceleratedBody_ )
        {
            partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody1,
                                         this,
                                         coefficientsParameter->getBlockIndices( ),
                                         std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
        else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
        {
            partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody2,
                                         this,
                                         coefficientsParameter->getBlockIndices( ),
                                         std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}

bool FullTwoBodySphericalHarmonicsGravityPartial::isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
{
    return integratedStateType == propagators::rotational_state &&
            ( stateReferencePoint.first == acceleratedBody_ || stateReferencePoint.first == acceleratingBody_ ) &&
            stateReferencePoint.second.empty( );
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtNonTranslationalStateOfAdditionalBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType,
        const bool addContribution )
{
    if( integratedStateType != propagators::rotational_state )
    {
        return;
    }

    const Eigen::Matrix< double, 3, 4 >& quaternionPartial =
            stateReferencePoint.first == acceleratedBody_ ? currentPartialWrtQuaternionOfBody1_ : currentPartialWrtQuaternionOfBody2_;
    if( addContribution )
    {
        partialMatrix.block( 0, 0, 3, 4 ) += quaternionPartial;
    }
    else
    {
        partialMatrix.block( 0, 0, 3, 4 ) -= quaternionPartial;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
