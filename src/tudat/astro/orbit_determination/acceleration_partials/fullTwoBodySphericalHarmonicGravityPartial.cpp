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
#include "tudat/astro/orbit_determination/estimatable_parameters/gravityFieldVariationParameters.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicCosineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/sphericalHarmonicSineCoefficients.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/tidalLoveNumber.h"
#include "tudat/basics/utilities.h"
#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/basic/sphericalHarmonicTransformations.h"
#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/basic/wignerDMatrices.h"

namespace tudat
{

namespace acceleration_partials
{

namespace detail
{

std::array< Eigen::Matrix3d, 4 > getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternionForFullTwoBodyTorque(
        const Eigen::Quaterniond& rotationFromInertialToBodyFixedFrame )
{
    std::vector< Eigen::Matrix3d > derivativeList( 4, Eigen::Matrix3d::Zero( ) );
    linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion(
            linear_algebra::convertQuaternionToVectorFormat( rotationFromInertialToBodyFixedFrame.inverse( ) ), derivativeList );

    std::array< Eigen::Matrix3d, 4 > derivativeArray;
    for( int i = 0; i < 4; i++ )
    {
        derivativeArray.at( i ) = derivativeList.at( i );
    }
    return derivativeArray;
}

Eigen::Matrix4d getLeftQuaternionMultiplicationMatrix( const Eigen::Vector4d& quaternion )
{
    Eigen::Matrix4d multiplicationMatrix = Eigen::Matrix4d::Zero( );
    multiplicationMatrix << quaternion( 0 ), -quaternion( 1 ), -quaternion( 2 ), -quaternion( 3 ), quaternion( 1 ), quaternion( 0 ),
            -quaternion( 3 ), quaternion( 2 ), quaternion( 2 ), quaternion( 3 ), quaternion( 0 ), -quaternion( 1 ), quaternion( 3 ),
            -quaternion( 2 ), quaternion( 1 ), quaternion( 0 );
    return multiplicationMatrix;
}

Eigen::Matrix4d getRightQuaternionMultiplicationMatrix( const Eigen::Vector4d& quaternion )
{
    Eigen::Matrix4d multiplicationMatrix = Eigen::Matrix4d::Zero( );
    multiplicationMatrix << quaternion( 0 ), -quaternion( 1 ), -quaternion( 2 ), -quaternion( 3 ), quaternion( 1 ), quaternion( 0 ),
            quaternion( 3 ), -quaternion( 2 ), quaternion( 2 ), -quaternion( 3 ), quaternion( 0 ), quaternion( 1 ), quaternion( 3 ),
            quaternion( 2 ), -quaternion( 1 ), quaternion( 0 );
    return multiplicationMatrix;
}

Eigen::MatrixXd computePartialOfQuaternionWrtRotationMatrixParameter(
        const Eigen::Quaterniond& rotationFromBodyFixedToInertial,
        const std::vector< Eigen::Matrix3d >& partialsOfRotationFromBodyFixedToInertial )
{
    const Eigen::Vector4d quaternionVector = linear_algebra::convertQuaternionToVectorFormat( rotationFromBodyFixedToInertial );
    std::vector< Eigen::Matrix3d > partialsOfRotationWrtQuaternion( 4, Eigen::Matrix3d::Zero( ) );
    linear_algebra::computePartialDerivativeOfRotationMatrixWrtQuaternion( quaternionVector, partialsOfRotationWrtQuaternion );

    Eigen::Matrix< double, 9, 3 > rotationWrtQuaternionVector = Eigen::Matrix< double, 9, 3 >::Zero( );
    for( int quaternionVectorIndex = 0; quaternionVectorIndex < 3; quaternionVectorIndex++ )
    {
        const Eigen::Matrix3d constrainedRotationPartial = partialsOfRotationWrtQuaternion.at( quaternionVectorIndex + 1 ) -
                quaternionVector( quaternionVectorIndex + 1 ) / quaternionVector( 0 ) * partialsOfRotationWrtQuaternion.at( 0 );
        rotationWrtQuaternionVector.col( quaternionVectorIndex ) =
                Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( constrainedRotationPartial.data( ) );
    }

    Eigen::MatrixXd quaternionWrtParameter = Eigen::MatrixXd::Zero( 4, partialsOfRotationFromBodyFixedToInertial.size( ) );
    for( unsigned int parameterIndex = 0; parameterIndex < partialsOfRotationFromBodyFixedToInertial.size( ); parameterIndex++ )
    {
        const Eigen::Matrix< double, 9, 1 > rotationPartial =
                Eigen::Map< const Eigen::Matrix< double, 9, 1 > >( partialsOfRotationFromBodyFixedToInertial.at( parameterIndex ).data( ) );
        const Eigen::Vector3d quaternionVectorWrtParameter = rotationWrtQuaternionVector.colPivHouseholderQr( ).solve( rotationPartial );
        quaternionWrtParameter.col( parameterIndex ).segment( 1, 3 ) = quaternionVectorWrtParameter;
        quaternionWrtParameter( 0, parameterIndex ) =
                -quaternionVector.segment( 1, 3 ).dot( quaternionVectorWrtParameter ) / quaternionVector( 0 );
    }
    return quaternionWrtParameter;
}

Eigen::Matrix< double, 3, 2 > computeCurrentBodyFixedBasisFunctionGradients(
        const Eigen::Vector3d& bodyFixedRelativePosition,
        const std::shared_ptr< basic_mathematics::SphericalHarmonicsCache >& sphericalHarmonicsCache,
        const double cosineOfLatitude,
        const double preMultiplier,
        const double equatorialRadiusRatioPower,
        const int totalDegree,
        const int totalOrder )
{
    const double legendrePolynomial = sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
    const double legendrePolynomialDerivative =
            sphericalHarmonicsCache->getLegendreCache( ).getLegendrePolynomialDerivative( totalDegree, totalOrder );

    Eigen::Matrix< double, 3, 2 > bodyFixedBasisFunctionGradients = Eigen::Matrix< double, 3, 2 >::Zero( );
    bodyFixedBasisFunctionGradients.col( 0 ) = coordinate_conversions::convertSphericalToCartesianGradient(
            basic_mathematics::computePotentialGradient( bodyFixedRelativePosition.norm( ),
                                                         equatorialRadiusRatioPower,
                                                         sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder ),
                                                         sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder ),
                                                         cosineOfLatitude,
                                                         preMultiplier,
                                                         totalDegree,
                                                         totalOrder,
                                                         1.0,
                                                         0.0,
                                                         legendrePolynomial,
                                                         legendrePolynomialDerivative ),
            bodyFixedRelativePosition );
    bodyFixedBasisFunctionGradients.col( 1 ) = coordinate_conversions::convertSphericalToCartesianGradient(
            basic_mathematics::computePotentialGradient( bodyFixedRelativePosition.norm( ),
                                                         equatorialRadiusRatioPower,
                                                         sphericalHarmonicsCache->getCosineOfMultipleLongitude( totalOrder ),
                                                         sphericalHarmonicsCache->getSineOfMultipleLongitude( totalOrder ),
                                                         cosineOfLatitude,
                                                         preMultiplier,
                                                         totalDegree,
                                                         totalOrder,
                                                         0.0,
                                                         1.0,
                                                         legendrePolynomial,
                                                         legendrePolynomialDerivative ),
            bodyFixedRelativePosition );
    return bodyFixedBasisFunctionGradients;
}

}  // namespace detail

//! Constructor.
/*!
 * Initializes caches required for analytical derivatives of the full two-body acceleration model
 * (Dirkx et al. (2019), Eqs. (47)-(49), used in Eq. (55)).
 */
FullTwoBodySphericalHarmonicsGravityPartial::FullTwoBodySphericalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > accelerationModel,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartialsOfBody1,
        const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartialsOfBody2,
        const std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >& tidalLoveNumberPartialInterfaces ):
    AccelerationPartial( acceleratedBody,
                         acceleratingBody,
                         accelerationModel,
                         basic_astrodynamics::full_two_body_spherical_harmonic_gravity ),
    accelerationModel_( accelerationModel ), effectiveMutualPotentialField_( accelerationModel_->getEffectiveMutualPotentialField( ) ),
    sphericalHarmonicsCache_(
            std::make_shared< basic_mathematics::SphericalHarmonicsCache >( *accelerationModel_->getSphericalHarmonicsCache( ) ) ),
    coefficientCombinationsToUse_( effectiveMutualPotentialField_->getCoefficientCombinationsToUse( ) ),
    rotationMatrixPartialsOfBody1_( rotationMatrixPartialsOfBody1 ), rotationMatrixPartialsOfBody2_( rotationMatrixPartialsOfBody2 ),
    tidalLoveNumberPartialInterfaces_( tidalLoveNumberPartialInterfaces )
{
    // Cache setup supports derivatives of the Dirkx et al. (2019), Eq. (49), expansion used in
    // Dirkx et al. (2019), Eq. (55), with effective coefficients defined through Eqs. (47)-(48).
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

        // Step 2: cache current rotations/positions and scalar factors used in Dirkx et al. (2019), Eq. (55), derivatives.
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

        // Step 4: update Jacobian of effective coefficients (Dirkx et al. (2019), Eqs. (47)-(48)) w.r.t transformed body-2 coefficients.
        effectiveMutualPotentialField_->computePartialsOfFullCoefficientsWrtTransformedCoefficients(
                currentEffectiveCoefficientsWrtTransformedBody2Coefficients_ );
        updateCurrentOrientationPartials( );
        for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberPartialInterfaces_.at( i )->update( currentTime );
            tidalLoveNumberPartialInterfaces_.at( i )->updateParameterPartials( );
        }

        currentPartialWrtVelocity_.setZero( );
        currentTime_ = currentTime;
    }
}

//! Update \partial a / \partial r in body-fixed and inertial frames.
/*!
 * Computes the Hessian of the effective potential summation (Dirkx et al. (2019), Eq. (49)) and maps it to Cartesian
 * coordinates, yielding the Jacobian entering derivatives of Dirkx et al. (2019), Eq. (55).
 */
void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentPositionPartial( )
{
    const double sineOfLatitude = std::sin( currentSphericalPosition_( 1 ) );
    const double cosineOfLatitude = std::cos( currentSphericalPosition_( 1 ) );
    const double preMultiplier = currentGravitationalParameter_ / currentDistance_;
    sphericalHarmonicsCache_->update( TUDAT_NAN, sineOfLatitude, currentSphericalPosition_( 2 ), TUDAT_NAN );

    Eigen::Matrix3d sphericalHessian = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d currentSphericalHessianContribution = Eigen::Matrix3d::Zero( );

    // Accumulate spherical Hessian contribution for each selected (l1,m1,l2,m2) term from Dirkx et al. (2019), Eq. (49).
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );

        const int totalDegree = degreeOfBody1 + degreeOfBody2;
        const double equatorialRadiusRatioPower = currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );

        // Expand to signed-order variants consistent with real effective coefficients (Dirkx et al. (2019), Eqs. (47)-(49)).
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
                // Dirkx et al. (2019), Eq. (49): per-term contribution to the spherical Hessian of the mutual potential.
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
    // Dirkx et al. (2019), Eq. (55): Cartesian Jacobian of acceleration from the gradient/Hessian of the mutual potential.
    currentPartialWrtPosition_ = currentRotationToInertialFrame_ * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;
}

//! Update \partial a / \partial C_eff and \partial a / \partial S_eff for all effective indices.
/*!
 * Evaluates basis-function gradients of Dirkx et al. (2019), Eq. (49), to obtain coefficient derivatives used in
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

                // Dirkx et al. (2019), Eq. (49): basis-gradient contribution for one signed (l,m) term.
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
    const Eigen::MatrixXd& cosineCoefficientsOfBody2 = effectiveMutualPotentialField_->getCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& sineCoefficientsOfBody2 = effectiveMutualPotentialField_->getSineCoefficientsOfBody2( );
    const Eigen::Vector4d quaternionVectorOfBody1 =
            linear_algebra::convertQuaternionToVectorFormat( accelerationModel_->getCurrentRotationFromInertialToBody1( ) );
    const Eigen::Matrix4d partialOfRelativeQuaternionWrtQuaternionOfBody2 =
            detail::getLeftQuaternionMultiplicationMatrix( quaternionVectorOfBody1 );

    basic_mathematics::computeDerivativeOfWignerDMatricesWrtQuaternion(
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
        Eigen::Vector3d coefficientContributionWrtQuaternionOfBody2 = Eigen::Vector3d::Zero( );
        for( int relativeQuaternionIndex = 0; relativeQuaternionIndex < 4; relativeQuaternionIndex++ )
        {
            coefficientContributionWrtQuaternionOfBody2 +=
                    partialOfMutualPotentialGradientWrtRelativeQuaternion.at( relativeQuaternionIndex ) *
                    partialOfRelativeQuaternionWrtQuaternionOfBody2( relativeQuaternionIndex, quaternionIndex );
        }

        currentPartialWrtQuaternionOfBody2_.col( quaternionIndex ) =
                currentRotationToInertialFrame_ * coefficientContributionWrtQuaternionOfBody2;
    }

    currentPartialWrtQuaternionOfBody1_ = computeCurrentPartialWrtQuaternionOfBody1( );
}

Eigen::Matrix< double, 3, 4 > FullTwoBodySphericalHarmonicsGravityPartial::computeCurrentPartialWrtQuaternionOfBody1( )
{
    const std::array< Eigen::Matrix3d, 4 > derivativeOfBody1RotationFromBodyFixedToInertialWrtQuaternion =
            detail::getDerivativeOfBodyFixedToInertialRotationMatrixWrtQuaternionForFullTwoBodyTorque(
                    accelerationModel_->getCurrentRotationFromInertialToBody1( ) );

    const Eigen::MatrixXd& cosineCoefficientsOfBody2 = effectiveMutualPotentialField_->getCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd& sineCoefficientsOfBody2 = effectiveMutualPotentialField_->getSineCoefficientsOfBody2( );
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
            -detail::getRightQuaternionMultiplicationMatrix( conjugatedQuaternionVectorOfBody2 ) *
            partialOfInertialToBodyQuaternionWrtBodyToInertialQuaternion;

    basic_mathematics::computeDerivativeOfWignerDMatricesWrtQuaternion(
            accelerationModel_->getCurrentRotationFromBody2ToBody1( ),
            effectiveMutualPotentialField_->getTransformationCache( )->getWignerDMatricesCache( ),
            derivativeOfWignerDMatricesWrtRelativeQuaternionScratch_ );

    const double sineOfLatitude = std::sin( currentSphericalPosition_( 1 ) );
    const double cosineOfLatitude = std::cos( currentSphericalPosition_( 1 ) );
    const double preMultiplier = currentGravitationalParameter_ / currentDistance_;
    sphericalHarmonicsCache_->update( TUDAT_NAN, sineOfLatitude, currentSphericalPosition_( 2 ), TUDAT_NAN );

    std::array< Eigen::MatrixXd, 4 > partialOfTransformedCosineCoefficientsBody2;
    std::array< Eigen::MatrixXd, 4 > partialOfTransformedSineCoefficientsBody2;
    for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
    {
        basic_mathematics::transformSphericalHarmonicCoefficientsWithWignerD(
                cosineCoefficientsOfBody2,
                sineCoefficientsOfBody2,
                derivativeOfWignerDMatricesWrtRelativeQuaternionScratch_.at( quaternionIndex ),
                partialOfTransformedCosineCoefficientsBody2.at( quaternionIndex ),
                partialOfTransformedSineCoefficientsBody2.at( quaternionIndex ),
                accelerationModel_->getAreCoefficientsNormalized( ) );
    }

    Eigen::Matrix< double, 3, 4 > partialWrtQuaternionOfBody1 = Eigen::Matrix< double, 3, 4 >::Zero( );
    Eigen::Matrix3d sphericalHessianContribution = Eigen::Matrix3d::Zero( );

    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get< 0 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody1 = std::get< 1 >( coefficientCombinationsToUse_.at( i ) );
        const int degreeOfBody2 = std::get< 2 >( coefficientCombinationsToUse_.at( i ) );
        const int orderOfBody2 = std::get< 3 >( coefficientCombinationsToUse_.at( i ) );
        if( degreeOfBody1 == 0 && orderOfBody1 == 0 )
        {
            // A body-1 point mass has no orientation-dependent acceleration contribution.
            continue;
        }

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
            const double effectiveCosineCoefficient = effectiveMutualPotentialField_->getEffectiveCosineCoefficient(
                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );
            const double effectiveSineCoefficient = effectiveMutualPotentialField_->getEffectiveSineCoefficient(
                    degreeOfBody1, signedOrderOfBody1, degreeOfBody2, signedOrderOfBody2 );

            const double legendrePolynomial =
                    sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
            const double legendrePolynomialDerivative =
                    sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomialDerivative( totalDegree, totalOrder );
            const Eigen::Vector3d sphericalGradientContribution =
                    basic_mathematics::computePotentialGradient( currentDistance_,
                                                                 equatorialRadiusRatioPower,
                                                                 sphericalHarmonicsCache_->getCosineOfMultipleLongitude( totalOrder ),
                                                                 sphericalHarmonicsCache_->getSineOfMultipleLongitude( totalOrder ),
                                                                 cosineOfLatitude,
                                                                 preMultiplier,
                                                                 totalDegree,
                                                                 totalOrder,
                                                                 effectiveCosineCoefficient,
                                                                 effectiveSineCoefficient,
                                                                 legendrePolynomial,
                                                                 legendrePolynomialDerivative );
            const Eigen::Vector3d bodyFixedGradientContribution = coordinate_conversions::convertSphericalToCartesianGradient(
                    sphericalGradientContribution, currentBodyFixedRelativePosition_ );

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
                    legendrePolynomial,
                    legendrePolynomialDerivative,
                    sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomialSecondDerivative( totalDegree, totalOrder ),
                    sphericalHessianContribution );

            const Eigen::Matrix3d sphericalToCartesianGradientMatrix =
                    coordinate_conversions::getSphericalToCartesianGradientMatrix( currentBodyFixedRelativePosition_ );
            Eigen::Matrix3d bodyFixedPartialWrtPositionContribution =
                    sphericalToCartesianGradientMatrix * sphericalHessianContribution * sphericalToCartesianGradientMatrix.transpose( );
            bodyFixedPartialWrtPositionContribution += coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
                    sphericalGradientContribution, currentBodyFixedRelativePosition_ );

            std::array< Eigen::Vector3d, 4 > coefficientContributionWrtRelativeQuaternion;
            for( int relativeQuaternionIndex = 0; relativeQuaternionIndex < 4; relativeQuaternionIndex++ )
            {
                const Eigen::Vector2d partialOfTransformedCoefficients =
                        ( Eigen::Vector2d( ) << partialOfTransformedCosineCoefficientsBody2.at( relativeQuaternionIndex )(
                                  degreeOfBody2, transformedOrderOfBody2 ),
                          partialOfTransformedSineCoefficientsBody2.at( relativeQuaternionIndex )( degreeOfBody2,
                                                                                                   transformedOrderOfBody2 ) )
                                .finished( );
                coefficientContributionWrtRelativeQuaternion.at( relativeQuaternionIndex ) =
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

            for( int quaternionIndex = 0; quaternionIndex < 4; quaternionIndex++ )
            {
                const Eigen::Vector3d partialOfBodyFixedPositionWrtQuaternionOfBody1 =
                        derivativeOfBody1RotationFromBodyFixedToInertialWrtQuaternion.at( quaternionIndex ).transpose( ) *
                        accelerationModel_->getCurrentRelativePosition( );

                Eigen::Vector3d coefficientContributionWrtQuaternionOfBody1 = Eigen::Vector3d::Zero( );
                for( int relativeQuaternionIndex = 0; relativeQuaternionIndex < 4; relativeQuaternionIndex++ )
                {
                    coefficientContributionWrtQuaternionOfBody1 +=
                            coefficientContributionWrtRelativeQuaternion.at( relativeQuaternionIndex ) *
                            partialOfRelativeQuaternionWrtQuaternionOfBody1( relativeQuaternionIndex, quaternionIndex );
                }

                partialWrtQuaternionOfBody1.col( quaternionIndex ) +=
                        derivativeOfBody1RotationFromBodyFixedToInertialWrtQuaternion.at( quaternionIndex ) *
                                bodyFixedGradientContribution +
                        currentRotationToInertialFrame_ *
                                ( bodyFixedPartialWrtPositionContribution * partialOfBodyFixedPositionWrtQuaternionOfBody1 +
                                  coefficientContributionWrtQuaternionOfBody1 );
            }
        }
    }

    return partialWrtQuaternionOfBody1;
}

//! Compute derivatives of transformed body-2 coefficients w.r.t. one original body-2 coefficient.
/*!
 * Uses the same spherical-harmonic rotation machinery as the model to obtain linearized transformed
 * coefficient mappings required for chain-rule terms in Dirkx et al. (2019), Eqs. (47)-(48).
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
    // Dirkx et al. (2019), Eqs. (47)-(48): these transformed-coefficient partials are the chain-rule inputs
    // for effective coefficients.
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

void FullTwoBodySphericalHarmonicsGravityPartial::wrtRotationModelParameter(
        Eigen::MatrixXd& partialMatrix,
        const bool wrtBody1,
        const estimatable_parameters::EstimatebleParametersEnum parameterType,
        const std::string& secondaryIdentifier )
{
    const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials =
            wrtBody1 ? rotationMatrixPartialsOfBody1_ : rotationMatrixPartialsOfBody2_;
    const std::shared_ptr< observation_partials::RotationMatrixPartial > rotationMatrixPartial =
            rotationMatrixPartials.at( std::make_pair( parameterType, secondaryIdentifier ) );
    const std::vector< Eigen::Matrix3d > currentRotationMatrixPartials =
            rotationMatrixPartial->calculatePartialOfRotationMatrixToBaseFrameWrParameter( currentTime_ );
    const Eigen::Quaterniond currentRotationFromBodyFixedToInertial =
            rotationMatrixPartial->getRotationModel( )->getRotationToBaseFrame( currentTime_ );

    const Eigen::MatrixXd partialOfQuaternionWrtParameter = detail::computePartialOfQuaternionWrtRotationMatrixParameter(
            currentRotationFromBodyFixedToInertial, currentRotationMatrixPartials );
    partialMatrix =
            ( wrtBody1 ? currentPartialWrtQuaternionOfBody1_ : currentPartialWrtQuaternionOfBody2_ ) * partialOfQuaternionWrtParameter;
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtPolynomialGravityFieldVariations(
        const bool wrtBody1,
        const std::vector< std::pair< int, int > >& cosineBlockIndices,
        const std::vector< std::pair< int, int > >& sineBlockIndices,
        const std::vector< std::vector< std::pair< int, int > > > powersPerCosineBlockIndex,
        const std::vector< std::vector< std::pair< int, int > > > powersPerSineBlockIndex,
        const double referenceEpoch,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::MatrixXd cosineCoefficientPartials;
    Eigen::MatrixXd sineCoefficientPartials;
    if( wrtBody1 )
    {
        wrtCosineCoefficientBlockOfBody1( cosineBlockIndices, cosineCoefficientPartials );
        wrtSineCoefficientBlockOfBody1( sineBlockIndices, sineCoefficientPartials );
    }
    else
    {
        wrtCosineCoefficientBlockOfBody2( cosineBlockIndices, cosineCoefficientPartials );
        wrtSineCoefficientBlockOfBody2( sineBlockIndices, sineCoefficientPartials );
    }

    partialMatrix.setZero( );

    int numberOfCosineParameters = 0;
    for( unsigned int i = 0; i < powersPerCosineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < powersPerCosineBlockIndex.at( i ).size( ); j++ )
        {
            partialMatrix.block( 0, powersPerCosineBlockIndex.at( i ).at( j ).first, 3, 1 ) +=
                    cosineCoefficientPartials.block( 0, i, 3, 1 ) *
                    std::pow( currentTime_ - referenceEpoch, powersPerCosineBlockIndex.at( i ).at( j ).second );
            numberOfCosineParameters++;
        }
    }

    for( unsigned int i = 0; i < powersPerSineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < powersPerSineBlockIndex.at( i ).size( ); j++ )
        {
            partialMatrix.block( 0, powersPerSineBlockIndex.at( i ).at( j ).first + numberOfCosineParameters, 3, 1 ) +=
                    sineCoefficientPartials.block( 0, i, 3, 1 ) *
                    std::pow( currentTime_ - referenceEpoch, powersPerSineBlockIndex.at( i ).at( j ).second );
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtPeriodicGravityFieldVariations(
        const bool wrtBody1,
        const std::vector< std::pair< int, int > >& cosineBlockIndices,
        const std::vector< std::pair< int, int > >& sineBlockIndices,
        const std::vector< std::vector< std::pair< int, int > > > periodsPerCosineBlockIndex,
        const std::vector< std::vector< std::pair< int, int > > > periodsPerSineBlockIndex,
        const std::vector< double >& frequencies,
        const double referenceEpoch,
        Eigen::MatrixXd& partialMatrix )
{
    Eigen::MatrixXd cosineCoefficientPartials;
    Eigen::MatrixXd sineCoefficientPartials;
    if( wrtBody1 )
    {
        wrtCosineCoefficientBlockOfBody1( cosineBlockIndices, cosineCoefficientPartials );
        wrtSineCoefficientBlockOfBody1( sineBlockIndices, sineCoefficientPartials );
    }
    else
    {
        wrtCosineCoefficientBlockOfBody2( cosineBlockIndices, cosineCoefficientPartials );
        wrtSineCoefficientBlockOfBody2( sineBlockIndices, sineCoefficientPartials );
    }

    partialMatrix.setZero( );

    int numberOfCosineParameters = 0;
    for( unsigned int i = 0; i < periodsPerCosineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < periodsPerCosineBlockIndex.at( i ).size( ); j++ )
        {
            const double frequency = frequencies.at( periodsPerCosineBlockIndex.at( i ).at( j ).second );
            partialMatrix.block( 0, 2 * periodsPerCosineBlockIndex.at( i ).at( j ).first, 3, 1 ) +=
                    cosineCoefficientPartials.block( 0, i, 3, 1 ) * std::cos( frequency * ( currentTime_ - referenceEpoch ) );
            partialMatrix.block( 0, 2 * periodsPerCosineBlockIndex.at( i ).at( j ).first + 1, 3, 1 ) +=
                    cosineCoefficientPartials.block( 0, i, 3, 1 ) * std::sin( frequency * ( currentTime_ - referenceEpoch ) );
            numberOfCosineParameters++;
        }
    }

    for( unsigned int i = 0; i < periodsPerSineBlockIndex.size( ); i++ )
    {
        for( unsigned int j = 0; j < periodsPerSineBlockIndex.at( i ).size( ); j++ )
        {
            const double frequency = frequencies.at( periodsPerSineBlockIndex.at( i ).at( j ).second );
            partialMatrix.block( 0, 2 * ( periodsPerSineBlockIndex.at( i ).at( j ).first + numberOfCosineParameters ), 3, 1 ) +=
                    sineCoefficientPartials.block( 0, i, 3, 1 ) * std::cos( frequency * ( currentTime_ - referenceEpoch ) );
            partialMatrix.block( 0, 2 * ( periodsPerSineBlockIndex.at( i ).at( j ).first + numberOfCosineParameters ) + 1, 3, 1 ) +=
                    sineCoefficientPartials.block( 0, i, 3, 1 ) * std::sin( frequency * ( currentTime_ - referenceEpoch ) );
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtTidalLoveNumber(
        const bool wrtBody1,
        const std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunctions,
        const int degree,
        const std::vector< int >& orders,
        const bool sumOrders,
        const int parameterSize,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, parameterSize );

    const std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > coefficientPartialsPerOrder = coefficientPartialFunctions( );
    if( coefficientPartialsPerOrder.size( ) != orders.size( ) )
    {
        throw std::runtime_error( "Error when computing full two-body acceleration Love-number partial, inconsistent order count." );
    }

    std::vector< std::pair< int, int > > singleCoefficientIndex( 1 );
    Eigen::MatrixXd partialWrtCosineCoefficient;
    Eigen::MatrixXd partialWrtSineCoefficient;
    for( unsigned int i = 0; i < orders.size( ); i++ )
    {
        singleCoefficientIndex[ 0 ] = std::make_pair( degree, orders.at( i ) );

        if( wrtBody1 )
        {
            wrtCosineCoefficientBlockOfBody1( singleCoefficientIndex, partialWrtCosineCoefficient );
            wrtSineCoefficientBlockOfBody1( singleCoefficientIndex, partialWrtSineCoefficient );
        }
        else
        {
            wrtCosineCoefficientBlockOfBody2( singleCoefficientIndex, partialWrtCosineCoefficient );
            wrtSineCoefficientBlockOfBody2( singleCoefficientIndex, partialWrtSineCoefficient );
        }

        const int singleOrderPartialSize = coefficientPartialsPerOrder.at( i ).cols( );
        const int startColumn = sumOrders ? 0 : static_cast< int >( i ) * singleOrderPartialSize;
        partialMatrix.block( 0, startColumn, 3, singleOrderPartialSize ) +=
                partialWrtCosineCoefficient * coefficientPartialsPerOrder.at( i ).block( 0, 0, 1, singleOrderPartialSize ) +
                partialWrtSineCoefficient * coefficientPartialsPerOrder.at( i ).block( 1, 0, 1, singleOrderPartialSize );
    }
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
    else if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
    {
        const bool isBody1Parameter = parameter->getParameterName( ).second.first == acceleratedBody_;
        const bool isBody2Parameter = parameter->getParameterName( ).second.first == acceleratingBody_;
        if( isBody1Parameter || isBody2Parameter )
        {
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials =
                    isBody1Parameter ? rotationMatrixPartialsOfBody1_ : rotationMatrixPartialsOfBody2_;
            if( rotationMatrixPartials.count(
                        std::make_pair( parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) ) ) == 0 )
            {
                throw std::runtime_error( "Error, missing full two-body acceleration rotation matrix partial for parameter " +
                                          std::to_string( parameter->getParameterName( ).first ) + " of " +
                                          parameter->getParameterName( ).second.first );
            }
            partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                                         this,
                                         std::placeholders::_1,
                                         isBody1Parameter,
                                         parameter->getParameterName( ).first,
                                         parameter->getSecondaryIdentifier( ) );
            numberOfRows = parameter->getParameterSize( );
        }
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

    if( estimatable_parameters::isParameterRotationMatrixProperty( parameter->getParameterName( ).first ) )
    {
        const bool isBody1Parameter = parameter->getParameterName( ).second.first == acceleratedBody_;
        const bool isBody2Parameter = parameter->getParameterName( ).second.first == acceleratingBody_;
        if( isBody1Parameter || isBody2Parameter )
        {
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials =
                    isBody1Parameter ? rotationMatrixPartialsOfBody1_ : rotationMatrixPartialsOfBody2_;
            if( rotationMatrixPartials.count(
                        std::make_pair( parameter->getParameterName( ).first, parameter->getSecondaryIdentifier( ) ) ) == 0 )
            {
                throw std::runtime_error( "Error, missing full two-body acceleration rotation matrix partial for parameter " +
                                          std::to_string( parameter->getParameterName( ).first ) + " of " +
                                          parameter->getParameterName( ).second.first );
            }
            partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtRotationModelParameter,
                                         this,
                                         std::placeholders::_1,
                                         isBody1Parameter,
                                         parameter->getParameterName( ).first,
                                         parameter->getSecondaryIdentifier( ) );
            numberOfRows = parameter->getParameterSize( );
        }
    }
    else if( estimatable_parameters::isParameterNonTidalGravityFieldVariationProperty( parameter->getParameterName( ).first ) )
    {
        const bool isBody1Parameter = parameter->getParameterName( ).second.first == acceleratedBody_;
        const bool isBody2Parameter = parameter->getParameterName( ).second.first == acceleratingBody_;
        if( isBody1Parameter || isBody2Parameter )
        {
            switch( parameter->getParameterName( ).first )
            {
                case estimatable_parameters::polynomial_gravity_field_variation_amplitudes: {
                    std::shared_ptr< estimatable_parameters::PolynomialGravityFieldVariationsParameters > polynomialVariationParameter =
                            std::dynamic_pointer_cast< estimatable_parameters::PolynomialGravityFieldVariationsParameters >( parameter );
                    if( polynomialVariationParameter == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when creating full two-body acceleration partial w.r.t. polynomial gravity variations, cast "
                                "failed." );
                    }
                    const std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex =
                            polynomialVariationParameter->getIndexAndPowerPerCosineBlockIndex( );
                    const std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex =
                            polynomialVariationParameter->getIndexAndPowerPerSineBlockIndex( );
                    partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtPolynomialGravityFieldVariations,
                                                 this,
                                                 isBody1Parameter,
                                                 utilities::createVectorFromMapKeys( indexAndPowerPerCosineBlockIndex ),
                                                 utilities::createVectorFromMapKeys( indexAndPowerPerSineBlockIndex ),
                                                 utilities::createVectorFromMapValues( indexAndPowerPerCosineBlockIndex ),
                                                 utilities::createVectorFromMapValues( indexAndPowerPerSineBlockIndex ),
                                                 polynomialVariationParameter->getPolynomialVariationModel( )->getReferenceEpoch( ),
                                                 std::placeholders::_1 );
                    numberOfRows = parameter->getParameterSize( );
                    break;
                }
                case estimatable_parameters::periodic_gravity_field_variation_amplitudes: {
                    std::shared_ptr< estimatable_parameters::PeriodicGravityFieldVariationsParameters > periodicVariationParameter =
                            std::dynamic_pointer_cast< estimatable_parameters::PeriodicGravityFieldVariationsParameters >( parameter );
                    if( periodicVariationParameter == nullptr )
                    {
                        throw std::runtime_error(
                                "Error when creating full two-body acceleration partial w.r.t. periodic gravity variations, cast failed." );
                    }
                    const std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPeriodPerCosineBlockIndex =
                            periodicVariationParameter->getIndexAndPowerPerCosineBlockIndex( );
                    const std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPeriodPerSineBlockIndex =
                            periodicVariationParameter->getIndexAndPowerPerSineBlockIndex( );
                    partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtPeriodicGravityFieldVariations,
                                                 this,
                                                 isBody1Parameter,
                                                 utilities::createVectorFromMapKeys( indexAndPeriodPerCosineBlockIndex ),
                                                 utilities::createVectorFromMapKeys( indexAndPeriodPerSineBlockIndex ),
                                                 utilities::createVectorFromMapValues( indexAndPeriodPerCosineBlockIndex ),
                                                 utilities::createVectorFromMapValues( indexAndPeriodPerSineBlockIndex ),
                                                 periodicVariationParameter->getPeriodicVariationModel( )->getFrequencies( ),
                                                 periodicVariationParameter->getPeriodicVariationModel( )->getReferenceEpoch( ),
                                                 std::placeholders::_1 );
                    numberOfRows = parameter->getParameterSize( );
                    break;
                }
                default:
                    break;
            }
        }
    }
    else if( estimatable_parameters::isParameterTidalProperty( parameter->getParameterName( ).first ) )
    {
        const bool isBody1Parameter = parameter->getParameterName( ).second.first == acceleratedBody_;
        const bool isBody2Parameter = parameter->getParameterName( ).second.first == acceleratingBody_;
        if( isBody1Parameter || isBody2Parameter )
        {
            if( parameter->getParameterName( ).first == estimatable_parameters::mode_coupled_tidal_love_numbers )
            {
                throw std::runtime_error(
                        "Error, full two-body acceleration partials w.r.t. mode-coupled tidal Love numbers are not "
                        "implemented." );
            }

            std::shared_ptr< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > > tidalLoveNumber =
                    std::dynamic_pointer_cast< estimatable_parameters::TidalLoveNumber< Eigen::VectorXd > >( parameter );
            if( tidalLoveNumber == nullptr )
            {
                throw std::runtime_error(
                        "Error when getting full two-body acceleration partial w.r.t. tidal Love number, parameter cast failed." );
            }

            int maximumDegree = 0;
            int maximumOrder = 0;
            for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
            {
                maximumDegree = std::max( maximumDegree,
                                          static_cast< int >( isBody1Parameter ? std::get< 0 >( coefficientCombinationsToUse_.at( i ) )
                                                                               : std::get< 2 >( coefficientCombinationsToUse_.at( i ) ) ) );
                maximumOrder = std::max( maximumOrder,
                                         static_cast< int >( isBody1Parameter ? std::get< 1 >( coefficientCombinationsToUse_.at( i ) )
                                                                              : std::get< 3 >( coefficientCombinationsToUse_.at( i ) ) ) );
            }

            std::pair< int, std::pair< int, int > > currentTidalPartialOutput;
            std::vector< std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > > coefficientPartialFunctions;
            std::pair< int, int > selectedDegreeAndOrder = std::make_pair( -1, -1 );
            for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
            {
                currentTidalPartialOutput =
                        tidalLoveNumberPartialInterfaces_.at( i )->setParameterPartialFunction( parameter, maximumDegree, maximumOrder );
                if( currentTidalPartialOutput.first > 0 )
                {
                    if( numberOfRows == 0 )
                    {
                        numberOfRows = currentTidalPartialOutput.first;
                        selectedDegreeAndOrder = currentTidalPartialOutput.second;
                    }
                    else
                    {
                        if( numberOfRows != currentTidalPartialOutput.first )
                        {
                            throw std::runtime_error(
                                    "Error when getting full two-body acceleration partial w.r.t. tidal Love number, inconsistent "
                                    "parameter sizes found." );
                        }
                        if( selectedDegreeAndOrder != currentTidalPartialOutput.second )
                        {
                            throw std::runtime_error(
                                    "Error when getting full two-body acceleration partial w.r.t. tidal Love number, inconsistent "
                                    "degree/order metadata found." );
                        }
                    }

                    std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > currentInterface =
                            tidalLoveNumberPartialInterfaces_.at( i );
                    std::pair< int, int > currentDegreeAndOrder = currentTidalPartialOutput.second;
                    coefficientPartialFunctions.push_back( [ currentInterface, parameter, currentDegreeAndOrder ]( ) {
                        return std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >(
                                currentInterface->getCurrentVectorParameterPartial( parameter, currentDegreeAndOrder ) );
                    } );
                }
            }

            if( coefficientPartialFunctions.size( ) > 0 )
            {
                std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunction =
                        [ coefficientPartialFunctions ]( ) {
                            std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > totalCoefficientPartials(
                                    coefficientPartialFunctions.front( )( ) );
                            for( unsigned int functionIndex = 1; functionIndex < coefficientPartialFunctions.size( ); functionIndex++ )
                            {
                                const std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > > currentCoefficientPartials(
                                        coefficientPartialFunctions.at( functionIndex )( ) );
                                if( currentCoefficientPartials.size( ) != totalCoefficientPartials.size( ) )
                                {
                                    throw std::runtime_error(
                                            "Error when summing full two-body tidal coefficient partials, inconsistent number of orders "
                                            "found." );
                                }
                                for( unsigned int orderIndex = 0; orderIndex < totalCoefficientPartials.size( ); orderIndex++ )
                                {
                                    if( currentCoefficientPartials.at( orderIndex ).rows( ) !=
                                                totalCoefficientPartials.at( orderIndex ).rows( ) ||
                                        currentCoefficientPartials.at( orderIndex ).cols( ) !=
                                                totalCoefficientPartials.at( orderIndex ).cols( ) )
                                    {
                                        throw std::runtime_error(
                                                "Error when summing full two-body tidal coefficient partials, inconsistent matrix sizes "
                                                "found." );
                                    }
                                    totalCoefficientPartials.at( orderIndex ) += currentCoefficientPartials.at( orderIndex );
                                }
                            }
                            return totalCoefficientPartials;
                        };
                partialFunction = std::bind( &FullTwoBodySphericalHarmonicsGravityPartial::wrtTidalLoveNumber,
                                             this,
                                             isBody1Parameter,
                                             coefficientPartialFunction,
                                             tidalLoveNumber->getDegree( ),
                                             tidalLoveNumber->getOrders( ),
                                             tidalLoveNumber->getSumOrders( ),
                                             parameter->getParameterSize( ),
                                             std::placeholders::_1 );
            }
        }
    }
    else if( parameter->getParameterName( ).first == estimatable_parameters::spherical_harmonics_cosine_coefficient_block )
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

    const bool wrtBody1 = stateReferencePoint.first == acceleratedBody_;
    const Eigen::Matrix< double, 3, 4 >& quaternionPartial =
            wrtBody1 ? currentPartialWrtQuaternionOfBody1_ : currentPartialWrtQuaternionOfBody2_;
    Eigen::MatrixXd rotationalStatePartial = Eigen::MatrixXd::Zero( 3, 7 );
    rotationalStatePartial.block( 0, 0, 3, 4 ) = quaternionPartial;

    if( addContribution )
    {
        partialMatrix.block( 0, 0, 3, 7 ) += rotationalStatePartial;
    }
    else
    {
        partialMatrix.block( 0, 0, 3, 7 ) -= rotationalStatePartial;
    }
}

}  // namespace acceleration_partials

}  // namespace tudat
