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

namespace tudat
{

namespace acceleration_partials
{

FullTwoBodySphericalHarmonicsGravityPartial::FullTwoBodySphericalHarmonicsGravityPartial(
        const std::string& acceleratedBody,
        const std::string& acceleratingBody,
        const std::shared_ptr< gravitation::FullTwoBodySphericalHarmonicAcceleration > accelerationModel ):
    AccelerationPartial(
            acceleratedBody,
            acceleratingBody,
            accelerationModel,
            basic_astrodynamics::full_two_body_spherical_harmonic_gravity ),
    accelerationModel_( accelerationModel ),
    effectiveMutualPotentialField_( accelerationModel_->getEffectiveMutualPotentialField( ) ),
    sphericalHarmonicsCache_( std::make_shared< basic_mathematics::SphericalHarmonicsCache >(
            *accelerationModel_->getSphericalHarmonicsCache( ) ) ),
    coefficientCombinationsToUse_( effectiveMutualPotentialField_->getCoefficientCombinationsToUse( ) )
{
    sphericalHarmonicsCache_->getLegendreCache( ).setComputeSecondDerivatives( 1 );

    const int numberOfEffectiveCoefficients = effectiveMutualPotentialField_->getTotalVectorSize( );
    currentPartialsWrtEffectiveCoefficients_.resize(
            numberOfEffectiveCoefficients, Eigen::Matrix< double, 3, 2 >::Zero( ) );
    currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.resize( numberOfEffectiveCoefficients, Eigen::Matrix2d::Zero( ) );

    effectiveIndicesForCoefficientCombinations_.resize( coefficientCombinationsToUse_.size( ) );
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        effectiveIndicesForCoefficientCombinations_.at( i ) = effectiveMutualPotentialField_->getEffectiveIndex(
                std::get<0>(coefficientCombinationsToUse_.at( i )),
                std::get<1>(coefficientCombinationsToUse_.at( i )),
                std::get<2>(coefficientCombinationsToUse_.at( i )),
                std::get<3>(coefficientCombinationsToUse_.at( i )) );
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::update( const double currentTime )
{
    if( !( currentTime_ == currentTime ) )
    {
        accelerationModel_->updateMembers( currentTime );

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

        updateCurrentPositionPartial( );
        updateCurrentPartialsWrtEffectiveCoefficients( );

        effectiveMutualPotentialField_->computePartialsOfFullCoefficientsWrtTransformedCoefficients(
                currentEffectiveCoefficientsWrtTransformedBody2Coefficients_ );

        currentPartialWrtVelocity_.setZero( );
        currentTime_ = currentTime;
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentPositionPartial( )
{
    const double sineOfLatitude = std::sin( currentSphericalPosition_( 1 ) );
    const double cosineOfLatitude = std::cos( currentSphericalPosition_( 1 ) );
    const double preMultiplier = currentGravitationalParameter_ / currentDistance_;
    sphericalHarmonicsCache_->update( TUDAT_NAN, sineOfLatitude, currentSphericalPosition_( 2 ), TUDAT_NAN );

    Eigen::Matrix3d sphericalHessian = Eigen::Matrix3d::Zero( );
    Eigen::Matrix3d currentSphericalHessianContribution = Eigen::Matrix3d::Zero( );

    int totalOrder = 0;
    bool computeTerm = false;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( i ));
        const int orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( i ));
        const int degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( i ));
        const int orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( i ));

        const int totalDegree = degreeOfBody1 + degreeOfBody2;
        const double equatorialRadiusRatioPower = currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );
        const double effectiveCosineCoefficient =
                effectiveMutualPotentialField_->getEffectiveCosineCoefficient(
                        degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );
        const double effectiveSineCoefficient =
                effectiveMutualPotentialField_->getEffectiveSineCoefficient(
                        degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

        for( unsigned int j = 0; j < 4; j++ )
        {
            switch( j )
            {
            case 0:
                totalOrder = std::abs( orderOfBody1 + orderOfBody2 );
                computeTerm = true;
                break;
            case 1:
                totalOrder = std::abs( -orderOfBody1 + orderOfBody2 );
                computeTerm = ( orderOfBody1 != 0 );
                break;
            case 2:
                totalOrder = std::abs( orderOfBody1 - orderOfBody2 );
                computeTerm = ( orderOfBody2 != 0 );
                break;
            case 3:
                totalOrder = std::abs( -orderOfBody1 - orderOfBody2 );
                computeTerm = !( orderOfBody1 == 0 || orderOfBody2 == 0 );
                break;
            default:
                computeTerm = false;
                break;
            }

            if( computeTerm )
            {
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
                sphericalHessian += currentSphericalHessianContribution;
            }
        }
    }

    const Eigen::Matrix3d sphericalToCartesianGradientMatrix =
            coordinate_conversions::getSphericalToCartesianGradientMatrix( currentBodyFixedRelativePosition_ );
    currentBodyFixedPartialWrtPosition_ =
            sphericalToCartesianGradientMatrix * sphericalHessian * sphericalToCartesianGradientMatrix.transpose( );
    currentBodyFixedPartialWrtPosition_ += coordinate_conversions::getDerivativeOfSphericalToCartesianGradient(
            currentBodyFixedSphericalGradient_,
            currentBodyFixedRelativePosition_ );
    currentPartialWrtPosition_ =
            currentRotationToInertialFrame_ * currentBodyFixedPartialWrtPosition_ * currentRotationToBodyFixedFrame_;
}

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

    int totalOrder = 0;
    bool computeTerm = false;
    for( unsigned int i = 0; i < coefficientCombinationsToUse_.size( ); i++ )
    {
        const int degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( i ));
        const int orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( i ));
        const int degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( i ));
        const int orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( i ));

        const int totalDegree = degreeOfBody1 + degreeOfBody2;
        const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( i );
        const double equatorialRadiusRatioPower = currentRadius1Powers_.at( degreeOfBody1 ) * currentRadius2Powers_.at( degreeOfBody2 );

        for( unsigned int j = 0; j < 4; j++ )
        {
            switch( j )
            {
            case 0:
                totalOrder = std::abs( orderOfBody1 + orderOfBody2 );
                computeTerm = true;
                break;
            case 1:
                totalOrder = std::abs( -orderOfBody1 + orderOfBody2 );
                computeTerm = ( orderOfBody1 != 0 );
                break;
            case 2:
                totalOrder = std::abs( orderOfBody1 - orderOfBody2 );
                computeTerm = ( orderOfBody2 != 0 );
                break;
            case 3:
                totalOrder = std::abs( -orderOfBody1 - orderOfBody2 );
                computeTerm = !( orderOfBody1 == 0 || orderOfBody2 == 0 );
                break;
            default:
                computeTerm = false;
                break;
            }

            if( computeTerm )
            {
                const double legendrePolynomial =
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomial( totalDegree, totalOrder );
                const double legendrePolynomialDerivative =
                        sphericalHarmonicsCache_->getLegendreCache( ).getLegendrePolynomialDerivative( totalDegree, totalOrder );

                currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ).block( 0, 0, 3, 1 ) +=
                        basic_mathematics::computePotentialGradient(
                                currentDistance_,
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
                        basic_mathematics::computePotentialGradient(
                                currentDistance_,
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
                            currentPartialsWrtEffectiveCoefficients_.at( i ).block( 0, 0, 3, 1 ),
                            currentBodyFixedRelativePosition_ );
            currentBodyFixedPartialsWrtEffectiveCoefficients.block( 0, 1, 3, 1 ) =
                    coordinate_conversions::convertSphericalToCartesianGradient(
                            currentPartialsWrtEffectiveCoefficients_.at( i ).block( 0, 1, 3, 1 ),
                            currentBodyFixedRelativePosition_ );
            currentPartialsWrtEffectiveCoefficients_.at( i ) =
                    currentRotationToInertialFrame_ * currentBodyFixedPartialsWrtEffectiveCoefficients;
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::updateCurrentTransformedBody2CoefficientPartials(
        const int degree,
        const int order,
        const bool wrtCosineCoefficient,
        Eigen::MatrixXd& transformedCosinePartials,
        Eigen::MatrixXd& transformedSinePartials )
{
    const Eigen::MatrixXd transformedCosineCoefficients = effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd transformedSineCoefficients = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );

    Eigen::MatrixXd basisCosineCoefficients = Eigen::MatrixXd::Zero(
            transformedCosineCoefficients.rows( ), transformedCosineCoefficients.cols( ) );
    Eigen::MatrixXd basisSineCoefficients = Eigen::MatrixXd::Zero(
            transformedSineCoefficients.rows( ), transformedSineCoefficients.cols( ) );

    if( degree < basisCosineCoefficients.rows( ) && order < basisCosineCoefficients.cols( ) )
    {
        if( wrtCosineCoefficient )
        {
            basisCosineCoefficients( degree, order ) = 1.0;
        }
        else
        {
            basisSineCoefficients( degree, order ) = 1.0;
        }
    }

    effectiveMutualPotentialField_->getTransformationCache( )->transformCoefficientsAtDegree(
            basisCosineCoefficients,
            basisSineCoefficients,
            transformedCosinePartials,
            transformedSinePartials,
            accelerationModel_->getAreCoefficientsNormalized( ) );
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody1(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::MatrixXd transformedCosineCoefficientsOfBody2 = effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd transformedSineCoefficientsOfBody2 = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( j ));
            const int orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( j ));
            const int degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( j ));
            const int orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( j ));

            if( blockIndices.at( i ).first == degreeOfBody1 && blockIndices.at( i ).second == orderOfBody1 )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const double coefficientMultiplier =
                        effectiveMutualPotentialField_->getMultiplier(
                                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

                if( degreeOfBody2 < transformedCosineCoefficientsOfBody2.rows( ) &&
                    orderOfBody2 < transformedCosineCoefficientsOfBody2.cols( ) )
                {
                    Eigen::Vector2d effectiveCoefficientsPartial;
                    effectiveCoefficientsPartial( 0 ) = coefficientMultiplier *
                            transformedCosineCoefficientsOfBody2( degreeOfBody2, orderOfBody2 );
                    effectiveCoefficientsPartial( 1 ) = coefficientMultiplier *
                            transformedSineCoefficientsOfBody2( degreeOfBody2, orderOfBody2 );
                    partialMatrix.block( 0, i, 3, 1 ) +=
                            currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) * effectiveCoefficientsPartial;
                }
            }
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody1(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    const Eigen::MatrixXd transformedCosineCoefficientsOfBody2 = effectiveMutualPotentialField_->getTransformedCosineCoefficientsOfBody2( );
    const Eigen::MatrixXd transformedSineCoefficientsOfBody2 = effectiveMutualPotentialField_->getTransformedSineCoefficientsOfBody2( );

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody1 = std::get<0>(coefficientCombinationsToUse_.at( j ));
            const int orderOfBody1 = std::get<1>(coefficientCombinationsToUse_.at( j ));
            const int degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( j ));
            const int orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( j ));

            if( blockIndices.at( i ).first == degreeOfBody1 && blockIndices.at( i ).second == orderOfBody1 )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const double coefficientMultiplier =
                        effectiveMutualPotentialField_->getMultiplier(
                                degreeOfBody1, orderOfBody1, degreeOfBody2, orderOfBody2 );

                if( degreeOfBody2 < transformedCosineCoefficientsOfBody2.rows( ) &&
                    orderOfBody2 < transformedCosineCoefficientsOfBody2.cols( ) )
                {
                    Eigen::Vector2d effectiveCoefficientsPartial;
                    effectiveCoefficientsPartial( 0 ) = -coefficientMultiplier *
                            transformedSineCoefficientsOfBody2( degreeOfBody2, orderOfBody2 );
                    effectiveCoefficientsPartial( 1 ) = coefficientMultiplier *
                            transformedCosineCoefficientsOfBody2( degreeOfBody2, orderOfBody2 );
                    partialMatrix.block( 0, i, 3, 1 ) +=
                            currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) * effectiveCoefficientsPartial;
                }
            }
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody2(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    Eigen::MatrixXd transformedCosinePartials;
    Eigen::MatrixXd transformedSinePartials;

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        const int degree = blockIndices.at( i ).first;
        const int order = blockIndices.at( i ).second;

        updateCurrentTransformedBody2CoefficientPartials(
                degree, order, true, transformedCosinePartials, transformedSinePartials );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( j ));
            const int orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( j ));

            if( degreeOfBody2 == degree && degreeOfBody2 < transformedCosinePartials.rows( ) &&
                orderOfBody2 < transformedCosinePartials.cols( ) )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const Eigen::Vector2d transformedCoefficientPartials =
                        ( Eigen::Vector2d( ) << transformedCosinePartials( degreeOfBody2, orderOfBody2 ),
                          transformedSinePartials( degreeOfBody2, orderOfBody2 ) )
                                .finished( );
                partialMatrix.block( 0, i, 3, 1 ) += currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) *
                        ( currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.at( effectiveIndex ) *
                          transformedCoefficientPartials );
            }
        }
    }
}

void FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody2(
        const std::vector< std::pair< int, int > >& blockIndices,
        Eigen::MatrixXd& partialMatrix )
{
    partialMatrix = Eigen::MatrixXd::Zero( 3, blockIndices.size( ) );

    Eigen::MatrixXd transformedCosinePartials;
    Eigen::MatrixXd transformedSinePartials;

    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        const int degree = blockIndices.at( i ).first;
        const int order = blockIndices.at( i ).second;

        updateCurrentTransformedBody2CoefficientPartials(
                degree, order, false, transformedCosinePartials, transformedSinePartials );

        for( unsigned int j = 0; j < coefficientCombinationsToUse_.size( ); j++ )
        {
            const int degreeOfBody2 = std::get<2>(coefficientCombinationsToUse_.at( j ));
            const int orderOfBody2 = std::get<3>(coefficientCombinationsToUse_.at( j ));

            if( degreeOfBody2 == degree && degreeOfBody2 < transformedCosinePartials.rows( ) &&
                orderOfBody2 < transformedCosinePartials.cols( ) )
            {
                const int effectiveIndex = effectiveIndicesForCoefficientCombinations_.at( j );
                const Eigen::Vector2d transformedCoefficientPartials =
                        ( Eigen::Vector2d( ) << transformedCosinePartials( degreeOfBody2, orderOfBody2 ),
                          transformedSinePartials( degreeOfBody2, orderOfBody2 ) )
                                .finished( );
                partialMatrix.block( 0, i, 3, 1 ) += currentPartialsWrtEffectiveCoefficients_.at( effectiveIndex ) *
                        ( currentEffectiveCoefficientsWrtTransformedBody2Coefficients_.at( effectiveIndex ) *
                          transformedCoefficientPartials );
            }
        }
    }
}

std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
FullTwoBodySphericalHarmonicsGravityPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::function< void( Eigen::MatrixXd& ) > partialFunction;
    return std::make_pair( partialFunction, 0 );
}

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
            throw std::runtime_error(
                    "Error when creating mutual extended-body SH partial w.r.t. cosine coefficients, cast failed." );
        }

        if( parameter->getParameterName( ).second.first == acceleratedBody_ )
        {
            partialFunction = std::bind(
                    &FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody1,
                    this,
                    coefficientsParameter->getBlockIndices( ),
                    std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
        else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
        {
            partialFunction = std::bind(
                    &FullTwoBodySphericalHarmonicsGravityPartial::wrtCosineCoefficientBlockOfBody2,
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
            throw std::runtime_error(
                    "Error when creating mutual extended-body SH partial w.r.t. sine coefficients, cast failed." );
        }

        if( parameter->getParameterName( ).second.first == acceleratedBody_ )
        {
            partialFunction = std::bind(
                    &FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody1,
                    this,
                    coefficientsParameter->getBlockIndices( ),
                    std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
        else if( parameter->getParameterName( ).second.first == acceleratingBody_ )
        {
            partialFunction = std::bind(
                    &FullTwoBodySphericalHarmonicsGravityPartial::wrtSineCoefficientBlockOfBody2,
                    this,
                    coefficientsParameter->getBlockIndices( ),
                    std::placeholders::_1 );
            numberOfRows = coefficientsParameter->getParameterSize( );
        }
    }

    return std::make_pair( partialFunction, numberOfRows );
}

}  // namespace acceleration_partials

}  // namespace tudat
