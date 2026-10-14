/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          Garmier, R., Barriot, J.-P. (2001), Celestial Mechanics and Dynamical Astronomy, 79(4), 235-275.
 *          Hu, X., Jekeli, C. (2015), Journal of Geodesy, 89(2), 159-177.
 *          See tudat/math/basic/ellipsoidalHarmonics.h for the full citations.
 */

#include "tudat/astro/gravitation/ellipsoidalHarmonicsGravityField.h"

#include <cmath>
#include <stdexcept>

#include "tudat/basics/basicTypedefs.h"

namespace tudat
{

namespace gravitation
{

using namespace tudat::basic_mathematics;

EllipsoidalHarmonicsGravityCache::EllipsoidalHarmonicsGravityCache( const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                                    const double firstFocalParameter,
                                                                    const double secondFocalParameter,
                                                                    const int quadratureOrder,
                                                                    const std::shared_ptr< EllipsoidalHarmonicRadialTable > radialTable ):
    lameFunctionCoefficients_( lameFunctionCoefficients ),
    legendreCache_( lameFunctionCoefficients.getMaximumDegree( ), lameFunctionCoefficients.getMaximumDegree( ), true ),
    legendreCacheWithSecondDerivatives_( lameFunctionCoefficients.getMaximumDegree( ), lameFunctionCoefficients.getMaximumDegree( ), true ),
    legendreCacheWithDerivatives_( lameFunctionCoefficients.getMaximumDegree( ), lameFunctionCoefficients.getMaximumDegree( ), true ),
    firstFocalParameter_( firstFocalParameter ), secondFocalParameter_( secondFocalParameter ), radialTable_( radialTable )
{
    const int maximumDegree = lameFunctionCoefficients.getMaximumDegree( );
    radialFunctions_ = createEllipsoidalHarmonicCoefficients( maximumDegree );
    radialFunctionDerivatives_ = createEllipsoidalHarmonicCoefficients( maximumDegree );
    radialFunctionSecondDerivatives_ = createEllipsoidalHarmonicCoefficients( maximumDegree );
    forwardRadialFunctions_ = createEllipsoidalHarmonicCoefficients( maximumDegree );
    forwardRadialFunctionDerivatives_ = createEllipsoidalHarmonicCoefficients( maximumDegree );
    backwardRadialFunctions_ = createEllipsoidalHarmonicCoefficients( maximumDegree );
    backwardRadialFunctionDerivatives_ = createEllipsoidalHarmonicCoefficients( maximumDegree );

    // The Legendre derivatives are singular on the long axis of the body, so the cache used for
    // the potential leaves them out and only the one used for the gradient computes them.
    legendreCache_.setComputeFirstDerivatives( false );
    legendreCache_.setComputeSecondDerivatives( false );
    legendreCacheWithDerivatives_.setComputeFirstDerivatives( true );
    legendreCacheWithDerivatives_.setComputeSecondDerivatives( false );
    legendreCacheWithSecondDerivatives_.setComputeFirstDerivatives( true );
    legendreCacheWithSecondDerivatives_.setComputeSecondDerivatives( true );

    computeGaussLegendreNodesAndWeights( quadratureOrder, quadratureNodes_, quadratureWeights_ );

    // Starting from NaN makes the first update always run.
    currentBodyFixedPosition_ = ( Eigen::Vector3d( ) << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    positionOfCachedDerivatives_ = currentBodyFixedPosition_;
    positionOfCachedSecondDerivatives_ = currentBodyFixedPosition_;
    ellipsoidalCoordinates_ = Eigen::Vector3d::Zero( );
}

void EllipsoidalHarmonicsGravityCache::setRadialTable( const std::shared_ptr< EllipsoidalHarmonicRadialTable > radialTable )
{
    radialTable_ = radialTable;

    // Anything cached was computed the other way, so force a recomputation.
    currentBodyFixedPosition_ = ( Eigen::Vector3d( ) << TUDAT_NAN, TUDAT_NAN, TUDAT_NAN ).finished( );
    positionOfCachedDerivatives_ = currentBodyFixedPosition_;
    positionOfCachedSecondDerivatives_ = currentBodyFixedPosition_;
}

void EllipsoidalHarmonicsGravityCache::evaluateRadialFunctions( const double rho,
                                                                EllipsoidalHarmonicCoefficients& values,
                                                                EllipsoidalHarmonicCoefficients& derivatives ) const
{
    if( radialTable_ != nullptr && radialTable_->isRadiusInRange( rho ) )
    {
        radialTable_->evaluateWithDerivatives( rho, values, derivatives );
    }
    else
    {
        computeRadialFunctionsAndDerivatives( lameFunctionCoefficients_,
                                              rho,
                                              firstFocalParameter_,
                                              secondFocalParameter_,
                                              quadratureNodes_,
                                              quadratureWeights_,
                                              values,
                                              derivatives );
    }
}

void EllipsoidalHarmonicsGravityCache::update( const Eigen::Vector3d& currentBodyFixedPosition )
{
    if( currentBodyFixedPosition != currentBodyFixedPosition_ )
    {
        currentBodyFixedPosition_ = currentBodyFixedPosition;

        ellipsoidalCoordinates_ =
                convertCartesianToEllipsoidalCoordinates( currentBodyFixedPosition, firstFocalParameter_, secondFocalParameter_ );

        if( radialTable_ != nullptr && radialTable_->isRadiusInRange( ellipsoidalCoordinates_( 0 ) ) )
        {
            radialTable_->evaluate( ellipsoidalCoordinates_( 0 ), radialFunctions_ );
        }
        else
        {
            computeRadialFunctions( lameFunctionCoefficients_,
                                    ellipsoidalCoordinates_( 0 ),
                                    firstFocalParameter_,
                                    secondFocalParameter_,
                                    quadratureNodes_,
                                    quadratureWeights_,
                                    radialFunctions_ );
        }

        legendreCache_.update( std::sin( ellipsoidalCoordinates_( 1 ) ) );
    }
}

void EllipsoidalHarmonicsGravityCache::updateIncludingDerivatives( const Eigen::Vector3d& currentBodyFixedPosition )
{
    update( currentBodyFixedPosition );

    if( currentBodyFixedPosition != positionOfCachedDerivatives_ )
    {
        positionOfCachedDerivatives_ = currentBodyFixedPosition;

        // This recomputes the values as well, but it is the same work: the radial parts and
        // their derivatives share the one integral that dominates the cost.
        evaluateRadialFunctions( ellipsoidalCoordinates_( 0 ), radialFunctions_, radialFunctionDerivatives_ );

        legendreCacheWithDerivatives_.update( std::sin( ellipsoidalCoordinates_( 1 ) ) );
    }
}

void EllipsoidalHarmonicsGravityCache::updateIncludingSecondDerivatives( const Eigen::Vector3d& currentBodyFixedPosition )
{
    update( currentBodyFixedPosition );

    if( currentBodyFixedPosition != positionOfCachedSecondDerivatives_ )
    {
        positionOfCachedSecondDerivatives_ = currentBodyFixedPosition;

        // The first derivatives are computed here too, so they count as cached. Both Legendre
        // caches have to be refreshed for that to be true, otherwise a later request for the
        // gradient would skip its own update and read values that were never filled in.
        positionOfCachedDerivatives_ = currentBodyFixedPosition;

        // The first derivatives come from the same route as the gradient uses, and the second
        // ones by differencing them, so a Hessian never leaves behind values that disagree with
        // what the gradient would have produced at the same position.
        const double rho = ellipsoidalCoordinates_( 0 );
        evaluateRadialFunctions( rho, radialFunctions_, radialFunctionDerivatives_ );

        const double step = rho * 1.0E-6;
        evaluateRadialFunctions( rho + step, forwardRadialFunctions_, forwardRadialFunctionDerivatives_ );
        evaluateRadialFunctions( rho - step, backwardRadialFunctions_, backwardRadialFunctionDerivatives_ );

        for( int degree = 0; degree <= lameFunctionCoefficients_.getMaximumDegree( ); degree++ )
        {
            for( int groupIndex = 0; groupIndex < 4; groupIndex++ )
            {
                const EllipsoidalHarmonicTermGroup group = static_cast< EllipsoidalHarmonicTermGroup >( groupIndex );
                radialFunctionSecondDerivatives_.getGroup( degree, group ) =
                        ( forwardRadialFunctionDerivatives_.getGroup( degree, group ) -
                          backwardRadialFunctionDerivatives_.getGroup( degree, group ) ) /
                        ( 2.0 * step );
            }
        }

        legendreCacheWithSecondDerivatives_.update( std::sin( ellipsoidalCoordinates_( 1 ) ) );
        legendreCacheWithDerivatives_.update( std::sin( ellipsoidalCoordinates_( 1 ) ) );
    }
}

EllipsoidalHarmonicsGravityField::EllipsoidalHarmonicsGravityField( const Eigen::Vector3d& referenceEllipsoidSemiAxes,
                                                                    const EllipsoidalHarmonicCoefficients& coefficients,
                                                                    const EllipsoidalHarmonicMatrices& lameFunctionCoefficients,
                                                                    const EllipsoidalHarmonicMatrices& transformFactors,
                                                                    const int quadratureOrder,
                                                                    const int radialTableGridPoints,
                                                                    const double radialTableLargestRadius,
                                                                    const std::string& fixedReferenceFrame,
                                                                    const std::function< void( ) > updateInertiaTensor ):
    GravityFieldModel( TUDAT_NAN, updateInertiaTensor ), referenceEllipsoidSemiAxes_( referenceEllipsoidSemiAxes ),
    firstFocalParameter_( std::sqrt( referenceEllipsoidSemiAxes( 0 ) * referenceEllipsoidSemiAxes( 0 ) -
                                     referenceEllipsoidSemiAxes( 1 ) * referenceEllipsoidSemiAxes( 1 ) ) ),
    secondFocalParameter_( std::sqrt( referenceEllipsoidSemiAxes( 0 ) * referenceEllipsoidSemiAxes( 0 ) -
                                      referenceEllipsoidSemiAxes( 2 ) * referenceEllipsoidSemiAxes( 2 ) ) ),
    maximumDegree_( coefficients.getMaximumDegree( ) ), quadratureOrder_( quadratureOrder ),
    lameFunctionCoefficients_( lameFunctionCoefficients ), transformFactors_( transformFactors ),
    fixedReferenceFrame_( fixedReferenceFrame )
{
    // The semi-axes must be strictly ordered, otherwise the focal parameters are not real and
    // the coordinates this model uses are not defined.
    if( !( referenceEllipsoidSemiAxes( 0 ) > referenceEllipsoidSemiAxes( 1 ) &&
           referenceEllipsoidSemiAxes( 1 ) > referenceEllipsoidSemiAxes( 2 ) && referenceEllipsoidSemiAxes( 2 ) > 0.0 ) )
    {
        throw std::runtime_error(
                "Error when creating ellipsoidal harmonics gravity field, the reference ellipsoid "
                "semi-axes must be given longest first and strictly decreasing, but they are " +
                std::to_string( referenceEllipsoidSemiAxes( 0 ) ) + ", " + std::to_string( referenceEllipsoidSemiAxes( 1 ) ) + ", " +
                std::to_string( referenceEllipsoidSemiAxes( 2 ) ) + "." );
    }

    if( lameFunctionCoefficients.getMaximumDegree( ) != maximumDegree_ || transformFactors.getMaximumDegree( ) != maximumDegree_ )
    {
        throw std::runtime_error(
                "Error when creating ellipsoidal harmonics gravity field, the coefficients, the Lame "
                "polynomials and the transform factors must all have the same maximum degree, but "
                "they have " +
                std::to_string( maximumDegree_ ) + ", " + std::to_string( lameFunctionCoefficients.getMaximumDegree( ) ) + " and " +
                std::to_string( transformFactors.getMaximumDegree( ) ) + "." );
    }

    gravityCache_ = std::make_shared< EllipsoidalHarmonicsGravityCache >(
            lameFunctionCoefficients_, firstFocalParameter_, secondFocalParameter_, quadratureOrder_ );

    // The coefficients are supplied normalised to the surface of the reference ellipsoid, so the
    // radial parts are evaluated there once and divided out.
    Eigen::VectorXd quadratureNodes;
    Eigen::VectorXd quadratureWeights;
    computeGaussLegendreNodesAndWeights( quadratureOrder_, quadratureNodes, quadratureWeights );

    EllipsoidalHarmonicCoefficients radialFunctionsAtSurface = createEllipsoidalHarmonicCoefficients( maximumDegree_ );
    computeRadialFunctions( lameFunctionCoefficients_,
                            referenceEllipsoidSemiAxes_( 0 ),
                            firstFocalParameter_,
                            secondFocalParameter_,
                            quadratureNodes,
                            quadratureWeights,
                            radialFunctionsAtSurface );

    coefficients_ = normaliseCoefficientsToReferenceSurface( coefficients, radialFunctionsAtSurface );

    // Far from the body only the degree zero term survives, and its radial part tends to one over
    // the ellipsoidal radius. The normalised degree zero coefficient is therefore exactly the
    // gravitational parameter.
    resetGravitationalParameter( coefficients_.kTerms.at( 0 )( 0 ) );

    // Tabulating the radial parts against the ellipsoidal radius is what makes repeated
    // evaluation cheap; outside the tabulated range they are integrated exactly instead.
    if( radialTableGridPoints > 0 )
    {
        const double largestRadius =
                ( radialTableLargestRadius > 0.0 ) ? radialTableLargestRadius : 50.0 * referenceEllipsoidSemiAxes_( 0 );

        radialTable_ = std::make_shared< EllipsoidalHarmonicRadialTable >( lameFunctionCoefficients_,
                                                                           referenceEllipsoidSemiAxes_( 0 ),
                                                                           largestRadius,
                                                                           firstFocalParameter_,
                                                                           secondFocalParameter_,
                                                                           quadratureOrder_,
                                                                           radialTableGridPoints );
        gravityCache_->setRadialTable( radialTable_ );
    }
}

void EllipsoidalHarmonicsGravityField::checkPositionIsOutsideReferenceEllipsoid( const Eigen::Vector3d& bodyFixedPosition,
                                                                                 const double ellipsoidalRadius ) const
{
    if( ellipsoidalRadius < referenceEllipsoidSemiAxes_( 0 ) )
    {
        throw std::runtime_error( "Error when evaluating ellipsoidal harmonics gravity field, the position (" +
                                  std::to_string( bodyFixedPosition( 0 ) ) + ", " + std::to_string( bodyFixedPosition( 1 ) ) + ", " +
                                  std::to_string( bodyFixedPosition( 2 ) ) +
                                  ") lies inside the reference ellipsoid, where the expansion does not converge. Its "
                                  "ellipsoidal radius is " +
                                  std::to_string( ellipsoidalRadius ) + " and the reference ellipsoid reaches " +
                                  std::to_string( referenceEllipsoidSemiAxes_( 0 ) ) + "." );
    }
}

double EllipsoidalHarmonicsGravityField::getGravitationalPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    gravityCache_->update( bodyFixedPosition );
    checkPositionIsOutsideReferenceEllipsoid( bodyFixedPosition, gravityCache_->getEllipsoidalRadius( ) );

    return sumEllipsoidalHarmonicTerms( coefficients_,
                                        gravityCache_->getRadialFunctions( ),
                                        transformFactors_,
                                        gravityCache_->getLegendreCache( ),
                                        gravityCache_->getAngleAboutLongAxis( ) );
}

Eigen::Vector3d EllipsoidalHarmonicsGravityField::computeGradientByDifferencingPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    // Central differences, with a step chosen to balance the truncation error against the
    // rounding error of the potential.
    const double step = bodyFixedPosition.norm( ) * 1.0E-5;

    Eigen::Vector3d gradient;
    for( int i = 0; i < 3; i++ )
    {
        Eigen::Vector3d forwardPosition = bodyFixedPosition;
        Eigen::Vector3d backwardPosition = bodyFixedPosition;
        forwardPosition( i ) += step;
        backwardPosition( i ) -= step;

        gradient( i ) = ( getGravitationalPotential( forwardPosition ) - getGravitationalPotential( backwardPosition ) ) / ( 2.0 * step );
    }

    return gradient;
}

Eigen::Vector3d EllipsoidalHarmonicsGravityField::getGradientOfPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    // The potential cache works everywhere outside the reference ellipsoid, so it is safe to
    // consult before deciding which route to take.
    gravityCache_->update( bodyFixedPosition );
    checkPositionIsOutsideReferenceEllipsoid( bodyFixedPosition, gravityCache_->getEllipsoidalRadius( ) );

    // Close to the long axis the ellipsoidal coordinates are degenerate and the analytic
    // gradient loses accuracy, so the potential is differenced instead.
    if( std::fabs( std::cos( gravityCache_->getEllipsoidalCoordinates( )( 1 ) ) ) < longAxisToleranceForGradient_ )
    {
        return computeGradientByDifferencingPotential( bodyFixedPosition );
    }

    gravityCache_->updateIncludingDerivatives( bodyFixedPosition );
    checkPositionIsOutsideReferenceEllipsoid( bodyFixedPosition, gravityCache_->getEllipsoidalRadius( ) );

    return computeEllipsoidalHarmonicPotentialGradient( coefficients_,
                                                        gravityCache_->getRadialFunctions( ),
                                                        gravityCache_->getRadialFunctionDerivatives( ),
                                                        transformFactors_,
                                                        gravityCache_->getLegendreCacheWithDerivatives( ),
                                                        gravityCache_->getEllipsoidalCoordinates( ),
                                                        firstFocalParameter_,
                                                        secondFocalParameter_ );
}

Eigen::Matrix3d EllipsoidalHarmonicsGravityField::computeHessianByDifferencingGradient( const Eigen::Vector3d& bodyFixedPosition )
{
    const double step = bodyFixedPosition.norm( ) * 1.0E-5;

    Eigen::Matrix3d hessian;
    for( int i = 0; i < 3; i++ )
    {
        Eigen::Vector3d forwardPosition = bodyFixedPosition;
        Eigen::Vector3d backwardPosition = bodyFixedPosition;
        forwardPosition( i ) += step;
        backwardPosition( i ) -= step;

        hessian.col( i ) = ( getGradientOfPotential( forwardPosition ) - getGradientOfPotential( backwardPosition ) ) / ( 2.0 * step );
    }

    // Differencing column by column leaves the result very slightly asymmetric, which the
    // Hessian of a potential cannot be.
    return 0.5 * ( hessian + hessian.transpose( ) );
}

Eigen::Matrix3d EllipsoidalHarmonicsGravityField::getHessianOfPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    gravityCache_->update( bodyFixedPosition );
    checkPositionIsOutsideReferenceEllipsoid( bodyFixedPosition, gravityCache_->getEllipsoidalRadius( ) );

    // Same degeneracy as for the gradient, and worse: the second derivatives of the Legendre
    // functions grow faster still towards the long axis.
    if( std::fabs( std::cos( gravityCache_->getEllipsoidalCoordinates( )( 1 ) ) ) < longAxisToleranceForGradient_ )
    {
        return computeHessianByDifferencingGradient( bodyFixedPosition );
    }

    gravityCache_->updateIncludingSecondDerivatives( bodyFixedPosition );

    return computeEllipsoidalHarmonicPotentialHessian( coefficients_,
                                                       gravityCache_->getRadialFunctions( ),
                                                       gravityCache_->getRadialFunctionDerivatives( ),
                                                       gravityCache_->getRadialFunctionSecondDerivatives( ),
                                                       transformFactors_,
                                                       gravityCache_->getLegendreCacheWithSecondDerivatives( ),
                                                       gravityCache_->getEllipsoidalCoordinates( ),
                                                       firstFocalParameter_,
                                                       secondFocalParameter_ );
}

double EllipsoidalHarmonicsGravityField::getLaplacianOfPotential( const Eigen::Vector3d& bodyFixedPosition )
{
    return getHessianOfPotential( bodyFixedPosition ).trace( );
}

}  // namespace gravitation

}  // namespace tudat
