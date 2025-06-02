/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/math/basic/basicMathematicsFunctions.h"
#include "tudat/math/basic/coordinateConversions.h"
#include "tudat/astro/basic_astro/stateVectorIndices.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicPartialFunctions.h"

namespace tudat
{

namespace acceleration_partials
{

using namespace orbital_element_conversions;

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
void computePotentialSphericalHessian( const double distance,
                                       const double radiusPowerTerm,
                                       const double cosineOfOrderLongitude,
                                       const double sineOfOrderLongitude,
                                       const double cosineOfLatitude,
                                       const double sineOfLatitude,
                                       const double preMultiplier,
                                       const int degree,
                                       const int order,
                                       const double cosineHarmonicCoefficient,
                                       const double sineHarmonicCoefficient,
                                       const double legendrePolynomial,
                                       const double legendrePolynomialDerivative,
                                       const double legendrePolynomialSecondDerivative,
                                       Eigen::Matrix3d& sphericalHessian )
{
    const double degreePlusOne = static_cast<double>(degree + 1);
    const double degreePlusTwo = degreePlusOne + 1.0;
    const double orderAsDouble = static_cast<double>(order);
    const double inverseDistance = 1.0 / distance;
    const double inverseDistanceSquared = inverseDistance * inverseDistance;

    const double combinedHarmonicSum = cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude;
    const double combinedHarmonicDifference = sineHarmonicCoefficient * cosineOfOrderLongitude - cosineHarmonicCoefficient * sineOfOrderLongitude;

    const double legendreTimesHarmonicSum = legendrePolynomial * combinedHarmonicSum;
    const double firstDerivativeTimesHarmonicSum = legendrePolynomialDerivative * combinedHarmonicSum;
    const double secondDerivativeTimesHarmonicSum = legendrePolynomialSecondDerivative * combinedHarmonicSum;
    const double legendreTimesHarmonicDifference = legendrePolynomial * combinedHarmonicDifference;
    const double firstDerivativeTimesHarmonicDifference = legendrePolynomialDerivative * combinedHarmonicDifference;

    sphericalHessian(0, 0) = degreePlusOne * degreePlusTwo * inverseDistanceSquared * legendreTimesHarmonicSum;
    sphericalHessian(1, 0) = -degreePlusOne * inverseDistance * cosineOfLatitude * firstDerivativeTimesHarmonicSum;
    sphericalHessian(0, 1) = sphericalHessian(1, 0);
    sphericalHessian(2, 0) = -degreePlusOne * orderAsDouble * inverseDistance * legendreTimesHarmonicDifference;
    sphericalHessian(0, 2) = sphericalHessian(2, 0);
    sphericalHessian(1, 1) = (cosineOfLatitude * cosineOfLatitude * secondDerivativeTimesHarmonicSum - sineOfLatitude * firstDerivativeTimesHarmonicSum);
    sphericalHessian(2, 1) = orderAsDouble * cosineOfLatitude * firstDerivativeTimesHarmonicDifference;
    sphericalHessian(1, 2) = sphericalHessian(2, 1);
    sphericalHessian(2, 2) = -orderAsDouble * orderAsDouble * legendreTimesHarmonicSum;

    sphericalHessian *= preMultiplier * radiusPowerTerm;
}

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
void computePotentialSphericalHessian( const Eigen::Vector3d& sphericalPosition,
                                       const double referenceRadius,
                                       const double preMultiplier,
                                       const int degree,
                                       const int order,
                                       const double cosineHarmonicCoefficient,
                                       const double sineHarmonicCoefficient,
                                       const double legendrePolynomial,
                                       const double legendrePolynomialDerivative,
                                       const double legendrePolynomialSecondDerivative,
                                       Eigen::Matrix3d& sphericalHessian )
{
    computePotentialSphericalHessian( sphericalPosition( radiusIndex ),
                                      basic_mathematics::raiseToIntegerPower( referenceRadius / sphericalPosition( radiusIndex ),
                                                                              static_cast< double >( degree ) + 1.0 ),
                                      std::cos( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                                      std::sin( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                                      std::cos( sphericalPosition( latitudeIndex ) ),
                                      std::sin( sphericalPosition( latitudeIndex ) ),
                                      preMultiplier,
                                      degree,
                                      order,
                                      cosineHarmonicCoefficient,
                                      sineHarmonicCoefficient,
                                      legendrePolynomial,
                                      legendrePolynomialDerivative,
                                      legendrePolynomialSecondDerivative,
                                      sphericalHessian );
}

//! Function to compute the spherical Hessian of a single term of a spherical harmonic potential
void computePotentialSphericalHessian( const Eigen::Vector3d& sphericalPosition,
                                       const double preMultiplier,
                                       const int degree,
                                       const int order,
                                       const double cosineHarmonicCoefficient,
                                       const double sineHarmonicCoefficient,
                                       const basic_mathematics::SphericalHarmonicsCache& sphericalHarmonicsCache,
                                       Eigen::Matrix3d& sphericalHessian,
                                       const bool checkSphericalHarmonicsConsistency )
{
    if( checkSphericalHarmonicsConsistency )
    {
        computePotentialSphericalHessian( sphericalPosition( 0 ),
                                          sphericalHarmonicsCache.getReferenceRadiusRatioPowers( degree + 1 ),
                                          sphericalHarmonicsCache.getCosineOfMultipleLongitude( order ),
                                          sphericalHarmonicsCache.getSineOfMultipleLongitude( order ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getCurrentPolynomialParameterComplement( ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getCurrentPolynomialParameter( ),
                                          preMultiplier,
                                          degree,
                                          order,
                                          cosineHarmonicCoefficient,
                                          sineHarmonicCoefficient,
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomial( degree, order ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomialDerivative( degree, order ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomialSecondDerivative( degree, order ),
                                          sphericalHessian );
    }
    else
    {
        computePotentialSphericalHessian( sphericalPosition( 0 ),
                                          sphericalHarmonicsCache.getReferenceRadiusRatioPowers( degree + 1 ),
                                          sphericalHarmonicsCache.getCosineOfMultipleLongitude( order ),
                                          sphericalHarmonicsCache.getSineOfMultipleLongitude( order ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getCurrentPolynomialParameterComplement( ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getCurrentPolynomialParameter( ),
                                          preMultiplier,
                                          degree,
                                          order,
                                          cosineHarmonicCoefficient,
                                          sineHarmonicCoefficient,
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomialWithoutCheck( degree, order ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomialDerivativeWithoutCheck( degree, order ),
                                          sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomialSecondDerivativeWithoutCheck( degree, order ),
                                          sphericalHessian );
    }
}


//! Calculate partial of spherical harmonic acceleration w.r.t. a set of cosine coefficients
void calculateSphericalHarmonicGravityWrtCCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const basic_mathematics::SphericalHarmonicsCache& sphericalHarmonicsCache,
        const std::vector< std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix,
        const int maximumAccelerationDegree,
        const int maximumAccelerationOrder )
{
    double preMultiplier = gravitionalParameter / referenceRadius;
    const basic_mathematics::LegendreCache& legendreCache = sphericalHarmonicsCache.getLegendreCacheConst( );

    int degree, order;
    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        degree = blockIndices.at( i ).first;
        order = blockIndices.at( i ).second;

        // Calculate and set partial of current degree and order.
        if( degree <= maximumAccelerationDegree && order <= maximumAccelerationOrder )
        {
            partialsMatrix.block( 0, i, 3, 1 ) = basic_mathematics::computePotentialGradient(
                    sphericalPosition( radiusIndex ),
                    sphericalHarmonicsCache.getReferenceRadiusRatioPowers( degree + 1 ),
                    sphericalHarmonicsCache.getCosineOfMultipleLongitude( order ),
                    sphericalHarmonicsCache.getSineOfMultipleLongitude( order ),
                    legendreCache.getCurrentPolynomialParameterComplement( ),
                    preMultiplier,
                    degree,
                    order,
                    1.0,
                    0.0,
                    legendreCache.getLegendrePolynomial( degree, order ),
                    legendreCache.getLegendrePolynomialDerivative( degree, order ) );
        }
        else
        {
            partialsMatrix.block( 0, i, 3, 1 ).setZero( );
        }
    }

    // Transform partials to Cartesian position and integration frame.
    partialsMatrix = bodyFixedToIntegrationFrame * sphericalToCartesianGradientMatrix * partialsMatrix;
}

//! Calculate partial of spherical harmonic acceleration w.r.t. a set of sine coefficients
void calculateSphericalHarmonicGravityWrtSCoefficients(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double gravitionalParameter,
        const basic_mathematics::SphericalHarmonicsCache& sphericalHarmonicsCache,
        const std::vector< std::pair< int, int > >& blockIndices,
        const Eigen::Matrix3d& sphericalToCartesianGradientMatrix,
        const Eigen::Matrix3d& bodyFixedToIntegrationFrame,
        Eigen::MatrixXd& partialsMatrix,
        const int maximumAccelerationDegree,
        const int maximumAccelerationOrder )
{
    double preMultiplier = gravitionalParameter / referenceRadius;
    const basic_mathematics::LegendreCache& legendreCache = sphericalHarmonicsCache.getLegendreCacheConst( );

    int degree, order;
    for( unsigned int i = 0; i < blockIndices.size( ); i++ )
    {
        degree = blockIndices.at( i ).first;
        order = blockIndices.at( i ).second;

        // Calculate and set partial of current degree and order.
        if( degree <= maximumAccelerationDegree && order <= maximumAccelerationOrder )
        {
            partialsMatrix.block( 0, i, 3, 1 ) = basic_mathematics::computePotentialGradient(
                    sphericalPosition( radiusIndex ),
                    sphericalHarmonicsCache.getReferenceRadiusRatioPowers( degree + 1 ),
                    sphericalHarmonicsCache.getCosineOfMultipleLongitude( order ),
                    sphericalHarmonicsCache.getSineOfMultipleLongitude( order ),
                    legendreCache.getCurrentPolynomialParameterComplement( ),
                    preMultiplier,
                    degree,
                    order,
                    0.0,
                    1.0,
                    legendreCache.getLegendrePolynomial( degree, order ),
                    legendreCache.getLegendrePolynomialDerivative( degree, order ) );
        }
        else
        {
            partialsMatrix.block( 0, i, 3, 1 ).setZero( );
        }
    }

    // Transform partials to Cartesian position and integration frame.
    partialsMatrix = bodyFixedToIntegrationFrame * sphericalToCartesianGradientMatrix * partialsMatrix;
}

}  // namespace acceleration_partials

}  // namespace tudat
