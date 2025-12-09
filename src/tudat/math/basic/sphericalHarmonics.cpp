/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include <cmath>

#include <Eigen/Core>

#include "tudat/math/basic/sphericalHarmonics.h"
#include "tudat/math/basic/basicMathematicsFunctions.h"

namespace tudat
{
namespace basic_mathematics
{

//! Update maximum degree and order of cache
void SphericalHarmonicsCache::resetMaximumDegreeAndOrder( const int maximumDegree, const int maximumOrder )
{
    maximumDegree_ = maximumDegree;
    maximumOrder_ = maximumOrder;

    if( maximumOrder_ > maximumDegree_ )
    {
        maximumOrder_ = maximumDegree_;
    }
    legendreCache_.resetMaximumDegreeAndOrder( maximumDegree_, maximumOrder_ );

    sinesOfLongitude_.resize( maximumOrder_ + 1 );
    cosinesOfLongitude_.resize( maximumOrder_ + 1 );
    referenceRadiusRatioPowers_.resize( maximumDegree_ + 2 );
}

//! Compute the gradient of a single term of a spherical harmonics potential field.
Eigen::Vector3d computePotentialGradientFromLatitudeDerivative( const double distance,
                                                                const double radiusPowerTerm,
                                                                const double cosineOfOrderLongitude,
                                                                const double sineOfOrderLongitude,
                                                                const double preMultiplier,
                                                                const int degree,
                                                                const int order,
                                                                const double cosineHarmonicCoefficient,
                                                                const double sineHarmonicCoefficient,
                                                                const double legendrePolynomial,
                                                                const double legendrePolynomialDerivativeWrLatitude )
{
    const double harmonicSum = cosineHarmonicCoefficient * cosineOfOrderLongitude + sineHarmonicCoefficient * sineOfOrderLongitude;
    const double harmonicDifference = sineHarmonicCoefficient * cosineOfOrderLongitude - cosineHarmonicCoefficient * sineOfOrderLongitude;
    const double factor = preMultiplier * radiusPowerTerm;
    const double degreePlusOne = static_cast< double >( degree ) + 1.0;
    const double orderDouble = static_cast< double >( order );

    Eigen::Vector3d gradient;
    gradient.x( ) = -1.0 / distance * degreePlusOne * legendrePolynomial * harmonicSum;
    gradient.y( ) = legendrePolynomialDerivativeWrLatitude * harmonicSum;
    gradient.z( ) = orderDouble * legendrePolynomial * harmonicDifference;
    return factor * gradient;
}

//! Compute the gradient of a single term of a spherical harmonics potential field.
Eigen::Vector3d computePotentialGradient( const double distance,
                                          const double radiusPowerTerm,
                                          const double cosineOfOrderLongitude,
                                          const double sineOfOrderLongitude,
                                          const double cosineOfLatitude,
                                          const double preMultiplier,
                                          const int degree,
                                          const int order,
                                          const double cosineHarmonicCoefficient,
                                          const double sineHarmonicCoefficient,
                                          const double legendrePolynomial,
                                          const double legendrePolynomialDerivative )
{
    return computePotentialGradientFromLatitudeDerivative( distance,
                                                           radiusPowerTerm,
                                                           cosineOfOrderLongitude,
                                                           sineOfOrderLongitude,
                                                           preMultiplier,
                                                           degree,
                                                           order,
                                                           cosineHarmonicCoefficient,
                                                           sineHarmonicCoefficient,
                                                           legendrePolynomial,
                                                           legendrePolynomialDerivative * cosineOfLatitude );
}

//! Compute the gradient of a single term of a spherical harmonics potential field at the poles.
Eigen::Vector3d computePotentialGradientNearPole( const Eigen::Vector3d& sphericalPosition,
                                                  const double preMultiplier,
                                                  const int degree,
                                                  const int order,
                                                  const double cosineHarmonicCoefficient,
                                                  const double sineHarmonicCoefficient,
                                                  const SphericalHarmonicsCache& sphericalHarmonicsCache )
{
    double legendrePolynomialDerivativeWithRespectToColatitude = 0.0;

    if( order == 0 )
    {
        const double incrementedOrderLegendrePolynomial =
                sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomial( degree, 1 );
        legendrePolynomialDerivativeWithRespectToColatitude =
                -std::sqrt( static_cast< double >( degree ) * ( static_cast< double >( degree ) + 1.0 ) / 2.0 ) *
                incrementedOrderLegendrePolynomial;
    }
    else if( order == degree )
    {
        const double decrementedOrderLegendrePolynomial =
                sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomial( degree, degree - 1 );
        legendrePolynomialDerivativeWithRespectToColatitude =
                std::sqrt( static_cast< double >( degree ) / 2.0 ) * decrementedOrderLegendrePolynomial;
        if( degree == 1 )
        {
            legendrePolynomialDerivativeWithRespectToColatitude *= sqrt( 2.0 );
        }
    }
    else
    {
        const double decrementedOrderLegendrePolynomial =
                sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomial( degree, order - 1 );
        const double incrementedOrderLegendrePolynomial =
                sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomial( degree, order + 1 );
        double persistenceTerm = TUDAT_NAN;
        double decayTerm = TUDAT_NAN;

        if( order == 1 )
        {
            persistenceTerm = std::sqrt( 2.0 * static_cast< double >( degree ) * ( static_cast< double >( degree ) + 1.0 ) );
            decayTerm = std::sqrt( static_cast< double >( degree - 1 ) * static_cast< double >( degree + 2 ) );
        }
        else
        {
            persistenceTerm = std::sqrt( static_cast< double >( degree + order ) * static_cast< double >( degree - order + 1 ) );
            decayTerm = std::sqrt( static_cast< double >( degree + order + 1 ) * static_cast< double >( degree - order ) );
        }

        legendrePolynomialDerivativeWithRespectToColatitude =
                0.5 * ( persistenceTerm * decrementedOrderLegendrePolynomial - decayTerm * incrementedOrderLegendrePolynomial );
    }

    const double distance = sphericalPosition( radiusIndex );
    const double radiusPowerTerm = sphericalHarmonicsCache.getReferenceRadiusRatioPowers( degree + 1 );
    const double cosineOfOrderLongitude = sphericalHarmonicsCache.getCosineOfMultipleLongitude( order );
    const double sineOfOrderLongitude = sphericalHarmonicsCache.getSineOfMultipleLongitude( order );
    const double legendrePolynomial = sphericalHarmonicsCache.getLegendreCacheConst( ).getLegendrePolynomial( degree, order );

    return computePotentialGradientFromLatitudeDerivative( distance,
                                                           radiusPowerTerm,
                                                           cosineOfOrderLongitude,
                                                           sineOfOrderLongitude,
                                                           preMultiplier,
                                                           degree,
                                                           order,
                                                           cosineHarmonicCoefficient,
                                                           sineHarmonicCoefficient,
                                                           legendrePolynomial,
                                                           -legendrePolynomialDerivativeWithRespectToColatitude );
}

//! Compute the gradient of a single term of a spherical harmonics potential field.
Eigen::Vector3d computePotentialGradient( const Eigen::Vector3d& sphericalPosition,
                                          const double referenceRadius,
                                          const double preMultiplier,
                                          const int degree,
                                          const int order,
                                          const double cosineHarmonicCoefficient,
                                          const double sineHarmonicCoefficient,
                                          const double legendrePolynomial,
                                          const double legendrePolynomialDerivative )
{
    return computePotentialGradient( sphericalPosition( radiusIndex ),
                                     basic_mathematics::raiseToIntegerPower( referenceRadius / sphericalPosition( radiusIndex ),
                                                                             static_cast< double >( degree ) + 1.0 ),
                                     std::cos( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                                     std::sin( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                                     std::cos( sphericalPosition( latitudeIndex ) ),
                                     preMultiplier,
                                     degree,
                                     order,
                                     cosineHarmonicCoefficient,
                                     sineHarmonicCoefficient,
                                     legendrePolynomial,
                                     legendrePolynomialDerivative );
}

//! Compute the gradient of a single term of a spherical harmonics potential field.
Eigen::Vector3d computePotentialGradient( const Eigen::Vector3d& sphericalPosition,
                                          const double preMultiplier,
                                          const int degree,
                                          const int order,
                                          const double cosineHarmonicCoefficient,
                                          const double sineHarmonicCoefficient,
                                          const double legendrePolynomial,
                                          const double legendrePolynomialDerivative,
                                          const SphericalHarmonicsCache& sphericalHarmonicsCache )
{
    return computePotentialGradient( sphericalPosition( radiusIndex ),
                                     sphericalHarmonicsCache.getReferenceRadiusRatioPowers( degree + 1 ),
                                     sphericalHarmonicsCache.getCosineOfMultipleLongitude( order ),
                                     sphericalHarmonicsCache.getSineOfMultipleLongitude( order ),
                                     sphericalHarmonicsCache.getLegendreCacheConst( ).getCurrentPolynomialParameterComplement( ),
                                     preMultiplier,
                                     degree,
                                     order,
                                     cosineHarmonicCoefficient,
                                     sineHarmonicCoefficient,
                                     legendrePolynomial,
                                     legendrePolynomialDerivative );
}

}  // namespace basic_mathematics
}  // namespace tudat
