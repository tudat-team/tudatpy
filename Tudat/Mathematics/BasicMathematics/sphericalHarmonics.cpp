/*    Copyright (c) 2010-2018, Delft University of Technology
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

#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

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
    legendreCache_->resetMaximumDegreeAndOrder( maximumDegree_, maximumOrder_ );

    sinesOfLongitude_.resize( maximumOrder_ + 1 );
    cosinesOfLongitude_.resize( maximumOrder_ + 1 );
    referenceRadiusRatioPowers_.resize( maximumDegree_ + 2 );
}


//! Compute the gradient of a single term of a spherical harmonics potential field.
Eigen::Vector3d computePotentialGradient(
        const double distance,
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
    // Return result.
    return ( Eigen::Vector3d( ) <<
             - preMultiplier / distance
             * radiusPowerTerm
             * ( static_cast< double >( degree ) + 1.0 ) * legendrePolynomial
             * ( cosineHarmonicCoefficient * cosineOfOrderLongitude
                 + sineHarmonicCoefficient * sineOfOrderLongitude ),
             preMultiplier * radiusPowerTerm
             * legendrePolynomialDerivative * cosineOfLatitude * (
                 cosineHarmonicCoefficient * cosineOfOrderLongitude
                 + sineHarmonicCoefficient * sineOfOrderLongitude ),
             preMultiplier * radiusPowerTerm
             * static_cast< double >( order ) * legendrePolynomial
             * ( sineHarmonicCoefficient * cosineOfOrderLongitude
                 - cosineHarmonicCoefficient * sineOfOrderLongitude ) ).finished( );
}

//! Compute the gradient of a single term of a spherical harmonics potential field.
Eigen::Vector3d computePotentialGradient(
        const Eigen::Vector3d& sphericalPosition,
        const double referenceRadius,
        const double preMultiplier,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        const double legendrePolynomial,
        const double legendrePolynomialDerivative )
{
    return computePotentialGradient(
                sphericalPosition( radiusIndex ),
                basic_mathematics::raiseToIntegerPower
                ( referenceRadius / sphericalPosition( radiusIndex ), static_cast< double >( degree ) + 1.0 ),
                std::cos( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                std::sin( static_cast< double >( order ) * sphericalPosition( longitudeIndex ) ),
                std::cos( sphericalPosition( latitudeIndex ) ), preMultiplier, degree, order,
                cosineHarmonicCoefficient, sineHarmonicCoefficient, legendrePolynomial,legendrePolynomialDerivative );
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
                                          const std::shared_ptr< SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    return computePotentialGradient(
                sphericalPosition( radiusIndex ),
                sphericalHarmonicsCache->getReferenceRadiusRatioPowers( degree + 1 ),
                sphericalHarmonicsCache->getCosineOfMultipleLongitude( order ),
                sphericalHarmonicsCache->getSineOfMultipleLongitude( order ),
                sphericalHarmonicsCache->getLegendreCache( )->getCurrentPolynomialParameterComplement( ),
                preMultiplier, degree, order,
                cosineHarmonicCoefficient, sineHarmonicCoefficient, legendrePolynomial,legendrePolynomialDerivative );
}

} // namespace basic_mathematics
} // namespace tudat
