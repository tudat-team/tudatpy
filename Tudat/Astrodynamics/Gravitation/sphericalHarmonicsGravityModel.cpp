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
#include <limits>
#include <stdexcept>

#include <boost/math/constants/constants.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/centralJ2J3GravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"
#include "Tudat/Mathematics/BasicMathematics/sphericalHarmonics.h"

namespace tudat
{
namespace gravitation
{

//! Compute gravitational acceleration due to multiple spherical harmonics terms, defined using geodesy-normalization.
Eigen::Vector3d computeGeodesyNormalizedGravitationalAccelerationSum(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameter,
        const double equatorialRadius,
        const Eigen::MatrixXd& cosineHarmonicCoefficients,
        const Eigen::MatrixXd& sineHarmonicCoefficients,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    // Set highest degree and order.
    const int highestDegree = cosineHarmonicCoefficients.rows( );
    const int highestOrder = cosineHarmonicCoefficients.cols( );

    // Declare spherical position vector.
    Eigen::Vector3d sphericalpositionOfBodySubjectToAcceleration = coordinate_conversions::
            convertCartesianToSpherical( positionOfBodySubjectToAcceleration );
    sphericalpositionOfBodySubjectToAcceleration( 1 ) = mathematical_constants::PI / 2.0 -
            sphericalpositionOfBodySubjectToAcceleration( 1 );

    double sineOfAngle = std::sin( sphericalpositionOfBodySubjectToAcceleration( 1 ) );
    sphericalHarmonicsCache->update( sphericalpositionOfBodySubjectToAcceleration( 0 ),
                                     sineOfAngle,
                                     sphericalpositionOfBodySubjectToAcceleration( 2 ),
                                     equatorialRadius );

    boost::shared_ptr< basic_mathematics::LegendreCache > legendreCacheReference =
            sphericalHarmonicsCache->getLegendreCache( );


    // Compute gradient premultiplier.
    const double preMultiplier = gravitationalParameter / equatorialRadius;

    // Initialize gradient vector.
    Eigen::Vector3d sphericalGradient = Eigen::Vector3d::Zero( );

    // Loop through all degrees.
    for ( int degree = 0; degree < highestDegree; degree++ )
    {
        // Loop through all orders.
        for ( int order = 0; ( order <= degree ) && ( order < highestOrder ); order++ )
        {
            // Compute geodesy-normalized Legendre polynomials.
            const double legendrePolynomial = legendreCacheReference->getLegendrePolynomial( degree, order );

            // Compute geodesy-normalized Legendre polynomial derivative.
            const double legendrePolynomialDerivative = legendreCacheReference->getLegendrePolynomialDerivative(
                        degree, order );

            // Compute the potential gradient of a single spherical harmonic term.
            sphericalGradient += basic_mathematics::computePotentialGradient(
                        sphericalpositionOfBodySubjectToAcceleration,
                        preMultiplier,
                        degree,
                        order,
                        cosineHarmonicCoefficients( degree, order ),
                        sineHarmonicCoefficients( degree, order ),
                        legendrePolynomial,
                        legendrePolynomialDerivative, sphericalHarmonicsCache );
        }
    }


    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector) and
    // return the resulting acceleration vector.
    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

//! Compute gravitational acceleration due to single spherical harmonics term.
Eigen::Vector3d computeSingleGeodesyNormalizedGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double gravitationalParameter,
        const double equatorialRadius,
        const int degree,
        const int order,
        const double cosineHarmonicCoefficient,
        const double sineHarmonicCoefficient,
        boost::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicsCache )
{
    // Declare spherical position vector.
    Eigen::Vector3d sphericalpositionOfBodySubjectToAcceleration = coordinate_conversions::
            convertCartesianToSpherical( positionOfBodySubjectToAcceleration );
    sphericalpositionOfBodySubjectToAcceleration( 1 ) = mathematical_constants::PI / 2.0 -
            sphericalpositionOfBodySubjectToAcceleration( 1 );


    double sineOfAngle = std::sin( sphericalpositionOfBodySubjectToAcceleration( 1 ) );
    sphericalHarmonicsCache->update( sphericalpositionOfBodySubjectToAcceleration( 0 ),
                                     sineOfAngle,
                                     sphericalpositionOfBodySubjectToAcceleration( 2 ),
                                     equatorialRadius );

    // Compute gradient premultiplier.
    const double preMultiplier = gravitationalParameter / equatorialRadius;

    // Compute geodesy-normalized Legendre polynomials.
    const double legendrePolynomial = sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomial( degree, order );

    // Compute geodesy-normalized Legendre polynomial derivative.
    const double legendrePolynomialDerivative =
            sphericalHarmonicsCache->getLegendreCache( )->getLegendrePolynomialDerivative( degree, order );

    // Compute the potential gradient of a single spherical harmonic term.
    Eigen::Vector3d sphericalGradient = basic_mathematics::computePotentialGradient(
                sphericalpositionOfBodySubjectToAcceleration,
                preMultiplier,
                degree,
                order,
                cosineHarmonicCoefficient,
                sineHarmonicCoefficient,
                legendrePolynomial,
                legendrePolynomialDerivative, sphericalHarmonicsCache );

    // Convert from spherical gradient to Cartesian gradient (which equals acceleration vector),
    // and return resulting acceleration vector.
    return coordinate_conversions::convertSphericalToCartesianGradient(
                sphericalGradient, positionOfBodySubjectToAcceleration );
}

} // namespace gravitation
} // namespace tudat
