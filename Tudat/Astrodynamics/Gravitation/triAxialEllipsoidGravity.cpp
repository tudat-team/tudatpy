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

#include <boost/math/special_functions/factorials.hpp>

#include "Tudat/Mathematics/BasicMathematics/legendrePolynomials.h"
#include "Tudat/Astrodynamics/Gravitation/triAxialEllipsoidGravity.h"

namespace tudat
{

namespace gravitation
{

//! Function to calculate (non-normalized) cosine spherical harmonic coefficient for a homogeneous
//! triaxial ellipsoid
double calculateCosineTermForTriaxialEllipsoidSphericalHarmonicGravity(
        const double aSquaredMinusCSquared, const double bSquaredMinusCSquared,
        const double referenceRadius, const int degree, const int order )
{
    // Initialize coefficient to zero
    double powerSeries = 0.0;

    // Only non-zero terms are for even degree and order
    if( !( ( degree % 2 != 0 || order % 2 != 0 ) ) )
    {
        using boost::math::factorial;

        // Calculate indices for algorithm
        int l = degree / 2;
        int m = order / 2;
        int maximumIndex = std::floor( ( l - m ) / 2 );

        // Evaluate single term in power series of final equation of Boyce (1997)
        for( int i = 0; i <= maximumIndex; i++ )
        {
            powerSeries +=
                    ( std::pow( ( aSquaredMinusCSquared - bSquaredMinusCSquared ) /
                                2.0, m + 2 * i ) *
                      std::pow( -( aSquaredMinusCSquared + bSquaredMinusCSquared ) /
                                2.0, l - m - 2 * i ) ) /
                    ( std::pow( 2.0, m + 2 * i ) * factorial< double >( l - m - 2 * i ) *
                      factorial< double >( m + i ) * factorial< double >( i ) );
        }

        // Calculate multiplier of power series in final equation of Boyce (1997)
        double multiplier = 3.0 / std::pow( referenceRadius, 2 * l ) *
                ( factorial< double >( l ) * factorial< double >( 2 * l - 2 * m ) ) /
                ( ( 2.0 * static_cast< double >( l ) + 3.0 ) * factorial< double >( 2 * l + 1 ) );
        if( order != 0 )
        {
            multiplier *= 2.0;
        }

        // Complete calculation of coefficient
        powerSeries *= multiplier;
    }
    return powerSeries;
}

//! Function to calculate triaxial ellipsoid reference radius
double calculateTriAxialEllipsoidReferenceRadius(
        const double axisA, const double axisB, const double axisC )
{
    return std::sqrt( 3.0 / ( 1.0 / ( axisA * axisA ) + 1.0 / ( axisB * axisB ) +
                              1.0 / ( axisC * axisC ) ) );
}

//! Function to calculate triaxial ellipsoid volume
double calculateTriAxialEllipsoidVolume(
        const double axisA, const double axisB, const double axisC )
{
    return 4.0 / 3.0 * mathematical_constants::PI * axisA * axisB * axisC;
}

//! Function to calculate (non-normalized) cosine spherical harmonic coefficients for a
//! homogeneous triaxial ellipsoid
Eigen::MatrixXd createTriAxialEllipsoidSphericalHarmonicCosineCoefficients(
        const double axisA, const double axisB, const double axisC,
        const int maximumDegree, const int maximumOrder )
{
    // Initialize vector to zeros
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd::Zero(
                maximumDegree + 1, maximumOrder + 1 );

    // Pre-calculate for single coefficient calculations.
    double aSquaredMinusCSquared = axisA * axisA - axisC * axisC;
    double bSquaredMinusCSquared = axisB * axisB - axisC * axisC;
    double referenceRadius =  calculateTriAxialEllipsoidReferenceRadius( axisA, axisB, axisC );

    // Iterate over all requested degrees.
    for( int i = 0; i <= maximumDegree; i++ )
    {
        // Iterate over all requested orders.
        for( int j = 0; ( ( j <= maximumOrder ) && ( j <= i ) ); j++ )
        {
            // Only even degree and order terms are non-zero
            if( ( i % 2 == 0 ) && ( j % 2 == 0 ) )
            {
                // Calculate coefficient at single degree and order.
                cosineCoefficients( i, j ) =
                        calculateCosineTermForTriaxialEllipsoidSphericalHarmonicGravity(
                            aSquaredMinusCSquared, bSquaredMinusCSquared, referenceRadius, i, j );
            }
        }
    }

    return cosineCoefficients;
}

//! Function to calculate (non-normalized) cosine and sine spherical harmonic coefficients for a
//! homogeneous triaxial ellipsoid
std::pair< Eigen::MatrixXd, Eigen::MatrixXd > createTriAxialEllipsoidSphericalHarmonicCoefficients(
        const double axisA, const double axisB, const double axisC,
        const int maximumDegree, const int maximumOrder )
{
    // Calculate cosine coefficients and add sine matrix (all zeroes)
    return std::make_pair( createTriAxialEllipsoidSphericalHarmonicCosineCoefficients(
                               axisA, axisB, axisC, maximumDegree, maximumOrder ),
                           Eigen::MatrixXd::Zero( maximumDegree + 1, maximumOrder + 1 ) );
}

//! Function to calculate (normalized) cosine and sine spherical harmonic coefficients for a
//! homogeneous triaxial ellipsoid
std::pair< Eigen::MatrixXd, Eigen::MatrixXd >
createTriAxialEllipsoidNormalizedSphericalHarmonicCoefficients(
        const double axisA, const double axisB, const double axisC,
        const int maximumDegree, const int maximumOrder )
{
    // Calculate non-normalized coefficients
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > unNormalizedCoefficients =
            createTriAxialEllipsoidSphericalHarmonicCoefficients(
                axisA, axisB, axisC, maximumDegree, maximumOrder );

    // Retrieve cosine coefficients
    Eigen::MatrixXd normalizedCosineCoefficients = unNormalizedCoefficients.first;

    // Iterate over all degrees and orders and normalized coefficients
    for( int i = 2; i < normalizedCosineCoefficients.rows( ); i++ )
    {
        for( int j = 0; ( ( j < normalizedCosineCoefficients.cols( ) ) &&
                                   ( j <= i ) ); j++ )
        {
            normalizedCosineCoefficients( i, j ) = normalizedCosineCoefficients( i, j ) /
                    basic_mathematics::calculateLegendreGeodesyNormalizationFactor( i, j );
        }
    }

    // Return geodesy-normalized coefficients.
    return std::make_pair( normalizedCosineCoefficients, unNormalizedCoefficients.second );
}

}

}
