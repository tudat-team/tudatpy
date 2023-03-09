/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "tudat/astro/gravitation/ringGravityField.h"

namespace tudat
{

namespace gravitation
{

//! Computes the gravitational potential of a one-dimensional ring.
double computeRingGravitationalPotential(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double ringRadius,
        const double gravitationalParameter,
        const double gravitationalConstant )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensity = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius * gravitationalConstant );

    double m = 4.0 * ringRadius * r / std::pow( p, 2.0 );
    // m = k^2
    double k = std::sqrt( std::abs( m ) );

    return 4.0 * lineDensity * ringRadius * boost::math::ellint_1( k ) / p;
}

//! Computes the gravitational acceleration of a one-dimensional ring
Eigen::Vector3d computeRingGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double ringRadius,
        const double gravitationalParameter,
        const double gravitationalConstant,
        const bool ellipticIntegralSFromDAndB )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );
    double q = std::sqrt( std::pow( r - ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensity = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius * gravitationalConstant );

    double m = 4.0 * ringRadius * r / std::pow( p, 2.0 );
    // m = k^2
    double k = std::sqrt( std::abs( m ) );

    Eigen::Vector3d acceleration;

    // Compute complete elliptic integrals
    double ellipticIntegralK = boost::math::ellint_1( k );
    double ellipticIntegralE = boost::math::ellint_2( k );
    double ellipticIntegralB = boost::math::ellint_rf( 0.0, 1.0 - m, 1.0 ) - boost::math::ellint_rd( 0.0, 1.0 - m, 1.0 ) / 3.0;
    double ellipticIntegralS;

    if ( ellipticIntegralSFromDAndB )
    {
        ellipticIntegralS = ( boost::math::ellint_d( k ) - ellipticIntegralB ) / m;
    }
    else
    {
        ellipticIntegralS = ( ( 2.0 - m ) * ellipticIntegralK - 2.0 * ellipticIntegralE ) / std::pow( m, 2.0 );
    }

    double Ar = 8.0 * lineDensity * ringRadius / std::pow( p, 3.0 ) * (
            ( std::pow( r, 2.0 ) + std::pow( z, 2.0 ) + std::pow( ringRadius, 2.0 ) * ellipticIntegralB / std::pow( q, 2.0 ) +
            2.0 * ringRadius * ( r + ringRadius ) * ellipticIntegralS / std::pow( p, 2.0 ) )
            );

    acceleration( 0 ) = - Ar * x;
    acceleration( 1 ) = - Ar * y;
    acceleration( 2 ) = - 4.0 * lineDensity * ringRadius * ellipticIntegralE / ( p * std::pow( q, 2.0) ) * z;

    return acceleration;
}

} // namespace gravitation

} // namespace tudat
