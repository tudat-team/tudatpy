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

void RingGravityCache::update (const Eigen::Vector3d& currentBodyFixedPosition)
{
    if ( currentBodyFixedPosition != currentBodyFixedPosition_ )
    {
        currentBodyFixedPosition_ = currentBodyFixedPosition;

        double x = currentBodyFixedPosition_( 0 );
        double y = currentBodyFixedPosition_( 1 );
        double z = currentBodyFixedPosition_( 2 );

        double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

        double p = std::sqrt( std::pow( r + ringRadius_, 2.0 ) + std::pow( z, 2.0 ) );

        double m = 4.0 * ringRadius_ * r / std::pow( p, 2.0 );
        // m = k^2
        double k = std::sqrt( std::abs( m ) );

        // Compute elliptic integrals
        currentEllipticIntegralK_ = boost::math::ellint_1( k );
        currentEllipticIntegralE_ = boost::math::ellint_2( k );
        currentEllipticIntegralB_ = boost::math::ellint_rf( 0.0, 1.0 - m, 1.0 ) - boost::math::ellint_rd( 0.0, 1.0 - m, 1.0 ) / 3.0;

        if ( ellipticIntegralSFromDAndB_ )
        {
            currentEllipticIntegralS_ = ( boost::math::ellint_d( k ) - currentEllipticIntegralB_ ) / m;
        }
        else
        {
            currentEllipticIntegralS_ = ( ( 2.0 - m ) * currentEllipticIntegralK_ - 2.0 * currentEllipticIntegralE_ ) / std::pow( m, 2.0 );
        }
    }
}

//! Computes the gravitational potential of a one-dimensional ring.
double computeRingGravitationalPotential(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double ringRadius,
        const double gravitationalParameter,
        const double ellipticIntegralK )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensityTimesGravitationalConst = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius );

    return 4.0 * lineDensityTimesGravitationalConst * ringRadius * ellipticIntegralK / p;
}

//! Computes the gravitational acceleration of a one-dimensional ring
Eigen::Vector3d computeRingGravitationalAcceleration(
        const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
        const double ringRadius,
        const double gravitationalParameter,
        const double ellipticIntegralB,
        const double ellipticIntegralE,
        const double ellipticIntegralS )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );
    double q = std::sqrt( std::pow( r - ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensityTimesGravitationalConst = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius );

    Eigen::Vector3d acceleration;

    double Ar = 8.0 * lineDensityTimesGravitationalConst * ringRadius / std::pow( p, 3.0 ) * (
            ( std::pow( r, 2.0 ) + std::pow( z, 2.0 ) + std::pow( ringRadius, 2.0 ) * ellipticIntegralB / std::pow( q, 2.0 ) +
            2.0 * ringRadius * ( r + ringRadius ) * ellipticIntegralS / std::pow( p, 2.0 ) )
            );

    acceleration( 0 ) = - Ar * x;
    acceleration( 1 ) = - Ar * y;
    acceleration( 2 ) = - 4.0 * lineDensityTimesGravitationalConst * ringRadius * ellipticIntegralE / ( p * std::pow( q, 2.0) ) * z;

    return acceleration;
}

} // namespace gravitation

} // namespace tudat
