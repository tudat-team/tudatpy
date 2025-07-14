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

void RingGravityCache::update( const Eigen::Vector3d& currentBodyFixedPosition )
{
    if( currentBodyFixedPosition != currentBodyFixedPosition_ )
    {
        currentBodyFixedPosition_ = currentBodyFixedPosition;

        double x = currentBodyFixedPosition_( 0 );
        double y = currentBodyFixedPosition_( 1 );
        double z = currentBodyFixedPosition_( 2 );

        // Fukushima (2010), eq. 2
        double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

        // Fukushima (2010), eq. 3
        double p = std::sqrt( std::pow( r + ringRadius_, 2.0 ) + std::pow( z, 2.0 ) );

        // Fukushima (2010), eq. 5
        double m = 4.0 * ringRadius_ * r / std::pow( p, 2.0 );
        // m = k^2
        double k = std::sqrt( std::abs( m ) );

        // Compute elliptic integrals

        try
        {
            currentEllipticIntegralK_ = boost::math::ellint_1( k );
        }
        // If boost throws overflow error: occurs for r=ringRadius singularity
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error(
                    "Error when computing ring gravity, likely as a result of trying to compute it for a point in or "
                    "very close to the ring (singularity). Caught exception: \n" +
                    std::string( caughtException.what( ) ) );
        }

        currentEllipticIntegralE_ = boost::math::ellint_2( k );
        // BOOST: https://www.boost.org/doc/libs/1_81_0/libs/math/doc/html/math_toolkit/ellint/ellint_intro.html
        currentEllipticIntegralB_ = boost::math::ellint_rf( 0.0, 1.0 - m, 1.0 ) - boost::math::ellint_rd( 0.0, 1.0 - m, 1.0 ) / 3.0;

        // If far from 1/m singularity: compute S using elliptic integrals
        if( m > 0.01 )
        {
            if( ellipticIntegralSFromDAndB_ )
            {
                // Fukushima (2010), eq. 32
                currentEllipticIntegralS_ = ( boost::math::ellint_d( k ) - currentEllipticIntegralB_ ) / m;
            }
            else
            {
                // Fukushima (2010), eq. 31
                currentEllipticIntegralS_ =
                        ( ( 2.0 - m ) * currentEllipticIntegralK_ - 2.0 * currentEllipticIntegralE_ ) / std::pow( m, 2.0 );
            }
        }
        // If close to 1/m singularity, compute S via Taylor series, according to Fukushima (2010), sec. A.1
        else
        {
            // Taylor series coefficients for m in [0, 0.1]
            std::vector< double > taylorCoefficients = { 0.204012532440038310, 0.159513582234205843, 0.130422818255893004,
                                                         0.111687838140976463, 0.098925188226691425, 0.089815348807960028,
                                                         0.083084759300136632, 0.077987984857306626, 0.074062924745595950,
                                                         0.071009059783923539, 0.068623059119746445, 0.066762755430661757,
                                                         0.065325983044110253 };
            double m0 = 0.05;

            currentEllipticIntegralS_ = 0;

            for( unsigned int taylorOrder = 0; taylorOrder <= 12; ++taylorOrder )
            {
                currentEllipticIntegralS_ += taylorCoefficients.at( taylorOrder ) * std::pow( m - m0, taylorOrder );
            }
        }
    }
}

//! Computes the gravitational potential of a one-dimensional ring.
double computeRingGravitationalPotential( const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
                                          const double ringRadius,
                                          const double gravitationalParameter,
                                          const double ellipticIntegralK )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    // Fukushima (2010), eq. 2
    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    // Fukushima (2010), eq. 3
    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensityTimesGravitationalConst = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius );

    // Fukushima (2010), eq. 1
    return 4.0 * lineDensityTimesGravitationalConst * ringRadius * ellipticIntegralK / p;
}

//! Computes the gravitational acceleration of a one-dimensional ring
Eigen::Vector3d computeRingGravitationalAcceleration( const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
                                                      const double ringRadius,
                                                      const double gravitationalParameter,
                                                      const double ellipticIntegralB,
                                                      const double ellipticIntegralE,
                                                      const double ellipticIntegralS )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    // Fukushima (2010), eq. 2
    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    // Fukushima (2010), eq. 3
    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );
    // Fukushima (2010), eq. 15
    double q = std::sqrt( std::pow( r - ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensityTimesGravitationalConst = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius );

    Eigen::Vector3d acceleration;

    // Fukushima (2010), eq. 30
    double Ar = 8.0 * lineDensityTimesGravitationalConst * ringRadius / std::pow( p, 3.0 ) *
            ( ( std::pow( r, 2.0 ) + std::pow( z, 2.0 ) - std::pow( ringRadius, 2.0 ) ) * ellipticIntegralB / std::pow( q, 2.0 ) +
              2.0 * ringRadius * ( r + ringRadius ) * ellipticIntegralS / std::pow( p, 2.0 ) );

    // Fukushima (2010), eq. 16
    acceleration( 0 ) = -Ar * x;
    acceleration( 1 ) = -Ar * y;
    // Fukushima (2010), eq. 16, 18
    acceleration( 2 ) = -4.0 * lineDensityTimesGravitationalConst * ringRadius * ellipticIntegralE / ( p * std::pow( q, 2.0 ) ) * z;

    return acceleration;
}

Eigen::Matrix3d computeRingHessianOfGravitationalPotential( const Eigen::Vector3d& positionOfBodySubjectToAcceleration,
                                                            const double ringRadius,
                                                            const double gravitationalParameter,
                                                            const double ellipticIntegralB,
                                                            const double ellipticIntegralE,
                                                            const double ellipticIntegralS,
                                                            const double ellipticIntegralK )
{
    double x = positionOfBodySubjectToAcceleration( 0 );
    double y = positionOfBodySubjectToAcceleration( 1 );
    double z = positionOfBodySubjectToAcceleration( 2 );

    // Fukushima (2010), eq. 2
    double r = std::sqrt( std::pow( x, 2.0 ) + std::pow( y, 2.0 ) );

    // Fukushima (2010), eq. 3
    double p = std::sqrt( std::pow( r + ringRadius, 2.0 ) + std::pow( z, 2.0 ) );
    // Fukushima (2010), eq. 15
    double q = std::sqrt( std::pow( r - ringRadius, 2.0 ) + std::pow( z, 2.0 ) );

    double lineDensityTimesGravitationalConst = gravitationalParameter / ( 2.0 * mathematical_constants::PI * ringRadius );

    // Derivatives of acceleration, derived via sympy

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Version 1 of equations: derived from Eq. 30 of Fukushima (2010) -> currently commented
    // Not sure why it leads to larger errors wrt numerical partials than version 2... could it be some error in the
    //      derivation? or maybe numerical cancellation occurring somewhere?

    //    double dax_dx =-8*lineDensityTimesGravitationalConst*ringRadius*x*(2*ellipticIntegralB*x/std::pow(q, 2) +
    //            0.5*ellipticIntegralS*(4*ringRadius*x/(std::pow(p, 2)*r) - 8*ringRadius*x*(r + ringRadius)/std::pow(p, 4))*
    //            (std::pow(r, 2) - std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) +
    //            (-2*ellipticIntegralB*x*(r - ringRadius)*(std::pow(r, 2) - std::pow(ringRadius, 2) + std::pow(z, 2))/
    //            std::pow(q, 4) + 2*ellipticIntegralS*ringRadius*x/std::pow(p, 2) -
    //            4*ellipticIntegralS*ringRadius*x*std::pow(r + ringRadius, 2)/std::pow(p, 4))/r -
    //            2*ringRadius*(r + ringRadius)*(0.5*ellipticIntegralB*std::pow(p, 2) + 1.0*ellipticIntegralS*
    //            (-1.0*std::pow(p, 2) + 4.0*r*ringRadius))*(4*ringRadius*x/(std::pow(p, 2)*r) -
    //            8*ringRadius*x*(r + ringRadius)/std::pow(p, 4))/(std::pow(p, 2)*(-std::pow(p, 2) +
    //            4*r*ringRadius)))/std::pow(p, 3) - 8*lineDensityTimesGravitationalConst*ringRadius*(ellipticIntegralB*(std::pow(r, 2) -
    //            std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) + 2*ellipticIntegralS*ringRadius*(r +
    //            ringRadius)/std::pow(p, 2))/std::pow(p, 3) + 24*lineDensityTimesGravitationalConst*ringRadius*std::pow(x, 2)*
    //            (r + ringRadius)*(ellipticIntegralB*(std::pow(r, 2) - std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) +
    //            2*ellipticIntegralS*ringRadius*(r + ringRadius)/std::pow(p, 2))/(std::pow(p, 5)*r);
    //
    //    double dax_dy =-8*lineDensityTimesGravitationalConst*ringRadius*x*(2*ellipticIntegralB*y/std::pow(q, 2) +
    //            0.5*ellipticIntegralS*(4*ringRadius*y/(std::pow(p, 2)*r) - 8*ringRadius*y*(r + ringRadius)/std::pow(p, 4))*
    //            (std::pow(r, 2) - std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) +
    //            (-2*ellipticIntegralB*y*(r - ringRadius)*(std::pow(r, 2) - std::pow(ringRadius, 2) +
    //            std::pow(z, 2))/std::pow(q, 4) + 2*ellipticIntegralS*ringRadius*y/std::pow(p, 2) -
    //            4*ellipticIntegralS*ringRadius*y*std::pow(r + ringRadius, 2)/std::pow(p, 4))/r -
    //            2*ringRadius*(r + ringRadius)*(0.5*ellipticIntegralB*std::pow(p, 2) + 1.0*ellipticIntegralS*(-1.0*std::pow(p, 2) +
    //            4.0*r*ringRadius))*(4*ringRadius*y/(std::pow(p, 2)*r) - 8*ringRadius*y*(r + ringRadius)/std::pow(p, 4))/
    //            (std::pow(p, 2)*(-std::pow(p, 2) + 4*r*ringRadius)))/std::pow(p, 3) +
    //            24*lineDensityTimesGravitationalConst*ringRadius*x*y*(r + ringRadius)*(ellipticIntegralB*(std::pow(r, 2) -
    //            std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) + 2*ellipticIntegralS*ringRadius*(r +
    //            ringRadius)/std::pow(p, 2))/(std::pow(p, 5)*r);
    //
    //    double day_dy =-8*lineDensityTimesGravitationalConst*ringRadius*y*(2*ellipticIntegralB*y/std::pow(q, 2) +
    //            0.5*ellipticIntegralS*(4*ringRadius*y/(std::pow(p, 2)*r) - 8*ringRadius*y*(r + ringRadius)/std::pow(p, 4))*(std::pow(r, 2)
    //            - std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) + (-2*ellipticIntegralB*y*(r - ringRadius)*(std::pow(r, 2) -
    //            std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 4) + 2*ellipticIntegralS*ringRadius*y/std::pow(p, 2) -
    //            4*ellipticIntegralS*ringRadius*y*std::pow(r + ringRadius, 2)/std::pow(p, 4))/r - 2*ringRadius*
    //            (r + ringRadius)*(0.5*ellipticIntegralB*std::pow(p, 2) + 1.0*ellipticIntegralS*(-1.0*std::pow(p, 2) +
    //            4.0*r*ringRadius))*(4*ringRadius*y/(std::pow(p, 2)*r) - 8*ringRadius*y*(r + ringRadius)/std::pow(p, 4))/
    //            (std::pow(p, 2)*(-std::pow(p, 2) + 4*r*ringRadius)))/std::pow(p, 3) -
    //            8*lineDensityTimesGravitationalConst*ringRadius*(ellipticIntegralB*(std::pow(r, 2) - std::pow(ringRadius, 2) +
    //            std::pow(z, 2))/std::pow(q, 2) + 2*ellipticIntegralS*ringRadius*(r + ringRadius)/std::pow(p, 2))/
    //            std::pow(p, 3) + 24*lineDensityTimesGravitationalConst*ringRadius*std::pow(y, 2)*(r + ringRadius)*(
    //            ellipticIntegralB*(std::pow(r, 2) - std::pow(ringRadius, 2) + std::pow(z, 2))/std::pow(q, 2) +
    //            2*ellipticIntegralS*ringRadius*(r + ringRadius)/std::pow(p, 2))/(std::pow(p, 5)*r);
    //
    //    double daz_dx =
    //            (-0.5*lineDensityTimesGravitationalConst*p*z*(-ellipticIntegralK + ellipticIntegralE)*(4*ringRadius*x/(std::pow(p, 2)*r)
    //            - 8*ringRadius*x*(r + ringRadius)/std::pow(p, 4))/std::pow(q, 2) + 8*lineDensityTimesGravitationalConst*
    //            ringRadius*x*z*(r - ringRadius)*ellipticIntegralE/(p*std::pow(q, 4)) + 4*lineDensityTimesGravitationalConst*
    //            ringRadius*x*z*(r + ringRadius)*ellipticIntegralE/(std::pow(p, 3)*std::pow(q, 2)))/r;
    //
    //    double daz_dy =
    //            (-0.5*lineDensityTimesGravitationalConst*p*z*(-ellipticIntegralK + ellipticIntegralE)*(4*ringRadius*y/(std::pow(p, 2)*r)
    //            - 8*ringRadius*y*(r + ringRadius)/std::pow(p, 4))/std::pow(q, 2) + 8*lineDensityTimesGravitationalConst*
    //            ringRadius*y*z*(r - ringRadius)*ellipticIntegralE/(p*std::pow(q, 4)) +
    //            4*lineDensityTimesGravitationalConst*ringRadius*y*z*(r + ringRadius)*ellipticIntegralE/(std::pow(p, 3)*
    //            std::pow(q, 2)))/r;
    //
    //    double daz_dz =
    //            -4*lineDensityTimesGravitationalConst*ringRadius*ellipticIntegralE/(p*std::pow(q, 2)) +
    //            8*lineDensityTimesGravitationalConst*ringRadius*std::pow(z, 2)*ellipticIntegralE/(p*std::pow(q, 4)) +
    //            4.0*lineDensityTimesGravitationalConst*ringRadius*std::pow(z, 2)*(-ellipticIntegralK + ellipticIntegralE)/
    //            (std::pow(p, 3)*std::pow(q, 2)) + 4*lineDensityTimesGravitationalConst*ringRadius*std::pow(z, 2)*
    //            ellipticIntegralE/(std::pow(p, 3)*std::pow(q, 2));

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    // Version 2 of equations: derived from Eq. 18 of Fukushima (2010)
    // Analytical partials are closer to numerical ones than for version 1

    // Fukushima (2010), eq. 5
    double m = 4.0 * ringRadius * r / std::pow( p, 2.0 );

    double dax_dx = -2 * lineDensityTimesGravitationalConst * ringRadius * x *
                    ( 0.125 * std::pow( p, 2 ) * ( ellipticIntegralE - ellipticIntegralK * ( 1 - m ) ) *
                              ( 4 * ringRadius * x / ( std::pow( p, 2 ) * r ) -
                                8 * ringRadius * x * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                              ( std::pow( r, 3 ) * ringRadius * ( 1 - m ) ) +
                      ( -2 * ellipticIntegralE * x * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) *
                                ( r - ringRadius ) / std::pow( q, 4 ) +
                        0.125 * std::pow( p, 2 ) * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) *
                                ( ellipticIntegralE - ellipticIntegralK ) *
                                ( 4 * ringRadius * x / ( std::pow( p, 2 ) * r ) -
                                  8 * ringRadius * x * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                                ( std::pow( q, 2 ) * ringRadius ) ) /
                              r +
                      ( 2 * ellipticIntegralE * x * ( std::pow( ringRadius, 2 ) + std::pow( z, 2 ) ) / std::pow( q, 2 ) -
                        2 * ellipticIntegralK * x ) /
                              std::pow( r, 4 ) ) /
                    p -
            2 * lineDensityTimesGravitationalConst * ringRadius *
                    ( ellipticIntegralE * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) / std::pow( q, 2 ) +
                      ellipticIntegralK / std::pow( r, 2 ) ) /
                    p +
            2 * lineDensityTimesGravitationalConst * ringRadius * std::pow( x, 2 ) * ( r + ringRadius ) *
                    ( ellipticIntegralE * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) / std::pow( q, 2 ) +
                      ellipticIntegralK / std::pow( r, 2 ) ) /
                    ( std::pow( p, 3 ) * r );

    double dax_dy = -2 * lineDensityTimesGravitationalConst * ringRadius * x *
                    ( 0.125 * std::pow( p, 2 ) * ( ellipticIntegralE - ellipticIntegralK * ( 1 - m ) ) *
                              ( 4 * ringRadius * y / ( std::pow( p, 2 ) * r ) -
                                8 * ringRadius * y * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                              ( std::pow( r, 3 ) * ringRadius * ( 1 - m ) ) +
                      ( -2 * ellipticIntegralE * y * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) *
                                ( r - ringRadius ) / std::pow( q, 4 ) +
                        0.125 * std::pow( p, 2 ) * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) *
                                ( ellipticIntegralE - ellipticIntegralK ) *
                                ( 4 * ringRadius * y / ( std::pow( p, 2 ) * r ) -
                                  8 * ringRadius * y * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                                ( std::pow( q, 2 ) * ringRadius ) ) /
                              r +
                      ( 2 * ellipticIntegralE * y * ( std::pow( ringRadius, 2 ) + std::pow( z, 2 ) ) / std::pow( q, 2 ) -
                        2 * ellipticIntegralK * y ) /
                              std::pow( r, 4 ) ) /
                    p +
            2 * lineDensityTimesGravitationalConst * ringRadius * x * y * ( r + ringRadius ) *
                    ( ellipticIntegralE * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) / std::pow( q, 2 ) +
                      ellipticIntegralK / std::pow( r, 2 ) ) /
                    ( std::pow( p, 3 ) * r );

    double day_dy = -2 * lineDensityTimesGravitationalConst * ringRadius * y *
                    ( 0.125 * std::pow( p, 2 ) * ( ellipticIntegralE - ellipticIntegralK * ( 1 - m ) ) *
                              ( 4 * ringRadius * y / ( std::pow( p, 2 ) * r ) -
                                8 * ringRadius * y * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                              ( std::pow( r, 3 ) * ringRadius * ( 1 - m ) ) +
                      ( -2 * ellipticIntegralE * y * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) *
                                ( r - ringRadius ) / std::pow( q, 4 ) +
                        0.125 * std::pow( p, 2 ) * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) *
                                ( ellipticIntegralE - ellipticIntegralK ) *
                                ( 4 * ringRadius * y / ( std::pow( p, 2 ) * r ) -
                                  8 * ringRadius * y * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                                ( std::pow( q, 2 ) * ringRadius ) ) /
                              r +
                      ( 2 * ellipticIntegralE * y * ( std::pow( ringRadius, 2 ) + std::pow( z, 2 ) ) / std::pow( q, 2 ) -
                        2 * ellipticIntegralK * y ) /
                              std::pow( r, 4 ) ) /
                    p -
            2 * lineDensityTimesGravitationalConst * ringRadius *
                    ( ellipticIntegralE * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) / std::pow( q, 2 ) +
                      ellipticIntegralK / std::pow( r, 2 ) ) /
                    p +
            2 * lineDensityTimesGravitationalConst * ringRadius * std::pow( y, 2 ) * ( r + ringRadius ) *
                    ( ellipticIntegralE * ( 1 + ( -std::pow( ringRadius, 2 ) - std::pow( z, 2 ) ) / std::pow( r, 2 ) ) / std::pow( q, 2 ) +
                      ellipticIntegralK / std::pow( r, 2 ) ) /
                    ( std::pow( p, 3 ) * r );

    double daz_dx =
            ( 8 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius * x * z * ( r - ringRadius ) /
                      ( p * std::pow( q, 4 ) ) +
              4 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius * x * z * ( r + ringRadius ) /
                      ( std::pow( p, 3 ) * std::pow( q, 2 ) ) -
              0.5 * lineDensityTimesGravitationalConst * p * z * ( ellipticIntegralE - ellipticIntegralK ) *
                      ( 4 * ringRadius * x / ( std::pow( p, 2 ) * r ) - 8 * ringRadius * x * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                      std::pow( q, 2 ) ) /
            r;

    double daz_dy =
            ( 8 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius * y * z * ( r - ringRadius ) /
                      ( p * std::pow( q, 4 ) ) +
              4 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius * y * z * ( r + ringRadius ) /
                      ( std::pow( p, 3 ) * std::pow( q, 2 ) ) -
              0.5 * lineDensityTimesGravitationalConst * p * z * ( ellipticIntegralE - ellipticIntegralK ) *
                      ( 4 * ringRadius * y / ( std::pow( p, 2 ) * r ) - 8 * ringRadius * y * ( r + ringRadius ) / std::pow( p, 4 ) ) /
                      std::pow( q, 2 ) ) /
            r;

    double daz_dz = -4 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius / ( p * std::pow( q, 2 ) ) +
            8 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius * std::pow( z, 2 ) / ( p * std::pow( q, 4 ) ) +
            4 * ellipticIntegralE * lineDensityTimesGravitationalConst * ringRadius * std::pow( z, 2 ) /
                    ( std::pow( p, 3 ) * std::pow( q, 2 ) ) +
            4.0 * lineDensityTimesGravitationalConst * ringRadius * std::pow( z, 2 ) * ( ellipticIntegralE - ellipticIntegralK ) /
                    ( std::pow( p, 3 ) * std::pow( q, 2 ) );

    Eigen::Matrix3d hessianMatrix;
    hessianMatrix( 0, 0 ) = dax_dx;
    hessianMatrix( 1, 1 ) = day_dy;
    hessianMatrix( 2, 2 ) = daz_dz;
    hessianMatrix( 0, 1 ) = dax_dy;
    hessianMatrix( 1, 0 ) = hessianMatrix( 0, 1 );
    hessianMatrix( 0, 2 ) = daz_dx;
    hessianMatrix( 2, 0 ) = hessianMatrix( 0, 2 );
    hessianMatrix( 1, 2 ) = daz_dy;
    hessianMatrix( 2, 1 ) = hessianMatrix( 1, 2 );

    return hessianMatrix;
}

}  // namespace gravitation

}  // namespace tudat
