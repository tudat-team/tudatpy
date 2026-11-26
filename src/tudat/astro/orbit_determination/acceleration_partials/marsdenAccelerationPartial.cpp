/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/astro/orbit_determination/acceleration_partials/marsdenAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

namespace
{

//! Skew-symmetric (cross-product) matrix: skew( a ) * b = a x b.
Eigen::Matrix3d skew( const Eigen::Vector3d& a )
{
    Eigen::Matrix3d s;
    s << 0.0, -a( 2 ), a( 1 ), a( 2 ), 0.0, -a( 0 ), -a( 1 ), a( 0 ), 0.0;
    return s;
}

//! Gradient of g(r) w.r.t. the delayed position d = r - v * tau  [1/m].
Eigen::Vector3d computeGFunctionGradientWrtDelayedPosition( const Eigen::Vector6d& relativeState,
                                                            const double alpha,
                                                            const double r0,
                                                            const double m,
                                                            const double n,
                                                            const double k,
                                                            const double tau )
{
    const Eigen::Vector3d delayedPosition = relativeState.segment( 0, 3 ) - relativeState.segment( 3, 3 ) * tau;
    const double delayedDistance = delayedPosition.norm( );

    // Scaled distance x = |d| / ( AU * r0 ) and g( x ) = alpha * x^-m * ( 1 + x^n )^-k
    const double x = delayedDistance / physical_constants::ASTRONOMICAL_UNIT / r0;
    const double g = alpha * std::pow( x, -m ) * std::pow( 1.0 + std::pow( x, n ), -k );

    // dg/dx = g * ( -m / x - k * n * x^( n - 1 ) / ( 1 + x^n ) )
    const double dgdx = g * ( -m / x - k * n * std::pow( x, n - 1.0 ) / ( 1.0 + std::pow( x, n ) ) );

    // dx/dd = d / ( |d| * AU * r0 )
    return dgdx * delayedPosition / ( delayedDistance * physical_constants::ASTRONOMICAL_UNIT * r0 );
}

//! Unit vectors of the RSW frame and their derivatives w.r.t. position and velocity (dR/dv = 0).
struct RswDerivatives {
    Eigen::Vector3d R, S, W;
    Eigen::Matrix3d dRdr, dSdr, dWdr, dSdv, dWdv;
};

RswDerivatives computeRswDerivatives( const Eigen::Vector6d& relativeState )
{
    RswDerivatives d;
    const Eigen::Vector3d r = relativeState.segment( 0, 3 );
    const Eigen::Vector3d v = relativeState.segment( 3, 3 );
    const Eigen::Matrix3d identity = Eigen::Matrix3d::Identity( );

    const double rNorm = r.norm( );
    const Eigen::Vector3d h = r.cross( v );
    const double hNorm = h.norm( );

    d.R = r / rNorm;
    d.W = h / hNorm;
    d.S = d.W.cross( d.R );

    // R = r / |r|
    d.dRdr = ( identity - d.R * d.R.transpose( ) ) / rNorm;

    // W = h / |h|, h = r x v  ->  dh/dr = -skew( v ), dh/dv = skew( r )
    const Eigen::Matrix3d dWdh = ( identity - d.W * d.W.transpose( ) ) / hNorm;
    d.dWdr = dWdh * ( -skew( v ) );
    d.dWdv = dWdh * skew( r );

    // S = W x R  ->  dS = skew( W ) dR - skew( R ) dW
    d.dSdr = skew( d.W ) * d.dRdr - skew( d.R ) * d.dWdr;
    d.dSdv = -skew( d.R ) * d.dWdv;

    return d;
}

}  // namespace

Eigen::Matrix3d calculatePartialOfMarsdenAccelerationWrtPositionOfAcceleratedBody( const Eigen::Vector6d& relativeState,
                                                                                   const double A1,
                                                                                   const double A2,
                                                                                   const double A3,
                                                                                   const double alpha,
                                                                                   const double r0,
                                                                                   const double m,
                                                                                   const double n,
                                                                                   const double k,
                                                                                   const double tau )
{
    const RswDerivatives rsw = computeRswDerivatives( relativeState );

    const double delayedDistanceAu =
            ( relativeState.segment( 0, 3 ) - relativeState.segment( 3, 3 ) * tau ).norm( ) / physical_constants::ASTRONOMICAL_UNIT;
    const double rRatio = delayedDistanceAu / r0;
    const double g = alpha * std::pow( rRatio, -m ) * std::pow( 1.0 + std::pow( rRatio, n ), -k );

    const Eigen::Vector3d dgdPosition = computeGFunctionGradientWrtDelayedPosition( relativeState, alpha, r0, m, n, k, tau );
    const Eigen::Vector3d directionInInertialFrame = A1 * rsw.R + A2 * rsw.S + A3 * rsw.W;

    // a = ( A1 R + A2 S + A3 W ) g  ->  da/dr = g ( A1 dR + A2 dS + A3 dW ) + ( A1 R + A2 S + A3 W ) dg/dr^T
    return g * ( A1 * rsw.dRdr + A2 * rsw.dSdr + A3 * rsw.dWdr ) + directionInInertialFrame * dgdPosition.transpose( );
}

Eigen::Matrix3d calculatePartialOfMarsdenAccelerationWrtVelocityOfAcceleratedBody( const Eigen::Vector6d& relativeState,
                                                                                   const double A1,
                                                                                   const double A2,
                                                                                   const double A3,
                                                                                   const double alpha,
                                                                                   const double r0,
                                                                                   const double m,
                                                                                   const double n,
                                                                                   const double k,
                                                                                   const double tau )
{
    const RswDerivatives rsw = computeRswDerivatives( relativeState );

    const double delayedDistanceAu =
            ( relativeState.segment( 0, 3 ) - relativeState.segment( 3, 3 ) * tau ).norm( ) / physical_constants::ASTRONOMICAL_UNIT;
    const double rRatio = delayedDistanceAu / r0;
    const double g = alpha * std::pow( rRatio, -m ) * std::pow( 1.0 + std::pow( rRatio, n ), -k );

    const Eigen::Vector3d dgdDelayedPosition = computeGFunctionGradientWrtDelayedPosition( relativeState, alpha, r0, m, n, k, tau );
    const Eigen::Vector3d directionInInertialFrame = A1 * rsw.R + A2 * rsw.S + A3 * rsw.W;

    // R does not depend on v; d = r - v tau  ->  dg/dv = -tau * dg/dd
    return g * ( A2 * rsw.dSdv + A3 * rsw.dWdv ) + directionInInertialFrame * ( -tau * dgdDelayedPosition ).transpose( );
}

//! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
std::pair< std::function< void( Eigen::MatrixXd& ) >, int > MarsdenAccelerationPartial::getParameterPartialFunctionDerivedAcceleration(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
{
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

    // Check dependencies. NOTE: the parameter names below are placeholders and must match the enum you add.
    if( parameter->getParameterName( ).second.first == this->getAcceleratedBody( ) )
    {
        switch( parameter->getParameterName( ).first )
        {
            case estimatable_parameters::marsden_a1_coefficient:
                partialFunctionPair = std::make_pair( std::bind( &MarsdenAccelerationPartial::wrtA1, this, std::placeholders::_1 ), 1 );
                break;
            case estimatable_parameters::marsden_a2_coefficient:
                partialFunctionPair = std::make_pair( std::bind( &MarsdenAccelerationPartial::wrtA2, this, std::placeholders::_1 ), 1 );
                break;
            case estimatable_parameters::marsden_a3_coefficient:
                partialFunctionPair = std::make_pair( std::bind( &MarsdenAccelerationPartial::wrtA3, this, std::placeholders::_1 ), 1 );
                break;
            case estimatable_parameters::marsden_tau:
                partialFunctionPair = std::make_pair( std::bind( &MarsdenAccelerationPartial::wrtTau, this, std::placeholders::_1 ), 1 );
                break;
            default:
                partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
                break;
        }
    }
    else
    {
        partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
    }

    return partialFunctionPair;
}

//! Function to calculate the partial w.r.t. the time offset tau.
void MarsdenAccelerationPartial::wrtTau( Eigen::MatrixXd& tauPartial )
{
    const Eigen::Vector6d relativeState = marsdenAcceleration_->getCurrentRelativeState( );

    // d = r - v tau  ->  dd/dtau = -v  ->  dg/dtau = -v . dg/dd
    const Eigen::Vector3d dgdDelayedPosition = computeGFunctionGradientWrtDelayedPosition( relativeState,
                                                                                           marsdenAcceleration_->getAlpha( ),
                                                                                           marsdenAcceleration_->getR0( ),
                                                                                           marsdenAcceleration_->getM( ),
                                                                                           marsdenAcceleration_->getN( ),
                                                                                           marsdenAcceleration_->getK( ),
                                                                                           marsdenAcceleration_->getTau( ) );
    const double dgdTau = -relativeState.segment( 3, 3 ).dot( dgdDelayedPosition );

    // The RSW rotation is taken at time t, so only g(r) depends on tau.
    const Eigen::Vector3d directionInInertialFrame = marsdenAcceleration_->getCurrentRotationRswToInertial( ) *
            Eigen::Vector3d( marsdenAcceleration_->getA1( ), marsdenAcceleration_->getA2( ), marsdenAcceleration_->getA3( ) );

    tauPartial = directionInInertialFrame * dgdTau;
}

}  // namespace acceleration_partials

}  // namespace tudat