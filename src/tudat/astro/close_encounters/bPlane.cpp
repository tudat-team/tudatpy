/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <cmath>
#include <stdexcept>
#include <string>

#include <Eigen/Geometry>

#include "tudat/astro/close_encounters/bPlane.h"

namespace tudat
{
namespace close_encounters
{

namespace
{

//! Minimum norm of the cross product of the (unit) reference direction and the incoming asymptote.
const double MINIMUM_REFERENCE_DIRECTION_ANGLE_SINE = 1.0E-8;

//! Quantities of the osculating hyperbola from which the asymptotes and the B-plane are computed.
struct OsculatingHyperbola {
    Eigen::Vector3d angularMomentum;
    Eigen::Vector3d eccentricityVector;
    double asymptoticVelocity;
};

OsculatingHyperbola computeOsculatingHyperbola( const Eigen::Vector6d& relativeState, const double centralBodyGravitationalParameter )
{
    if( !( centralBodyGravitationalParameter > 0.0 ) )
    {
        throw std::runtime_error( "Error when computing B-plane: gravitational parameter must be positive, but is " +
                                  std::to_string( centralBodyGravitationalParameter ) );
    }

    const Eigen::Vector3d position = relativeState.segment< 3 >( 0 );
    const Eigen::Vector3d velocity = relativeState.segment< 3 >( 3 );
    const double distance = position.norm( );

    OsculatingHyperbola hyperbola;
    hyperbola.angularMomentum = position.cross( velocity );
    hyperbola.eccentricityVector = velocity.cross( hyperbola.angularMomentum ) / centralBodyGravitationalParameter - position / distance;

    // Twice the specific orbital energy, equal to v_inf^2 for a hyperbola
    const double squaredAsymptoticVelocity = velocity.squaredNorm( ) - 2.0 * centralBodyGravitationalParameter / distance;
    if( !( squaredAsymptoticVelocity > 0.0 ) )
    {
        throw std::runtime_error( "Error when computing B-plane: relative orbit is not hyperbolic (eccentricity " +
                                  std::to_string( hyperbola.eccentricityVector.norm( ) ) + ", specific orbital energy " +
                                  std::to_string( squaredAsymptoticVelocity / 2.0 ) + " m^2/s^2)" );
    }
    hyperbola.asymptoticVelocity = std::sqrt( squaredAsymptoticVelocity );

    return hyperbola;
}

//! Compute the direction of the incoming or outgoing asymptote of an osculating hyperbola.
/*!
 *  The incoming asymptote is S = P / e + sqrt( 1 - 1 / e^2 ) Q, with P = e_vec / e and Q = h_hat x P. Since
 *  sqrt( e^2 - 1 ) = h v_inf / mu, this is S = ( mu e_vec + v_inf h x e_vec ) / ( mu e^2 ). The outgoing asymptote
 *  follows by reversing the sign of the eccentricity vector term.
 */
Eigen::Vector3d computeAsymptoteDirection( const OsculatingHyperbola& hyperbola,
                                           const double centralBodyGravitationalParameter,
                                           const bool isIncoming )
{
    const double eccentricityVectorSign = isIncoming ? 1.0 : -1.0;
    const Eigen::Vector3d asymptoteDirection = eccentricityVectorSign * centralBodyGravitationalParameter * hyperbola.eccentricityVector +
            hyperbola.asymptoticVelocity * hyperbola.angularMomentum.cross( hyperbola.eccentricityVector );

    // Normalised explicitly (instead of dividing by mu e^2) to remove rounding errors
    return asymptoteDirection.normalized( );
}

}  // namespace

Eigen::Matrix3d BPlaneData::getRotationMatrixToBPlaneFrame( ) const
{
    Eigen::Matrix3d rotationToBPlaneFrame;
    rotationToBPlaneFrame.row( 0 ) = tUnitVector.transpose( );
    rotationToBPlaneFrame.row( 1 ) = rUnitVector.transpose( );
    rotationToBPlaneFrame.row( 2 ) = sUnitVector.transpose( );
    return rotationToBPlaneFrame;
}

Eigen::Vector3d computeIncomingAsymptoteDirection( const Eigen::Vector6d& relativeState, const double centralBodyGravitationalParameter )
{
    return computeAsymptoteDirection(
            computeOsculatingHyperbola( relativeState, centralBodyGravitationalParameter ), centralBodyGravitationalParameter, true );
}

Eigen::Vector3d computeOutgoingAsymptoteDirection( const Eigen::Vector6d& relativeState, const double centralBodyGravitationalParameter )
{
    return computeAsymptoteDirection(
            computeOsculatingHyperbola( relativeState, centralBodyGravitationalParameter ), centralBodyGravitationalParameter, false );
}

BPlaneData makeBPlane( const Eigen::Vector6d& relativeState,
                       const double centralBodyGravitationalParameter,
                       const Eigen::Vector3d& referenceDirection )
{
    const OsculatingHyperbola hyperbola = computeOsculatingHyperbola( relativeState, centralBodyGravitationalParameter );

    BPlaneData bPlane;
    bPlane.asymptoticVelocity = hyperbola.asymptoticVelocity;
    bPlane.sUnitVector = computeAsymptoteDirection( hyperbola, centralBodyGravitationalParameter, true );

    // B-plane axes, from the reference direction
    const double referenceDirectionNorm = referenceDirection.norm( );
    if( !( referenceDirectionNorm > 0.0 ) )
    {
        throw std::runtime_error( "Error when computing B-plane: reference direction has zero norm" );
    }
    const Eigen::Vector3d referenceCrossS = ( referenceDirection / referenceDirectionNorm ).cross( bPlane.sUnitVector );
    if( referenceCrossS.norm( ) < MINIMUM_REFERENCE_DIRECTION_ANGLE_SINE )
    {
        throw std::runtime_error(
                "Error when computing B-plane: reference direction is (almost) parallel to the incoming asymptote, "
                "B-plane axes are undefined" );
    }
    bPlane.tUnitVector = referenceCrossS.normalized( );
    bPlane.rUnitVector = bPlane.sUnitVector.cross( bPlane.tUnitVector );

    // The angular momentum along the incoming asymptote is h = v_inf B x S, so that B = S x h / v_inf
    bPlane.bVector = bPlane.sUnitVector.cross( hyperbola.angularMomentum ) / hyperbola.asymptoticVelocity;
    bPlane.impactParameter = bPlane.bVector.norm( );
    bPlane.bT = bPlane.bVector.dot( bPlane.tUnitVector );
    bPlane.bR = bPlane.bVector.dot( bPlane.rUnitVector );

    return bPlane;
}

double computeCaptureRadius( const double bodyRadius, const double centralBodyGravitationalParameter, const double asymptoticVelocity )
{
    if( !( bodyRadius > 0.0 ) )
    {
        throw std::runtime_error( "Error when computing capture radius: body radius must be positive, but is " +
                                  std::to_string( bodyRadius ) );
    }
    if( !( centralBodyGravitationalParameter > 0.0 ) )
    {
        throw std::runtime_error( "Error when computing capture radius: gravitational parameter must be positive, but is " +
                                  std::to_string( centralBodyGravitationalParameter ) );
    }
    if( !( asymptoticVelocity > 0.0 ) )
    {
        throw std::runtime_error( "Error when computing capture radius: asymptotic velocity must be positive, but is " +
                                  std::to_string( asymptoticVelocity ) );
    }
    return bodyRadius *
            std::sqrt( 1.0 + 2.0 * centralBodyGravitationalParameter / ( bodyRadius * asymptoticVelocity * asymptoticVelocity ) );
}

}  // namespace close_encounters

}  // namespace tudat
