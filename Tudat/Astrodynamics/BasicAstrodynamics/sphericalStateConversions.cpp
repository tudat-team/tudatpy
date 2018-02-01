/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <vector>
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/coordinateConversions.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/sphericalStateConversions.h"
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

namespace orbital_element_conversions
{

//! Calculate current heading angle.
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return std::atan2( velocityInVerticalFrame( 1 ), velocityInVerticalFrame( 0 ) );
}

//! Calculatre current flight path angle.
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame )
{
    return -std::asin( velocityInVerticalFrame( 2 ) / velocityInVerticalFrame.norm( ) );
}

//! Function to convert a Cartesian to a spherical orbital state.
Eigen::Vector6d convertCartesianToSphericalOrbitalState(
        const Eigen::Vector6d& bodyFixedCartesianState )
{
    Eigen::Vector6d sphericalOrbitalState;

    // Compute and set spherical position
    Eigen::Vector3d sphericalPosition = coordinate_conversions::convertCartesianToSpherical< double >(
                bodyFixedCartesianState.segment( 0, 3 ) );
    sphericalOrbitalState( radiusIndex ) = sphericalPosition( 0 );
    sphericalOrbitalState( latitudeIndex ) = mathematical_constants::PI / 2.0 - sphericalPosition( 1 );
    sphericalOrbitalState( longitudeIndex ) = sphericalPosition( 2 );

    // Convert velocity to vertical frame.
    Eigen::Vector3d verticalFrameVelocity =
            reference_frames::getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
                sphericalOrbitalState( longitudeIndex ), sphericalOrbitalState( latitudeIndex ) ) *
            bodyFixedCartesianState.segment( 3, 3 );

    // Set velocity in spherical orbital state.
    sphericalOrbitalState( speedIndex ) = verticalFrameVelocity.norm( );
    sphericalOrbitalState( flightPathIndex ) = calculateFlightPathAngle( verticalFrameVelocity );
    sphericalOrbitalState( headingAngleIndex ) = calculateHeadingAngle( verticalFrameVelocity );

    return sphericalOrbitalState;

}

//! Function to convert a spherical orbital to a Cartesian state.
Eigen::Vector6d convertSphericalOrbitalToCartesianState(
        const Eigen::Vector6d& sphericalOrbitalState )
{
    Eigen::Vector6d cartesianState;

    // Compute Cartesian position
    Eigen::Vector3d sphericalPosition = sphericalOrbitalState.segment( 0, 3 );
    sphericalPosition( 1 ) = mathematical_constants::PI / 2.0 - sphericalOrbitalState( 1 );
    cartesianState.segment( 0, 3 ) = coordinate_conversions::convertSphericalToCartesian(
                sphericalPosition );

    // Check whether flight path angle and heading angle are valid (corresponding velocity components are zero otherwise)
    bool isFlightPathAngleValid = ( sphericalOrbitalState( flightPathIndex ) == sphericalOrbitalState( flightPathIndex ) );
    bool isHeadingAngleValid = ( sphericalOrbitalState( headingAngleIndex ) == sphericalOrbitalState( headingAngleIndex ) );

    // Compute velocity in vertical frame.
    Eigen::Vector3d velocityInVerticalFrame = Eigen::Vector3d::Zero( );
    if( isFlightPathAngleValid && isHeadingAngleValid )
    {
        velocityInVerticalFrame( 0 ) = sphericalOrbitalState( speedIndex ) *
                std::cos( sphericalOrbitalState( flightPathIndex ) ) *
                std::cos( sphericalOrbitalState( headingAngleIndex ) );
        velocityInVerticalFrame( 1 ) = sphericalOrbitalState( speedIndex ) *
                std::cos( sphericalOrbitalState( flightPathIndex ) ) *
                std::sin( sphericalOrbitalState( headingAngleIndex ) );
    }

    if( isFlightPathAngleValid )
    {
        velocityInVerticalFrame( 2 ) = -sphericalOrbitalState( speedIndex ) *
                std::sin( sphericalOrbitalState( flightPathIndex ) );
    }

    // Set velocity in body-fixed frame.
    cartesianState.segment( 3, 3 ) = reference_frames::getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
                sphericalOrbitalState( longitudeIndex ), sphericalOrbitalState( latitudeIndex ) ) * velocityInVerticalFrame;

    return cartesianState;
}

} // namespace tudat

} // namespace orbital_element_conversions

