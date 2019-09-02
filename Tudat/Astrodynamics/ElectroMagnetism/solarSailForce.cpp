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

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/solarSailForce.h"
#include <iostream>
#include <math.h>
#include <vector>

namespace tudat
{
namespace electro_magnetism
{

//! Compute solar sail force using a non-ideal reflective model.
Eigen::Vector3d computeSolarSailForce(
        const double frontEmissivityCoefficient,
        const double backEmissivityCoefficient,
        const double frontLambertianCoefficient,
        const double backLambertianCoefficient,
        const double reflectivityCoefficient,
        const double specularReflectionCoefficient,
        const Eigen::Vector3d& normalisedVectorToSource,
        const Eigen::Vector3d& normalisedVelocityVector,
        const double radiationPressure,
        const double area,
        const double coneAngle,
        const double clockAngle)
{

    // Define constant A to simplify the formulae of the solar sail model.
    double A = frontLambertianCoefficient * reflectivityCoefficient * ( 1.0 - specularReflectionCoefficient )
            + ( 1.0 - reflectivityCoefficient ) * ( frontEmissivityCoefficient * frontLambertianCoefficient - backEmissivityCoefficient * backLambertianCoefficient)
            / ( frontEmissivityCoefficient + backEmissivityCoefficient );

    // Compute cosinus and sinus of the sail cone angle.
    double cosinusConeAngle = cos( coneAngle );
    double sinusConeAngle = sin( coneAngle );

    // Compute force magnitude.
    double forceMagnitude = radiationPressure * area * cosinusConeAngle
            * std::sqrt( std::pow( cosinusConeAngle * ( 1.0 + reflectivityCoefficient * specularReflectionCoefficient ) + A, 2 )
                         + std::pow( ( 1.0 - specularReflectionCoefficient * reflectivityCoefficient ) * sinusConeAngle , 2 ) );


    // Compute phi and theta angles, to later define the force direction.
    double phi = atan2( ( ( 1.0 - specularReflectionCoefficient * reflectivityCoefficient ) * cosinusConeAngle * sinusConeAngle ),
                        ( ( 1.0 + reflectivityCoefficient * specularReflectionCoefficient ) * std::pow( cosinusConeAngle, 2 )
                          + A * cosinusConeAngle ) );

    double theta = coneAngle - phi;


    // Compute the normalised direction vector of the force defined in the (r, theta, k) frame.
    Eigen::Vector3d forceDirectionLocalFrame = ( Eigen::Vector3d() << cos( theta ), sin( theta ) * sin( clockAngle ),
                                                 sin( theta ) * cos( clockAngle ) ).finished();

    // Compute normalised vector from source to target.
    Eigen::Vector3d normalisedVectorFromSource = - normalisedVectorToSource;


    // Compute the rotation matrix from local to inertial reference frame.
    Eigen::Matrix3d rotationMatrixFromLocalToInertial;

    rotationMatrixFromLocalToInertial(0,0) = normalisedVectorFromSource[0];
    rotationMatrixFromLocalToInertial(0,1) =
            normalisedVelocityVector[0] * ( std::pow( normalisedVectorFromSource[2], 2 ) + std::pow( normalisedVectorFromSource[1], 2 ) )
            - normalisedVectorFromSource[0] *
            ( normalisedVectorFromSource[1] * normalisedVelocityVector[1] + normalisedVectorFromSource[2] * normalisedVelocityVector[2] );
    rotationMatrixFromLocalToInertial(0,2) = normalisedVectorFromSource.cross( normalisedVelocityVector )[0];

    rotationMatrixFromLocalToInertial(1,0) = normalisedVectorFromSource[1];
    rotationMatrixFromLocalToInertial(1,1) =
            normalisedVelocityVector[1] *  ( std::pow( normalisedVectorFromSource[2], 2 ) + std::pow( normalisedVectorFromSource[0], 2 ) )
            - normalisedVectorFromSource[1] *
            ( normalisedVectorFromSource[0] * normalisedVelocityVector[0] + normalisedVectorFromSource[2] * normalisedVelocityVector[2] );
    rotationMatrixFromLocalToInertial(1,2) = normalisedVectorFromSource.cross( normalisedVelocityVector )[1];

    rotationMatrixFromLocalToInertial(2,0) = normalisedVectorFromSource[2];
    rotationMatrixFromLocalToInertial(2,1) =
            normalisedVelocityVector[2] * ( std::pow( normalisedVectorFromSource[1], 2 ) + std::pow( normalisedVectorFromSource[0], 2 ) )
            - normalisedVectorFromSource[2] *
            ( normalisedVectorFromSource[1] * normalisedVelocityVector[1] + normalisedVectorFromSource[0] * normalisedVelocityVector[0] );
    rotationMatrixFromLocalToInertial(2,2) = normalisedVectorFromSource.cross( normalisedVelocityVector )[2];

    // Compute force direction in inertial frame.
    Eigen::Vector3d forceDirection = rotationMatrixFromLocalToInertial * forceDirectionLocalFrame;

    // Return solar sail force in inertial frame.
    return forceMagnitude * forceDirection;
}

} // namespace electro_magnetism
} // namespace tudat
