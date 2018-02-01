/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *
 *
 */

#include <cmath>

#include <Eigen/Geometry>

#include "Tudat/Mathematics/BasicMathematics/rotationAboutArbitraryAxis.h"

namespace tudat
{
namespace basic_mathematics
{

//! Compute rotation of point about arbitrary axis
Eigen::Vector3d computeRotationOfPointAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfPoint )
{

    //Declare and initialize rotation matrix
    Eigen::Matrix3d rotationMatrix = Eigen::Matrix3d::Zero( );

    // Compute rotation matrix using AngleAxis object.
    rotationMatrix = Eigen::AngleAxisd( angleOfRotation, axisOfRotation.normalized( ) );

    // Compute initial of position of point with respect to origin of rotation.
    const Eigen::Vector3d initialPositionOfPointWithRespectToOriginOfRotation =
            initialPositionOfPoint - originOfRotation;

    // Compute rotation of point about axis of rotation with respect to origin of rotation.
    const Eigen::Vector3d rotatedPositionWithRespectToOriginOfRotation =
            rotationMatrix * initialPositionOfPointWithRespectToOriginOfRotation;

    //Return position with respect to the chosen arbitrary origin after rotation about
    //arbitrary axis.
    return rotatedPositionWithRespectToOriginOfRotation + originOfRotation;

}

//! Compute rotation of vector about arbitrary axis
Eigen::Vector3d computeRotationOfVectorAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfVectorTail,
        const Eigen::Vector3d& initialVector )
{

    // Compute rotation of the tail of vector. Resulted position is with respect to the chosen
    // arbitrary origin.
    Eigen::Vector3d rotatedPositionOfVectorTail =
            computeRotationOfPointAboutArbitraryAxis( originOfRotation, angleOfRotation,
                                                      axisOfRotation, initialPositionOfVectorTail );

    // Compute rotation of the head of vector. Resulted position is with respect to the chosen
    // arbitrary origin.
    Eigen::Vector3d rotatedPositionOfVectorHead =
            computeRotationOfPointAboutArbitraryAxis( originOfRotation, angleOfRotation,
                                                      axisOfRotation,
                                                      initialPositionOfVectorTail + initialVector );

    // Return rotated vector with respect to the chosen arbitrary origin
    return rotatedPositionOfVectorHead - rotatedPositionOfVectorTail;

}

} // namespace basic_mathematics
} // namespace tudat
