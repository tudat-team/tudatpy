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
 *    Murray, G. Rotation Matrices and Formulas java script, RotationMatrix.java available at
 *        https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas, 2011.
 *        last accessed: 20th January, 2014.
 *    Wikipedia. Quaternions and spatial rotation,
 *       http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation, last access: 4 December,
 *       2013.
 *
 */

#ifndef TUDAT_ROTATION_ABOUT_ARBITRARY_AXIS_H
#define TUDAT_ROTATION_ABOUT_ARBITRARY_AXIS_H

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{

//! Compute rotation of point about arbitrary axis.
/*!
 * Computes the rotation of a point about an arbitrary axis. This function computes the position
 * of the rotated point by computing the displacement vector between the origin of rotation and
 * the initalPosition of point and then applying a matrix as shown in section 6 of
 * http://inside.mines.edu/~gmurray/ArbitraryAxisRotation/ArbitraryAxisRotation.html.
 * A visualization of the computation being performed can be found at
 * http://twist-and-shout.appspot.com/. Units for originOfRotation, and initialPositionOfPoint must
 * be the same. All positions are assumed to be with respect to a common arbitrary origin. Note
 * that this arbitrary origin is not necessarily the origin of rotation.
 * \param originOfRotation Position of tail of axis of rotation vector with respect to a chosen
 *        arbitrary origin.
 * \param angleOfRotation Angle about which the point rotates with respect to the axis of rotation
 *        [rad].
 * \param axisOfRotation Unit vector pointed in the direction of rotation determined by the
 *        right-hand rule.
 * \param initialPositionOfPoint Position of point before rotation with respect to the same chosen
 *        arbitrary origin.
 * \return Position of point after rotation about the axis of rotation with respect to the same
 *         chosen arbitrary origin.
 */
Eigen::Vector3d computeRotationOfPointAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfPoint );

//! Compute rotation of vector about an arbitrary axis.
/*!
 * Computes rotation of vector about an arbitrary axis. This function uses
 * computeRotationOfPointAboutArbitraryAxis to rotate the position of the head and tail of the
 * vector and compute the resultant vector defined as the difference of the two rotated points.
 * Units for originOfRotation, initialPositionOfVectorTail and initialVector must be the same.
 * \param originOfRotation Position of tail of axis of rotation vector with respect to a chosen
 *        arbitrary origin.
 * \param angleOfRotation Angle over which the vector is rotated. A positive angle is determined
 *        using the right hand rule [rad].
 * \param axisOfRotation Unit vector pointing in the direction of rotation determined by the
 *        right-hand rule.
 * \param initialPositionOfVectorTail Initial position of tail of vector which will be rotated with
 *        respect to this same chosen arbitrary origin.
 * \param initialVector Vector to be rotated with respect to this same chosen arbitrary origin.
 * \return Vector after rotation about the arbitrary axis with respect to this same chosen
 *        arbitrary origin.
 * \sa computeRotationOfPointAboutArbitraryAxis.
 */
Eigen::Vector3d computeRotationOfVectorAboutArbitraryAxis(
        const Eigen::Vector3d& originOfRotation,
        const double angleOfRotation,
        const Eigen::Vector3d& axisOfRotation,
        const Eigen::Vector3d& initialPositionOfVectorTail,
        const Eigen::Vector3d& initialVector );

} // namespace basic_mathematics
} // namespace tudat

#endif // TUDAT_ROTATION_ABOUT_ARBITRARY_AXIS_H
