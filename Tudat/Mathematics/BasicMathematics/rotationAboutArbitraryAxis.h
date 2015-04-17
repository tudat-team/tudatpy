/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      131109    S. Hirsh          File Created.
 *      140220    T. Roegiers       Textual changes in Doxygen comments.
 *      140410    S. Hirsh          Minor revisions.
 *      140411    T. Roegiers       Textual changes to Doxygen comments.
 *
 *    References
 *    Murray, G. Rotation Matrices and Formulas java script, RotationMatrix.java available at
 *        https://sites.google.com/site/glennmurray/Home/rotation-matrices-and-formulas, 2011.
 *        last accessed: 20th January, 2014.
 *    Wikipedia. Quaternions and spatial rotation,
 *       http://en.wikipedia.org/wiki/Quaternions_and_spatial_rotation, last access: 4 December,
 *       2013.
 *
 *    Notes
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
