/*    Copyright (c) 2010-2013, Delft University of Technology
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
 *      110519    F.M. Engelen      Creation of code.
 *      110628    K. Kumar          Minor comment and layout changes; changed
 *                                  input arguments to pass-by-reference.
 *      110701    K. Kumar          Removed function definitions to be implemented in a future
 *                                  code-check; updated file path.
 *      110718    F.M. Engelen      Added Quaternion transformations and ItoE tranformation.
 *      110726    K. Kumar          Minor modifications.
 *      110809    F.M. Engelen      Replaced the normal planet fixed, with local vertical frame.
 *      110824    J. Leloux         Corrected doxygen documentation.
 *      120530    E.A.G. Heeren     Namespace update.
 *      120614    P. Musegaas       Corrected include guard.
 *
 *    References
 *      Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *
 *    Notes
 *
 */

#ifndef TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H
#define TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H

#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{
namespace reference_frames
{

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
/*!
 * Returns tranformation matrix from rotating planetocentric reference frame (R) to inertial
 * reference frame (I).
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Reference frame (R) to inertial reference frame (I) transformation matrix.
 */
Eigen::Matrix3d getRotatingPlanetocentricToInertialFrameTransformationMatrix(
        double angleFromXItoXR );

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
/*!
 * Returns tranformation quaternion from rotating planetocentric reference frame (R) to inertial
 * reference frame (I). It's an eigen library transformation and can be applied directly to a
 * vector Vector_new = Quaternion * Vector_old.
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Reference frame (R) to inertial reference frame (I) transformation quaternion.
 */
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
    double angleFromXItoXR );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
/*!
 * Returns transformation matrix from inertial referenceframe (I) to the rotating planetocentric
 * reference frame (R).
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Inertial (I) to planetocentric reference frame (R) transformation matrix.
 */
Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix( double angleFromXItoXR );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion quaternion.
/*!
 * Returns transformation quaternion from inertial referenceframe (I) to the rotating
 * planetocentric reference frame (R). It's an eigen library transformation and can be applied
 * directly to a vector Vector_new = Quaternion * Vector_old.
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Inertial (I) to planetocentric reference frame (R) transformation quaternion.
 */
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        double angleFromXItoXR );

//! Create a Quaterniond rotation state object from four quaternion values in vector 4d.
/*!
 * Creates a Quaterniond rotation state object from four quaternion values.
 * This function is not related to any specific rotation matrix, but can be used for general
 * purposes. It's an eigen library transformation and can be applied
 * directly to a vector Vector_new = Quaternion * Vector_old.
 * Note that is also possible to create a quaternion object directly from a Vector4d,
 * but Eigen will rearrange the order of the coefficients ([q2 q3 q4 q1]). This function
 * retreives the individual entries of the Vector4d an uses them as four doubles as input
 * arguments for the constructor of the Quateriond. From the Eigen code documentation:
 * \warning Note the order of the arguments: the real \a w coefficient first,
 *              while internally the coefficients are stored in the following order:
 *              [\c x, \c y, \c z, \c w]
 * \param vectorWithQuaternion A vector containing the quaternions of the rotation state
 * \return Transformation quaternion.
 */
Eigen::Quaterniond getQuaternionObjectFromQuaternionValues(
        const Eigen::Vector4d& vectorWithQuaternion );

//! Get Aerodynamic (airspeed-based) (AA) to body reference frame (B) tranformation matrix.
/*!
 * Returns transformation matrix from Aerodynamic (airspeed-based) (AA) to body reference
 * frame (B).
 * \param angleOfAttack Angle of attack [rad].
 * \param angleOfSideslip Angle of sideslip [rad].
 * \return Transformation matrix.
 */
Eigen::Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
        double angleOfAttack, double angleOfSideslip );

//! Get Aerodynamic (airspeed-based) (AA) to body reference frame (B) tranformation quaternion.
/*!
 * Returns transformation quaternion from Aerodynamic (airspeed-based) (AA) to body reference
 * frame (B).
 * \param angleOfAttack Angle of attack [rad].
 * \param angleOfSideslip Angle of sideslip [rad].
 * \return Transformation quaternion.
 */
Eigen::Quaterniond getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
        double angleOfAttack, double angleOfSideslip );

//! Get transformation quaternion from Planetocentric (R) to the Local vertical (V) frame.
/*!
 * Returns the frame transformation quaternion from the Planetocentric (R) to the Local vertical
 * (V) reference frame. The Z axis is alligned with the local gravity vector. Whether or not,
 * this is in the direction of the center of the central body, depedens which kind of latitude
 * is provided (geocentric, geodetic, gravitation latitude)
 * The X axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The planetocentric latitude.
 * \return Transformation quaternion from Planetocentric (R) to the local vertical (V) frame.
 */
Eigen::Quaterniond getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
    double longitude, double latitude );

//! Get transformation quaternion from local vertical (V) to the Planetocentric frame (R).
/*!
 * Returns the frame transformation quaternion from the local vertical (V) to the
 * Planetocentric frame (R)reference frame. The Z axis is alligned with the local gravity vector.
 * Whether or not, this is in the direction of the center of the central body, depedens which kind
 * of latitude is provided (geocentric, geodetic, gravitation latitude) The X-axis is directed to
 * the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The planetocentric latitude in [rad].
 * \return Transformation quaternion from local vertical (V) to the Planetocentric (R) frame.
 */
Eigen::Quaterniond getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
    double longitude, double latitude );

} // namespace reference_frames
} // namespace tudat

#endif // TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H
