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
 *      110519    F.M. Engelen      File created.
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
 *      130121    K. Kumar          Updated functions to be const-correct.
 *      130219    D. Dirkx          Migrated from personal code.
 *      130312    A. Ronse          Added V-T, TA-AA and AA-B transformations.
 *
 *
 *    References
 *      Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *      Seidelmann, P. K. (Ed.). (2005). Explanatory supplement to the astronomical almanac.
 *              Univ Science Books.
 *
 *    Notes
 *
 */

#ifndef TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H
#define TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H

#include <cmath>
#include <vector>

#include <boost/function.hpp>

#include <Eigen/Core>
#include <Eigen/Geometry>

namespace tudat
{
namespace reference_frames
{

//! Wrapper function to transform a vector to a different frame from a single rotation function.
/*!
 * Wrapper function to transform a vector to a different frame from a single rotation function.
 * \param originalVector Vector that is to be rotated to a new frame
 * \param rotation Function returning the current rotation to the new frame
 * \return Vector originalVector, transformed to new frame.
 */
Eigen::Vector3d transformVector(
        const Eigen::Vector3d& originalVector,
        const boost::function< Eigen::Quaterniond( ) > rotation );

//! Wrapper function to transform a vector to a different frame from a single transformation function.
/*!
 * Wrapper function to transform a vector to a different frame from a single transformation function.
 * \param originalVector Vector that is to be transformed to a new frame
 * \param transformationFunction Function transforming a vector to a new frame
 * \return Vector originalVector, transformed to new frame.
 */
Eigen::Vector3d transformVector(
        const boost::function< Eigen::Vector3d( ) > originalVector,
        const boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > transformationFunction );

//! Wrapper function to transform a vector to a different frame from a list of transformation function.
/*!
 * Wrapper function to transform a vector to a different frame from a list of transformation function.
 * \param originalVector Vector that is to be transformed to a new frame
 * \param rotationsList List of transformation function, each of which transforms a vector to a new frame. The functions
 * in this list are called in descending order.
 * \return Vector originalVector, transformed to new frame.
 */
Eigen::Vector3d transformVector(
        const Eigen::Vector3d& originalVector,
        const std::vector< boost::function< Eigen::Vector3d( const Eigen::Vector3d& ) > >& rotationsList );

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
        const double angleFromXItoXR );

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
    const double angleFromXItoXR );

//! Get rotation from planet-fixed to inertial frame.
/*!
 * Retuns rotation from planet-fixed to inertial frame, assuming that the equatorial plane is not
 * equal to x-y plane of inertial frame. Orientation of body-fixed frame is obtained from right
 * ascension and declination of body's pole and the location of the prime meridian
 * (Seidelmann et al. 2005).
 * \param declinationOfPole Declination of body's pole in inertial frame.
 * \param rightAscensionOfPole Right ascension of body's pole in inertial frame.
 * \param longitudeOfPrimeMeridian Longitude of body prime meridian.
 * \return Rotation quaternion computed.
 */
Eigen::Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
/*!
 * Returns transformation matrix from inertial referenceframe (I) to the rotating planetocentric
 * reference frame (R).
 * \param angleFromXItoXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Inertial (I) to planetocentric reference frame (R) transformation matrix.
 */
Eigen::Matrix3d getInertialToPlanetocentricFrameTransformationMatrix(
        const double angleFromXItoXR );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformation quaternion.
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
        const double angleFromXItoXR );

//! Get rotation from inertial to planet-fixed frame.
/*!
 * Returns rotation from inertial to planet-fixed frame, assuming that the equatorial plane is not
 * equal to x-y plane of inertial frame. Orientation of body-fixed frame is obtained from right
 * ascension and declination of body's pole and the location of the prime meridian
 * (Seidelmann et al. 2005).
 * \param declinationOfPole Declination of body's pole in inertial frame.
 * \param rightAscensionOfPole Right ascension of body's pole in inertial frame.
 * \param longitudeOfPrimeMeridian Longitude of body prime meridian.
 * \return Rotation quaternion computed.
 */
Eigen::Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        const double declinationOfPole,
        const double rightAscensionOfPole,
        const double longitudeOfPrimeMeridian );

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

//! Get transformation matrix from Planetocentric (R) to the Local vertical (V) frame.
/*!
 * Returns the frame transformation matrix from the Planetocentric (R) to the Local vertical
 * (V) reference frame. The Z-axis is aligned with the local gravity vector. Whether or not,
 * this is in the direction of the center of the central body, depends which kind of latitude
 * is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation matrix from Planetocentric (R) to the local vertical (V) frame.
 */
Eigen::Matrix3d getRotatingPlanetocentricToLocalVerticalFrameTransformationMatrix(
    const double longitude, const double latitude );

//! Get transformation quaternion from Planetocentric (R) to the local vertical (V) frame.
/*!
 * Returns the frame transformation quaternion from the Planetocentric (R) to the Local vertical
 * (V) reference frame. The Z-axis is aligned with the local gravity vector. Whether or not,
 * this is in the direction of the center of the central body, depends which kind of latitude
 * is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation quaternion from Planetocentric (R) to the local vertical (V) frame.
 */
Eigen::Quaterniond getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
    const double longitude, const double latitude );

//! Get transformation matrix from local vertical (V) to the Planetocentric frame (R).
/*!
 * Returns the frame transformation matrix from the local vertical (V) to the
 * Planetocentric frame (R) reference frame. The Z-axis is aligned with the local gravity vector.
 * Whether or not, this is in the direction of the center of the central body, depends which kind
 * of latitude is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation matrix from local vertical (V) to the Planetocentric (R) frame.
 */
Eigen::Matrix3d getLocalVerticalToRotatingPlanetocentricFrameTransformationMatrix(
    const double longitude, const double latitude );

//! Get transformation quaternion from local vertical (V) to the Planetocentric frame (R).
/*!
 * Returns the frame transformation quaternion from the local vertical (V) to the 
 * Planetocentric (R) reference frame. The Z-axis is aligned with the local gravity vector.
 * Whether or not, this is in the direction of the center of the central body, depends which kind
 * of latitude is provided (geocentric, geodetic, gravitation latitude).
 * The X-axis is directed to the north.
 * \param longitude The longitude in the planetocentric reference frame in [rad].
 * \param latitude The latitude in the planetocentric reference frame in [rad].
 * \return Transformation quaternion from local vertical (V) to the Planetocentric (R) frame.
 */
Eigen::Quaterniond getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
    const double longitude, const double latitude );

//! Get transformation matrix from the TA/TG to the V-frame.
/*!
 * Returns the frame transformation matrix from the trajectory frame (T) to the local vertical (V)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation matrix from trajectory frame (T) to the local vertical frame (V).
 */
Eigen::Matrix3d getTrajectoryToLocalVerticalFrameTransformationMatrix(
        const double flightPathAngle, const double headingAngle );

//! Get transformation quaternion from the TA/TG to the V-frame.
/*!
 * Returns the frame transformation quaternion from the trajectory (T) to the local vertical (V)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation quaternion from trajectory frame (T) to the local vertical frame (V).
 */
Eigen::Quaterniond getTrajectoryToLocalVerticalFrameTransformationQuaternion(
        const double flightPathAngle, const double headingAngle );

//! Get transformation matrix from the local V- to TA/TG-frame.
/*!
 * Returns the frame transformation matrix from the local vertical (V) to the trajectory (T)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation matrix from the local vertical (V) to the trajectory (T) frame.
 */
Eigen::Matrix3d getLocalVerticalFrameToTrajectoryTransformationMatrix(
        const double flightPathAngle, const double headingAngle );

//! Get transformation quaternion from V- to the TA/TG-frame.
/*!
 * Returns the transformation quaternion from the local vertical (V) to the trajectory (T)
 * reference frame. Depending on whether the flight-path angle and heading angle express the
 * velocity relative to the rotating planetocentric frame or relative to the local atmosphere,
 * the resulting rotation holds for the groundspeed-based trajectory (TG) or airspeed-based
 * trajectory (TA) frame respectively.
 * \param flightPathAngle The trajectory's flight-path angle, positive upwards, in [rad].
 * \param headingAngle The trajectory's heading angle with respect to the North in [rad].
 * \return Transformation quaternion from the local vertical (V) to the trajectory (T) frame.
 */
Eigen::Quaterniond getLocalVerticalFrameToTrajectoryTransformationQuaternion(
        const double flightPathAngle, const double headingAngle );

//! Get transformation matrix from the TA- to the AA-frame.
/*!
 * Returns the transformation matrix from the airspeed-based trajectory (TA) to the airspeed-based
 * aerodynamic frame (AA). These frames differ from each other only by the bank-angle, representing
 * one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation matrix from the TA- to the AA-frame.
 */
Eigen::Matrix3d getTrajectoryToAerodynamicFrameTransformationMatrix( const double bankAngle );

//! Get transformation quaternion from the TA- to the AA-frame.
/*!
 * Returns the transformation quaternion from the airspeed-based trajectory (TA) to the airspeed-
 * based aerodynamic frame (AA). These frames differ from each other only by the bank-angle,
 * representing one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation quaternion from the TA- to the AA-frame.
 */
Eigen::Quaterniond getTrajectoryToAerodynamicFrameTransformationQuaternion( 
    const double bankAngle );

//! Get transformation matrix from the AA- to the TA-frame.
/*!
 * Returns the transformation matrix from the airspeed-based aerodynamic (AA) to the airspeed-based
 * trajectory frame (TA). These frames differ from each other only by the bank-angle, representing
 * one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation matrix from the AA- to the TA-frame.
 */
Eigen::Matrix3d getAerodynamicToTrajectoryFrameTransformationMatrix( const double bankAngle );

//! Get transformation quaternion from the AA- to the TA-frame.
/*!
 * Returns the transformation quaternion from the airspeed-based aerodynamic (AA) to the airspeed-
 * based trajectory frame (TA). These frames differ from each other only by the bank-angle,
 * representing one Euler-rotation.
 * \param bankAngle The object's bank angle in [rad].
 * \return Transformation quaternion from the AA- to the TA-frame.
 */
Eigen::Quaterniond getAerodynamicToTrajectoryFrameTransformationQuaternion(
        const double bankAngle );

//! Get transformation matrix fom the B- to the AA-frame.
/*!
 * Returns the transformation matrix from the body-fixed (B) to the airspeed-based aerodynamic
 * frame (AA).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation matrix from the B- to the AA-frame.
 */
Eigen::Matrix3d getBodyToAirspeedBasedAerodynamicFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip );

//! Get transformation quaternion fom the B- to the AA-frame.
/*!
 * Returns the transformation quaternion from the body-fixed (B) to the airspeed-based aerodynamic
 * frame (AA).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation quaternion from the B- to the AA-frame.
 */
Eigen::Quaterniond getBodyToAirspeedBasedAerodynamicFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip );

//! Get transformation matrix fom the AA- to the B-frame.
/*!
 * Returns the transformation matrix from the airspeed-based aerodynamic (AA) to the body-fixed
 * frame (B).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation matrix from the AA- to the B-frame.
 */
Eigen::Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
        const double angleOfAttack, const double angleOfSideslip );

//! Get transformation quaternion fom the AA- to the B-frame.
/*!
 * Returns the transformation quaternion from the airspeed-based aerodynamic (AA) to the body-fixed
 * frame (B).
 * \param angleOfAttack The angle of attack in [rad].
 * \param angleOfSideslip The angle of sideslip in [rad].
 * \return Transformation quaternion from the AA- to the B-frame.
 */
Eigen::Quaterniond getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
        const double angleOfAttack, const double angleOfSideslip );

//! Calculate current heading angle.
/*!
 * Calculate heading angle from velocity in vertical (LVLH) frame.
 * \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 * \return Current heading angle.
 */
double calculateHeadingAngle( const Eigen::Vector3d& velocityInVerticalFrame );

//! Calculate current flight path angle. Angle is defined positive upwards.
/*!
 *  Calculate flight path angle from velocity in vertical (LVLH) frame.
 *  Angle is defined positive upwards.
 *  \param velocityInVerticalFrame Current Cartesian velocity in vertical frame.
 *  \return Current flight path angle.
 */
double calculateFlightPathAngle( const Eigen::Vector3d& velocityInVerticalFrame );

} // namespace reference_frames
} // namespace tudat

#endif // TUDAT_REFERENCE_FRAME_TRANSFORMATIONS_H
