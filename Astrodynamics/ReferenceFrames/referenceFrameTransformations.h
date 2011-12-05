/*! \file referenceFrameTransformations.h
 *    This file contains the definition of the reference frame transformations namespace included
 *    in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 7
 *    Check status      : Checked
 *
 *    Checker           : F. M. Engelen
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : F.M.Engelen@student.tudelft.nl
 *
 *    Checker           : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Date created      : 19 May, 2011
 *    Last modified     : 24 August, 2011
 *
 *    References
 *      Muller, J.A., et al. Flight Dynamics Lecture Notes, TU Delft, February 2007.
 *      Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110519    F.M. Engelen      First creation of code.
 *      110628    K. Kumar          Minor comment and layout changes; changed
 *                                  input arguments to pass-by-reference.
 *      110701    K. Kumar          Removed function definitions to be implemented in a future
 *                                  code-check; updated file path.
 *      110718    F.M. Engelen      Added Quaternion transformations and ItoE tranformation.
 *      110726    K. Kumar          Minor modifications.
 *      110809    F.M. Engelen      Replaced the normal planet fixed, with local vertical frame.
 *      110824    J. Leloux         Corrected doxygen documentation.
 */

#ifndef FRAMETRANSFORMATION_H
#define FRAMETRANSFORMATION_H

// Include statements.
#include <cmath>
#include "Mathematics/LinearAlgebra/linearAlgebra.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

// Using declarations.
using Eigen::AngleAxisd;
using Eigen::Quaterniond;

namespace reference_frame_transformations
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
Matrix3d getRotatingPlanetocentricToInertialFrameTransformationMatrix(
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
Quaterniond getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
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
Matrix3d getInertialToPlanetocentricFrameTransformationMatrix( double angleFromXItoXR );

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
Quaterniond getInertialToPlanetocentricFrameTransformationQuaternion(
        double angleFromXItoXR );

//! Create a Quaterniond rotation state object from four quaternion values in vector 4d.
/*!
 * Creates a Quaterniond rotation state object from four quaternion values.
 * This function is not related to any specific rotation matrix, but can be used for general
 * purposes. It's an eigen library transformation and can be applied
 * directly to a vector Vector_new = Quaternion * Vector_old.
 * Note that is also possible to create a quaternion object directly from a vector4d,
 * but Eigen will rearrange the order of the coefficients ([q2 q3 q4 q1]). This function
 * retreives the individual entries of the Vector4d an uses them as four doubles as input
 * arguments for the constructor of the Quateriond. From the Eigen code documentation:
 * \warning Note the order of the arguments: the real \a w coefficient first,
 *              while internally the coefficients are stored in the following order:
 *              [\c x, \c y, \c z, \c w]
 * \param vectorWithQuaternion A vector containing the quaternions of the rotation state
 * \return Transformation quaternion.
 */
Quaterniond getQuaternionObjectFromQuaternionValues( const Vector4d& vectorWithQuaternion );

//! Get Aerodynamic (airspeed-based) (AA) to body reference frame (B) tranformation matrix.
/*!
 * Returns transformation matrix from Aerodynamic (airspeed-based) (AA) to body reference
 * frame (B).
 * \param angleOfAttack Angle of attack [rad].
 * \param angleOfSideslip Angle of sideslip [rad].
 * \return Transformation matrix.
 */
Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix(
        double angleOfAttack, double angleOfSideslip );

//! Get Aerodynamic (airspeed-based) (AA) to body reference frame (B) tranformation quaternion.
/*!
 * Returns transformation quaternion from Aerodynamic (airspeed-based) (AA) to body reference
 * frame (B).
 * \param angleOfAttack Angle of attack [rad].
 * \param angleOfSideslip Angle of sideslip [rad].
 * \return Transformation quaternion.
 */
Quaterniond getAirspeedBasedAerodynamicToBodyFrameTransformationQuaternion(
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
Quaterniond getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
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
Quaterniond getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
    double longitude, double latitude );

}

}

#endif // FRAMETRANSFORMATION_H

// End of file.
