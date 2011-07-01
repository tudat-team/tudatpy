/*! \file referenceFrameTransformations.h
 *    This file contains the definition of the reference frame transformations
 *    namespace included in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 3
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
 *    Last modified     : 1 July, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
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
 *      110701    K. Kumar          Removed function definitions to be
 *                                  implemented in a future code-checkl;
 *                                  updated file path.
 */

#ifndef FRAMETRANSFORMATION_H
#define FRAMETRANSFORMATION_H

// Include statements.
#include <cmath>
#include "linearAlgebra.h"

// Using statments.
using Eigen::AngleAxisd;

namespace reference_frame_transformations
{

//! Get rotating planetocentric (R) to intertial (I) reference frame transformation matrix.
/*!
 * Returns tranformation matrix from rotating planetocentric reference frame (R) to inertial
 * reference frame (I).
 * \param  angleBetweenXIandXR Angle between X-axis of the planetocentric reference frame and
 *          X-axis of the inertial reference frame in [rad]. This angle is same as
 *          the rotational rate of the central body [rad/s] times the time from epoch [s].
 * \return Reference frame (R) to inertial reference frame (I) transformation matrix.
 */
Matrix3d getRotatingPlanetocentricToInertialFrameTransformation(
    const double& angleBetweenXIandXR );

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
/*!
  * Returns transformation matrix from inertial referenceframe (I) to the rotating planetocentric
  * reference frame (R).
  * \param  angleBetweenXIandXR Angle between X-axis of the planetocentric reference frame and
  *             X-axis of the inertial reference frame in [rad]. This angle is same as
  *             the rotational rate of the central body [rad/s] times the time from epoch [s].
  * \return Inertial (I) to planetocentric reference frame (R) transformation matrix.
  */
Matrix3d getInertialToPlanetocentricFrameTransformation( const double& angleBetweenXIandXR );

//! Get Aerodynamic (airspeed based) (AA) to body reference frame (B) tranformation matrix.
/*!
 * Returns transformation matrix from Aerodynamic (airspeed based) (AA) to body reference
 * frame (B).
 * \param angleOfAttack Angle of attack [rad].
 * \param angleOfSideslip Angle of sideslip [rad].
 * \return Transformation matrix.
 */
Matrix3d getAirspeedBasedAerodynamicToBodyFrameTransformation(
        const double& angleOfAttack, const double& angleOfSideslip );

}

#endif // FRAMETRANSFORMATION_H

// End of file.
