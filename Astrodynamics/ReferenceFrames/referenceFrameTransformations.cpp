/*! \file referenceFrameTransformations.cpp
 *    This file contains the implementation of the reference frame transformations
 *    namespace included in Tudat.
 *
 *    Path              : /Astrodynamics/ReferenceFrames/
 *    Version           : 6
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
 *    Last modified     : 26 July, 2011
 *
 *    References
 *      Muller, J.A., et al. Flight Dynamics Lecture Notes, TU Delft, February 2007.
 *      Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *
 *    Notes
 *    Because potential speed improvement it was chosen to use AngleAxisd and quaternions
 *    but to get things working, the rotation angle inputted in angleAxisd need to be inverted.
 *    In the future it might be better to change it to write out the complete transformation for
 *    clearity, or work with directional cosine matrices.
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
 *      110701    K. Kumar          Updated file path.
 *      110718    F.M. Engelen      Repaired incorrect sign for angle in R2I and I2R
 *                                  Added Quaternion transformations and ItoE tranformation.
 *      110726    K. Kumar          Minor modifications.
 *      110809    F.M. Engelen      Applied the minus one correction for angleAxisD,
 *                                  changed to local vertical frame.
 */

// Include statements.
#include "referenceFrameTransformations.h"

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
Matrix3d reference_frame_transformations::
getRotatingPlanetocentricToInertialFrameTransformationMatrix(
        const double& angleFromXItoXR )
{
    // Declare local variables.
    // Declare local matrix.
    Matrix3d localMatrix_;

    // Set local matrix.
    localMatrix_ = reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationMatrix( angleFromXItoXR );

    // Return transformation matrix.
    return localMatrix_.transpose();
}

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
Quaterniond reference_frame_transformations::
getRotatingPlanetocentricToInertialFrameTransformationQuaternion(
    const double& angleFromXItoXR )
{
    // Compute transformation quaternion
    // Note the sign change, because how angleAxisd is defined.
    AngleAxisd eigenRotationObject = AngleAxisd( -1.0 * -angleFromXItoXR, Vector3d::UnitZ( ) );
    Quaterniond frameTransformationQuaternion = Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
Matrix3d reference_frame_transformations::getInertialToPlanetocentricFrameTransformationMatrix(
        const double& angleFromXItoXR )
{
    // Compute rotation about Z-Axis.
     // Note the sign change, because how angleAxisd is defined.
    AngleAxisd eigenRotationObject = AngleAxisd( -1.0 * angleFromXItoXR, Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return eigenRotationObject.toRotationMatrix( );
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion quaternion.
Quaterniond
  reference_frame_transformations::getInertialToPlanetocentricFrameTransformationQuaternion(
        const double& angleFromXItoXR )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    AngleAxisd eigenRotationObject = AngleAxisd( -1.0 * angleFromXItoXR, Vector3d::UnitZ( ) );
    Quaterniond  frameTransformationQuaternion = Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Creates a Quaterniond rotation state object from four quaternion values in a Vector4d
Quaterniond reference_frame_transformations::
getQuaternionObjectFromQuaternionValues(
    const Vector4d& vectorWithQuaternion )
{
    // Set transformation quaternion.
    Quaterniond frameTransformationQuaternion = Quaterniond(
            vectorWithQuaternion( 0 ),
            vectorWithQuaternion( 1 ),
            vectorWithQuaternion( 2 ),
            vectorWithQuaternion( 3 ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get Aerodynamic (airspeed based) (AA) to body reference frame (B) tranformation matrix.
Matrix3d reference_frame_transformations::
getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( const double& angleOfAttack,
                                                            const double& angleOfSideslip )
{
    // Compute rotation by side-slip angle about Z-Axis, followed by rotation negative angle of
    // attack angle about Y-Axis.
    // Note the sign change, because how angleAxisd is defined.
    AngleAxisd eigenRotationObject = AngleAxisd( -1.0 * angleOfAttack, Vector3d::UnitY( ) ) *
                            AngleAxisd( -1.0 * -angleOfSideslip, Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return eigenRotationObject.toRotationMatrix( );
}

//! Get transformation quaternion from Planetocentric (R) to the Local vertical (V) frame.
Quaterniond reference_frame_transformations::
        getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
            const double& longitude, const double& latitude )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    AngleAxisd RotationAroundZaxis = AngleAxisd( -1.0 * longitude, Vector3d::UnitZ( ) );
    AngleAxisd RotationAroundYaxis = AngleAxisd(
            -1.0 * ( -latitude - M_PI / 2.0 ), Vector3d::UnitY( ) );
    Quaterniond frameTransformationQuaternion = Quaterniond(
            ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation quaternion from local vertical (V) to the Planetocentric frame (R).
Quaterniond reference_frame_transformations::
        getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
            const double& longitude, const double& latitude )
{
    // Compute transformation quaternion.
    // Note the sign change (-1.0), because how angleAxisd is defined.
    AngleAxisd RotationAroundZaxis = AngleAxisd( -1.0 * -longitude, Vector3d::UnitZ( ) );
    AngleAxisd RotationAroundYaxis =
            AngleAxisd( -1.0 * ( latitude + M_PI / 2.0 ), Vector3d::UnitY( ) );
    Quaterniond  frameTransformationQuaternion = Quaterniond(
            ( RotationAroundZaxis * RotationAroundYaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

// End of file.
