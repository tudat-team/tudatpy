/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *
 *    References
 */

// Temporary notes (move to class/function doxygen):
// Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
// 
// Because potential speed improvement it was chosen to use AngleAxisd and quaternions
// but to get things working, the rotation angle inputted in angleAxisd need to be inverted.
// In the future it might be better to change it to write out the complete transformation for
// clarity, or work with directional cosine matrices.
// 

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>
#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
Eigen::Matrix3d reference_frame_transformations::
getRotatingPlanetocentricToInertialFrameTransformationMatrix( double angleFromXItoXR )
{
    // Declare local variables.
    // Declare local matrix.
    Eigen::Matrix3d localMatrix_;

    // Set local matrix.
    localMatrix_ = reference_frame_transformations::
            getInertialToPlanetocentricFrameTransformationMatrix( angleFromXItoXR );

    // Return transformation matrix.
    return localMatrix_.transpose( );
}

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
Eigen::Quaterniond reference_frame_transformations::
getRotatingPlanetocentricToInertialFrameTransformationQuaternion( double angleFromXItoXR )
{
    // Compute transformation quaternion
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * -angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion matrix.
Eigen::Matrix3d reference_frame_transformations::
getInertialToPlanetocentricFrameTransformationMatrix( double angleFromXItoXR )
{
    // Compute rotation about Z-Axis.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return eigenRotationObject.toRotationMatrix( );
}

//! Get inertial (I) to rotating planetocentric (R) reference frame transformtion quaternion.
Eigen::Quaterniond reference_frame_transformations::
getInertialToPlanetocentricFrameTransformationQuaternion( double angleFromXItoXR )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd eigenRotationObject = Eigen::AngleAxisd( -1.0 * angleFromXItoXR,
                                                               Eigen::Vector3d::UnitZ( ) );

    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond( eigenRotationObject );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Create a Quaterniond rotation state object from four quaternion values in a Vector4d
Eigen::Quaterniond reference_frame_transformations::
getQuaternionObjectFromQuaternionValues( const Eigen::Vector4d& vectorWithQuaternion )
{
    // Set transformation quaternion.
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                vectorWithQuaternion( 0 ), vectorWithQuaternion( 1 ),
                vectorWithQuaternion( 2 ), vectorWithQuaternion( 3 ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get Aerodynamic (airspeed-based) (AA) to body reference frame (B) tranformation matrix.
Eigen::Matrix3d reference_frame_transformations::
getAirspeedBasedAerodynamicToBodyFrameTransformationMatrix( double angleOfAttack,
                                                            double angleOfSideslip )
{
    // Declare local variables.
    // Declare local transformation matrix.
    Eigen::Matrix3d transformationMatrix_;

    // Compute rotation by side-slip angle about Z-Axis, followed by rotation negative angle of
    // attack angle about Y-Axis.
    // Note the sign change, because how angleAxisd is defined.
    transformationMatrix_ = Eigen::AngleAxisd( -1.0 * angleOfAttack, Eigen::Vector3d::UnitY( ) )
            * Eigen::AngleAxisd( -1.0 * -angleOfSideslip, Eigen::Vector3d::UnitZ( ) );

    // Return transformation matrix.
    return transformationMatrix_;
}

//! Get transformation quaternion from Planetocentric (R) to the Local vertical (V) frame.
Eigen::Quaterniond reference_frame_transformations::
getRotatingPlanetocentricToLocalVerticalFrameTransformationQuaternion(
        double longitude, double latitude )
{
    // Compute transformation quaternion.
    // Note the sign change, because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * longitude, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundYaxis = Eigen::AngleAxisd(
                -1.0 * ( -latitude - mathematics::PI / 2.0 ), Eigen::Vector3d::UnitY( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundYaxis * RotationAroundZaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

//! Get transformation quaternion from local vertical (V) to the Planetocentric frame (R).
Eigen::Quaterniond reference_frame_transformations::
getLocalVerticalToRotatingPlanetocentricFrameTransformationQuaternion(
    double longitude, double latitude )
{
    // Compute transformation quaternion.
    // Note the sign change (-1.0), because how angleAxisd is defined.
    Eigen::AngleAxisd RotationAroundZaxis = Eigen::AngleAxisd(
                -1.0 * -longitude, Eigen::Vector3d::UnitZ( ) );
    Eigen::AngleAxisd RotationAroundYaxis =
            Eigen::AngleAxisd( -1.0 * ( latitude + mathematics::PI / 2.0 ),
                               Eigen::Vector3d::UnitY( ) );
    Eigen::Quaterniond frameTransformationQuaternion = Eigen::Quaterniond(
                ( RotationAroundZaxis * RotationAroundYaxis ) );

    // Return transformation quaternion.
    return frameTransformationQuaternion;
}

} // namespace tudat
