/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      110701    K. Kumar          Updated file path.
 *      110718    F.M. Engelen      Repaired incorrect sign for angle in R2I and I2R
 *                                  Added Quaternion transformations and ItoE tranformation.
 *      110726    K. Kumar          Minor modifications.
 *      110809    F.M. Engelen      Applied the minus one correction for angleAxisD,
 *                                  changed to local vertical frame.
 *      120530    E.A.G. Heeren     Namespace update.
 *
 *    References
 *      Mooij, E. The Motion of a vehicle in a Planetary Atmosphere, TU Delft, 1997.
 *
 *    Because potential speed improvement it was chosen to use AngleAxisd and quaternions
 *    but to get things working, the rotation angle inputted in angleAxisd need to be inverted.
 *    In the future it might be better to change it to write out the complete transformation for
 *    clarity, or work with directional cosine matrices.
 *
 */

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/ReferenceFrames/referenceFrameTransformations.h"

namespace tudat
{
namespace reference_frames
{
//! Get rotating planetocentric (R) to inertial (I) reference frame transformation matrix.
Eigen::Matrix3d
getRotatingPlanetocentricToInertialFrameTransformationMatrix( double angleFromXItoXR )
{
    // Declare local variables.
    // Declare local matrix.
    Eigen::Matrix3d localMatrix_;

    // Set local matrix.
    localMatrix_ = reference_frames::
            getInertialToPlanetocentricFrameTransformationMatrix( angleFromXItoXR );

    // Return transformation matrix.
    return localMatrix_.transpose( );
}

//! Get rotating planetocentric (R) to inertial (I) reference frame transformation quaternion.
Eigen::Quaterniond
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
Eigen::Matrix3d
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
Eigen::Quaterniond
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
Eigen::Quaterniond
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
Eigen::Matrix3d
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
Eigen::Quaterniond
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
Eigen::Quaterniond
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

} // namespace reference_frames
} // namespace tudat
