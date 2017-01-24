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
 *      130219    D. Dirkx          Migrated from personal code.
 *      130227    R.C.A. Boon       Changed header indentation, minor textual changes, updated
 *                                  commenting.
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

namespace tudat
{
namespace ephemerides
{

//! Get rotation quaternion to target frame from base frame.
Eigen::Quaterniond SimpleRotationalEphemeris::getRotationToTargetFrame(
        const double secondsSinceEpoch )
{
    // Determine number of seconds since initial rotational state, as set by constructor.
    double inputSecondsSinceEpoch = secondsSinceEpoch;

    // Determine rotation angle compared to initial rotational state.
    double rotationAngle = basic_mathematics::computeModulo(
                ( inputSecondsSinceEpoch - initialSecondsSinceEpoch_ ) * rotationRate_,
                2.0 * mathematical_constants::PI );

    // Calculate and return rotation to base frame.
    return reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
                rotationAngle ) * initialRotationToTargetFrame_;
}

//! Function to calculate the derivative of the rotation matrix from target frame to original frame.
Eigen::Matrix3d SimpleRotationalEphemeris::getDerivativeOfRotationToTargetFrame(
        const double secondsSinceEpoch )
{
    // Determine number of seconds since initial rotational state, as set by constructor.
    double inputSecondsSinceEpoch = secondsSinceEpoch;

    // Determine rotation angle compared to initial rotational state.
    double rotationAngle = basic_mathematics::computeModulo(
                ( inputSecondsSinceEpoch - initialSecondsSinceEpoch_ ) * rotationRate_,
                2.0 * mathematical_constants::PI );

    // Calculate derivative of rotation matrix.
    return rotationRate_ * auxiliaryMatrix_ * tudat::reference_frames::
            getInertialToPlanetocentricFrameTransformationQuaternion( rotationAngle )
            * Eigen::Matrix3d( initialRotationToTargetFrame_ );
}

//! Function to reset the right ascension and declination of body's north pole.
void SimpleRotationalEphemeris::resetInitialPoleRightAscensionAndDeclination(
        const double rightAscension, const double declination )
{
    // Recalculate initial rotation quaternion
    initialRotationToTargetFrame_ =
        reference_frames::getInertialToPlanetocentricFrameTransformationQuaternion(
            declination, rightAscension, initialEulerAngles_.z( ) );

    // Reset angles in vector of Euler angles.
    initialEulerAngles_.x( ) = rightAscension;
    initialEulerAngles_.y( ) = declination;
}

} // namespace tudat
} // namespace ephemerides
