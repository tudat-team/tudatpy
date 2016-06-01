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
 *      120716    A. Ronse          File created.
 *      120718    D. Dirkx          Splitted up functions into Coriolis, centripetal and Euler.
 *
 *    References
 *      J.S. Torok, Analytical Mechanics, Wiley-Interscience, 2000.
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/ReferenceFrames/apparentAccelerationModel.h"

namespace tudat
{
namespace reference_frames
{

//! Compute apparent acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeApparentAcceleration(
        const Eigen::Vector3d& accelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& angularAccelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame,
        const Eigen::Vector3d& velocityOfBodyInNonInertialReferenceFrame )
{
    return -accelerationOfNonInertialReferenceFrame +
            computeCentripetalAcceleration( angularVelocityOfNonInertialReferenceFrame,
                                            positionOfBodyInNonInertialReferenceFrame ) +
            computeCoriolisAcceleration( angularVelocityOfNonInertialReferenceFrame,
                                         velocityOfBodyInNonInertialReferenceFrame ) +
            computeEulerAcceleration( angularAccelerationOfNonInertialReferenceFrame,
                                      positionOfBodyInNonInertialReferenceFrame );
}

//! Compute centripetal acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeCentripetalAcceleration(
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame )
{
    return -angularVelocityOfNonInertialReferenceFrame.cross(
                ( angularVelocityOfNonInertialReferenceFrame.cross(
                      positionOfBodyInNonInertialReferenceFrame ) ) );
}

//! Compute Coriolis acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeCoriolisAcceleration(
        const Eigen::Vector3d& angularVelocityOfNonInertialReferenceFrame,
        const Eigen::Vector3d& velocityOfBodyInNonInertialReferenceFrame )
{
    return -2.0 * ( angularVelocityOfNonInertialReferenceFrame.cross(
                       velocityOfBodyInNonInertialReferenceFrame ) );
}

//! Compute Euler acceleration due to non-inertiality of reference frame.
Eigen::Vector3d computeEulerAcceleration(
        const Eigen::Vector3d& angularAccelerationOfNonInertialReferenceFrame,
        const Eigen::Vector3d& positionOfBodyInNonInertialReferenceFrame )
{
    return -angularAccelerationOfNonInertialReferenceFrame.cross(
                positionOfBodyInNonInertialReferenceFrame );
}

//! Update member variables used by apparent acceleration model.
void ApparentAccelerationModel::updateMembers( const double currentTime )
{
    if( !( this->currentTime_ == currentTime ) )
    {
        currentAccelerationOfNonInertialReferenceFrame_ =
                accelerationOfNonInertialReferenceFrameFunction_( );
        currentAngularVelocityOfNonInertialReferenceFrame_ =
                angularVelocityOfNonInertialReferenceFrameFunction_( );
        currentAngularAccelerationOfNonInertialReferenceFrame_ =
                angularAccelerationOfNonInertialReferenceFrameFunction_( );
        currentPositionOfBodyInNonInertialReferenceFrame_ =
                positionOfBodyInNonInertialReferenceFrameFunction_( );
        currentVelocityOfBodyInNonInertialReferenceFrame_ =
                velocityOfBodyInNonInertialReferenceFrameFunction_( );
    }
}

} // namespace reference_frames
} // namespace tudat
