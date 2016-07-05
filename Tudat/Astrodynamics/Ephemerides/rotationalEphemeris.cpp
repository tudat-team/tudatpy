/*    Copyright (c) 2010-2015, Delft University of Technology
   *   All rights reserved.
   *
   *   Redistribution and use in source and binary forms, with or without modification, are
   *   permitted provided that the following conditions are met:
   *     - Redistributions of source code must retain the above copyright notice, this list of
   *       conditions and the following disclaimer.
   *     - Redistributions in binary form must reproduce the above copyright notice, this list of
   *       conditions and the following disclaimer in the documentation and/or other materials
   *       provided with the distribution.
   *     - Neither the name of the Delft University of Technology nor the names of its contributors
   *       may be used to endorse or promote products derived from this software without specific
   *       prior written permission.
   *
   *   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
   *   OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
   *   MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
   *   COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
   *   EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
   *   GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
   *   AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
   *   NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
   *   OF THE POSSIBILITY OF SUCH DAMAGE.
   *
   *   Changelog
   *     YYMMDD    Author            Comment
   *     130219    D. Dirkx          Migrated from personal code.
   *     130227    R.C.A. Boon       Changed include guard, improved commenting.
   *
   *   References
   *
   *   Notes
   *
   */

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

namespace tudat
{
namespace ephemerides
{

//! Function to calculate the rotational velocity vector of frame B w.r.t frame A.
Eigen::Vector3d getRotationalVelocityVectorInBaseFrameFromMatrices(
        const Eigen::Matrix3d& rotationToTargetFrame,
        const Eigen::Matrix3d& rotationMatrixToGlobalFrameDerivative )
{
    Eigen::Matrix3d crossProductMatrix =
            rotationMatrixToGlobalFrameDerivative * rotationToTargetFrame;
    return ( Eigen::Vector3d( ) << crossProductMatrix( 2, 1 ),
             crossProductMatrix( 0, 2 ), crossProductMatrix( 1, 0 ) ).finished( );

}

//! Function to calculate the time derivative of rotation matrix from frame A to frame B.
Eigen::Matrix3d getDerivativeOfRotationMatrixToFrame(
        const Eigen::Matrix3d& rotationToTargetFrame,
        const Eigen::Vector3d& rotationalVelocityVectorOfTargetFrameInBaseFrame )
{
    return linear_algebra::getCrossProductMatrix(
                -rotationToTargetFrame * rotationalVelocityVectorOfTargetFrameInBaseFrame ) *
            rotationToTargetFrame;
}

} // namespace tudat
} // namespace ephemerides

