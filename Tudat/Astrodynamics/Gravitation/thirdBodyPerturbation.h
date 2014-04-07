/*    Copyright (c) 2010-2014, Delft University of Technology
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
 *      130307    R.C.A. Boon       File created.
 *      130308    D. Dirkx          Modified to add variable central body position.
 *
 *    References
 *      Wakker, K.F. Astrodynamics I, Delft University of Technology, 2010.
 *      Montebruck O, Gill E. Satellite Orbits, Corrected Third Printing, Springer, 2005.
 *
 *    Notes
 *
 */

#ifndef TUDAT_THIRD_BODY_PERTURBATION_H
#define TUDAT_THIRD_BODY_PERTURBATION_H

#include <Eigen/Core>

namespace tudat
{
namespace gravitation
{

//! Compute perturbing acceleration by third body.
/*!
 * Computes the perturbing acceleration on a point mass in orbit about a central body (point mass),
 * caused by a third body (point mass). This acceleration is expressed with respect to a
 * pseudo-inertial reference frame centered at the central body body, i.e, a non-rotating frame
 * aligned with an inertial reference frame. This acceleration may be summed with other third body
 * accelerations in order to get the total gravitational perturbation acceleration of N bodies
 * (Wakker, 2010; Montenbruck & Gill, 2005).
 * \param gravitationalParameterOfPerturbingBody The gravitational parameter of the perturbing body,
 *      i.e., the third body.                                                             [m^3/s^2]
 * \param positionOfPerturbingBody The position of the third body, the body that
 *      causes the perturbations, in Cartesian coordinates.
 *          positionOfPerturbingBody( 0 ) = x-position                                          [m]
 *          positionOfPerturbingBody( 1 ) = y-position                                          [m]
 *          positionOfPerturbingBody( 2 ) = z-position                                          [m]
 * \param positionOfAffectedBody The position of the second body, the body that
 *      experiences the perturbation, in Cartesian coordinates.
 *          positionOfAffectedBody( 0 ) = x-position                                            [m]
 *          positionOfAffectedBody( 1 ) = y-position                                            [m]
 *          positionOfAffectedBody( 2 ) = z-position                                            [m]
 * \param positionOfCentralBody The position of the central body, the body w.r.t. which the
 *      acceleration is calculated the perturbation, in Cartesian coordinates (default=origin).
 *          positionOfAffectedBody( 0 ) = x-position                                            [m]
 *          positionOfAffectedBody( 1 ) = y-position                                            [m]
 *          positionOfAffectedBody( 2 ) = z-position                                            [m]
 * \return perturbingAcceleration The perturbing acceleration in Cartesian components, that gives
 *      the difference of the total acceleration with the central two-body acceleration.
 *          perturbingAcceleration( 0 ) = x-acceleration                                    [m/s^2]
 *          perturbingAcceleration( 1 ) = y-acceleration                                    [m/s^2]
 *          perturbingAcceleration( 2 ) = z-acceleration                                    [m/s^2]
 */
Eigen::Vector3d computeThirdBodyPerturbingAcceleration(
        const double gravitationalParameterOfPerturbingBody,
        const Eigen::Vector3d& positionOfPerturbingBody,
        const Eigen::Vector3d& positionOfAffectedBody,
        const Eigen::Vector3d& positionOfCentralBody = Eigen::Vector3d::Zero( ) );

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_THIRD_BODY_PERTURBATION_H
