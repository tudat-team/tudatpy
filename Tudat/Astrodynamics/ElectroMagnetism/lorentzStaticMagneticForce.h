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
 *      131109	  S. Hirsh			File Created. 
 *
 *    References
 *		
 *
 *    Notes
 *
 */

#ifndef TUDAT_LORENTZ_STATIC_MAGNETIC_FORCE_H
#define TUDAT_LORENTZ_STATIC_MAGNETIC_FORCE_H

#include <Eigen/Core>

namespace tudat
{
namespace electro_magnetism
{

//! Compute Lorentz force due to static magnetic field
/*!
 * Computes Lorentz force due to static magnetic field on a particle in an
 * inertial reference frame.
 * Assumes the particle is point-like.The following equation is used to calculate the force:
 * \f[
 *		\bar{f} = q * \bar{v} x \bar{B}
 * \f] 
 * where \f$f\f$ is the force on the particle, \f$q\f$ is the charge of the accelerating
 * particle, \f$\bar{v}\f$ is the velocity of the accelerating particle, and \f$\bar{B}\f$
 * is the local magnetic field. Note: This force does not take into account
 * the Coulomb force between the particle and the source.
 * \param velocityOfBodySubjectToAcceleration Velocity of body which is being accelerated
 *			by the Lorentz force.                                                             [m/s]
 * \param chargeOfBodySubjectToAcceleration Charge of body which is being accelerated by
 *			Lorentz force. 																	  [C]
 * \param localMagneticField local magnetic field at position of body subject to acceleration [T·m]
 * \return Lorentz force due to static magnetic field 				                          [N]
 */

Eigen::Vector3d computeLorentzForceDueToStaticMagneticField(
        const Eigen::Vector3d& velocityOfBodySubjectToAcceleration,
        const Eigen::Vector3d& localMagneticField,
        const double chargeOfBodySubjectToAcceleration );

} // namespace electro_magnetism
} // namespace tudat

#endif // TUDAT_LORENTZ_STATIC_MAGNETIC_FORCE_H
