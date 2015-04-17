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
 *      110829    L. van der Ham    File created.
 *      120221    L. van der Ham    Update of code and namespace.
 *      120227    K. Kumar          Updated code to new Tudat N Commandments; added enum to
 *                                  replace "magic numbers"; renamed file.
 *      120306    K. Kumar          Added Doxygen equation.
 *      120307    K. Kumar          Moved file.
 *      120321    K. Kumar          Moved enum for state element indices to
 *                                  stateDerivativeCircularRestrictedThreeBodyProblem.h.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Notes
 *
 */

#ifndef TUDAT_JACOBI_ENERGY_H
#define TUDAT_JACOBI_ENERGY_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Gravitation/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute Jacobi energy.
/*!
 * Returns the value of the Jacobi energy in normalized units of the CRTBP, given in the following
 * equation (Wakker, 2007):
 * \f[
 *      C_{J} = x^{2} + y^{2} + 2 * \left(\frac{1-\mu}{r_{1}} + \frac{\mu}{r_{2}}\right) - V^{2},
 * \f]
 * where \f$C_{J}\f$ is the Jacobi energy, \f$x\f$ and \f$y\f$ are Cartesian position coordinates
 * of the third body, \f$mu\f$ is the gravitational parameter of the secondary, \f$r_{1}\f$ and
 * \f$r_{2}\f$ are the distances from the third body to the primaries, and \f$V\f$ is the magnitude
 * of the velocity of the third body. This equation assumes that the system has been normalized
 * such that the mean motion of the primiaries, and their mutual distance is equal to 1.
 * \param massParameter Mass parameter of three-body system.
 * \param stateInNormalizedUnits State of third body in normalized units.
 * \return Jacobi energy in normalized units.
 */
double computeJacobiEnergy( const double massParameter,
                            const Eigen::VectorXd& stateInNormalizedUnits );

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat

#endif // TUDAT_JACOBI_ENERGY_H
