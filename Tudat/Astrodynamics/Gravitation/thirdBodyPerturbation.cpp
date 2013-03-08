/*    Copyright (c) 2010-2013, Delft University of Technology
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

#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"

namespace tudat
{
namespace gravitation
{

//! Compute perturbing acceleration by third body
Eigen::Vector3d computeThirdBodyPerturbingAcceleration(
        const double gravitationalParameterOfPerturbingBody,
        const Eigen::Vector3d& positionOfPerturbingBody,
        const Eigen::Vector3d& positionOfAffectedBody,
        const Eigen::Vector3d& positionOfCentralBody )
// Using chapter 4 of (Wakker, 2010).
{
    // Return acceleration.
    return computeGravitationalAcceleration( positionOfAffectedBody,
                                             gravitationalParameterOfPerturbingBody,
                                             positionOfPerturbingBody ) -
            computeGravitationalAcceleration( positionOfCentralBody,
                                              gravitationalParameterOfPerturbingBody,
                                              positionOfPerturbingBody );
}

} // namespace gravitation
} // namespace tudat
