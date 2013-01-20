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
 *      110617    F.M. Engelen      Creation of code.
 *      110822    D. Dirkx          Removed pointer to double member, minor changes.
 *      120324    K. Kumar          Added missing include-statements; used TUDAT_NAN in class
 *                                  constructor initialization list instead of -0.0; minor
 *                                  corrections in Doxygen comments; updated input parameters for
 *                                  free functions; added astrodynamics namespace layer.
 *      121020    D. Dirkx          Update to new acceleration model architecture.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_AERODYNAMIC_FORCE_H
#define TUDAT_AERODYNAMIC_FORCE_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include <TudatCore/Basics/utilityMacros.h>
#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic force in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic force. It takes primitive types as arguments to
 * perform the calculations. Therefore, these quantities (dynamicPressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param aerodynamicCoefficients. Aerodynamic coefficients in right-handed reference frame.
 * \return Resultant aerodynamic force, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::Vector3d computeAerodynamicForce( const double dynamicPressure,
                                         const double referenceArea,
                                         const Eigen::Vector3d& aerodynamicCoefficients )
{
    return dynamicPressure * referenceArea * aerodynamicCoefficients;
}

//! Compute the aerodynamic force in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic force. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and coefficients are retrieved.
 * \return Resultant aerodynamic force, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::Vector3d computeAerodynamicForce(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface )
{
    return computeAerodynamicForce( dynamicPressure,
                                    coefficientInterface->getReferenceArea( ),
                                    coefficientInterface->getCurrentForceCoefficients( ) );
}

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_FORCE_H
