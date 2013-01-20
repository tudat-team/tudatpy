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
 *      120322    D. Dirkx          Creation of code.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_AERODYNAMIC_ROTATIONAL_ACCELERATION_H
#define TUDAT_AERODYNAMIC_ROTATIONAL_ACCELERATION_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicMoment.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic rotational acceleration in same reference frame as coefficients.
/*!
 * This function calculates the aerodynamic rotational acceleration. It takes primitive types
 * as arguments to perform the calculations. Therefor, these quantities (dynamic pressure,
 * reference area and aerodynamic coefficients) have to computed before passing them to this
 * function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param referenceLength Reference length of the aerodynamic coefficients. Note that this
 *          reference length is used for all three independent directions.
 * \param aerodynamicCoefficients. Aerodynamic rotational acceleration coefficients in
 *          right-handed reference frame.
 * \return Resultant aerodynamic rotational acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given, but with opposite sign. i.e., a positive drag
 *          coefficient will give a negative force in -x direction (in the aerodynamic frame).
 */
Eigen::Vector3d computeAerodynamicRotationalAcceleration(
        const double dynamicPressure, const double referenceArea, const double referenceLength,
        const Eigen::Vector3d& momentCoefficients, const double vehicleMass )
{
    return computeAerodynamicMoment( dynamicPressure, referenceArea,
                                     referenceLength, momentCoefficients ) / vehicleMass;
}

//! Compute the aerodynamic rotational acceleration in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic rotational acceleration. It takes the dynamic pressure
 *  and an aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area and length are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and length, and the rotational acceleration coefficients are retrieved.
 * \return Resultant aerodynamic rotational acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::Vector3d computeAerodynamicRotationalAcceleration(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface,
        const double vehicleMass )
{
    return computeAerodynamicMoment( dynamicPressure, coefficientInterface ) / vehicleMass;
}

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ROTATIONAL_ACCELERATION_H
