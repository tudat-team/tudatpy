/*    Copyright (c) 2010-2012 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120322    D. Dirkx          First creation of code.
 *
 *    References
 *
 */

#ifndef TUDAT_AERODYNAMIC_ROTATIONAL_ACCELERATION_H
#define TUDAT_AERODYNAMIC_ROTATIONAL_ACCELERATION_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicMoment.h"

namespace tudat
{
namespace astrodynamics
{
namespace rotational_acceleration_models
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
Eigen::MatrixXd computeAerodynamicRotationalAcceleration(
        const double dynamicPressure, const double referenceArea, const double referenceLength,
        const Eigen::Vector3d& momentCoefficients, const double vehicleMass )
{
    return moment_models::computeAerodynamicMoment( dynamicPressure, referenceArea,
                                                    referenceLength, momentCoefficients )
            / vehicleMass;
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
Eigen::MatrixXd computeAerodynamicRotationalAcceleration(
        const double dynamicPressure, AerodynamicCoefficientInterface& coefficientInterface,
        const double vehicleMass )
{
    return moment_models::computeAerodynamicMoment( dynamicPressure, coefficientInterface )
            / vehicleMass;
}

} // namespace rotational_acceleration_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ROTATIONAL_ACCELERATION_H
