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
 *      110617    D. Dirkx          First creation of code.
 *      120324    K. Kumar          Minor Doxygen comment corrections, added astrodynamics
 *                                  namespace layer; added missing Eigen include-statement.
 *
 *    References
 *
 */

#ifndef TUDAT_AERODYNAMIC_ACCELERATION_H
#define TUDAT_AERODYNAMIC_ACCELERATION_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicForce.h"

namespace tudat
{
namespace astrodynamics
{
namespace acceleration_models
{

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
/*!
 * This function computes the aerodynamic acceleration. It takes primitive types as arguments to
 * perform the calculations. Therefore, these quantities (dynamicPressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the acceleration flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param aerodynamicCoefficients. Aerodynamic coefficients in right-handed reference frame.
 * \param vehicleMass Mass of vehicle undergoing acceleration.
 * \return Resultant aerodynamic acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::VectorXd computeAerodynamicAcceleration( const double dynamicPressure,
                                                const double referenceArea,
                                                const Eigen::Vector3d& aerodynamicCoefficients,
                                                const double vehicleMass )
{
    return force_models::computeAerodynamicForce( dynamicPressure,referenceArea,
                                                  aerodynamicCoefficients ) / vehicleMass;
}

//! Compute the aerodynamic acceleration in same reference frame as input coefficients.
/*!
 * This function computes the aerodynamic acceleration. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the acceleration flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and coefficients are retrieved.
 * \param vehicleMass Mass of vehicle undergoing acceleration.
 * \return Resultant aerodynamic acceleration, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::MatrixXd computeAerodynamicAcceleration(
        const double dynamicPressure, AerodynamicCoefficientInterface& coefficientInterface,
        const double vehicleMass )
{
    return force_models::computeAerodynamicForce( dynamicPressure, coefficientInterface )
            / vehicleMass;
}

} // namespace acceleration_models
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_ACCELERATION_H
