/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_AERODYNAMIC_MOMENT_H
#define TUDAT_AERODYNAMIC_MOMENT_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicCoefficientInterface.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic moment in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic moment. It takes primitive types as arguments to
 * perform the calculations. Therefor, these quantities (dynamic pressure, reference area and
 * aerodynamic coefficients) have to computed before passing them to this function.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param referenceArea Reference area of the aerodynamic coefficients.
 * \param referenceLength Reference length of the aerodynamic coefficients. Note that this
 *          reference length is used for all three independent directions.
 * \param momentCoefficients Aerodynamic moment coefficients in right-handed reference frame.
 * \return Resultant aerodynamic moment, given in reference frame in which the
 *          aerodynamic coefficients were given, but with opposite sign. i.e., a positive drag
 *          coefficient will give a negative force in -x direction (in the aerodynamic frame).
 */
Eigen::Vector3d computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const double referenceLength,
                                          const Eigen::Vector3d& momentCoefficients );

//! Compute the aerodynamic moment in same reference frame as input coefficients.
/*!
 * This function calculates the aerodynamic moment. It takes the dynamic pressure and an
 * aerodynamic coefficient interface as input. The coefficient interface has to have been
 * updated with current vehicle conditions before being passed to this function. Aerodynamic
 * coefficients and reference area and length are then retrieved from it.
 * \param dynamicPressure Dynamic pressure at which the body undergoing the force flies.
 * \param coefficientInterface AerodynamicCoefficientInterface class from which reference area
 *          and length, and the moment coefficients are retrieved.
 * \return Resultant aerodynamic moment, given in reference frame in which the
 *          aerodynamic coefficients were given.
 */
Eigen::Vector3d computeAerodynamicMoment(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface );

} // namespace aerodynamics
} // namespace tudat

#endif // TUDAT_AERODYNAMIC_MOMENT_H
