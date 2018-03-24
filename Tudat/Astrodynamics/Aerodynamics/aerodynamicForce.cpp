/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicForce.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic force in same reference frame as input coefficients.
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
