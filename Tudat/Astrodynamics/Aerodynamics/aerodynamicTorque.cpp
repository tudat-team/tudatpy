/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicTorque.h"

namespace tudat
{
namespace aerodynamics
{

//! Compute the aerodynamic moment in same reference frame as input coefficients.
Eigen::Vector3d computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const double referenceLength,
                                          const Eigen::Vector3d& momentCoefficients )
{
    return dynamicPressure * referenceArea * referenceLength * momentCoefficients;
}

//! Compute the aerodynamic moment in same reference frame as input coefficients.
Eigen::Vector3d computeAerodynamicMoment( const double dynamicPressure, const double referenceArea,
                                          const Eigen::Vector3d& referenceLengths,
                                          const Eigen::Vector3d& momentCoefficients )
{
    return dynamicPressure * referenceArea * referenceLengths.cwiseProduct( momentCoefficients );
}

//! Calculates the aerodynamic moment in same reference frame as input coefficients.
Eigen::Vector3d computeAerodynamicMoment(
        const double dynamicPressure,
        AerodynamicCoefficientInterfacePointer coefficientInterface )
{
    return computeAerodynamicMoment( dynamicPressure,
                                     coefficientInterface->getReferenceArea( ),
                                     coefficientInterface->getReferenceLength( ),
                                     coefficientInterface->getCurrentMomentCoefficients( ) );
}


} // namespace aerodynamics
} // namespace tudat
