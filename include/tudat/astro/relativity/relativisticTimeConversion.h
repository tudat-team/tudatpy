/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RELATIVISTIVTIMECONVERSION_H
#define TUDAT_RELATIVISTIVTIMECONVERSION_H

#include <Eigen/Core>
#include <vector>

#include "tudat/basics/basicTypedefs.h"
namespace tudat
{

namespace relativity
{

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential
/*!
 * Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, from a speed and scalar potential, allowing
 * for the possibility of EP violation (nominally none)
 * \param relativeSpeed Speed of object for which proper time is to be computed, w.r.t. frame of which coordinate time is used
 * \param gravitationalScalarPotential Scalar potential at location of object for which proper time is to be computed
 * \param equivalencePrincipleLpiViolationParameter Violation parameter of equivalence principle (default to GR value of 0.0
 * \return Proper time rate minus one (d tau/dt - 1.0)
 */
double calculateFirstCentralBodyProperTimeRateDifference(
        const double computationPointSpeed,
        const double gravitationalScalarPotential,
        const double equivalencePrincipleLpiViolationParameter = 0.0 );

//! Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, for a 1/c^2 potential from a static mass monopole
/*!
 * Function to compute proper-time rate w.r.t. coordinate time, minus 1.0, for a 1/c^2 potential from a static mass monopole
 * located at the origin of the reference frame of which the coordinate time is used.,
 * allowing for the possibility of EP violation (nominally none)
 * \param computationPointState Cartesian state of point where proper time rate is to be computed
 * \param centralBodyGravitationalParameter Newtonian gravitational potnetial of central body.
 * \param equivalencePrincipleLpiViolationParameter Violation parameter of equivalence principle (default to GR value of 0.0
 * \return Proper time rate minus one (d tau/dt - 1.0)
 */
double calculateFirstCentralBodyProperTimeRateDifference(
        const Eigen::Vector6d& computationPointState,
        const std::vector< Eigen::Vector6d >& perturbedStates,
        const std::vector< double >& centralBodyGravitationalParameters,
        const double equivalencePrincipleLpiViolationParameter = 0.0 );

}

}

#endif // TUDAT_RELATIVISTIVTIMECONVERSION_H
