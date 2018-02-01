/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#ifndef TUDAT_JACOBI_ENERGY_H
#define TUDAT_JACOBI_ENERGY_H

#include <Eigen/Core>

namespace tudat
{
namespace gravitation
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

} // namespace gravitation
} // namespace tudat

#endif // TUDAT_JACOBI_ENERGY_H
