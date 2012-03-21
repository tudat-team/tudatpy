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
 *      110829    L. van der Ham    First creation of code.
 *      120221    L. van der Ham    Update of code and namespace.
 *      120227    K. Kumar          Updated code to new Tudat N Commandments; added enum to
 *                                  replace "magic numbers"; renamed file.
 *      120306    K. Kumar          Added Doxygen equation.
 *      120307    K. Kumar          Moved file.
 *      120321    K. Kumar          Moved enum for state element indices to
 *                                  stateDerivativeCircularRestrictedThreeBodyProblem.h.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#ifndef TUDAT_JACOBI_ENERGY_H
#define TUDAT_JACOBI_ENERGY_H

#include <Eigen/Core>

#include <Tudat/Astrodynamics/Gravitation/stateDerivativeCircularRestrictedThreeBodyProblem.h>

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
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

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_JACOBI_ENERGY_H
