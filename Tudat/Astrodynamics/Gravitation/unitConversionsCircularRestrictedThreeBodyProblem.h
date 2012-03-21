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
 *      110717    L. van der Ham    First creation of code.
 *      111108    K. Kumar          Simplified conversion functions; changed "normalized" to
 *                                  "dimensionless".
 *      120307    K. Kumar          Updated Doxygen comments; made functions const-correct;
 *                                  updated to satisfy Tudat protocols; moved file.
 *
 *    References
 *        Wakker, K.F.,"Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Position dimension-scale is distance between the two massive bodies in the CRTBP.
 *    Time dimension-scale is based on orbital period of 2*pi.
 *
 */

#ifndef TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <Eigen/Core>

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Convert dimensionless Cartesian state to dimensional units.
/*!
 * Convert dimensionless Cartesian state to dimensional units.
 * \param dimensionlessCartesianState Dimensionless Cartesian state vector.
 *          Order is important!
 *          dimensionlessCartesianState( 0 ) = x-position coordinate,                           [m]
 *          dimensionlessCartesianState( 1 ) = y-position coordinate,                           [m]
 *          dimensionlessCartesianState( 2 ) = z-position coordinate,                           [m]
 *          dimensionlessCartesianState( 3 ) = x-velocity coordinate,                         [m/s]
 *          dimensionlessCartesianState( 4 ) = y-velocity coordinate,                         [m/s]
 *          dimensionlessCartesianState( 5 ) = z-velocity coordinate.                         [m/s]
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body
 *          [kg m^3 s^-2].
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body
 *          [kg m^3 s^-2].
 * \param distanceBetweenPrimaries Distance between primaries.
 * \return Dimensional Cartesian state.
 */
Eigen::VectorXd convertDimensionlessCartesianStateToDimensionalUnits(
        const Eigen::VectorXd& dimensionlessCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

//! Convert dimensionless time to dimensional units.
/*!
 * Convert dimensionless time to dimensional units.
 * \param timeInNormalizedUnits Time in normalized units.
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body [kg m^3 s^-2].
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body
 *           [kg m^3 s^-2].
 * \param distanceBetweenPrimaries Distance between primaries.
 * \return Dimensional time [s].
 */
double convertDimensionlessTimeToDimensionalTime(
        const double timeInDimensionlessUnits, const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat

#endif // TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
