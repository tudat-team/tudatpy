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
 *        Wakker, K.F.,"Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 *    Notes
 *      Position dimension-scale is distance between the two massive bodies in the CRTBP.
 *      Time dimension-scale is based on orbital period of 2*pi.
 *
 */

#ifndef TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
#define TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
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
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body.   [m^3 s^-2]
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body.
 *                                                                                       [m^3 s^-2]
 * \param distanceBetweenPrimaries Distance between primaries.                                  [m]
 * \return Dimensional Cartesian state.
 */
Eigen::VectorXd convertDimensionlessCartesianStateToDimensionalUnits(
        const Eigen::Vector6d& dimensionlessCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

//! Convert dimensionless time to dimensional units.
/*!
 * Convert dimensionless time to dimensional units.
 * \param timeInDimensionlessUnits Time in normalized units.
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body.   [m^3 s^-2]
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body.
 *                                                                                       [m^3 s^-2]
 * \param distanceBetweenPrimaries Distance between primaries.                                  [m]
 * \return Dimensional time.                                                                    [s]
 */
double convertDimensionlessTimeToDimensionalTime(
        const double timeInDimensionlessUnits, const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat

#endif // TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
