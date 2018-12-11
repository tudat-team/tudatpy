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


//! Convert dimensional Cartesian state to dimensionless state.
/*!
 * Convert dimensional Cartesian state to dimensionless state.
 * \param dimensionalCartesianState Dimensional Cartesian state vector.                       [m]
 *          Order is important!
 *          dimensionalCartesianState( 0 ) = x-position coordinate,                           [m]
 *          dimensionalCartesianState( 1 ) = y-position coordinate,                           [m]
 *          dimensionalCartesianState( 2 ) = z-position coordinate,                           [m]
 *          dimensionalCartesianState( 3 ) = x-velocity coordinate,                         [m/s]
 *          dimensionalCartesianState( 4 ) = y-velocity coordinate,                         [m/s]
 *          dimensionalCartesianState( 5 ) = z-velocity coordinate.                         [m/s]
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body.   [m^3 s^-2]
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body.
 *                                                                                       [m^3 s^-2]
 * \param distanceBetweenPrimaries Distance between primaries.
 * \return Dimensionless Cartesian state.
 */
Eigen::VectorXd convertDimensionalCartesianStateToDimensionlessState(
        const Eigen::Vector6d& dimensionalCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries );

//! Convert dimensional time to dimensionless time.
/*!
 * Convert dimensional time to dimensionless time.
 * \param timeDimensional Time in dimensional units                                             [s]
 * \param gravitationalParameterOfPrimaryBody Gravitational parameter of primary body.   [m^3 s^-2]
 * \param gravitationalParameterOfSecondaryBody Gravitational parameter of secondary body.
 *                                                                                       [m^3 s^-2]
 * \param distanceBetweenPrimaries Distance between primaries.                                  [m]
 * \return Adimensional time.
 */
double convertDimensionalTimeToDimensionlessTime(
        const double dimensionalTime,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries);


//!Function to transform normalized co-rotating coordinates into cartesian ones
Eigen::Vector6d convertCorotatingNormalizedToCartesianCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& normalizedState,
        const double normalizedTime );

//! Function to transform cartesian coordinates into co-rotating normalized ones
Eigen::Vector6d convertCartesianToCorotatingNormalizedCoordinates(
        const double gravitationalParameterPrimary,
        const double gravitationalParameterSecondary,
        const double distancePrimarySecondary,
        const Eigen::Vector6d& cartesianState,
        const double time );

} // namespace circular_restricted_three_body_problem
} // namespace tudat

#endif // TUDAT_UNIT_CONVERSIONS_CIRCULAR_RESTRICTED_THREE_BODY_PROBLEM_H
