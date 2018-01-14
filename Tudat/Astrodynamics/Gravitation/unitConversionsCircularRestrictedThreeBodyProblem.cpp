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

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/unitConversionsCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Convert dimensionless Cartesian state to dimensional units.
Eigen::VectorXd convertDimensionlessCartesianStateToDimensionalUnits(
        const Eigen::Vector6d &dimensionlessCartesianState,
        const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    // Declare dimensional Cartesian state.
    Eigen::Vector6d dimensionalCartesianState;

    // Convert position to dimensional units.
    dimensionalCartesianState.segment( 0, 3 )
            = dimensionlessCartesianState.segment( 0, 3 ) * distanceBetweenPrimaries;

    // Convert velocity to dimensional units.
    dimensionalCartesianState.segment( 3, 3 )
            = dimensionlessCartesianState.segment( 3, 3 )
            * std::sqrt( ( gravitationalParameterOfPrimaryBody
                           + gravitationalParameterOfSecondaryBody ) / distanceBetweenPrimaries );

    // Return state in dimensional units.
    return dimensionalCartesianState;
}

//! Convert dimensionless time to dimensional units.
double convertDimensionlessTimeToDimensionalTime(
        const double timeInDimensionlessUnits, const double gravitationalParameterOfPrimaryBody,
        const double gravitationalParameterOfSecondaryBody,
        const double distanceBetweenPrimaries )
{
    return timeInDimensionlessUnits * std::sqrt( std::pow( distanceBetweenPrimaries, 3.0 )
                                                 / ( gravitationalParameterOfPrimaryBody
                                                     + gravitationalParameterOfSecondaryBody ) );
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat
