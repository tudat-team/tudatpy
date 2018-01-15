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

#include <cmath>

#include "Tudat/Basics/utilityMacros.h"

#include "Tudat/Astrodynamics/Propagators/stateDerivativeCircularRestrictedThreeBodyProblem.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{
namespace propagators
{

//! Compute state derivative.
Eigen::Vector6d StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative(
        const double time, const Eigen::Vector6d& cartesianState )
{
    using namespace orbital_element_conversions;

    TUDAT_UNUSED_PARAMETER( time );

    // Compute distance to primary body.
    const double xCoordinateToPrimaryBodySquared =
            ( cartesianState( xCartesianPositionIndex ) + massParameter )
            * ( cartesianState( xCartesianPositionIndex ) + massParameter );

    const double yCoordinateSquared = cartesianState( yCartesianPositionIndex )
            * cartesianState( yCartesianPositionIndex );

    const double zCoordinateSquared = cartesianState( zCartesianPositionIndex )
            * cartesianState( zCartesianPositionIndex );

    const double normDistanceToPrimaryBodyCubed = pow(
                xCoordinateToPrimaryBodySquared + yCoordinateSquared + zCoordinateSquared, 1.5 );

    // Compute distance to secondary body.
    const double xCoordinateSecondaryBodySquared =
            ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) )
            * ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) );

    double normDistanceToSecondaryBodyCubed = pow(
                xCoordinateSecondaryBodySquared + yCoordinateSquared + zCoordinateSquared, 1.5 );

    // Compute derivative of state.
    Eigen::VectorXd stateDerivative( 6 );

    stateDerivative.segment( xCartesianPositionIndex, 3 ) = cartesianState.segment( xCartesianVelocityIndex, 3 );

    stateDerivative( 3 ) = cartesianState( xCartesianPositionIndex )
            - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
            * ( cartesianState( xCartesianPositionIndex ) + massParameter )
            - ( massParameter / normDistanceToSecondaryBodyCubed )
            * ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) )
            + 2.0 * cartesianState( yCartesianVelocityIndex );
    stateDerivative( 4 ) = cartesianState( yCartesianPositionIndex )
            * ( 1.0 - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                - ( massParameter / normDistanceToSecondaryBodyCubed ) )
            - 2.0 * cartesianState( xCartesianVelocityIndex );
    stateDerivative( 5 ) = -cartesianState( zCartesianPositionIndex )
            * ( ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                + ( massParameter / normDistanceToSecondaryBodyCubed ) );

    // Return computed state derivative.
    return stateDerivative;
}

} // namespace propagators

} // namespace tudat
