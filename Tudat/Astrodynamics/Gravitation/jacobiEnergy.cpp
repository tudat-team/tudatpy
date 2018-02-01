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

#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/stateVectorIndices.h"

namespace tudat
{
namespace gravitation
{

//! Compute Jacobi energy.
double computeJacobiEnergy( const double massParameter,
                            const Eigen::VectorXd& stateInNormalizedUnits )
{
    using namespace orbital_element_conversions;

    // Compute squared distances for efficient computation.
    const double xCoordinateToPrimaryBodySquared =
            ( massParameter + stateInNormalizedUnits( xCartesianPositionIndex ) )
            * ( massParameter + stateInNormalizedUnits( xCartesianPositionIndex ) );

    const double yCoordinateSquared = stateInNormalizedUnits( yCartesianPositionIndex )
            * stateInNormalizedUnits( yCartesianPositionIndex );

    const double zCoordinateSquared = stateInNormalizedUnits( zCartesianPositionIndex )
            * stateInNormalizedUnits( zCartesianPositionIndex );

    const double xCoordinateToSecondaryBodySquared =
            ( 1.0 - massParameter - stateInNormalizedUnits( xCartesianPositionIndex ) )
              * ( 1.0 - massParameter - stateInNormalizedUnits( xCartesianPositionIndex ) );

    const double distanceToPrimaryBody = std::sqrt(
                xCoordinateToPrimaryBodySquared + yCoordinateSquared + zCoordinateSquared );

    const double distanceToSecondaryBody = std::sqrt(
                xCoordinateToSecondaryBodySquared + yCoordinateSquared + zCoordinateSquared );

    return stateInNormalizedUnits( xCartesianPositionIndex ) * stateInNormalizedUnits( xCartesianPositionIndex )
            + yCoordinateSquared
            + 2.0 * ( 1.0 - massParameter ) / distanceToPrimaryBody
            + 2.0 * massParameter / distanceToSecondaryBody
            - stateInNormalizedUnits.segment( xCartesianVelocityIndex, 3 ).squaredNorm( );
}

} // namespace gravitation
} // namespace tudat
