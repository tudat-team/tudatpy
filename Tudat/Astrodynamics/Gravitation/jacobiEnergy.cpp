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
 *      120306    K. Kumar          Removed erroneous Doxygen comments.
 *      120307    K. Kumar          Moved file.
 *      120426    K. Kumar          Updated code to compute Jacobi energy more efficiently.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute Jacobi energy.
double computeJacobiEnergy( const double massParameter,
                            const Eigen::VectorXd& stateInNormalizedUnits )
{
    // Compute squared distances for efficient computation.
    const double xCoordinateToPrimaryBodySquared =
            ( massParameter + stateInNormalizedUnits( xPositionIndex ) )
            * ( massParameter + stateInNormalizedUnits( xPositionIndex ) );

    const double yCoordinateSquared = stateInNormalizedUnits( yPositionIndex )
            * stateInNormalizedUnits( yPositionIndex );

    const double zCoordinateSquared = stateInNormalizedUnits( zPositionIndex )
            * stateInNormalizedUnits( zPositionIndex );

    const double xCoordinateToSecondaryBodySquared =
            ( 1.0 - massParameter - stateInNormalizedUnits( xPositionIndex ) )
              * ( 1.0 - massParameter - stateInNormalizedUnits( xPositionIndex ) );

    const double distanceToPrimaryBody = std::sqrt(
                xCoordinateToPrimaryBodySquared + yCoordinateSquared + zCoordinateSquared );

    const double distanceToSecondaryBody = std::sqrt(
                xCoordinateToSecondaryBodySquared + yCoordinateSquared + zCoordinateSquared );

    return stateInNormalizedUnits( xPositionIndex ) * stateInNormalizedUnits( xPositionIndex )
            + yCoordinateSquared
            + 2.0 * ( 1.0 - massParameter ) / distanceToPrimaryBody
            + 2.0 * massParameter / distanceToSecondaryBody
            - stateInNormalizedUnits.segment( xVelocityIndex, 3 ).squaredNorm( );
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat
