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
 *      110905    L. van der Ham    First creation of code.
 *      110922    L. van der Ham    Added computation of Jacobian energy constant.
 *      111107    L. van der Ham    Removed option planar.
 *      111110    K. Kumar          Minor comment and naming modifications; error removed from
 *                                  z-direction equation.
 *      120221    L. van der Ham    Removed computation Jacobian energy.
 *      120309    K. Kumar          Updated code to latest Tudat standards; updated
 *                                  computeStateDerivative() function.
 *      120426    K. Kumar          Updated code to compute state derivative more efficiently.
 *
 *    References
 *        Wakker, K.F., "Astrodynamics I, AE4-874", Delft University of Technology, 2007.
 *
 */

#include <cmath>

#include <TudatCore/Basics/utilityMacros.h>

#include "Tudat/Astrodynamics/Gravitation/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace astrodynamics
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute state derivative.
Eigen::VectorXd StateDerivativeCircularRestrictedThreeBodyProblem::
computeStateDerivative( const double time, const Eigen::VectorXd& cartesianState )
{
    TUDAT_UNUSED_PARAMETER( time );

    // Compute distance to primary body.
    const double xCoordinateToPrimaryBodySquared =
            ( cartesianState( xPositionIndex ) + massParameter )
            * ( cartesianState( xPositionIndex ) + massParameter );

    const double yCoordinateSquared = cartesianState( yPositionIndex )
            * cartesianState( yPositionIndex );

    const double zCoordinateSquared = cartesianState( zPositionIndex )
            * cartesianState( zPositionIndex );

    const double normDistanceToPrimaryBodyCubed = pow(
                xCoordinateToPrimaryBodySquared + yCoordinateSquared + zCoordinateSquared, 1.5 );

    // Compute distance to secondary body.
    const double xCoordinateSecondaryBodySquared =
            ( cartesianState( xPositionIndex ) - ( 1.0 - massParameter ) )
            * ( cartesianState( xPositionIndex ) - ( 1.0 - massParameter ) );

    double normDistanceToSecondaryBodyCubed = pow(
                xCoordinateSecondaryBodySquared + yCoordinateSquared + zCoordinateSquared, 1.5 );

    // Compute derivative of state.
    Eigen::VectorXd stateDerivative( 6 );

    stateDerivative.segment( xPositionIndex, 3 ) = cartesianState.segment( xVelocityIndex, 3 );

    stateDerivative( xAccelerationIndex ) = cartesianState( xPositionIndex )
            - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
            * ( cartesianState( xPositionIndex ) + massParameter )
            - ( massParameter / normDistanceToSecondaryBodyCubed )
            * ( cartesianState( xPositionIndex ) - ( 1.0 - massParameter ) )
            + 2.0 * cartesianState( yVelocityIndex );
    stateDerivative( yAccelerationIndex ) = cartesianState( yPositionIndex )
            * ( 1.0 - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                - ( massParameter / normDistanceToSecondaryBodyCubed ) )
            - 2.0 * cartesianState( xVelocityIndex );
    stateDerivative( zAccelerationIndex ) = -cartesianState( zPositionIndex )
            * ( ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                + ( massParameter / normDistanceToSecondaryBodyCubed ) );

    // Return computed state derivative.
    return stateDerivative;
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace astrodynamics
} // namespace tudat
