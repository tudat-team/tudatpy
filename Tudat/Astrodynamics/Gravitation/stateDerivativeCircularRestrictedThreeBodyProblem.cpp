/*    Copyright (c) 2010-2015, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110905    L. van der Ham    File created.
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
 *    Notes
 *
 */

#include <cmath>

#include "Tudat/Basics/utilityMacros.h"

#include "Tudat/Astrodynamics/Gravitation/stateDerivativeCircularRestrictedThreeBodyProblem.h"

namespace tudat
{
namespace gravitation
{
namespace circular_restricted_three_body_problem
{

//! Compute state derivative.
Eigen::Vector6d StateDerivativeCircularRestrictedThreeBodyProblem::computeStateDerivative(
        const double time, const Eigen::Vector6d& cartesianState )
{
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

    stateDerivative( xAccelerationIndex ) = cartesianState( xCartesianPositionIndex )
            - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
            * ( cartesianState( xCartesianPositionIndex ) + massParameter )
            - ( massParameter / normDistanceToSecondaryBodyCubed )
            * ( cartesianState( xCartesianPositionIndex ) - ( 1.0 - massParameter ) )
            + 2.0 * cartesianState( yCartesianVelocityIndex );
    stateDerivative( yAccelerationIndex ) = cartesianState( yCartesianPositionIndex )
            * ( 1.0 - ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                - ( massParameter / normDistanceToSecondaryBodyCubed ) )
            - 2.0 * cartesianState( xCartesianVelocityIndex );
    stateDerivative( zAccelerationIndex ) = -cartesianState( zCartesianPositionIndex )
            * ( ( ( 1.0 - massParameter ) / normDistanceToPrimaryBodyCubed )
                + ( massParameter / normDistanceToSecondaryBodyCubed ) );

    // Return computed state derivative.
    return stateDerivative;
}

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat
