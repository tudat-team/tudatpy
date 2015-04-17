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
 *      110829    L. van der Ham    File created.
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
 *    Notes
 *
 */

#include <cmath>

#include "Tudat/Astrodynamics/Gravitation/jacobiEnergy.h"

namespace tudat
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

} // namespace circular_restricted_three_body_problem
} // namespace gravitation
} // namespace tudat
