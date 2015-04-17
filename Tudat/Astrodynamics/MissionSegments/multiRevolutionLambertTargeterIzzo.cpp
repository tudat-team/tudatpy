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
 *      121213    R.C.A. Boon       Creation of code.
 *      130117    R.C.A. Boon       Added solved boolean flag to applicable member functions (with
 *                                  help from S. Billemont). Moved constructor to header file.
 *                                  Moved debugged getMaximumNumberOfRevolutions to source file.
 *      130211    R.C.A. Boon       Added hasSolution flag to root finder.
 *      120227    S. Billemont      Removed hasSolution in favor of an exception.
 *      130325    R.C.A. Boon       Removed superfluous sanity check of number of revolutions in
 *                                  execute() function, fixed bug in computation of maximumNumberOf-
 *                                  Revolutions.
 *
 *    References
 *      PyKEP toolbox, Dario Izzo, ESA Advanced Concepts Team.
 *      Richard H. An Introduction to the Mathematics and Methods of Astrodynamics, Revised
 *          Edition.
 *      Battin, AIAA Education Series.
 *
 *    Notes
 *
 */

#include <cmath>

#include <boost/format.hpp>
#include <boost/math/special_functions.hpp> // for asinh and acosh

#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/Astrodynamics/MissionSegments/multiRevolutionLambertTargeterIzzo.h"
#include "Tudat/Mathematics/BasicMathematics/convergenceException.h"

namespace tudat
{
namespace mission_segments
{

//! Compute solution for N revolutions and branch.
void MultiRevolutionLambertTargeterIzzo::computeForRevolutionsAndBranch(
        const int aNumberOfRevolutions, const bool aIsRightBranch )
{
    // Adjust parameters for new solution
    numberOfRevolutions = aNumberOfRevolutions;
    isRightBranch = aIsRightBranch;

    // Check whether number of revolutions is possible
    sanityCheckNumberOfRevolutions( );

    // Execute problem solving for new solution
    execute( );
}

//! Get maximum number of revolutions calculated.
int MultiRevolutionLambertTargeterIzzo::getMaximumNumberOfRevolutions( )
{
    if ( !solved )
    {
        transformDimensions( );
        sanityCheckNumberOfRevolutions( );
    }

    return maximumNumberOfRevolutions;
}

//! Sanity check number of revolutions.
void MultiRevolutionLambertTargeterIzzo::sanityCheckNumberOfRevolutions( )
{
    // If not yet defined, calculate number of revolutions possible.
    if ( maximumNumberOfRevolutions == NO_MAXIMUM_REVOLUTIONS )
    {
        // Temporarily store specified number, as numberOfRevolutions is needed to calculate max
        // (this is a tricky way to work, but on the other hand this makes this approach decidedly
        // different from PyKEP routines and it also happens only once per object).
        int copyOfOriginalNumberOfRevolutions = numberOfRevolutions;

        // Calculate first guess of maximum, by dividing the time of flight of the minimum energy
        // ellipse by the normalized time of flight.
        numberOfRevolutions = static_cast< int >(
                    normalizedTimeOfFlight / (
                        mathematical_constants::PI / 2.0
                        * std::sqrt( 2.0 * normalizedSemiPerimeter
                                     * normalizedSemiPerimeter
                                     * normalizedSemiPerimeter ) ) );

        // If the current guess for the maximum is non-zero, then additional analysis is required to
        // determine the correct maximum.
        if( numberOfRevolutions != 0)
        {
            // The following try-block is meant to check whether the solution converges or not. If
            // the current guess for the maximum number of revolutions is correct, then the problem
            // will converge. If it does not, an exception will be thrown stating that it did not
            // converge. Catching this exception allows to decrease the guess only when the
            // exception occurs, and not under other circumstances.
            try
            {
                // Compute root (no further information is required)
                computeRootTimeOfFlight();
            }
            catch( basic_mathematics::ConvergenceException )
            {
                // If the rootfinder did not converge, then the current guess is wrong and needs to
                // be decreased
                numberOfRevolutions--;
            }
        }
        // No further analysis is needed of the current guess is equal to zero.

        // Maximum is now found.
        maximumNumberOfRevolutions = numberOfRevolutions;

        // Reinstating original number of revolutions specified.
        numberOfRevolutions = copyOfOriginalNumberOfRevolutions;
    }

    // Default: compare maximum with specified number of revolutions.
    // If specified is larger than maximum, no solution is possible.
    if ( numberOfRevolutions > maximumNumberOfRevolutions )
    {
        // Throw exception.
        BOOST_THROW_EXCEPTION( std::runtime_error( 
            ( boost::format(
                "Number of revolutions specified in Lambert problem is larger than possible.\n"
                "Specified number of revolutions %d while the maximum is %d" 
            ) % numberOfRevolutions % maximumNumberOfRevolutions).str( )
        ) );
    }
    // Else, nothing wrong.
}

//! Execute solving procedure (for multiple revolutions).
void MultiRevolutionLambertTargeterIzzo::execute( )
{
    // Sanity checks.
    sanityCheckTimeOfFlight( );
    sanityCheckGravitationalParameter( );

    // Transform dimensions.
    transformDimensions( );

    /*// Sanity check for number of revolutions (must be after dimension removal).
    sanityCheckNumberOfRevolutions( );*/

    if ( numberOfRevolutions == 0 )
    {
        // call base class function that works on zero revolutions.
        ZeroRevolutionLambertTargeterIzzo::execute( );
    }
    else
    {
        // Solve multi-rev root.
        double multipleRevolutionXParameter = computeRootTimeOfFlight( );

        // Reconstruct velocities.
        computeVelocities( multipleRevolutionXParameter );
    }

    solved = true;
}

//! Compute time-of-flight using Lagrange's equation (for multiple revolutions).
double MultiRevolutionLambertTargeterIzzo::computeTimeOfFlight( const double xParameter )
{
    // Determine semi-major axis.
    const double semiMajorAxis = normalizedMinimumEnergySemiMajorAxis
            / ( 1.0 - xParameter * xParameter );

    // If x < 1, the solution is an ellipse.
    if ( xParameter < 1.0 )
    {
        // Alpha parameter in Lagrange's equation (no explanation available).
        const double alphaParameter = 2.0 * std::acos( xParameter );

        // Beta parameter in Lagrange's equation (no explanation available).
        double betaParameter;

        // If long transfer arc.
        if ( isLongway )
        {
            betaParameter = -2.0 * std::asin(
                        std::sqrt( ( normalizedSemiPerimeter - normalizedChord )
                                   / ( 2.0 * semiMajorAxis ) ) );
        }
        // Otherwise short transfer arc.
        else
        {
            betaParameter = 2.0 * std::asin(
                        std::sqrt( ( normalizedSemiPerimeter - normalizedChord )
                                   / ( 2.0 * semiMajorAxis ) ) );
        }

        // Time-of-flight according to Lagrange including multiple revolutions.
        const double timeOfFlight = semiMajorAxis * std::sqrt( semiMajorAxis ) *
                ( ( alphaParameter - std::sin( alphaParameter ) )
                  - ( betaParameter - std::sin( betaParameter ) )
                  + 2.0 * mathematical_constants::PI
                  * numberOfRevolutions );

        return timeOfFlight;
    }
    // Otherwise it is a hyperbola.
    else
    {
        // Alpha parameter in Lagrange's equation (no explanation available).
        const double alphaParameter = 2.0 * boost::math::acosh( xParameter );

        // Beta parameter in Lagrange's equation (no explanation available).
        double betaParameter;

        // If long transfer arc.
        if ( isLongway )
        {
            betaParameter = -2.0 * boost::math::asinh( std::sqrt( ( normalizedSemiPerimeter
                                                                    - normalizedChord )
                                                                  / ( -2.0 * semiMajorAxis ) ) );
        }
        // Otherwise short transfer arc
        else
        {
            betaParameter = 2.0 * boost::math::asinh( std::sqrt( ( normalizedSemiPerimeter
                                                                   - normalizedChord )
                                                                 / ( -2.0 * semiMajorAxis ) ) );
        }

        // Time-of-flight according to Lagrange.
        const double timeOfFlightLagrange = -semiMajorAxis * std::sqrt( -semiMajorAxis ) *
                ( ( std::sinh( alphaParameter ) - alphaParameter )
                  - ( std::sinh( betaParameter ) - betaParameter ) );

        return timeOfFlightLagrange;
    }
}

//! Solve the time of flight equation for x (for multiple revolutions).
double MultiRevolutionLambertTargeterIzzo::computeRootTimeOfFlight( )
{
    using mathematical_constants::PI;

    // Define initial guesses for abcissae (x) and ordinates (y).
    double x1, x2;

    if ( isRightBranch )
    { // right branch solution.
        x1 = std::tan( .7234 * PI / 2.0 );
        x2 = std::tan( .5234 * PI / 2.0 );
    }
    else
    { // left branch solution.
        x1 = std::tan( -.5234 * PI / 2.0 );
        x2 = std::tan( -.2234 * PI / 2.0 );
    }

    double y1 = computeTimeOfFlight( std::atan( x1 ) * 2.0 / PI ) - normalizedTimeOfFlight;

    double y2 = computeTimeOfFlight( std::atan( x2 ) * 2.0 / PI ) - normalizedTimeOfFlight;

    // Declare and initialize root-finding parameters.
    double rootFindingError = 1.0, xNew = 0.0, yNew = 0.0;
    int iterator = 0;

    // Root-finding loop.
    while ( ( rootFindingError > convergenceTolerance ) && ( y1 != y2 )
            && ( iterator < maximumNumberOfIterations ) )
    {
        // Update iterator.
        iterator++;

        // Compute new x-value.
        xNew = ( x1 * y2 - y1 * x2 ) / ( y2 - y1 );

        // Compute corresponding y-value.
        yNew = computeTimeOfFlight( std::atan( xNew ) * 2.0 / PI ) - normalizedTimeOfFlight;

        // Update abcissae and ordinates.
        x1 = x2;
        y1 = y2;
        x2 = xNew;
        y2 = yNew;

        // Compute root-finding error.
        rootFindingError = std::fabs( x1 - xNew );
    }

    // Verify that root-finder has converged.
    if ( iterator == maximumNumberOfIterations )
    {
        BOOST_THROW_EXCEPTION( basic_mathematics::ConvergenceException( 
            ( boost::format(
                "Multi-Revolution Lambert targeter failed to converge to a solution.\n"
                "Reached the maximum number of iterations: %d"
            ) % maximumNumberOfIterations).str( )
        ) );
    }

    // Revert to x parameter.
    double xParameter = std::atan( xNew ) * 2.0 / PI;
    return xParameter;
}

// Add compute maximum number of revolutions routine?

} // namespace mission_segments
} // namespace tudat
