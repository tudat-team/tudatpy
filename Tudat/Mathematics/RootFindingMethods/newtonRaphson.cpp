/*    Copyright (c) 2010-2012, Delft University of Technology
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
 *      101111    E. Iorfida        Creation of code.
 *      101118    E. Iorfida        Added algorithm implementation.
 *      101121    E. Iorfida        Modified algorithm and added Doxygen comments.
 *      101216    E. Iorfida        Modified punctuation.
 *      110111    E. Iorfida        Deleted useless lines, and modified punctuation.
 *      110111    K. Kumar          Updated based on changes to .h file; changed do-while loop to
 *                                  for loop to prevent infinite loop; added "End of file."
 *      110119    K. Kumar          Updated code to work with adaptor and abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required; changed filename.
 *      110124    K. Kumar          Added cerr statement for when loop does not converge.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120712    P. Musegaas       Changed absolute tolerance into a safe variant of relative
 *                                  tolerance.
 *
 *    References
 *
 */

#include <cmath>

#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"

namespace tudat
{

//! Set adaptor class for Newton-Raphson.
void NewtonRaphson::setNewtonRaphsonAdaptor( NewtonRaphsonBase* pointerToNewtonRaphsonBase )
{
    pointerToNewtonRaphsonBase_ = pointerToNewtonRaphsonBase;
    pointerToGlobalFunction_ = NULL;
    pointerToGlobalFirstDerivativeFunction_ = NULL;
}

//! Execute Newton-Raphson method.
void NewtonRaphson::execute( )
{
    // Algorithm implementation.
    // First step.
    nextValueOfRoot_ = initialGuessOfRoot_;

    for ( unsigned int i = 0; i < maximumNumberOfIterations_; i++ )
    {
        // Update current value of root.
        currentValueOfRoot_ = nextValueOfRoot_;

        // Compute next value of root using the following algorithm:
        // x_n+1 = x_n - F(x_n)/F'(x_n).

        // Check if necessary global mathematical functions are set.
        if ( pointerToGlobalFunction_ && pointerToGlobalFirstDerivativeFunction_ )
        {
            // Use pointers to global functions.
            nextValueOfRoot_ = currentValueOfRoot_
                    - ( *pointerToGlobalFunction_ )( currentValueOfRoot_ )
                    / ( *pointerToGlobalFirstDerivativeFunction_ )( currentValueOfRoot_ );
        }

        // Else check if class with member mathematical functions is set.
        else if ( pointerToNewtonRaphsonBase_ )
        {
            // Use pointers to member functions.
            nextValueOfRoot_ = currentValueOfRoot_
                    - pointerToNewtonRaphsonBase_->computeFunction( currentValueOfRoot_ )
                    /  pointerToNewtonRaphsonBase_->computeFirstDerivativeFunction(
                        currentValueOfRoot_ );
        }

        // Check if difference between successive iterations of the  root-finding method satisfies
        // the set tolerance and break from loop if it is satisfied.
        if ( std::fabs( ( nextValueOfRoot_ - currentValueOfRoot_ ) / currentValueOfRoot_ ) <=
             relativeTolerance_ )
        {
            break;
        }

        // Check if the root-finding method converges to zero. This check is required because a
        // relative tolerance is used.
        if ( std::fabs( nextValueOfRoot_ ) < zeroRepresentation_ )
        {
            // Set the root to 0.0 and return.
            nextValueOfRoot_ = 0.0;
            break;
        }

        // If the end of the loop is reached before the tolerance is satisfied, print cerr
        // statement to indicate to the user that there method did not converge.
        if ( i == maximumNumberOfIterations_ - 1 )
        {
            // Reset values of root to initialized values.
            currentValueOfRoot_ = -0.0;
            nextValueOfRoot_ = -0.0;

            std::cerr << "The Newton-Raphson algorithm did not converge after "
                      << maximumNumberOfIterations_ << " iterations." << std::endl;
        }
    }
}

//! Overload ostream to print class information.
std::ostream& operator<<( std::ostream& stream, NewtonRaphson& newtonRaphson )
{
    stream << "This is a NewtonRaphson object" << std::endl;
    stream << "The maximum number of iterations is set to: "
           << newtonRaphson.getMaximumNumberOfIterations( )
           << "The tolerance is set to: " << newtonRaphson.getRelativeTolerance( )
           << "The initial guess of root is set to: " << newtonRaphson.initialGuessOfRoot_
           << "The computed root is: " << newtonRaphson.getComputedRootOfFunction( ) << std::endl;

    // Return stream.
    return stream;
}

} // namespace tudat
