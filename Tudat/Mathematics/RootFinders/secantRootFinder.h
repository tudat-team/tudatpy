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
 *      120208    S. Billemont      Creation of code.
 *      140219    E. Brandon        Adapted to current Tudat root-finder structure.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2007.
 *      Weisstein, Eric W. "Secant Method." From MathWorld -- A Wolfram Web Resource.
 *          http://mathworld.wolfram.com/SecantMethod.html, retrieved on 19/02/2014.
 *
 *    Notes
 *
 */

#ifndef TUDAT_SECANT_ROOT_FINDER_H
#define TUDAT_SECANT_ROOT_FINDER_H

#include <cmath>

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace root_finders
{

//! Secant root-finding method.
/*!
 * The secant root-finding method uses a succession of roots of secant lines to better aspproximate
 * the root. The method requires a valid function (no derivatives are required), and two initial
 * values ($x_0$ and $x_1$), which should ideally be chosen to lie close to the root, and where
 * $x_1$ is the most accurate guess. It is not required that the initial values bracket the root.
 * Note that this root finder does not always converge.
 *
 * The algorithm is similar to the Newton-Raphson root-finder, but uses the slope between the two
 * points as the value for the function derivative [Weisstein, 2014]:
 *
 * \f[
 *  x_{n+1} = x_n - F\left(x_n\right) \frac{\left(x_n - x_{n-1}\right)}
 *            {F\left(x_n\right) - F\left(x_{n-1}\right)}
 * \f]
 *
 * Defined shorthand notations:
 *  SecantRootFinderCore< double >   =>   SecantRootFinder
 *
 * \tparam DataType Data type used to represent floating-point values.
 */
template< typename DataType = double >
class SecantRootFinderCore : public RootFinderCore< DataType >
{
public:

    //! Usefull type definition for the function pointer (from base class).
    typedef typename RootFinderCore< DataType >::FunctionPointer FunctionPointer;

    //! Usefull type definition for the termination function (from base class).
    typedef typename RootFinderCore< DataType >::TerminationFunction TerminationFunction;

    //! Constructor taking the general termination function and the least accurate initial guess.
    /*!
     * Constructor of the Secant root-finder, taking the general termination function (function
     * determining whether to terminate the root-finding process), and the least accurate value of
     * the initial guess ($x_0$). (Note that although it is prefered that the least accurate guess
     * is passed to the constructor, this is not required, as the algorithm will determine the most
     * accurate guess and switches the values if necessary.)
     *
     * \param terminationFunction The function specifying the termination conditions of the
     *  root-finding process. \sa RootFinderCore::terminationFunction
     * \param initialGuessOfRootOne First point used to initiate the Secant root-finder algorithm.
     *          (Default is 0.5)
     */
    SecantRootFinderCore( TerminationFunction terminationFunction,
                          const DataType initialGuessOfRootOne = 0.5 )
        : RootFinderCore< DataType >( terminationFunction ),
          initialGuessOfRootOne_( initialGuessOfRootOne )
    { }

    //! Constructor taking typical convergence criteria and the least accurate initial guess.
    /*!
     * Constructor of the Secant root-finder, taking the maximum number of iterations, the relative
     * tolerance for the independent variable, and the least accurate value of the initial guess
     * ($x_0$). (Note that although it is prefered that the least accurate guess is passed to the
     * constructor, this is not required, as the algorithm will determine the most accurate guess
     * and switches the values if necessary.) If desired, a custom convergence function can be
     * provided to the alternative constructor.
     *
     * \param relativeXTolerance Relative difference between the root solution of two subsequent
     *          solutions below which convergence is reached.
     * \param maxIterations Maximum number of iterations after which the root finder is
     *          terminated, i.e. convergence is assumed.
     * \param initialGuessOfRootOne First point used to initiate the Secant root-finder algorithm.
     *          (Default is 0.5)
     */
    SecantRootFinderCore( const double relativeXTolerance, const unsigned int maxIterations,
                          const DataType initialGuessOfRootOne = 0.5 );

    //! Default destructor.
    ~SecantRootFinderCore( ) { }

    //! Find a root of the function provided as input.
    /*!
     * Finds a root of the function provided as input, using the termination function set by the
     * constructor. This method takes the second initial value ($x_1$), which is the most accurate
     * guess. (Note that although it is prefered that the most accurate guess is passed to this
     * function, this is not required, as the algorithm will determine the most accurate guess
     * and switches the values if necessary.)
     *
     * \param rootFunction Function to find the root of.
     * \param initialGuess The second, most accurate, initial guess of the root.
     * \return Root of the rootFunction that is found.
     *
     * \throws ConvergenceExeption If the solution does not converge to a root value.
     */
    DataType execute( const FunctionPointer rootFunction, const DataType initialGuess )
    {
        // Set the root function.
        this->rootFunction = rootFunction;

        // Start at the two initial values that are used in the algorithm, and compute the function
        // values.
        DataType lastRootValue          = TUDAT_NAN;
        DataType currentRootValue       = initialGuessOfRootOne_;
        DataType nextRootValue          = initialGuess;
        DataType lastFunctionValue      = TUDAT_NAN;
        DataType currentFunctionValue   = this->rootFunction->evaluate( currentRootValue );
        DataType nextFunctionValue      = this->rootFunction->evaluate( nextRootValue );

        // Check if the next root value is the most accurate guess. If not, switch the values.
        if( std::fabs( currentFunctionValue ) < std::fabs( nextFunctionValue ) )
        {
            // Switch the root value.
            DataType temporaryRoot  = currentRootValue;
            currentRootValue        = nextRootValue;
            nextRootValue           = temporaryRoot;

            // Switch the function value
            DataType temporaryValue = currentFunctionValue;
            currentFunctionValue    = nextFunctionValue;
            nextFunctionValue       = temporaryValue;
        }

        // Loop counter.
        unsigned int counter = 1;

        // Loop until we have a solution with sufficient accuracy.
        do
        {
            // Save the old values.
            lastRootValue           = currentRootValue;
            lastFunctionValue       = currentFunctionValue;
            currentRootValue        = nextRootValue;
            currentFunctionValue    = nextFunctionValue;

            // Compute next value of root using the following algorithm (see class documentation):
            nextRootValue           = currentRootValue - currentFunctionValue
                                      * ( currentRootValue - lastRootValue )
                                      / ( currentFunctionValue - lastFunctionValue );
            nextFunctionValue       = this->rootFunction->evaluate( nextRootValue );

            // Update the counter.
            counter++;
        }
        while( !this->terminationFunction( nextRootValue, currentRootValue, nextFunctionValue,
                                           currentFunctionValue, counter ) );

        return nextRootValue;

    }

    //! Set a new value for the first point used in the Secant algorithm.
    /**
     * Sets a new value for least accurate value of the initial guess ($x_0$) used in the Secant
     * algorithm. (Note that although it is prefered that the least accurate guess is passed to the
     * constructor, this is not required, as the algorithm will determine the most accurate guess
     * and switches the values if necessary.)
     *
     * \param newInitialGuessOfRootOne A new value of the first point used to initiate the Secant
     *          root-finder algorithm.
     */
    void setNewInitialGuess( DataType newInitialGuessOfRootOne )
    {
        this->initialGuessOfRootOne_ = newInitialGuessOfRootOne;
    }

protected:

private:

    //! Initial value of the first point (least accurate guess) in the Secant algorithm.
    double initialGuessOfRootOne_;

};

//! Constructor taking typical convergence criteria and the least accurate initial guess.
template< typename DataType >
SecantRootFinderCore< DataType >::SecantRootFinderCore( const double relativeXTolerance,
                                                        const unsigned int maxIterations,
                                                        const DataType initialGuessOfRootOne )
    : RootFinderCore< DataType >(
          boost::bind( &termination_conditions::RootRelativeToleranceTerminationCondition::
                       checkTerminationCondition, boost::make_shared<
                       termination_conditions::RootRelativeToleranceTerminationCondition >(
                           relativeXTolerance, maxIterations ), _1, _2, _3, _4, _5 ) ),
      initialGuessOfRootOne_( initialGuessOfRootOne )
{ }

// Some handy typedefs.
typedef SecantRootFinderCore< double > SecantRootFinder;
typedef boost::shared_ptr< SecantRootFinder > SecantRootFinderPointer;

} // namespace root_finders
} // namespace tudat

#endif // TUDAT_SECANT_ROOT_FINDER_H
