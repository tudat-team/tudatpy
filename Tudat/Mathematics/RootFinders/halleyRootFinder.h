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
 *      150417    D. Dirkx          Made modifications for templated root finding.
 *
 *    References
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2007.
 *      Weisstein, Eric W. "Halley's Method." From MathWorld -- A Wolfram Web Resource.
 *          http://mathworld.wolfram.com/HalleysMethod.html, retrieved on 19/02/2014.
 *
 *    Notes
 *
 */

#ifndef TUDAT_HALLEY_ROOT_FINDER_H
#define TUDAT_HALLEY_ROOT_FINDER_H

#include <boost/bind.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace root_finders
{

//! Halley's root-finding method.
/*!
 * Halley's root-finding method, also known as the tangent hyperbolas method or Halley's rational
 * formula. It requires a function and its first two derivatives. To start the algorithm, also an
 * initial guess of the root location is needed.
 *
 * It only makes sense to use Halley's method when it is easy to compute the second derivative
 * \f$f''(x_i)\f$, often from other pieces of functions that are already being used in the
 * computation of \f$f(x_i)\f$ and \f$f'(x_i)\f$. Otherwise, Newton-Rapshon with an increased number
 * of iterations is preferable. [Press et al., 2007, p.463]
 *
 * The iterative scheme is given by: [Weisstein, 2014]
 *
 * \f[
 *  x_{n+1} = x_n-\frac{2 F\left(x_n\right) F'\left(x_n\right)}
 *                                 {2 F'\left(x_n\right){}^2-F\left(x_n\right) F''\left(x_n\right)}
 * \f]
 *
 * Defined shorthand notations:
 *  HalleyRootFinderCore< double >   =>   HalleyRootFinder
 *
 * \tparam DataType Data type used to represent floating-point values.
 */
template< typename DataType = double >
class HalleyRootFinderCore : public RootFinderCore< DataType >
{
public:

    //! Usefull type definition for the function pointer (from base class).
    typedef typename RootFinderCore< DataType >::FunctionPointer FunctionPointer;

    //! Usefull type definition for the termination function (from base class).
    typedef typename RootFinderCore< DataType >::TerminationFunction TerminationFunction;

    //! Constructor taking the general termination function.
    /*!
     * Constructor of Halley's root-finding method, taking the termination function (function
     * determining whether to terminate the root-finding process).
     *
     * \param terminationFunction The function specifying the termination conditions of the
     *          root-finding process. \sa RootFinderCore::terminationFunction
     */
    HalleyRootFinderCore( TerminationFunction terminationFunction )
        : RootFinderCore< DataType >( terminationFunction )
    { }

    //! Constructor taking typical convergence criteria.
    /*!
     * Constructor of Halley's root-finding method, taking the maximum number of iterations and the
     * relative tolerance for the independent variable. If desired, a custom convergence function
     * can provided to the alternative constructor.

     *  \param relativeXTolerance Relative difference between the root solution of two subsequent
     *          solutions below which convergence is reached.
     *  \param maxIterations Maximum number of iterations after which the root finder is
     *          terminated, i.e. convergence is assumed.
     */
    HalleyRootFinderCore( const DataType relativeXTolerance, const unsigned int maxIterations )
        : RootFinderCore< DataType >(
              boost::bind(
                  &termination_conditions::RootRelativeToleranceTerminationCondition< DataType >::
                  checkTerminationCondition, boost::make_shared<
                  termination_conditions::RootRelativeToleranceTerminationCondition< DataType > >(
                      relativeXTolerance, maxIterations ), _1, _2, _3, _4, _5 ) )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~HalleyRootFinderCore( ){ }

    //! Find a root of the function provided as input.
    /*!
     * Find a root of the function provided as input, using the termination function set by the
     * constructor. Not only the function but also the first two derivatives are required.
     *
     * \param rootFunction Function to find the root of.
     * \param initialGuess The initial guess of the root.
     * \return Root of the rootFunction that is found.
     *
     * \throws ConvergenceExeption If the solution does not converge to a root value.
     */
    DataType execute( const FunctionPointer rootFunction, const DataType initialGuess )
    {
        // Set the root function.
        this->rootFunction = rootFunction;

        // Start at initial guess, and compute the function value, its first and second derivative.
        DataType currentRootValue             = TUDAT_NAN;
        DataType nextRootValue                = initialGuess;
        DataType currentFunctionValue         = TUDAT_NAN;
        DataType nextFunctionValue            = this->rootFunction->evaluate( nextRootValue );
        DataType currentFirstDerivativeValue  = TUDAT_NAN;
        DataType nextFirstDerivativeValue     = this->rootFunction->
                computeDerivative( 1, nextRootValue );
        DataType currentSecondDerivativeValue = TUDAT_NAN;
        DataType nextSecondDerivativeValue    = this->rootFunction->
                computeDerivative( 2, nextRootValue );

        // Loop counter.
        unsigned int counter = 1;

        // Loop until we have a solution with sufficient accuracy.
        do
        {
            // Save the old values.
            currentRootValue             = nextRootValue;
            currentFunctionValue         = nextFunctionValue;
            currentFirstDerivativeValue  = nextFirstDerivativeValue;
            currentSecondDerivativeValue = nextSecondDerivativeValue;

            // Compute next value of root using the following algorithm (see class documentation):
            nextRootValue               = currentRootValue
                    - ( ( 2.0 * currentFunctionValue * currentFirstDerivativeValue )
                        / ( 2.0 * currentFirstDerivativeValue * currentFirstDerivativeValue
                            - currentFunctionValue * currentSecondDerivativeValue ) );
            nextFunctionValue           = this->rootFunction->evaluate( nextRootValue );
            nextFirstDerivativeValue    = this->rootFunction->computeDerivative( 1, nextRootValue );
            nextSecondDerivativeValue   = this->rootFunction->computeDerivative( 2, nextRootValue );

            // Update the counter.
            counter++;
        }
        while( !this->terminationFunction( nextRootValue, currentRootValue, nextFunctionValue,
                                           currentFunctionValue, counter ) );

        return nextRootValue;
    }

protected:

private:

};

// Some handy typedefs.
typedef HalleyRootFinderCore< double > HalleyRootFinder;
typedef boost::shared_ptr< HalleyRootFinder > HalleyRootFinderPointer;

} // namespace root_finders
} // namespace tudat

#endif // TUDAT_HALLEY_ROOT_FINDER_H
