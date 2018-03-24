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
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2002.
 *      Press W.H., et al. Numerical Recipes in C++: The Art of Scientific Computing. Cambridge
 *          University Press, February 2007.
 *      Weisstein, Eric W. "Halley's Method." From MathWorld -- A Wolfram Web Resource.
 *          http://mathworld.wolfram.com/HalleysMethod.html, retrieved on 19/02/2014.
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
