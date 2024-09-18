/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_NEWTON_RAPHSON_H
#define TUDAT_NEWTON_RAPHSON_H





#include <memory>

#include "tudat/math/root_finders/rootFinder.h"
#include "tudat/math/root_finders/terminationConditions.h"
#include "tudat/math/basic/convergenceException.h"

namespace tudat
{
namespace root_finders
{

//! Newton-Raphson rootfinder.
/*!
 * Rootfinder using Newton-Raphson's Method. It requires a function, and its first derivative.
 * To start the process, also one initial guess of the root value is required.
 *
 * The iterative scheme is given by:
 *
 * \f[
 *  x_{n+1} = x_n - \frac{F\left(x_n\right)}{F'\left(x_n\right)}
 * \f]
 *
 * Defined shorthand notations:
 *  NewtonRaphsonCore< double, double >   =>   NewtonRaphson
 *
 * \tparam DataType Data type used to represent floating-point values.
 */
template< typename DataType = double >
class NewtonRaphson : public RootFinder< DataType >
{
public:

    //! Useful type definition for the function pointer (from base class)
    typedef typename RootFinder< DataType >::FunctionPointer FunctionPointer;

    //! Useful type definition for the termination function (from base class)
    typedef typename RootFinder< DataType >::TerminationFunction TerminationFunction;

    //! This is the constructor taking the general termination function
    /*!
     *  This is the constructor taking the general termination function
     *  \param terminationFunction The function specifying the termination conditions of the
     *  root-finding process \sa RootFinderCore::terminationFunction
     */
    NewtonRaphson( TerminationFunction terminationFunction )
        : RootFinder< DataType >( terminationFunction )
    { }

    //! Constructor taking typical convergence criteria.
    /*!
     *  Constructor taking maximum number of iterations and relative tolerance for independent
     *  variable. If desired, a custom convergence function can provided to the alternative
     *  constructor
     *  \param relativeIndependentVariableTolerance Relative difference between the root solution of two subsequent
     *  solutions below which convergence is reached.
     *  \param maxIterations Maximum number of iterations after which the root finder is
     *  terminated, i.e. convergence is assumed
     */
    NewtonRaphson( const DataType relativeIndependentVariableTolerance, const unsigned int maxIterations )
        : RootFinder< DataType >(
              createTerminationCondition( relativeIndependentVariableTolerance, TUDAT_NAN, TUDAT_NAN, maxIterations ) )
    {}

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~NewtonRaphson( ) { }

    //! Find a root of the function provided as input.
    /*!
     *  Find a root of the function provided as input, using the termination function set by the
     *  constructor
     * \param rootFunction Function to find root of.
     * \param initialGuess The initial guess of the root.
     * \throws ConvergenceExeption If the solution does not converge to a root value.
     * \return Root of the rootFunction that is found
     */
    DataType execute( const FunctionPointer rootFunction, const DataType initialGuess )
    {
        // Set the root function.
        this->rootFunction = rootFunction;

        // Start at initial guess, and compute the function value and its first derivative.
        DataType currentRootValue       = TUDAT_NAN;
        DataType nextRootValue          = initialGuess;
        DataType currentFunctionValue   = TUDAT_NAN;
        DataType nextFunctionValue      = this->rootFunction->evaluate( nextRootValue );
        DataType currentDerivativeValue = TUDAT_NAN;
        DataType nextDerivativeValue    = this->rootFunction->computeDerivative( 1, nextRootValue );

        if( !( nextFunctionValue == mathematical_constants::getFloatingInteger< DataType >( 0 ) ) )
        {
            // Loop counter.
            unsigned int counter = 1;

            // Loop until we have a solution with sufficient accuracy.
            do
            {
                // Save the old values.
                currentRootValue       = nextRootValue;
                currentFunctionValue   = nextFunctionValue;
                currentDerivativeValue = nextDerivativeValue;

                // Compute next value of root using the following algorithm (see class documentation):
                nextRootValue          = currentRootValue -
                        currentFunctionValue / currentDerivativeValue;
                nextFunctionValue      = this->rootFunction->evaluate( nextRootValue );
                nextDerivativeValue    = this->rootFunction->computeDerivative( 1, nextRootValue );

                // Update the counter.
                counter++;
            }
            while( nextFunctionValue != mathematical_constants::getFloatingInteger< DataType >( 0 ) &&
                   !this->terminationFunction_( nextRootValue, currentRootValue, nextFunctionValue,
                                                currentFunctionValue, counter ) );
        }

        return nextRootValue;
    }

protected:

private:

};

} // namespace root_finders
} // namespace tudat

#endif // TUDAT_NEWTON_RAPHSON_H
