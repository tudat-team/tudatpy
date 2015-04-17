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
 *      101111    E. Iorfida        File created.
 *      101116    E. Iorfida        Added setFunction and execute.
 *      101121    E. Iorfida        Added Doxygen comments.
 *      110111    E. Iorfida        Deleted useless lines, and modified punctuation.
 *      110111    K. Kumar          Changed variable and function names to be more descriptive;
 *                                  added "End of file."
 *      110114    K. Kumar          Removed circular code dependency.
 *      110119    K. Kumar          Updated code to work with adaptor and abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required; changed filename.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120208    S. Billemont      Move to new root_finders codebase.
 *      120402    T. Secretin       Code-check.
 *      120726    S. Billemont      Restructuring. Implemented new termination conditions.
 *      120810    P. Musegaas       Code-check.
 *      140312    E. Brandon        Removed unneccesary include statement, minor correction Doxygen
 *                                  comments.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_NEWTON_RAPHSON_H
#define TUDAT_NEWTON_RAPHSON_H

#include <boost/bind.hpp>
#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"
#include "Tudat/Mathematics/BasicMathematics/convergenceException.h"

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
class NewtonRaphsonCore : public RootFinderCore< DataType >
{
public:

    //! Usefull type definition for the function pointer (from base class)
    typedef typename RootFinderCore< DataType >::FunctionPointer FunctionPointer;

    //! Usefull type definition for the termination function (from base class)
    typedef typename RootFinderCore< DataType >::TerminationFunction TerminationFunction;

    //! This is the constructor taking the general termination function
    /*!
     *  This is the constructor taking the general termination function
     *  \param terminationFunction The function specifying the termination conditions of the
     *  root-finding process \sa RootFinderCore::terminationFunction
     */
    NewtonRaphsonCore( TerminationFunction terminationFunction )
        : RootFinderCore< DataType >( terminationFunction )
    { }

    //! Constructor taking typical convergence criteria.
    /*!
     *  Constructor taking maximum number of iterations and relative tolerance for independent
     *  variable. If desired, a custom convergence function can provided to the alternative
     *  constructor
     *  \param relativeXTolerance Relative difference between the root solution of two subsequent
     *  solutions below which convergence is reached.
     *  \param maxIterations Maximum number of iterations after which the root finder is
     *  terminated, i.e. convergence is assumed
     */
    NewtonRaphsonCore( const double relativeXTolerance, const unsigned int maxIterations );
	
    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~NewtonRaphsonCore( ) { }

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
        DataType nextDerivativeValue    = this->rootFunction->
                computeDerivative( 1, nextRootValue );

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
        while( !this->terminationFunction( nextRootValue, currentRootValue, nextFunctionValue,
                                           currentFunctionValue, counter ) );
		
        return nextRootValue;
    }

protected:

private:
};

//! Constructor taking typical convergence criteria.
template< typename DataType >
NewtonRaphsonCore< DataType >::NewtonRaphsonCore( const double relativeXTolerance,
                                                  const unsigned int maxIterations )
    : RootFinderCore< DataType >(
          boost::bind( &termination_conditions::RootRelativeToleranceTerminationCondition::
                       checkTerminationCondition, boost::make_shared<
                       termination_conditions::RootRelativeToleranceTerminationCondition >(
                           relativeXTolerance, maxIterations ), _1, _2, _3, _4, _5 ) )
{ }

// Some handy typedefs.
typedef NewtonRaphsonCore< > NewtonRaphson;
typedef boost::shared_ptr< NewtonRaphson > NewtonRaphsonPointer;

} // namespace root_finders
} // namespace tudat

#endif // TUDAT_NEWTON_RAPHSON_H
