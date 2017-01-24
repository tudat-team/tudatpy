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
 *      Weisstein, Eric W. "Bisection." From MathWorld -- A Wolfram Web Resource.
 *          http://mathworld.wolfram.com/Bisection.html, retrieved on 19/02/2014.
 *
 *    Notes
 *
 */

#ifndef TUDAT_BISECTION_H
#define TUDAT_BISECTION_H

#include <boost/bind.hpp>
#include <boost/exception/all.hpp>
#include <boost/format.hpp>
#include <boost/format/free_funcs.hpp>
#include <boost/make_shared.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Basics/utilityMacros.h"
#include "Tudat/Mathematics/RootFinders/rootFinder.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace root_finders
{

//! Bisection root-finding method.
/*!
 * The bisection root-finding method, a basic and robust root-finder that will always find a root
 * given that a root exists, the function is continuous on the interval, and that it is bracketed
 * by the lower and upper bound (required). For this method only the function of which the zero is
 * sought is required, and no derivatives. It is recommended to use this method for validation only,
 * as it is relatively slow.
 *
 * It works by repeatedly shrinking the [lowerbound, upperbound] interval until the root has been
 * found with sufficient accuracy. The shrinking is done by dividing the interval in half and
 * evaluating in which interval the root lies. This sub-interval is then kept for the next
 * iteration. The process is continued until the interval is sufficiently small.
 * [Press et al., 2007]
 *
 * Defined shorthand notations:
 *  BisectionCore< double >   =>   Bisection
 *
 * \tparam DataType Data type used to represent floating-point values.
 */
template< typename DataType = double >
class BisectionCore : public RootFinderCore< DataType >
{
public:

    //! Useful type definition for the function pointer (from base class).
    typedef typename RootFinderCore< DataType >::FunctionPointer FunctionPointer;

    //! Useful type definition for the termination function (from base class).
    typedef typename RootFinderCore< DataType >::TerminationFunction TerminationFunction;

    //! Constructor taking a general termination function and the bracket of the solution.
    /*!
     * Constructor of the Bisection root-finder, taking the termination function (function
     * determining whether to terminate the root-finding process) and the search interval with an
     * upper and lower bound. It is required that the function values at upper and lower bound have
     * an opposite sign. The default interval is [-1, 1].
     *
     * \param terminationFunction The function specifying the termination conditions of the
     *          root-finding process. \sa RootFinderCore::terminationFunction
     * \param lowerBound Lower bound of the interval containing a root. (Default is -1.0).
     * \param upperBound Upper bound of the interval containing a root. (Default is 1.0).
     */
    BisectionCore( const TerminationFunction terminationFunction,
                   const DataType lowerBound = -1.0,
                   const DataType upperBound = 1.0 )
        : RootFinderCore< DataType >( terminationFunction ),
          lowerBound_( lowerBound ),
          upperBound_( upperBound )
    { }

    //! Constructor taking typical convergence criteria and the bracket of the solution.
    /*!
     * Constructor of the Bisection root-finder, taking the maximum number of iterations, the
     * relative tolerance for the independent variable, and the search interval with upper and
     * lower bound. It is required that the function values at upper and lower bound have an
     * opposite sign. The default interval is [-1, 1]. If desired, a custom convergence function
     * can be provided to the alternative constructor.
     *
     * \param relativeXTolerance Relative difference between the root solution of two subsequent
     *          solutions below which convergence is reached.
     * \param maxIterations Maximum number of iterations after which the root finder is
     *          terminated, i.e. convergence is assumed.
     * \param lowerBound Lower bound of the interval containing a root. (Default is -1.0).
     * \param upperBound Upper bound of the interval containing a root. (Default is 1.0).
     */
    BisectionCore( const DataType relativeXTolerance, const unsigned int maxIterations,
                   const DataType lowerBound = -1.0, const DataType upperBound = 1.0 ):
        RootFinderCore< DataType >(
            boost::bind(
                &termination_conditions::RootRelativeToleranceTerminationCondition< DataType >::
                checkTerminationCondition, boost::make_shared<
                termination_conditions::RootRelativeToleranceTerminationCondition< DataType > >(
                    relativeXTolerance, maxIterations ), _1, _2, _3, _4, _5 ) ),
        lowerBound_( lowerBound ),
        upperBound_( upperBound )
    { }

    //! Default destructor.
    ~BisectionCore( ) { }

    //! Find a root of the function provided as input.
    /*!
     * Find a root of the function provided as input, using the termination function set by the
     * constructor. (Note that the initial guess is not used, but is a requirement of the
     * root-finder architecture.)
     *
     * \param rootFunction Function to find the root of.
     * \param initialGuess The initial guess of the root. (Not used, default is 0.0).
     * \return Root of the rootFunction that is found.
     *
     * \throws ConvergenceExeption If the solution does not converge to a root value.
     * \throws std::runtime_error If the interval does not bracket the root.
     */
    DataType execute( const FunctionPointer rootFunction, const DataType initialGuess = 0.0 )
    {
        // The value of the initialGuess is not used.
        TUDAT_UNUSED_PARAMETER( initialGuess );

        // Set the root function.
        this->rootFunction = rootFunction;

        // Initialize previous values.
        DataType previousRootValue = TUDAT_NAN;
        DataType previousRootFunctionValue = TUDAT_NAN;

        // Duplicate the interval and use this duplicate to shrink the interval.
        // Compute the midpoint of the interval, and take this as the current root.
        DataType currentLowerBound = lowerBound_;
        DataType currentUpperBound = upperBound_;
        DataType rootValue = ( currentLowerBound + currentUpperBound ) / 2.0;

        // Find the corresponding function values at the important interval points (lower bound,
        // upper bound and midpoint).
        DataType currentLowerBoundFunctionValue = this->rootFunction->evaluate( currentLowerBound );
        DataType currentUpperBoundFunctionValue = this->rootFunction->evaluate( currentUpperBound );
        DataType rootFunctionValue = this->rootFunction->evaluate( rootValue );

        // Validate that upperbound and lowerbound function values have different signs
        // (requirement).
        if( currentLowerBoundFunctionValue * currentUpperBoundFunctionValue > 0.0 )
        {
            boost::throw_exception( boost::enable_error_info( std::runtime_error(
                                                                  boost::str( boost::format(
                                                                                  "The Bisection algorithm requires that the values at the upper "
                                                                                  "and lower bounds have a different sign." ) ) ) ) );
        }

        // Loop counter.
        unsigned int counter = 1;

        // Loop until we have a solution with sufficient accuracy.
        do
        {
            // Save old values.
            previousRootValue = rootValue;
            previousRootFunctionValue = rootFunctionValue;

            // Check which subinterval to keep, by maintaining endpoints with opposite function
            // value signs.
            if( rootFunctionValue * currentLowerBoundFunctionValue < 0.0 )
            {
                // Different sign, hence the upper bound is replaced.
                currentUpperBound = rootValue;
                currentUpperBoundFunctionValue = rootFunctionValue;
            }
            else
            {
                // Same sign, hence the lower bound is replaced.
                currentLowerBound = rootValue;
                currentLowerBoundFunctionValue = rootFunctionValue;
            }

            // Compute the new midpoint of the interval and its function value.
            rootValue = ( currentLowerBound + currentUpperBound ) / 2.0;
            rootFunctionValue = this->rootFunction->evaluate( rootValue );

            // Sanity check.
            if( currentLowerBoundFunctionValue * currentUpperBoundFunctionValue > 0.0 )
            {
                boost::throw_exception( boost::enable_error_info( std::runtime_error(
                                                                      boost::str( boost::format(
                                                                                      "The Bisection algorithm requires that the values at the upper "
                                                                                      "and lower bounds have a different sign, error during iteration." ) ) ) ) );
            }

            counter++;
        }
        while( !this->terminationFunction( rootValue, previousRootValue, rootFunctionValue,
                                           previousRootFunctionValue, counter ) );

        return rootValue;

    }

    //! Reset the bracket of the solution.
    /*!
     * Resets the search interval with an upper and lower bound. It is required that the function
     * values at upper and lower bound have an opposite sign.
     *
     * \param lowerBound Lower bound of the interval containing a root.
     * \param upperBound Upper bound of the interval containing a root.
     */
    void resetBoundaries( const DataType lowerBound, const DataType upperBound )
    {
        this->lowerBound_ = lowerBound;
        this->upperBound_ = upperBound;
    }

protected:

private:

    //! Lower bound of the bracket containing the solution.
    double lowerBound_;

    //! Upper bound of the bracket containing the solution.
    double upperBound_;

};

// Some handy typedefs.
typedef BisectionCore< double > Bisection;
typedef boost::shared_ptr< Bisection > BisectionPointer;

} // namespace root_finders
} // namespace tudat

#endif // TUDAT_BISECTION_H
