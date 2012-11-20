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
 *      101116    E. Iorfida        Added set-get functions.
 *      101121    E. Iorfida        Added Doxygen comments.
 *      110111    E. Iorfida        Deleted useless lines, and modified punctuation.
 *      110111    K. Kumar          Renamed functions and variables to be more accurate and
 *                                  descriptive; added "End of file." comment; added example of
 *                                  pointer-to-member-function for LambertTargeter.
 *      110119    K. Kumar          Updated code to work with adaptor and  abstract base
 *                                  implementation so that pointer-to-member functions are not
 *                                  required.
 *      110124    E. Iorfida        Added set/get functions for maximum number of iterations.
 *      110810    J. Leloux         Corrected doxygen documentation.
 *      110905    S. Billemont      Reorganized includes.
 *                                  Moved (con/de)structors and getter/setters to header.
 *      120712    P. Musegaas       Changed absolute tolerance into a safe variant of relative
 *                                  tolerance.
 *      120208    S. Billemont      Move to new root_finders codebase.
 *      120402    T. Secretin       Code-check.
 *      120726    S. Billemont      Restructuring. Implemented new termination conditions.
 *      120810    P. Musegaas       Code-check, various edits.
 *
 *    References
 *
 *    Notes
 */

#ifndef TUDAT_ROOT_FINDER_H
#define TUDAT_ROOT_FINDER_H

#include <boost/function.hpp>
#include <boost/shared_ptr.hpp>

#include <TudatCore/Mathematics/BasicMathematics/mathematicalConstants.h>

#include "Tudat/Mathematics/BasicMathematics/function.h"
#include "Tudat/Mathematics/RootFinders/terminationConditions.h"

namespace tudat
{
namespace root_finders
{

//! Root-finder class.
/*!
 * A description of root-finding algorithms. These algorihms take a Function, and find a root
 * using passed along data, such as starting point or a solution bracket.
 * \tparam DataType Data type used to represent floating-point values.
 */
template< typename DataType = double >
class RootFinderCore
{
public:

    //! Definition of the function whose root we have to determine.
    typedef boost::shared_ptr< basic_mathematics::Function< DataType, DataType > > FunctionPointer;
    typedef boost::function< bool( DataType, DataType,
                                   DataType, DataType, unsigned int ) > TerminationFunction;
	
    //! Default constructor.
    /*!
     * \param rootFunction_ The function whose root is to be determined.
     * \param terminationFunction_ The function specifying the termination conditions of the
     *                                      root-finding process.
     */
    RootFinderCore( TerminationFunction terminationFunction_ )
        : terminationFunction( terminationFunction_ )
    { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~RootFinderCore( ) { }

    //! Get the function subject to the rootfinding algorithm.
    /*!
     * Returns the function subject to the rootfinding algorithm.
     * \return rootFunction Pointer to the function subject to the rootfinding algorithm.
     */
    const FunctionPointer getFunction( )
    {
        return rootFunction;
    }

    //! Find a root of the set function.
    /*!
     * Try to find the root of the current function, using a specified technique.
     * \param aRootFunction Function to find root of.
     * \param initialGuess The initial guess of the root.
     * \throws ConvergenceExeption If the solution does not converge to a root value.
     */
    virtual DataType execute( const FunctionPointer aRootFunction,
                              const DataType initialGuess ) = 0;

protected:

    //! Function to find the root of.
    /*!
     * The root-finder tries to find a root of this function.
     */
    FunctionPointer rootFunction;
	
    //! Function specifying the termination conditions.
    /*!
     * The rootfinder will continue improving the solution of the root until the termination
     * function returns true.
     */
    TerminationFunction terminationFunction;

private:
};

//! Typedef for a root-finder with double data type.
typedef RootFinderCore< > RootFinder;

//! Typedef for a shared-pointer to a root-finder with double data type.
typedef boost::shared_ptr< RootFinder > RootFinderPointer;

} // namespace root_finders
} // namespace tudat

#endif // TUDAT_ROOT_FINDER_H
