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

#ifndef TUDAT_ROOT_FINDER_H
#define TUDAT_ROOT_FINDER_H

#include <functional>
#include <memory>

#include "tudat/math/basic/mathematicalConstants.h"

#include "tudat/math/basic/function.h"
#include "tudat/math/root_finders/terminationConditions.h"

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

    //! Typedef of the function whose root we have to determine.
    typedef std::shared_ptr< basic_mathematics::Function< DataType, DataType > > FunctionPointer;

    //! Typedef of the function determining whether to terminate the root finding (i.e convergence)
    typedef std::function< bool( DataType, DataType,
                                   DataType, DataType, unsigned int ) > TerminationFunction;
        
    //! Constructor taking custom termination function.
    /*!
     *  Constructor taking custom termination function.
     *  \param terminationFunction The function specifying the termination conditions of the
     *  root-finding process \sa RootFinderCore::terminationFunction_
     */
    RootFinderCore( TerminationFunction terminationFunction )
        : terminationFunction_( terminationFunction )
    { }

    RootFinderCore( std::shared_ptr< TerminationCondition< DataType > > terminationCondition )
        : terminationFunction_(
              std::bind( &TerminationCondition< DataType >::checkTerminationCondition, terminationCondition,
                         std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                         std::placeholders::_4, std::placeholders::_5 ) )
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

    //! Find a root of the function provided as input.
    /*!
     *  Find a root of the function provided as input, using the termination function set by the
     *  constructor and the technique specified in the derived class
     * \param aRootFunction Function to find root of.
     * \param initialGuess The initial guess of the root.
     * \throws ConvergenceExeption If the solution does not converge to a root value.
     * \return Root of the rootFunction that is found
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
     * function returns true. The five input variables of the TerminationFunction typedef are:
     * current root value; previous root value; current function value; previous function value;
     * number of iterations. Its output is true if the algorithm is to terminate, false otherwise.
     */
    TerminationFunction terminationFunction_;

private:

};

//! Typedef for a root-finder with double data type.
typedef RootFinderCore< > RootFinder;

//! Typedef for a shared-pointer to a root-finder with double data type.
typedef std::shared_ptr< RootFinder > RootFinderPointer;

} // namespace root_finders

} // namespace tudat

#endif // TUDAT_ROOT_FINDER_H
