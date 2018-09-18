/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_CREATE_ROOT_FINDERS_H
#define TUDAT_CREATE_ROOT_FINDERS_H

#include <Tudat/Mathematics/RootFinders/rootFinder.h>
#include <Tudat/Mathematics/RootFinders/bisection.h>
#include <Tudat/Mathematics/RootFinders/halleyRootFinder.h>
#include <Tudat/Mathematics/RootFinders/newtonRaphson.h>
#include <Tudat/Mathematics/RootFinders/secantRootFinder.h>

namespace tudat
{

namespace root_finders
{

//! Enum defininf types of root finder
enum RootFinderType
{
    bisection_root_finder,
    halley_root_finder,
    newton_raphson_root_finder,
    secant_root_finder
};

//! Class to define settings for a root finder
class RootFinderSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param rootFinderType Type of root finder to be used
     * \param terminationTolerance Convergence tolerance to be used for independent variable
     * \param maximumNumberOfIterations Maximum number of iterations to be used by root finder.
     */
    RootFinderSettings( const RootFinderType rootFinderType,
                        const double terminationTolerance, const unsigned int maximumNumberOfIterations  ):
        rootFinderType_( rootFinderType ), terminationTolerance_( terminationTolerance ),
        maximumNumberOfIterations_( maximumNumberOfIterations ){ }

    //! Destructor
    ~RootFinderSettings( ){ }

    //! Type of root finder to be used
    RootFinderType rootFinderType_;

    //! Convergence tolerance to be used for independent variable
    double terminationTolerance_;

    //! Maximum number of iterations to be used by root finder.
    unsigned int maximumNumberOfIterations_;

};

//! Function to determine whether a root finder requires any analytical derivatives
/*!
 * Function to determine whether a root finder requires any analytical derivatives, based on settings for root finder
 * \param rootFinderSettings Settings to be used to create a root finder
 * \return True if any analytical derivatives are required, false otherwise
 */
bool doesRootFinderRequireDerivatives( const std::shared_ptr< RootFinderSettings > rootFinderSettings );

//! Function to create a root finder
/*!
 * Function to create a root finde
 * \param rootFinderSettings Settings to be used to create a root finder
 * \param lowerBound Lower bound of search space of independent variables (default NaN; not needed by all root finders).
 * \param upperBound Lower bound of search space of independent variables (default NaN; not needed by all root finders).
 * \param previousGuess Initial guess for root position (default NaN; not needed by all root finders).
 */
template< typename DataType = double >
std::shared_ptr< RootFinderCore< DataType > > createRootFinder(
        const std::shared_ptr< RootFinderSettings > rootFinderSettings,
        const DataType lowerBound = TUDAT_NAN, const DataType upperBound = TUDAT_NAN,
        const DataType previousGuess = TUDAT_NAN )
{
    std::shared_ptr< RootFinderCore< DataType > > rootFinder;
    switch( rootFinderSettings->rootFinderType_ )
    {
    case bisection_root_finder:
        if( !( lowerBound == lowerBound ) || !( upperBound == upperBound ) )
        {
            throw std::runtime_error( "Error when making bisection root finder, lower/upped bound not provided" );
        }
        rootFinder = std::make_shared< BisectionCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_,
                    lowerBound, upperBound );
        break;
    case halley_root_finder:
        rootFinder = std::make_shared< HalleyRootFinderCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_ );
        break;
    case newton_raphson_root_finder:
        rootFinder = std::make_shared< NewtonRaphsonCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_ );
        break;
    case secant_root_finder:
        if( !( previousGuess == previousGuess ) )
        {
            throw std::runtime_error( "Error when making secant root finder, initial guess not provided" );
        }
        rootFinder = std::make_shared< SecantRootFinderCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_,
                    previousGuess );
        break;
    default:
        throw std::runtime_error( "Error when creating root finder, did not recognize root finder type" );
    }
    return rootFinder;
}

} // namespace root_finders

} // namespace tudat

#endif // TUDAT_CREATE_ROOT_FINDERS_H
