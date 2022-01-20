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

#ifndef TUDAT_CREATE_ROOT_FINDERS_H
#define TUDAT_CREATE_ROOT_FINDERS_H

#include <tudat/math/root_finders/rootFinder.h>
#include <tudat/math/root_finders/bisection.h>
#include <tudat/math/root_finders/halleyRootFinder.h>
#include <tudat/math/root_finders/newtonRaphson.h>
#include <tudat/math/root_finders/secantRootFinder.h>
#include <tudat/math/root_finders/terminationConditions.h>

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
                        const double relativeIndependentVariableTolerance,
                        const double absoluteIndependentVariableTolerance,
                        const double rootFunctionTolerance,
                        const unsigned int maximumNumberOfIterations,
                        const MaximumIterationHandling maximumIterationHandling = throw_exception ):
        rootFinderType_( rootFinderType ),
        relativeIndependentVariableTolerance_( relativeIndependentVariableTolerance ),
        absoluteIndependentVariableTolerance_( absoluteIndependentVariableTolerance ),
        rootFunctionTolerance_( rootFunctionTolerance ),
        maximumNumberOfIterations_( maximumNumberOfIterations ),
        maximumIterationHandling_( maximumIterationHandling )
    { }

    //! Destructor
    ~RootFinderSettings( ){ }

    //! Type of root finder to be used
    RootFinderType rootFinderType_;

    double relativeIndependentVariableTolerance_;

    double absoluteIndependentVariableTolerance_;

    double rootFunctionTolerance_;

    unsigned int maximumNumberOfIterations_;

    MaximumIterationHandling maximumIterationHandling_;

};

inline std::shared_ptr< RootFinderSettings > bisectionRootFinderSettings(
        const double relativeIndependentVariableTolerance = TUDAT_NAN,
        const double absoluteIndependentVariableTolerance = TUDAT_NAN,
        const double rootFunctionTolerance = TUDAT_NAN,
        const unsigned int maximumNumberOfIterations = 1000,
        const MaximumIterationHandling maximumIterationHandling = throw_exception )
{
    return std::make_shared< RootFinderSettings >(
                bisection_root_finder, relativeIndependentVariableTolerance, absoluteIndependentVariableTolerance,
                rootFunctionTolerance, maximumNumberOfIterations, maximumIterationHandling );
}

inline std::shared_ptr< RootFinderSettings > newtonRaphsonRootFinderSettings(
        const double relativeIndependentVariableTolerance = TUDAT_NAN,
        const double absoluteIndependentVariableTolerance = TUDAT_NAN,
        const double rootFunctionTolerance = TUDAT_NAN,
        const unsigned int maximumNumberOfIterations = 1000,
        const MaximumIterationHandling maximumIterationHandling = throw_exception )
{
    return std::make_shared< RootFinderSettings >(
                newton_raphson_root_finder, relativeIndependentVariableTolerance, absoluteIndependentVariableTolerance,
                rootFunctionTolerance, maximumNumberOfIterations, maximumIterationHandling );
}

inline std::shared_ptr< RootFinderSettings > halleyRootFinderSettings(
        const double relativeIndependentVariableTolerance = TUDAT_NAN,
        const double absoluteIndependentVariableTolerance = TUDAT_NAN,
        const double rootFunctionTolerance = TUDAT_NAN,
        const unsigned int maximumNumberOfIterations = 1000,
        const MaximumIterationHandling maximumIterationHandling = throw_exception )
{
    return std::make_shared< RootFinderSettings >(
                halley_root_finder, relativeIndependentVariableTolerance, absoluteIndependentVariableTolerance,
                rootFunctionTolerance, maximumNumberOfIterations, maximumIterationHandling );
}

inline std::shared_ptr< RootFinderSettings > secantRootFinderSettings(
        const double relativeIndependentVariableTolerance = TUDAT_NAN,
        const double absoluteIndependentVariableTolerance = TUDAT_NAN,
        const double rootFunctionTolerance = TUDAT_NAN,
        const unsigned int maximumNumberOfIterations = 1000,
        const MaximumIterationHandling maximumIterationHandling = throw_exception )
{
    return std::make_shared< RootFinderSettings >(
                secant_root_finder, relativeIndependentVariableTolerance, absoluteIndependentVariableTolerance,
                rootFunctionTolerance, maximumNumberOfIterations, maximumIterationHandling );
}


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
std::shared_ptr< RootFinder< DataType > > createRootFinder(
        const std::shared_ptr< RootFinderSettings > rootFinderSettings,
        const DataType lowerBound = TUDAT_NAN, const DataType upperBound = TUDAT_NAN,
        const DataType previousGuess = TUDAT_NAN )
{
    std::shared_ptr< TerminationCondition< DataType > > terminationCondition =
            createTerminationCondition< DataType >(
                rootFinderSettings->relativeIndependentVariableTolerance_,
                rootFinderSettings->absoluteIndependentVariableTolerance_,
                rootFinderSettings->rootFunctionTolerance_,
                rootFinderSettings->maximumNumberOfIterations_,
                rootFinderSettings->maximumIterationHandling_ );

    auto terminationFunction =
            std::bind( &TerminationCondition< DataType >::checkTerminationCondition, terminationCondition,
                       std::placeholders::_1, std::placeholders::_2, std::placeholders::_3,
                       std::placeholders::_4, std::placeholders::_5 );


    std::shared_ptr< RootFinder< DataType > > rootFinder;
    switch( rootFinderSettings->rootFinderType_ )
    {
    case bisection_root_finder:
        if( !( lowerBound == lowerBound ) || !( upperBound == upperBound ) )
        {
            throw std::runtime_error( "Error when making bisection root finder, lower/upped bound not provided" );
        }
        rootFinder = std::make_shared< Bisection< DataType > >(
                    terminationFunction, lowerBound, upperBound );
        break;
    case halley_root_finder:
        rootFinder = std::make_shared< HalleyRootFinder< DataType > >( terminationFunction );
        break;
    case newton_raphson_root_finder:
        rootFinder = std::make_shared< NewtonRaphson< DataType > >( terminationFunction );
        break;
    case secant_root_finder:
        if( !( previousGuess == previousGuess ) )
        {
            throw std::runtime_error( "Error when making secant root finder, initial guess not provided" );
        }
        rootFinder = std::make_shared< SecantRootFinder< DataType > >(
                    terminationFunction,
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
