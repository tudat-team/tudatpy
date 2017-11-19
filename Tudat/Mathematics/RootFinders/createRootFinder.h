/*    Copyright (c) 2010-2017, Delft University of Technology
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

enum RootFinderType
{
    bisection_root_finder,
    halley_root_finder,
    newton_raphson_root_finder,
    secant_root_finder
};

class RootFinderSettings
{
public:
    RootFinderSettings( const RootFinderType rootFinderType,
                        const double terminationTolerance, const unsigned int maximumNumberOfIterations  ):
        rootFinderType_( rootFinderType ), terminationTolerance_( terminationTolerance ),
        maximumNumberOfIterations_( maximumNumberOfIterations ){ }

    ~RootFinderSettings( ){ }

    RootFinderType rootFinderType_;

    double terminationTolerance_;

    unsigned int maximumNumberOfIterations_;

};

template< typename DataType = double >
boost::shared_ptr< RootFinderCore< DataType > > createRootFinder(
        const boost::shared_ptr< RootFinderSettings > rootFinderSettings,
        const DataType lowerBound = TUDAT_NAN, const DataType upperBound = TUDAT_NAN,
        const DataType previousGuess = TUDAT_NAN )
{
    boost::shared_ptr< RootFinderCore< DataType > > rootFinder;
    switch( rootFinderSettings->rootFinderType_ )
    {
    case bisection_root_finder:
        if( !( lowerBound == lowerBound ) || !( upperBound == upperBound ) )
        {
            throw std::runtime_error( "Error when making bisection root finder, lower/upped bound not provided" );
        }
        rootFinder = boost::make_shared< BisectionCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_,
                    lowerBound, upperBound );
        break;
    case halley_root_finder:
        rootFinder = boost::make_shared< HalleyRootFinderCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_ );
        break;
    case newton_raphson_root_finder:

        rootFinder = boost::make_shared< NewtonRaphsonCore< DataType > >(
                    rootFinderSettings->terminationTolerance_, rootFinderSettings->maximumNumberOfIterations_ );
        break;
    case secant_root_finder:
        if( !( previousGuess == previousGuess ) )
        {
            throw std::runtime_error( "Error when making secant root finder, initial guess not provided" );
        }
        rootFinder = boost::make_shared< SecantRootFinderCore< DataType > >(
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
