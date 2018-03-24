/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/Mathematics/RootFinders/createRootFinder.h"

namespace tudat
{

namespace root_finders
{

//! Function to determine whether a root finder requires any analytical derivatives
bool doesRootFinderRequireDerivatives( const boost::shared_ptr< RootFinderSettings > rootFinderSettings )
{
    bool rootFinderRequireDerivatives = -1;
    switch( rootFinderSettings->rootFinderType_ )
    {
    case bisection_root_finder:
        rootFinderRequireDerivatives = false;
        break;
    case halley_root_finder:
        rootFinderRequireDerivatives = true;
        break;
    case newton_raphson_root_finder:
        rootFinderRequireDerivatives = true;
        break;
    case secant_root_finder:
        rootFinderRequireDerivatives = false;
        break;
    default:
        throw std::runtime_error( "Error when getting root finder derivative needs, root finder type not found" );
    }
    return rootFinderRequireDerivatives;
}


}

}
