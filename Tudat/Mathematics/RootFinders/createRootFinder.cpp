
#include "Tudat/Mathematics/RootFinders/createRootFinder.h"

namespace tudat
{

namespace root_finders
{

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
