#include "Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h"

namespace tudat
{

namespace propagators
{

boost::shared_ptr< ephemerides::ReferenceFrameManager > createFrameManager(
        const simulation_setup::NamedBodyMap& bodyMap )
{
    // Get ephemerides from bodies
    std::map< std::string, boost::shared_ptr< ephemerides::Ephemeris > > ephemerides;
    for( simulation_setup::NamedBodyMap::const_iterator bodyIterator = bodyMap.begin( );
         bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        if( bodyIterator->second->getEphemeris( ) != NULL )
        {
            ephemerides[ bodyIterator->first ] = bodyIterator->second->getEphemeris( );
        }
    }
    return boost::make_shared< ephemerides::ReferenceFrameManager >(
                ephemerides );
}

}

}

