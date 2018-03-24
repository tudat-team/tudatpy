/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

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

