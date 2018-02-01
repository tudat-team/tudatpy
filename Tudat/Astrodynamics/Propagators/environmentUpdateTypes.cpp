/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <algorithm>
#include "Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h"

namespace tudat
{

namespace propagators
{

//! Function to extend existing list of required environment update types
void addEnvironmentUpdates( std::map< propagators::EnvironmentModelsToUpdate,
                            std::vector< std::string > >& environmentUpdateList,
                            const std::map< propagators::EnvironmentModelsToUpdate,
                            std::vector< std::string > > updatesToAdd )
{
    // Iterate over all environment update types.
    for( std::map< propagators::EnvironmentModelsToUpdate,
             std::vector< std::string > >::const_iterator
         environmentUpdateIterator = updatesToAdd.begin( );
         environmentUpdateIterator != updatesToAdd.end( ); environmentUpdateIterator++ )
    {
        bool addCurrentUpdate = 0;

        // Iterate over all updated bodies.
        for( unsigned int i = 0; i < environmentUpdateIterator->second.size( ); i++ )
        {
            addCurrentUpdate = 0;

            // Check if current update type exists.
            if( environmentUpdateList.count( environmentUpdateIterator->first ) == 0 )
            {
                addCurrentUpdate = 1;
            }
            // Check if current body exists for update type.
            else if( std::find( environmentUpdateList.at( environmentUpdateIterator->first ).begin( ),
                                environmentUpdateList.at( environmentUpdateIterator->first ).end( ),
                                environmentUpdateIterator->second.at( i ) ) ==
                     environmentUpdateList.at( environmentUpdateIterator->first ).end( ) )
            {
                addCurrentUpdate = 1;
            }

            // Add update type if required.
            if( addCurrentUpdate )
            {
                environmentUpdateList[ environmentUpdateIterator->first ].push_back(
                            environmentUpdateIterator->second.at( i ) );
            }
        }
    }
}


}

}

