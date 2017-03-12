/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>


#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"

namespace tudat
{

namespace observation_models
{

//! Function to create list of observation models sorted by observable type and link ends from list only sorted in link ends.
SortedObservationSettingsMap convertUnsortedToSortedObservationSettingsMap(
        const ObservationSettingsMap& unsortedObservationSettingsMap )
{
    SortedObservationSettingsMap sortedObservationSettingsMap;

    for( ObservationSettingsMap::const_iterator iterator = unsortedObservationSettingsMap.begin( );
         iterator != unsortedObservationSettingsMap.end( ); iterator++ )
    {
        sortedObservationSettingsMap[ iterator->second->observableType_ ][ iterator->first ] =
                iterator->second;
    }
    return sortedObservationSettingsMap;
}

} // namespace observation_models

} // namespace tudat

