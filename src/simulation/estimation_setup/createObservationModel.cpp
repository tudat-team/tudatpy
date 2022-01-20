/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>

#include <functional>
#include <boost/make_shared.hpp>


#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace observation_models
{


std::vector< LinkEnds > getObservationModelListLinkEnds(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelList )
{
    std::vector< LinkEnds > linkEndsList;
    for( unsigned int i = 0; i < observationModelList.size( ); i++ )
    {
        linkEndsList.push_back( observationModelList.at( i )->linkEnds_ );
    }
    return linkEndsList;
}

std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortObservationModelSettingsByType(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings )
{
    std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortedObservationModelSettings;
    for( unsigned int i = 0; i < observationModelSettings.size( ); i++ )
    {
        sortedObservationModelSettings[ observationModelSettings.at( i )->observableType_ ].push_back(
                observationModelSettings.at( i ) );
    }
    return sortedObservationModelSettings;
}



} // namespace observation_models

} // namespace tudat

