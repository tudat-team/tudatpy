/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEDIFFERENCEDOBSERVABLEPARTIALS_H
#define TUDAT_CREATEDIFFERENCEDOBSERVABLEPARTIALS_H

#include "tudat/astro/orbit_determination/observation_partials/differencedOneWayRangeRatePartial.h"
#include "tudat/simulation/estimation_setup/createDirectObservationPartials.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace observation_partials
{


template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > >
createDifferencedObservationPartial(
        const observation_models::ObservableType undifferencedObservableType,
        const std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial,
        const std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial )
{
    using namespace observation_models;

    std::shared_ptr< ObservationPartial< ObservationSize > > differencedPartial;
    switch( undifferencedObservableType )
    {
    case one_way_differenced_range:
    {
        if( std::dynamic_pointer_cast< DirectObservationPartial< ObservationSize > >( firstPartial ) == nullptr )
        {
            throw std::runtime_error( "Error when creating one-way differenced range partial; first input object type is incompatible" );
        }
        else if( std::dynamic_pointer_cast< DirectObservationPartial< ObservationSize > >( firstPartial )->getObservableType( ) != one_way_range )
        {
            throw std::runtime_error( "Error when creating one-way differenced range partial; first input observable type is incompatible" );
        }

        if( std::dynamic_pointer_cast< DirectObservationPartial< ObservationSize > >( secondPartial ) == nullptr )
        {
            throw std::runtime_error( "Error when creating one-way differenced range partial; second input object type is incompatible" );
        }
        else if( std::dynamic_pointer_cast< DirectObservationPartial< ObservationSize > >( secondPartial )->getObservableType( ) != one_way_range )
        {
            throw std::runtime_error( "Error when creating one-way differenced range partial; second input observable type is incompatible" );
        }
        differencedPartial = std::make_shared< DifferencedObservablePartial< ObservationSize > >(
                    firstPartial, secondPartial, &observation_models::getDifferencedOneWayRangeScalingFactor,
                    getUndifferencedTimeAndStateIndices( one_way_differenced_range ) );
    }
    default:
        throw std::runtime_error( "Error when creating differenced observable partial; observable " + getObservableName( undifferencedObservableType ) +
                                  " is not differenced. " );

    }
    return differencedPartial;
}



template< typename ParameterType, typename TimeType, int ObservationSize >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > createDifferencedObservablePartials(
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > observationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool useBiasPartials = true )
{
    using namespace observation_models;

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
    std::shared_ptr< PositionPartialScaling > > differencedPartialsAndScaling;

    LinkEnds linkEnds = observationModel->getLinkEnds( );
    ObservableType undifferencedObservableType = getUndifferencedObservableType(
                observationModel->getObservableType( ) );

    auto undifferencedObservationModels = getUndifferencedObservationModels( observationModel );

    std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType > > undifferencedObservationModelFirst =
            undifferencedObservationModels.first;
    std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType > > undifferencedObservationModelSecond =
            undifferencedObservationModels.second;

    std::map< LinkEnds, std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType> > > inputList;

    inputList[ undifferencedObservationModelFirst->getLinkEnds( ) ] = undifferencedObservationModelFirst;
    std::vector< std::shared_ptr< LightTimeCorrection > > lightTimeCorrectionsFirst =
            getLightTimeCorrectionsList( inputList ).begin( )->second.at( 0 );
    inputList.clear( );

    inputList[ undifferencedObservationModelSecond->getLinkEnds( ) ] = undifferencedObservationModelSecond;
    std::vector< std::shared_ptr< LightTimeCorrection > > lightTimeCorrectionsSecond =
            getLightTimeCorrectionsList( inputList ).begin( )->second.at( 0 );

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > >
            firstUndifferencedObservablePartials =
    createSingleLinkObservationPartials< ParameterType, ObservationSize, TimeType >(
            linkEnds, undifferencedObservableType, bodies, parametersToEstimate, lightTimeCorrectionsFirst, false, undifferencedObservationModelFirst );

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > >
            secondUndifferencedObservablePartials =
    createSingleLinkObservationPartials< ParameterType, ObservationSize, TimeType >(
            linkEnds, undifferencedObservableType, bodies, parametersToEstimate, lightTimeCorrectionsSecond, false, undifferencedObservationModelSecond );

    if( firstUndifferencedObservablePartials.first.size( ) != secondUndifferencedObservablePartials.first.size( ) )
    {
        throw std::runtime_error( "Error when making differenced observation partials of type " + getObservableName( observationModel->getObservableType( ) ) +
                                  ", size of constituent partials is not equal" );
    }

    int numberOfPartials = firstUndifferencedObservablePartials.first.size( );

    auto firstPartialsIterator = firstUndifferencedObservablePartials.first.begin( );
    auto secondPartialsIterator = secondUndifferencedObservablePartials.first.begin( );

    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > differencedObservationPartialList;

    // Iterate over all one-way range partials and create one-way range rate partial from them.
    for( int j = 0; j < numberOfPartials; j++ )
    {
        if( firstPartialsIterator->first !=
                secondPartialsIterator->first )
        {
            throw std::runtime_error(
                        "Error when making differenced observation partials of type " + getObservableName( observationModel->getObservableType( ) ) + " parameter indices did not match" );
        }
        if( firstPartialsIterator->second->getParameterIdentifier( ) !=
                secondPartialsIterator->second->getParameterIdentifier( ) )
        {
            throw std::runtime_error(
                        "Error when making differenced observation partials of type " + getObservableName( observationModel->getObservableType( ) ) + " parameters did not match" );
        }
        else
        {
            // Create range rate partial.
            differencedObservationPartialList[ firstPartialsIterator->first ] =
                    createDifferencedObservationPartial(
                            undifferencedObservableType,
                            firstPartialsIterator->second,
                            secondPartialsIterator->second );

//                    std::make_shared< DifferencedOneWayRangeRatePartial >(
//                        firstPartialsIterator->second->getParameterIdentifier( ),
//                        firstPartialsIterator->second,
//                        secondPartialsIterator->second );
        }

        // Increment range partial iterators.
        firstPartialsIterator++;
        secondPartialsIterator++;
    }

    differencedPartialsAndScaling = std::make_pair(
                differencedObservationPartialList, ObservationPartialScalingCreator< ObservationSize >::
                template createDifferencedPositionPartialScalingObject< ParameterType, TimeType >(
                    observationModel->getObservableType( ), firstUndifferencedObservablePartials.second, secondUndifferencedObservablePartials.second, bodies ) );
    return differencedPartialsAndScaling;


}

template< typename ParameterType, typename TimeType, int ObservationSize >
std::map< observation_models::LinkEnds,
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > > createDifferencedObservablePartialsList(
        const std::map< observation_models::LinkEnds,
        std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > > observationModelList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool useBiasPartials = true )
{
    std::map< observation_models::LinkEnds,
    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
    std::shared_ptr< PositionPartialScaling > > > partialsList;
    for( auto it : observationModelList )
    {
        partialsList[ it.first ] = createDifferencedObservablePartials(
                    it.second, bodies, parametersToEstimate, useBiasPartials );
    }

    return partialsList;
}

}

}
#endif // TUDAT_CREATEDIFFERENCEDOBSERVABLEPARTIALS_H
