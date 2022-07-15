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

#include "tudat/astro/orbit_determination/observation_partials/differencedObservationPartial.h"
#include "tudat/simulation/estimation_setup/createDirectObservationPartials.h"
#include "tudat/simulation/estimation_setup/createNWayRangePartials.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace observation_partials
{

template< int ObservationSize >
class DifferencedObservationPartialCreator
{
public:
    static std::shared_ptr< ObservationPartial< ObservationSize > > createDifferencedObservationPartial(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial,
            const std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial,
            const observation_models::LinkEnds& linkEnds );
};

template< >
class DifferencedObservationPartialCreator< 1 >
{
public:
    static std::shared_ptr< ObservationPartial< 1 > > createDifferencedObservationPartial(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< ObservationPartial< 1 > > firstPartial,
            const std::shared_ptr< ObservationPartial< 1 > > secondPartial,
            const observation_models::LinkEnds& linkEnds )
    {
        using namespace observation_models;

        std::shared_ptr< ObservationPartial< 1 > > differencedPartial;
        switch( differencedObservableType )
        {
        case one_way_differenced_range:
        {
            if( firstPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( firstPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating one-way differenced range partial; first input object type is incompatible" );
                }
                else if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( firstPartial )->getObservableType( ) != one_way_range )
                {
                    throw std::runtime_error( "Error when creating one-way differenced range partial; first input observable type is incompatible" );
                }
            }

            if( secondPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( secondPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating one-way differenced range partial; second input object type is incompatible" );
                }
                else if( std::dynamic_pointer_cast< DirectObservationPartial< 1 > >( secondPartial )->getObservableType( ) != one_way_range )
                {
                    throw std::runtime_error( "Error when creating one-way differenced range partial; second input observable type is incompatible" );
                }
            }
            differencedPartial = std::make_shared< DifferencedObservablePartial< 1 > >(
                        firstPartial, secondPartial, &observation_models::getDifferencedOneWayRangeScalingFactor,
                        getUndifferencedTimeAndStateIndices( one_way_differenced_range, linkEnds.size( ) ) );
            break;
        }
        case n_way_differenced_range:
        {
            if( firstPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< NWayRangePartial >( firstPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating n-way differenced range partial; first input object type is incompatible" );
                }
            }

            if( secondPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< NWayRangePartial >( secondPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating n-way differenced range partial; second input object type is incompatible" );
                }
            }

            differencedPartial = std::make_shared< DifferencedObservablePartial< 1 > >(
                        firstPartial, secondPartial, &observation_models::getDifferencedNWayRangeScalingFactor,
                        getUndifferencedTimeAndStateIndices( n_way_differenced_range, linkEnds.size( ) ) );
            break;
        }
        default:
            throw std::runtime_error( "Error when creating differenced observable partial (size 1); observable " + getObservableName( differencedObservableType ) +
                                      " is not differenced. " );

        }
        return differencedPartial;
    }
};

template< >
class DifferencedObservationPartialCreator< 2 >
{
public:
    static std::shared_ptr< ObservationPartial< 2 > > createDifferencedObservationPartial(
            const observation_models::ObservableType differencedObservableType,
            const std::shared_ptr< ObservationPartial< 2 > > firstPartial,
            const std::shared_ptr< ObservationPartial< 2 > > secondPartial,
            const observation_models::LinkEnds& linkEnds )
    {
        using namespace observation_models;

        std::shared_ptr< ObservationPartial< 2 > > differencedPartial;
        switch( differencedObservableType )
        {
        case relative_angular_position:
        {
            if( firstPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( firstPartial ) == nullptr )
                {
                    throw std::runtime_error( "Error when creating relative angular position partial; first input object type is incompatible" );
                }
                else if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( firstPartial )->getObservableType( ) != angular_position )
                {
                    throw std::runtime_error( "Error when creating relative angular position partial; first input observable type is incompatible" );
                }
            }

            if( secondPartial != nullptr )
            {
                if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( secondPartial ) == nullptr )
                {
                    std::cout<<secondPartial<<" "<<
                               std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( secondPartial )<<std::endl;
                    throw std::runtime_error( "Error when creating relative angular position partial; second input object type is incompatible" );
                }
                else if( std::dynamic_pointer_cast< DirectObservationPartial< 2 > >( secondPartial )->getObservableType( ) != angular_position )
                {
                    throw std::runtime_error( "Error when creating relative angular position partial; second input observable type is incompatible" );
                }
            }
            differencedPartial = std::make_shared< DifferencedObservablePartial< 2 > >(
                        firstPartial, secondPartial, [=]( const std::vector< double >&, const observation_models::LinkEndType ){ return 1.0; },
            getUndifferencedTimeAndStateIndices( relative_angular_position, linkEnds.size( ) ) );
            break;
        }
        default:
            throw std::runtime_error( "Error when creating differenced observable partial (size 2); observable " + getObservableName( differencedObservableType ) +
                                      " is not differenced. " );

        }
        return differencedPartial;
    }
};

template< int ObservationSize >
std::map< std::pair< int, int >, std::pair< std::shared_ptr< ObservationPartial< ObservationSize > >, std::shared_ptr< ObservationPartial< ObservationSize > > > >
mergeUndifferencedPartialContribution(
        const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& firstPartialList,
        const std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >& secondPartialList )
{
    std::map< std::pair< int, int >, std::pair< std::shared_ptr< ObservationPartial< ObservationSize > >, std::shared_ptr< ObservationPartial< ObservationSize > > > >
            mergedPartials;

    std::shared_ptr< ObservationPartial< ObservationSize > > firstPartial;
    std::shared_ptr< ObservationPartial< ObservationSize > > secondPartial;
    for( auto it : firstPartialList )
    {
        firstPartial = it.second;
        if( secondPartialList.count( it.first ) != 0 )
        {
            secondPartial = secondPartialList.at( it.first );
        }
        else
        {
            secondPartial = nullptr;
        }
        mergedPartials[ it.first ] = std::make_pair( firstPartial, secondPartial );
    }

    for( auto it : secondPartialList )
    {
        if( mergedPartials.count( it.first ) == 0 )
        {
            mergedPartials[ it.first ] = std::make_pair( nullptr, it.second );
        }
    }
    return mergedPartials;
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
    ObservableType differencedObservableType = observationModel->getObservableType( );
    ObservableType undifferencedObservableType = getUndifferencedObservableType(
                differencedObservableType );

    auto undifferencedObservationModels = UndifferencedObservationModelExtractor< ObservationSize >::extract( observationModel );

    std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType > > undifferencedObservationModelFirst =
            undifferencedObservationModels.first;
    std::shared_ptr< ObservationModel< ObservationSize, ParameterType, TimeType > > undifferencedObservationModelSecond =
            undifferencedObservationModels.second;

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > >
            firstUndifferencedObservablePartials =
            createSingleLinkObservationPartials< ParameterType, ObservationSize, TimeType >(
                undifferencedObservationModelFirst, bodies, parametersToEstimate, false );

    std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >, std::shared_ptr< PositionPartialScaling > >
            secondUndifferencedObservablePartials =
            createSingleLinkObservationPartials< ParameterType, ObservationSize, TimeType >(
                undifferencedObservationModelSecond, bodies, parametersToEstimate, false );

    std::map< std::pair< int, int >, std::pair< std::shared_ptr< ObservationPartial< ObservationSize > >, std::shared_ptr< ObservationPartial< ObservationSize > > > >
            mergedPartials = mergeUndifferencedPartialContribution(
                firstUndifferencedObservablePartials.first, secondUndifferencedObservablePartials.first );

    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > differencedObservationPartialList;

    // Iterate over all one-way range partials and create one-way range rate partial from them.
    for( auto it : mergedPartials )
    {
        // Create range rate partial.
        differencedObservationPartialList[ it.first] =
                DifferencedObservationPartialCreator< ObservationSize >::createDifferencedObservationPartial(
                    differencedObservableType,
                    it.second.first,
                    it.second.second,
                    linkEnds );
    }


    // Create bias partials
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =  parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        std::shared_ptr< ObservationPartial< ObservationSize > > currentDifferencedObservationPartial;
        if( isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first ) && useBiasPartials )
        {
            currentDifferencedObservationPartial = createObservationPartialWrtLinkProperty< ObservationSize >(
                        linkEnds, undifferencedObservableType, parameterIterator->second );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current doppler and current parameter)
        if( currentDifferencedObservationPartial != nullptr )
        {
            // Add partial to the list.
            std::pair< double, double > currentPair = std::pair< int, int >( parameterIterator->first,
                                                                             parameterIterator->second->getParameterSize( ) );
            differencedObservationPartialList[ currentPair ] = currentDifferencedObservationPartial;
        }
    }

    differencedPartialsAndScaling = std::make_pair(
                differencedObservationPartialList, ObservationPartialScalingCreator< ObservationSize >::
                template createDifferencedPositionPartialScalingObject< ParameterType, TimeType >(
                    differencedObservableType, firstUndifferencedObservablePartials.second,
                    secondUndifferencedObservablePartials.second, bodies ) );
    return differencedPartialsAndScaling;


}

}

}
#endif // TUDAT_CREATEDIFFERENCEDOBSERVABLEPARTIALS_H
