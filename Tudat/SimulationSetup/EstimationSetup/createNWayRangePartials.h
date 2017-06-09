#ifndef CREATENWAYRANGEPARTIALS_H
#define CREATENWAYRANGEPARTIALS_H

#include <vector>
#include <map>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/SimulationSetup/EstimationSetup/createOneWayRangePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/nWayRangePartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_partials
{

template< typename ParameterType >
std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > createNWayRangePartials(
        const observation_models::LinkEnds& nWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > lightTimeCorrections =
        std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > >( ) )

{
    SingleLinkObservationPartialList nWayRangePartialList;

    observation_models::LinkEnds currentLinkEnds;
    int numberOfLinkEnds = nWayRangeLinkEnds.size( );

    typedef std::map< int, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >
            OneWayRangePartialList;
    OneWayRangePartialList constituentOneWayRangePartials;

    std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > currentLightTimeCorrections;

    for( int i = 0; i < numberOfLinkEnds - 1; i++ )
    {
        currentLightTimeCorrections.clear( );
        if( lightTimeCorrections.size( ) > 0 )
        {
            currentLightTimeCorrections = lightTimeCorrections.at( i );
        }

        currentLinkEnds.clear( );
        currentLinkEnds[ observation_models::transmitter ] = nWayRangeLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i, numberOfLinkEnds ) );
        currentLinkEnds[ observation_models::receiver ] = nWayRangeLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i + 1, numberOfLinkEnds ) );

        constituentOneWayRangePartials[ i ] =
                createOneWayRangePartials( currentLinkEnds, bodyMap, parametersToEstimate, currentLightTimeCorrections );
    }

    std::map< int, boost::shared_ptr< OneWayRangeScaling > > oneWayRangeScalers;
    std::map< std::pair< int, int >, std::map< int, boost::shared_ptr< ObservationPartial< 1 > > > > sortedOneWayRangePartials;
    std::map< std::pair< int, int >, estimatable_parameters::EstimatebleParameterIdentifier > parameterIdList;

    for( OneWayRangePartialList::iterator oneWayPartialIterator = constituentOneWayRangePartials.begin( );
         oneWayPartialIterator != constituentOneWayRangePartials.end( ); oneWayPartialIterator++ )
    {
        oneWayRangeScalers[ oneWayPartialIterator->first ] = boost::dynamic_pointer_cast< OneWayRangeScaling >
                ( oneWayPartialIterator->second.second );

        for( SingleLinkObservationPartialList::iterator parameterIterator = oneWayPartialIterator->second.first.begin( );
             parameterIterator != oneWayPartialIterator->second.first.end( ); parameterIterator++ )
        {
            sortedOneWayRangePartials[ parameterIterator->first ][ oneWayPartialIterator->first ] = parameterIterator->second;

            if( parameterIdList.count( parameterIterator->first ) == 0 )
            {

                parameterIdList[ parameterIterator->first ] = parameterIterator->second->getParameterIdentifier( );

            }
            else if( parameterIdList.at( parameterIterator->first ) != parameterIterator->second->getParameterIdentifier( ) )
            {
                std::cerr<<"Error when making n way range partial, parameter indices are inconsistent"<<std::endl;
            }
        }
    }

    boost::shared_ptr< NWayRangeScaling > nWayRangeScaler = boost::make_shared< NWayRangeScaling >( oneWayRangeScalers );

    for( std::map< std::pair< int, int >, std::map< int, boost::shared_ptr< ObservationPartial< 1 > > > >::iterator sortedPartialIterator =
         sortedOneWayRangePartials.begin( ); sortedPartialIterator != sortedOneWayRangePartials.end( ); sortedPartialIterator++ )
    {
        nWayRangePartialList[ sortedPartialIterator->first ] = boost::make_shared< NWayRangePartial >(
                    nWayRangeScaler, sortedPartialIterator->second, parameterIdList.at( sortedPartialIterator->first ),
                    numberOfLinkEnds );
    }

    return std::make_pair( nWayRangePartialList, nWayRangeScaler );
}

template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >
createNWayRangePartials(
        const std::vector< observation_models::LinkEnds >& linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds, std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds, std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > partialMap;

    std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            currentLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) );
            if( currentLightTimeCorrections.size( ) != linkEnds.at( i ).size( ) - 1 )
            {
                std::cerr<<"Error when making n-way range partials, found light time correction partials for "<<currentLightTimeCorrections.size( )
                           <<" links, with "<<linkEnds.size( )<<" link ends"<<std::endl;
            }
        }
        else
        {
            currentLightTimeCorrections.clear( );
        }
        partialMap[ linkEnds.at( i ) ] = createNWayRangePartials< ParameterType >(
                    linkEnds.at( i ), bodyMap, parametersToEstimate, currentLightTimeCorrections );

    }

    return partialMap;
}

}

}

#endif // CREATENWAYRANGEPARTIALS_H
