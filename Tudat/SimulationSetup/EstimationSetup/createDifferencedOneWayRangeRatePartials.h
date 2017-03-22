#ifndef CREATEDIFFERENCEDONEWAYRANGERATEPARTIALS_H
#define CREATEDIFFERENCEDONEWAYRANGERATEPARTIALS_H

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createDifferencedOneWayRangeRatePartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createOneWayRangePartials.h"

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/differencedOneWayRangeRatePartial.h"

namespace tudat
{

namespace observation_partials
{

typedef std::map< observation_models::LinkEnds, std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > > PerLinkEndPerLightTimeSolutionCorrections;



typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;
void removeObservationBiasPartialFromList( SingleLinkObservationPartialList& observationPartialList );


std::pair< PerLinkEndPerLightTimeSolutionCorrections, PerLinkEndPerLightTimeSolutionCorrections > splitOneWayRangeRateLightTimeCorrectionsBetweenArcs(
        const PerLinkEndPerLightTimeSolutionCorrections& combinedCorrections );

template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > createDifferencedOneWayRangeRatePartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const PerLinkEndPerLightTimeSolutionCorrections& lightTimeCorrections =
        PerLinkEndPerLightTimeSolutionCorrections( ) )
{
    using namespace observation_partials;

    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > rangeRatePartials;


    std::pair< PerLinkEndPerLightTimeSolutionCorrections, PerLinkEndPerLightTimeSolutionCorrections > splitLightTimeCorrections =
            splitOneWayRangeRateLightTimeCorrectionsBetweenArcs( lightTimeCorrections );

    // Create one way range partials for link at start of arc.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > arcStartPartials
            = createOneWayRangePartials( linkEnds, bodyMap, parametersToEstimate, splitLightTimeCorrections.first );

    // Create one way range partials for link at end of arc.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList , boost::shared_ptr< PositionPartialScaling > > > arcEndPartials
            = createOneWayRangePartials( linkEnds, bodyMap, parametersToEstimate, splitLightTimeCorrections.second );

    for( std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >::iterator
         linkEndIterator = arcStartPartials.begin( ); linkEndIterator != arcStartPartials.end( ); linkEndIterator++ )
    {
        removeObservationBiasPartialFromList( linkEndIterator->second.first );
    }

    for( std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >::iterator
         linkEndIterator = arcEndPartials.begin( ); linkEndIterator != arcEndPartials.end( ); linkEndIterator++ )
    {
        removeObservationBiasPartialFromList( linkEndIterator->second.first );
    }


    if( arcStartPartials.size( ) != arcEndPartials.size( ) )
    {
        std::cerr<<"Error when making differenced one way range rate partials, arc start and end partial set size is not consistent"<<std::endl;
    }

    // Create iterators over created partial links.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >::iterator
            arcStartIterator = arcStartPartials.begin( );
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >::iterator
            arcEndIterator = arcEndPartials.begin( );

    // Pre-create iterators over one-way range partials for the single sets of link ends at the start and end of the arc.
    SingleLinkObservationPartialList::iterator currentLinkEndArcStartIterator;
    SingleLinkObservationPartialList::iterator currentLinkEndArcEndIterator;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < arcStartPartials.size( ); i++ )
    {
        SingleLinkObservationPartialList currentRangeRatePartialList;

        if( arcEndIterator->first != arcStartIterator->first )
        {
            std::cerr<<"Error when making differenced one way range rate partials, arc start and end partial link ends not consistent"<<std::endl;
        }

        // Set iterators over set of partials for current link ends.
        currentLinkEndArcStartIterator = arcStartIterator->second.first.begin( );
        currentLinkEndArcEndIterator = arcEndIterator->second.first.begin( );

        // Iterate over all one-way range partials and create one-way range rate partial from them.
        for( unsigned int j = 0; j < arcStartIterator->second.first.size( ); j++ )
        {
            if( currentLinkEndArcStartIterator->second->getParameterIdentifier( ) !=
                    currentLinkEndArcEndIterator->second->getParameterIdentifier( ) )
            {
                std::cerr<<"Error when making differenced one way range partial, parameters did not match"<<std::endl;
            }
            else
            {
                // Create range rate partial.
                currentRangeRatePartialList[ currentLinkEndArcStartIterator->first ] =
                        boost::make_shared< DifferencedOneWayRangeRatePartial >(
                            currentLinkEndArcStartIterator->second->getParameterIdentifier( ),
                            currentLinkEndArcStartIterator->second,
                            currentLinkEndArcEndIterator->second );
            }

            // Increment range partial iterators.
            currentLinkEndArcStartIterator++;
            currentLinkEndArcEndIterator++;
        }

        // Set range partials for current link end.
        rangeRatePartials[ arcStartIterator->first ] = std::make_pair(
                    currentRangeRatePartialList, boost::make_shared< OneWayRangeRateScaling >(
                        boost::dynamic_pointer_cast< OneWayRangeScaling >( arcStartIterator->second.second ),
                        boost::dynamic_pointer_cast< OneWayRangeScaling >( arcEndIterator->second.second ) ) );

        // Increment iterators to next link ends.
        arcStartIterator++;
        arcEndIterator++;
    }

    // Return partials
    return rangeRatePartials;
}

}

}
#endif // CREATEDIFFERENCEDONEWAYRANGERATEPARTIALS_H
