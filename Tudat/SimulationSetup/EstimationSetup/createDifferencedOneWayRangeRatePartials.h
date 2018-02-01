/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEDIFFERENCEDONEWAYRANGERATEPARTIALS_H
#define TUDAT_CREATEDIFFERENCEDONEWAYRANGERATEPARTIALS_H

#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createDifferencedOneWayRangeRatePartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createOneWayRangePartials.h"

#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/differencedOneWayRangeRatePartial.h"

namespace tudat
{

namespace observation_partials
{

//! Typedef for list of light time corrections for a list of link ends
typedef std::map< observation_models::LinkEnds,
std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > >
PerLinkEndPerLightTimeSolutionCorrections;

//! Typedef for single observation partial list, with key parameter start index and size
typedef std::map< std::pair< int, int >, boost::shared_ptr< ObservationPartial< 1 > > > SingleLinkObservationPartialList;

//! Function to split the total list of light-time corrections for one-way differenced range rate into list for either arc.
/*!
 * Function to split the total list of light-time corrections for one-way differenced range rate into list for either arc.
 * \param combinedCorrections Total list of light-time corrections for one-way differenced range rate, per link end. First
 * vector in map entry should always have size two: entry 0 is for arc start, entry 1 for arc end.
 * \return List of light-time corrections for two arcs of one-way differenced range rate.
 */
std::pair< PerLinkEndPerLightTimeSolutionCorrections, PerLinkEndPerLightTimeSolutionCorrections >
splitOneWayRangeRateLightTimeCorrectionsBetweenArcs(
        const PerLinkEndPerLightTimeSolutionCorrections& combinedCorrections );

//! Function to generate one-way differenced range partials for all parameters that are to be estimated, for all sets of link ends
/*!
 *  Function to generate one-way differenced range partials for all parameters that are to be estimated, for all sets of link ends
 *  The one-way differenced range partials are generated per set of link ends, using partials for the range of the arc start and
 *  arc end links The set of parameters and bodies that are to be estimated, as well as the set of link ends (each of which must
 *  contain a transmitter and receiever linkEndType) that are to be used.
 *  \param linkEnds Vector of all link ends for which  partials are to be calculated
 *  \param bodyMap List of all bodies that consitute the environment
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary one-way differenced range partials of a single
 *  link end, and scaling, object, used for scaling the position partial members of all partial in link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds,
std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >
createDifferencedOneWayRangeRatePartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::NamedBodyMap bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const PerLinkEndPerLightTimeSolutionCorrections& lightTimeCorrections =
        PerLinkEndPerLightTimeSolutionCorrections( ) )
{
    using namespace observation_partials;


    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > rangeRatePartials;

    // Retrieve separate light-time corrections for arc start and arc end.
    std::pair< PerLinkEndPerLightTimeSolutionCorrections, PerLinkEndPerLightTimeSolutionCorrections > splitLightTimeCorrections =
            splitOneWayRangeRateLightTimeCorrectionsBetweenArcs( lightTimeCorrections );

    // Create one way range partials for link at start of arc.
    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > arcStartPartials
            = createOneWayRangePartials( linkEnds, bodyMap, parametersToEstimate, splitLightTimeCorrections.first );

    // Create one way range partials for link at end of arc.
    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList , boost::shared_ptr< PositionPartialScaling > > > arcEndPartials
            = createOneWayRangePartials( linkEnds, bodyMap, parametersToEstimate, splitLightTimeCorrections.second );

    // Check output consistency
    if( arcStartPartials.size( ) != arcEndPartials.size( ) )
    {
        throw std::runtime_error(
                    "Error when making differenced one way range rate partials, arc start and end partial set size is not consistent" );
    }

    // Create iterators over created partial links.
    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >::iterator
            arcStartIterator = arcStartPartials.begin( );
    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >::iterator
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
            throw std::runtime_error(
                        "Error when making differenced one way range rate partials, arc start and end partial link ends not consistent" );
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
                throw std::runtime_error(
                            "Error when making differenced one way range partial, parameters did not match" );
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
#endif // TUDAT_CREATEDIFFERENCEDONEWAYRANGERATEPARTIALS_H
