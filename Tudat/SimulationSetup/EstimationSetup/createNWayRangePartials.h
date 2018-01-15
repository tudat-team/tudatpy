/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATENWAYRANGEPARTIALS_H
#define TUDAT_CREATENWAYRANGEPARTIALS_H

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

//! Function to generate n-way range partials and associated scaler for single link ends.
/*!
 *  Function to generate n-way range partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends set.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  The n-way range partials are built from one-way range partials of the constituent links
 *  \param nWayRangeLinkEnds Link ends (transmitter and receiever) for which n-way range partials are to be calculated
 *  (i.e. for which n-way range observations are to be processed).
 *  \param bodyMap List of all bodies, for creating n-way range partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default). First vector entry is
 *  index of link in n-way link ends, second vector is list of light-time corrections.
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary n-way range partials of a single link end, and NWayRangeScaling, object, used for
 *  scaling the position partial members of all NWayRangePartials in link end.
 */
template< typename ParameterType >
std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > createNWayRangePartials(
        const observation_models::LinkEnds& nWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > lightTimeCorrections =
        std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > >( ) )

{
    // Define return partial list
    SingleLinkObservationPartialList nWayRangePartialList;

    // Define list of constituent one-way partials.
    typedef std::map< int, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >
            OneWayRangePartialList;
    OneWayRangePartialList constituentOneWayRangePartials;

    // Getr number of link ends,
    observation_models::LinkEnds currentLinkEnds;
    int numberOfLinkEnds = nWayRangeLinkEnds.size( );

    // Iterate over all links in the n-way range observable
    std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > currentLightTimeCorrections;
    for( int i = 0; i < numberOfLinkEnds - 1; i++ )
    {
        currentLightTimeCorrections.clear( );
        if( lightTimeCorrections.size( ) > 0 )
        {
            currentLightTimeCorrections = lightTimeCorrections.at( i );
        }

        // Define links for current one-way range link
        currentLinkEnds.clear( );
        currentLinkEnds[ observation_models::transmitter ] = nWayRangeLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i, numberOfLinkEnds ) );
        currentLinkEnds[ observation_models::receiver ] = nWayRangeLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i + 1, numberOfLinkEnds ) );

        // Create onw-way range partials for current link
        constituentOneWayRangePartials[ i ] =
                createOneWayRangePartials( currentLinkEnds, bodyMap, parametersToEstimate, currentLightTimeCorrections, false );
    }

    // Retrieve sorted (by parameter index and link index) one-way range partials and (by link index) opne-way range partials
    std::map< int, boost::shared_ptr< OneWayRangeScaling > > oneWayRangeScalers;
    std::map< std::pair< int, int >, std::map< int, boost::shared_ptr< ObservationPartial< 1 > > > > sortedOneWayRangePartials;
    std::map< std::pair< int, int >, estimatable_parameters::EstimatebleParameterIdentifier > parameterIdList;
    for( OneWayRangePartialList::iterator oneWayPartialIterator = constituentOneWayRangePartials.begin( );
         oneWayPartialIterator != constituentOneWayRangePartials.end( ); oneWayPartialIterator++ )
    {
        // Retrieve one-way range paritals
        oneWayRangeScalers[ oneWayPartialIterator->first ] = boost::dynamic_pointer_cast< OneWayRangeScaling >
                ( oneWayPartialIterator->second.second );

        // Iterate over all one-way range partials of current link
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
                throw std::runtime_error( "Error when making n way range partial, parameter indices are inconsistent" );
            }
        }
    }


    // Create n-way range scaling object
    boost::shared_ptr< NWayRangeScaling > nWayRangeScaler = boost::make_shared< NWayRangeScaling >(
                oneWayRangeScalers, nWayRangeLinkEnds.size( ) );

    // Create n-way range partial object
    for( std::map< std::pair< int, int >, std::map< int, boost::shared_ptr< ObservationPartial< 1 > > > >::iterator sortedPartialIterator =
         sortedOneWayRangePartials.begin( ); sortedPartialIterator != sortedOneWayRangePartials.end( ); sortedPartialIterator++ )
    {
        nWayRangePartialList[ sortedPartialIterator->first ] = boost::make_shared< NWayRangePartial >(
                    nWayRangeScaler, sortedPartialIterator->second, parameterIdList.at( sortedPartialIterator->first ),
                    numberOfLinkEnds );
    }

    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =  parametersToEstimate->getVectorParameters( );
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        boost::shared_ptr< ObservationPartial< 1 > > currentNWayRangePartial;
        if( isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentNWayRangePartial = createObservationPartialWrtLinkProperty< 1 >(
                        nWayRangeLinkEnds, observation_models::n_way_range, parameterIterator->second );
        }

        // Check if partial is non-null
        if( currentNWayRangePartial != NULL )
        {
            // Add partial to the list.
            std::pair< double, double > currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            nWayRangePartialList[ currentPair ] = currentNWayRangePartial;
        }
    }

    return std::make_pair( nWayRangePartialList, nWayRangeScaler );
}

//! Function to generate n-way range partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate n-way range partials for all parameters that are to be estimated, for all sets of link ends.
 *  The n-way range partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  The n-way range partials are built from one-way range partials of the constituent links
 *  \param linkEnds List of all n-way link ends sets with observation models for which partials are to be created
 *  \param bodyMap List of all bodies, for creating n-way range partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default). First vector entry is
 *  index of link in n-way link ends, second vector is list of light-time corrections.
 *  \return Map of SingleLinkObservationPartialList, representing all necessary n-way range partials of a single link end,
 *  and NWayRangeScaling, object, used for scaling the position partial members of all NWayRangePartials in link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > >
createNWayRangePartials(
        const std::vector< observation_models::LinkEnds >& linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > > partialMap;
    std::vector< std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;

    // Iterate over all sets of link ends, and  create associated n-way range partials
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Retrieve light-time corrections
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            currentLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) );
            if( currentLightTimeCorrections.size( ) != linkEnds.at( i ).size( ) - 1 )
            {
                throw std::runtime_error(
                            "Error when making n-way range partials, found light time correction partials for " +
                            std::to_string( currentLightTimeCorrections.size( ) ) +
                            " links, with " + std::to_string( linkEnds.size( ) ) + " link ends" );
            }
        }
        else
        {
            currentLightTimeCorrections.clear( );
        }

        // Create n-way range partials for current LinkEnds
        partialMap[ linkEnds.at( i ) ] = createNWayRangePartials< ParameterType >(
                    linkEnds.at( i ), bodyMap, parametersToEstimate, currentLightTimeCorrections );

    }

    return partialMap;
}

}

}

#endif // TUDAT_CREATENWAYRANGEPARTIALS_H
