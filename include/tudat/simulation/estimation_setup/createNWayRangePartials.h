/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/interpolator.h"

#include "tudat/simulation/estimation_setup/createDirectObservationPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/nWayRangePartial.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/simulation/estimation_setup/createObsevationBiasPartial.h"

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
 *  \param bodies List of all bodies, for creating n-way range partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default). First vector entry is
 *  index of link in n-way link ends, second vector is list of light-time corrections.
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary n-way range partials of a single link end, and NWayRangeScaling, object, used for
 *  scaling the position partial members of all NWayRangePartials in link end.
 */
template< typename ParameterType, typename TimeType >
std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > createNWayRangePartials(
        const std::shared_ptr< observation_models::ObservationModel< 1, ParameterType, TimeType > > observationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable )
{
    using namespace observation_models;

    std::shared_ptr< observation_models::NWayRangeObservationModel< ParameterType, TimeType > >
            nWayRangeObservationModel =
            std::dynamic_pointer_cast< observation_models::NWayRangeObservationModel< ParameterType, TimeType > >(
                observationModel );
    if( nWayRangeObservationModel == nullptr )
    {
        throw std::runtime_error( "Error when creating n-way range partials; input observation model is not n-way range" );
    }

    observation_models::LinkEnds nWayRangeLinkEnds = nWayRangeObservationModel->getLinkEnds( );


    // Define return partial list
    SingleLinkObservationPartialList nWayRangePartialList;

    // Define list of constituent one-way partials.
    typedef std::map< int, std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > >
            OneWayRangePartialList;
    OneWayRangePartialList constituentOneWayRangePartials;

    // Getr number of link ends,
    observation_models::LinkEnds currentLinkEnds;
    int numberOfLinkEnds = nWayRangeLinkEnds.size( );

    // Iterate over all links in the n-way range observable
    std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > currentLightTimeCorrections;
    for( int i = 0; i < numberOfLinkEnds - 1; i++ )
    {
        currentLightTimeCorrections.clear( );
        currentLightTimeCorrections = nWayRangeObservationModel->getLightTimeCalculators( ).at( i )->getLightTimeCorrection( );

        // Define links for current one-way range link
        currentLinkEnds.clear( );
        currentLinkEnds[ observation_models::transmitter ] = nWayRangeLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i, numberOfLinkEnds ) );
        currentLinkEnds[ observation_models::receiver ] = nWayRangeLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i + 1, numberOfLinkEnds ) );

        // Create onw-way range partials for current link
        constituentOneWayRangePartials[ i ] =
                createSingleLinkObservationPartials< ParameterType, 1, TimeType >
                ( std::make_shared< OneWayRangeObservationModel< ParameterType, TimeType > >(
                      currentLinkEnds, nWayRangeObservationModel->getLightTimeCalculators( ).at( i ) ),
                  bodies, parametersToEstimate, false, true );
    }

    // Retrieve sorted (by parameter index and link index) one-way range partials and (by link index) opne-way range partials
    std::map< int, std::shared_ptr< OneWayRangeScaling > > oneWayRangeScalers;
    std::map< std::pair< int, int >, std::map< int, std::shared_ptr< ObservationPartial< 1 > > > > sortedOneWayRangePartials;
    std::map< std::pair< int, int >, estimatable_parameters::EstimatebleParameterIdentifier > parameterIdList;
    for( OneWayRangePartialList::iterator oneWayPartialIterator = constituentOneWayRangePartials.begin( );
         oneWayPartialIterator != constituentOneWayRangePartials.end( ); oneWayPartialIterator++ )
    {
        // Retrieve one-way range paritals
        oneWayRangeScalers[ oneWayPartialIterator->first ] = std::dynamic_pointer_cast< OneWayRangeScaling >
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

//            else if( parameterIdList.at( parameterIterator->first ) != parameterIterator->second->getParameterIdentifier( ) )
//            {
//                throw std::runtime_error( "Error when making n way range partial, parameter indices are inconsistent" );
//            }
        }
    }


    // Create n-way range scaling object
    std::shared_ptr< NWayRangeScaling > nWayRangeScaler = std::make_shared< NWayRangeScaling >(
                oneWayRangeScalers, nWayRangeLinkEnds.size( ) );

    // Create n-way range partial object
    for( std::map< std::pair< int, int >, std::map< int, std::shared_ptr< ObservationPartial< 1 > > > >::iterator sortedPartialIterator =
         sortedOneWayRangePartials.begin( ); sortedPartialIterator != sortedOneWayRangePartials.end( ); sortedPartialIterator++ )
    {
        nWayRangePartialList[ sortedPartialIterator->first ] = std::make_shared< NWayRangePartial >(
                    nWayRangeScaler, sortedPartialIterator->second, parameterIdList.at( sortedPartialIterator->first ),
                    numberOfLinkEnds );
    }

    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =  parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        std::shared_ptr< ObservationPartial< 1 > > currentNWayRangePartial;
        if( isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first ) &&
            !isParameterObservationLinkTimeProperty( parameterIterator->second->getParameterName( ).first ) &&
        !isPartialForDifferencedObservable )
        {
            currentNWayRangePartial = createObservationPartialWrtLinkProperty< 1 >(
                        nWayRangeLinkEnds, observation_models::n_way_range, parameterIterator->second, bodies );
        }

        // Check if partial is non-nullptr
        if( currentNWayRangePartial != nullptr )
        {
            // Add partial to the list.
            std::pair< double, double > currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            nWayRangePartialList[ currentPair ] = currentNWayRangePartial;
        }
    }

    return std::make_pair( nWayRangePartialList, nWayRangeScaler );
}


}

}

#endif // TUDAT_CREATENWAYRANGEPARTIALS_H
