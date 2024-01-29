/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEONEWAYDOPPLERPARTIALS_H
#define TUDAT_CREATEONEWAYDOPPLERPARTIALS_H


#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/interpolator.h"

#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/simulation/estimation_setup/createLightTimeCalculator.h"
#include "tudat/simulation/estimation_setup/createDirectObservationPartials.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/twoWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/twoWayDopplerObservationModel.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"
#include "tudat/simulation/estimation_setup/createObsevationBiasPartial.h"

namespace tudat
{

namespace observation_partials
{

//! Function to generate tow-way Doppler partials and associated scaler for single link ends.
/*!
 *  Function to generate tow-way Doppler partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends set.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  The tow-way Doppler partials are built from one-way Doppler partials of the constituent links
 *  \param twoWayDopplerLinkEnds Link ends (transmitter and receiever) for which tow-way Doppler partials are to be calculated
 *  (i.e. for which tow-way Doppler observations are to be processed).
 *  \param twoWayObservationModel Observation model for two-way Doppler for which partials are to be created
 *  \param bodies List of all bodies, for creating tow-way Doppler partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default). First vector entry is
 *  index of link in 2-way link ends (up and downlink), second vector is list of light-time corrections.
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary two-way Doppler partials of a single link end, and TwoWayDopplerScaling, object, used for
 *  scaling the position partial members of all TwoWayDopplerPartials in link end.
 */
template< typename ParameterType, typename TimeType >
std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > createTwoWayDopplerPartials(
        const std::shared_ptr< observation_models::ObservationModel< 1, ParameterType, TimeType > > observationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable = false )
{    
    std::shared_ptr< observation_models::TwoWayDopplerObservationModel< ParameterType, TimeType > >  twoWayObservationModel =
            std::dynamic_pointer_cast< observation_models::TwoWayDopplerObservationModel< ParameterType, TimeType > >(
                observationModel );
    observation_models::LinkEnds twoWayDopplerLinkEnds = twoWayObservationModel->getLinkEnds( );


    // Define list of constituent one-way partials.
    typedef std::vector< std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > >
            OneWayPartialList;
    OneWayPartialList constituentOneWayDopplerPartials;
    OneWayPartialList constituentOneWayRangePartials;

    int numberOfLinkEnds = 3;
    observation_models::LinkEnds currentLinkEnds;

    // Iterate over all links in the two-way Doppler observable
    std::vector< std::function< double( const double, const observation_models::LinkEndType ) > > oneWayDopplerModels;
    std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > currentLightTimeCorrections;
    for( int i = 0; i < numberOfLinkEnds - 1; i++ )
    {
        if( i == 0 )
        {
            currentLightTimeCorrections =
                    twoWayObservationModel->getUplinkDopplerCalculator( )->getLightTimeCalculator( )->getLightTimeCorrection( );
        }
        else
        {
            currentLightTimeCorrections =
                    twoWayObservationModel->getDownlinkDopplerCalculator( )->getLightTimeCalculator( )->getLightTimeCorrection( );
        }

        // Define links for current one-way Doppler link
        currentLinkEnds.clear( );
        currentLinkEnds[ observation_models::transmitter ] = twoWayDopplerLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i, 3 ) );
        currentLinkEnds[ observation_models::receiver ] = twoWayDopplerLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i + 1, 3 ) );


        std::shared_ptr< observation_models::OneWayDopplerObservationModel< ParameterType, TimeType > > currentDopplerModel;
        if( i == 0 )
        {
            oneWayDopplerModels.push_back(
                        observation_models::getSizeOneObservationFunctionAtDoublePrecisionFromObservationModel<
                        ParameterType, TimeType >( twoWayObservationModel->getUplinkDopplerCalculator( ) ) );
            currentDopplerModel = twoWayObservationModel->getUplinkDopplerCalculator( );

        }
        else if( i == 1 )
        {
            oneWayDopplerModels.push_back(
                        observation_models::getSizeOneObservationFunctionAtDoublePrecisionFromObservationModel<
                        ParameterType, TimeType >( twoWayObservationModel->getDownlinkDopplerCalculator( ) ) );
            currentDopplerModel = twoWayObservationModel->getDownlinkDopplerCalculator( );

        }


        // Create one-way Doppler partials for current link
        constituentOneWayDopplerPartials.push_back(
                    createSingleLinkObservationPartials< ParameterType, 1, TimeType >(
                        currentDopplerModel, bodies, parametersToEstimate, false ) );

        constituentOneWayRangePartials.push_back(
                    createSingleLinkObservationPartials< ParameterType, 1, TimeType >(
                        std::make_shared< observation_models::OneWayRangeObservationModel< ParameterType, TimeType > >(
                            currentLinkEnds, currentDopplerModel->getLightTimeCalculator( ) ), bodies,
                        parametersToEstimate, false ) );
    }

    // Retrieve sorted (by parameter index and link index) one-way range partials and (by link index) opne-way range partials
    std::vector< std::shared_ptr< OneWayDopplerScaling > > oneWayDopplerScalers;
    std::vector< std::shared_ptr< OneWayRangeScaling > > oneWayRangeScalings;

    std::map< std::pair< int, int >, std::map< int, std::shared_ptr< ObservationPartial< 1 > > > > sortedOneWayDopplerPartials;
    std::map< std::pair< int, int >, std::map< int, std::shared_ptr< ObservationPartial< 1 > > > > sortedOneWayRangePartials;

    std::map< std::pair< int, int >, estimatable_parameters::EstimatebleParameterIdentifier > parameterIdList;
    for( unsigned int i = 0; i < constituentOneWayDopplerPartials.size( ); i++ )
    {
        // Retrieve one-way Doppler paritals
        oneWayDopplerScalers.push_back( std::dynamic_pointer_cast< OneWayDopplerScaling >
                                        ( constituentOneWayDopplerPartials.at( i ).second ) );
        oneWayRangeScalings.push_back( std::dynamic_pointer_cast< OneWayRangeScaling >
                                       ( constituentOneWayRangePartials.at( i ).second ) );

        // Iterate over all one-way Doppler partials of current link
        for( SingleLinkObservationPartialList::iterator parameterIterator = constituentOneWayDopplerPartials.at( i ).first.begin( );
             parameterIterator != constituentOneWayDopplerPartials.at( i ).first.end( ); parameterIterator++ )
        {
            sortedOneWayDopplerPartials[ parameterIterator->first ][ i ] = parameterIterator->second;
            if( parameterIdList.count( parameterIterator->first ) == 0 )
            {
                parameterIdList[ parameterIterator->first ] = parameterIterator->second->getParameterIdentifier( );
            }
            else if( parameterIdList.at( parameterIterator->first ) != parameterIterator->second->getParameterIdentifier( ) )
            {
                throw std::runtime_error( "Error when making two way Doppler partial, parameter indices in Doppler are inconsistent" );
            }
        }

        // Iterate over all one-way range partials of current link
        for( SingleLinkObservationPartialList::iterator parameterIterator = constituentOneWayRangePartials.at( i ).first.begin( );
             parameterIterator != constituentOneWayRangePartials.at( i ).first.end( ); parameterIterator++ )
        {
            sortedOneWayRangePartials[ parameterIterator->first ][ i ] = parameterIterator->second;
            if( parameterIdList.count( parameterIterator->first ) == 0 )
            {
                parameterIdList[ parameterIterator->first ] = parameterIterator->second->getParameterIdentifier( );
            }
            else if( parameterIdList.at( parameterIterator->first ) != parameterIterator->second->getParameterIdentifier( ) )
            {
                throw std::runtime_error( "Error when making two way Doppler partial, parameter indices in range are inconsistent" );
            }
        }
    }

    // Create two-way Doppler scaling object
    std::shared_ptr< TwoWayDopplerScaling > twoWayDopplerScaler = std::make_shared< TwoWayDopplerScaling >(
                oneWayDopplerScalers, oneWayRangeScalings, oneWayDopplerModels,
                twoWayObservationModel->getNormalizeWithSpeedOfLight( ) ? 1.0 : physical_constants::SPEED_OF_LIGHT );

    // Define return partial list
    SingleLinkObservationPartialList twoWayDopplerPartialList;

    // Create two-way Doppler partial objects
    for( std::map< std::pair< int, int >, std::map< int, std::shared_ptr< ObservationPartial< 1 > > > >::iterator sortedPartialIterator =
         sortedOneWayDopplerPartials.begin( ); sortedPartialIterator != sortedOneWayDopplerPartials.end( ); sortedPartialIterator++ )
    {
        std::map< int, std::shared_ptr< ObservationPartial< 1 > > > currentDopplerPartialList = sortedPartialIterator->second;
        std::map< int, std::shared_ptr< ObservationPartial< 1 > > > currentRangePartialList;
        if( sortedOneWayRangePartials.count( sortedPartialIterator->first ) != 0 )
        {
            currentRangePartialList = sortedOneWayRangePartials.at( sortedPartialIterator->first );
        }

        twoWayDopplerPartialList[ sortedPartialIterator->first ] = std::make_shared< TwoWayDopplerPartial >(
                    twoWayDopplerScaler, currentDopplerPartialList, currentRangePartialList,
                    parameterIdList.at( sortedPartialIterator->first ), numberOfLinkEnds );
    }

    // Create two-way Doppler partial objects that only have range dependencies.
    for( std::map< std::pair< int, int >, std::map< int, std::shared_ptr< ObservationPartial< 1 > > > >::iterator sortedPartialIterator =
         sortedOneWayRangePartials.begin( ); sortedPartialIterator != sortedOneWayRangePartials.end( ); sortedPartialIterator++ )
    {
        if( sortedOneWayDopplerPartials.count( sortedPartialIterator->first ) == 0 )
        {
            std::map< int, std::shared_ptr< ObservationPartial< 1 > > > currentDopplerPartialList;
            std::map< int, std::shared_ptr< ObservationPartial< 1 > > > currentRangePartialList = sortedPartialIterator->second;

            twoWayDopplerPartialList[ sortedPartialIterator->first ] = std::make_shared< TwoWayDopplerPartial >(
                        twoWayDopplerScaler, currentDopplerPartialList, currentRangePartialList,
                        parameterIdList.at( sortedPartialIterator->first ), numberOfLinkEnds );
        }
    }

    // Create two-way Doppler partials
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =  parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        std::shared_ptr< ObservationPartial< 1 > > currentTwoWayDopplerPartial;
        if( isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first ) )
        {
            currentTwoWayDopplerPartial = createObservationPartialWrtLinkProperty< 1 >(
                        twoWayDopplerLinkEnds, observation_models::two_way_doppler, parameterIterator->second, bodies,
                        isPartialForDifferencedObservable );
        }


        // Check if partial is non-nullptr (i.e. whether dependency exists between current doppler and current parameter)
        if( currentTwoWayDopplerPartial != nullptr )
        {
            // Add partial to the list.
            std::pair< double, double > currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            twoWayDopplerPartialList[ currentPair ] = currentTwoWayDopplerPartial;
        }
    }

    return std::make_pair( twoWayDopplerPartialList, twoWayDopplerScaler );
}





}

}

#endif // TUDAT_CREATEONEWAYDOPPLERPARTIALS_H

