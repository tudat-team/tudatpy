/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEONEWAYRANGEPARTIALS_H
#define TUDAT_CREATEONEWAYRANGEPARTIALS_H


#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"
#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayRangePartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"

namespace tudat
{

namespace observation_partials
{

//! Function to generate one-way range partial wrt a single  parameter.
/*!
 *  Function to generate one-way range partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param oneWayRangeLinkEnds Link ends (transmitter and receiever) for which one-way range partials are to be calculated
 *  (i.e. for which  one-way range observations are to be processed).
 *  \param bodyMap List of all bodies, for creating one-way range partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param oneWayRangeScaler Object scale position partials to one-way range partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way range partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< ObservationPartial< 1 > > createOneWayRangePartialWrtParameter(
        const observation_models::LinkEnds oneWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< OneWayRangeScaling > oneWayRangeScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    std::shared_ptr< ObservationPartial< 1 > > oneWayRangePartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( oneWayRangeLinkEnds, bodyMap, parameterToEstimate );

        // Create one-range partials if any position partials are created (i.e. if any dependency exists).
        std::shared_ptr< OneWayRangePartial > testOneWayRangePartial  = std::make_shared< OneWayRangePartial >(
                    oneWayRangeScaler, positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );
        if( positionPartials.size( ) > 0 || testOneWayRangePartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) > 0 )
        {
            oneWayRangePartial = testOneWayRangePartial;
        }
    }

    // Return range partial object (nullptr if no dependency exists).
    return oneWayRangePartial;
}

//! Function to generate one-way range partial wrt a position of a body.
/*!
 *  Function to generate one-way range partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param oneWayRangeLinkEnds Link ends (transmitter and receiever) for which one-way range partials are to be calculated
 *  (i.e. for which one-way range observations are to be processed).
 *  \param bodyMap List of all bodies, for creating one-way range partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param oneWayRangeScaler Object scale position partials to one-way range partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way range partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< OneWayRangePartial > createOneWayRangePartialWrtBodyPosition(
        const observation_models::LinkEnds oneWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< OneWayRangeScaling > oneWayRangeScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) );

//! Function to generate one-way range partial wrt a rotational state of a body.
/*!
 *  Function to generate one-way range partial wrt a rotational state of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param oneWayRangeLinkEnds Link ends (transmitter and receiever) for which one-way range partials are to be calculated
 *  (i.e. for which one-way range observations are to be processed).
 *  \param bodyMap List of all bodies, for creating one-way range partial.
 *  \param bodyToEstimate Name of body wrt rotational state of which a partial is to be created.
 *  \param oneWayRangeScaler Object scale position partials to one-way range partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way range partial object wrt a current rotational state of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< OneWayRangePartial > createOneWayRangePartialWrtBodyRotationalState(
        const observation_models::LinkEnds oneWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const std::shared_ptr< OneWayRangeScaling > oneWayRangeScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects  );

//! Function to generate one-way range partials and associated scaler for single link ends.
/*!
 *  Function to generate one-way range partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends set.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  \param oneWayRangeLinkEnds Link ends (transmitter and receiever) for which one-way range partials are to be calculated
 *  (i.e. for which one-way range observations are to be processed).
 *  \param bodyMap List of all bodies, for creating one-way range partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary one-way range partials of a single link end, and OneWayRangeScaling, object, used for
 *  scaling the position partial members of all OneWayRangePartials in link end.
 */
template< typename ParameterType >
std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > createOneWayRangePartials(
        const observation_models::LinkEnds oneWayRangeLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrections =
        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >( ),
        const bool useBiasPartials = true )
{

    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        lightTimeCorrectionPartialObjects = observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections );
    }

    // Create scaling object, to be used for all one-way range partials in current link end.
    std::shared_ptr< OneWayRangeScaling > oneWayRangeScaling = std::make_shared< OneWayRangeScaling >( );

    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    SingleLinkObservationPartialList rangePartials;

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        std::string acceleratedBody;
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state ||
                initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;

            // Create position one-way range partial for current body
            std::shared_ptr< OneWayRangePartial > currentRangePartial = createOneWayRangePartialWrtBodyPosition(
                        oneWayRangeLinkEnds, bodyMap, acceleratedBody, oneWayRangeScaling, lightTimeCorrectionPartialObjects );

            // Check if partial is non-null (i.e. whether dependency exists between current range and current body)
            if( currentRangePartial != nullptr )
            {
                // Add partial to the list.
                currentPair = std::pair< int, int >( currentIndex, 6 );
                rangePartials[ currentPair ] = currentRangePartial;
            }

            // Increment current index by size of body initial state (6).
            currentIndex += 6;
        }        
        else if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_rotational_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;

            // Create position one-way range partial for current body
            std::shared_ptr< OneWayRangePartial > currentRangePartial = createOneWayRangePartialWrtBodyRotationalState(
                        oneWayRangeLinkEnds, bodyMap, acceleratedBody, oneWayRangeScaling, lightTimeCorrectionPartialObjects );

            // Check if partial is non-null (i.e. whether dependency exists between current range and current body)
            if( currentRangePartial != nullptr )
            {
                // Add partial to the list.
                currentPair = std::pair< int, int >( currentIndex, 7 );
                rangePartials[ currentPair ] = currentRangePartial;
            }

            // Increment current index by size of body initial state (6).
            currentIndex += 7;
        }
        else
        {
            throw std::runtime_error( "Error when making one way range partials, could not identify parameter " +
                       std::to_string( initialDynamicalParameters.at( i )->getParameterName( ).first ) );
        }
    }

    // Iterate over all double parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator =
         doubleParametersToEstimate.begin( ); parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position one-way range partial for current parameter
        std::shared_ptr< ObservationPartial< 1 > > currentRangePartial = createOneWayRangePartialWrtParameter(
                    oneWayRangeLinkEnds, bodyMap, parameterIterator->second,
                    oneWayRangeScaling, lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current range and current parameter)
        if( currentRangePartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            rangePartials[ currentPair ] = currentRangePartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =
            parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position one-way range partial for current parameter
        std::shared_ptr< ObservationPartial< 1 > > currentRangePartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentRangePartial = createOneWayRangePartialWrtParameter(
                    oneWayRangeLinkEnds, bodyMap, parameterIterator->second, oneWayRangeScaling,
                    lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentRangePartial = createObservationPartialWrtLinkProperty< 1 >(
                        oneWayRangeLinkEnds, observation_models::one_way_range, parameterIterator->second, useBiasPartials );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current range and current parameter)
        if( currentRangePartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            rangePartials[ currentPair ] = currentRangePartial;
        }
    }

    // Return complete set of partials and scaling object.
    return std::make_pair( rangePartials, oneWayRangeScaling );
}

//! Function to generate one-way range partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate one-way range partials for all parameters that are to be estimated, for all sets of link ends.
 *  The one-way range partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  \param linkEnds Vector of all link ends for which one-way range partials are to be calculated (i.e. for which one-way
 *  range observations are  to be processed).
 *  \param bodyMap List of all bodies, for creating one-way range partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Map of SingleLinkObservationPartialList, representing all necessary one-way range partials of a single link end,
 *  and OneWayRangeScaling, object, used for scaling the position partial members of all OneWayRangePartials in link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
std::shared_ptr< PositionPartialScaling > > > createOneWayRangePartials(
        const std::vector< observation_models::LinkEnds >& linkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ),
        const bool useBiasPartials = true )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
            std::shared_ptr< PositionPartialScaling > > > rangePartials;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Check if required link end types are present
        if( ( linkEnds[ i ].count( observation_models::receiver ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter ) == 0 ) )
        {
            throw std::runtime_error( "Error when making 1-way range partials, did not find both receiver and transmitter in link ends" );

        }

        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > singleLinkLightTimeCorrections;
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            if( lightTimeCorrections.at( linkEnds.at( i ) ).size( ) > 1 )
            {
                std::cerr << "Error when making one way range partials, light time corrections for "
                          << lightTimeCorrections.at( linkEnds.at( i ) ).size( ) << " links found" << std::endl;
            }
            else if( lightTimeCorrections.at( linkEnds.at( i ) ).size( ) == 1 )
            {
                singleLinkLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) ).at( 0 );
            }
        }

        // Create range partials for current link ends
        rangePartials[ linkEnds[ i ] ] = createOneWayRangePartials< ParameterType >(
                    linkEnds[ i ], bodyMap, parametersToEstimate, singleLinkLightTimeCorrections, useBiasPartials );
    }

    // Return complete set of link ends.
    return rangePartials;
}

}

}

#endif // TUDAT_CREATEONEWAYRANGEPARTIALS_H

