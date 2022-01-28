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
#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/twoWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/twoWayDopplerObservationModel.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"

namespace tudat
{

namespace observation_partials
{

//! Function to generate one-way doppler partial wrt a single  parameter.
/*!
 *  Function to generate one-way doppler partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param oneWayDopplerLinkEnds Link ends (transmitter and receiever) for which one-way doppler partials are to be calculated
 *  (i.e. for which  one-way doppler observations are to be processed).
 *  \param bodies List of all bodies, for creating one-way doppler partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param oneWayDopplerScaler Object scale position partials to one-way doppler partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way doppler partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< ObservationPartial< 1 > > createOneWayDopplerPartialWrtParameter(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    if( std::dynamic_pointer_cast< OneWayDopplerScaling >( oneWayDopplerScaler ) == nullptr )
    {
        throw std::runtime_error( "Error, expected one-way doppler scaling when making one-way doppler partial" );
    }
    std::shared_ptr< ObservationPartial< 1 > > oneWayDopplerPartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( oneWayDopplerLinkEnds, bodies, parameterToEstimate );

        // Create one-doppler partials if any position partials are created (i.e. if any dependency exists).
        std::shared_ptr< OneWayDopplerPartial > testOneWayDopplerPartial  = std::make_shared< OneWayDopplerPartial >(
                    std::dynamic_pointer_cast< OneWayDopplerScaling >( oneWayDopplerScaler ),
                    positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );

        if( positionPartials.size( ) > 0 || testOneWayDopplerPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) > 0
                || oneWayDopplerScaler->getProperTimeParameterDependencySize( parameterToEstimate->getParameterName( ) ) > 0 )
        {
            oneWayDopplerPartial = testOneWayDopplerPartial;
        }
    }

    // Return doppler partial object (nullptr if no dependency exists).
    return oneWayDopplerPartial;
}

//! Function to generate one-way doppler partial wrt a position of a body.
/*!
 *  Function to generate one-way doppler partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param oneWayDopplerLinkEnds Link ends (transmitter and receiever) for which one-way doppler partials are to be calculated
 *  (i.e. for which one-way doppler observations are to be processed).
 *  \param bodies List of all bodies, for creating one-way doppler partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param oneWayDopplerScaler Object scale position partials to one-way doppler partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way doppler partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< OneWayDopplerPartial > createOneWayDopplerPartialWrtBodyState(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< PositionPartialScaling > oneWayDopplerScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) );

//! Function to create an object that computes the scaling of the state partials to obtain proper time rate partials
/*!
 * Function to create an object that computes the scaling of the state partials to obtain proper time rate partials. A single
 * scaling object is used for a single link end of the one-way Doppler partials
 * \param dopplerProperTimeInterface Object that is used to computed proper-time rate in one-way Doppler modelkkl
 * \param oneWayDopplerLinkEnds Link ends of observable
 * \param linkEndAtWhichPartialIsComputed Link end for which proper-time partials are to be created
 * \return Scaling object for proper-time rate partials
 */
std::shared_ptr< OneWayDopplerProperTimeComponentScaling > createDopplerProperTimePartials(
        const std::shared_ptr< observation_models::DopplerProperTimeRateInterface > dopplerProperTimeInterface,
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const observation_models::LinkEndType linkEndAtWhichPartialIsComputed  );

//! Function to generate one-way doppler partials and associated scaler for single link end.
/*!
 *  Function to generate one-way doppler partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  \param oneWayDopplerLinkEnds Link ends (transmitter and receiever) for which one-way doppler partials are to be calculated
 *  (i.e. for which one-way doppler observations are to be processed).
 *  \param bodies List of all bodies, for creating one-way doppler partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param transmitterDopplerProperTimeInterface Proper time rate calculator for transmitter
 *  \param receiverDopplerProperTimeInterface Proper time rate calculator for receiver
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary one-way doppler partials of a single link end, and OneWayDopplerScaling, object, used for
 *  scaling the position partial members of all OneWayDopplerPartials in link end.
 */
template< typename ParameterType >
std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > createOneWayDopplerPartials(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::shared_ptr< observation_models::DopplerProperTimeRateInterface > transmitterDopplerProperTimeInterface,
        const std::shared_ptr< observation_models::DopplerProperTimeRateInterface > receiverDopplerProperTimeInterface,
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrections =
        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >( ),
        const bool useBiasPartials = true )
{
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        lightTimeCorrectionPartialObjects = observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections );
    }

    std::function< Eigen::Vector6d( const double )> transmitterNumericalStateDerivativeFunction =
            std::bind( &numerical_derivatives::computeCentralDifferenceFromFunction< Eigen::Vector6d, double >,
                         simulation_setup::getLinkEndCompleteEphemerisFunction< double, double >(
                             oneWayDopplerLinkEnds.at( observation_models::transmitter ), bodies ), std::placeholders::_1, 100.0,
                         numerical_derivatives::order8 );
    std::function< Eigen::Vector6d( const double )> receiverNumericalStateDerivativeFunction =
            std::bind( numerical_derivatives::computeCentralDifferenceFromFunction< Eigen::Vector6d, double >,
                         simulation_setup::getLinkEndCompleteEphemerisFunction< double, double >(
                             oneWayDopplerLinkEnds.at( observation_models::receiver ), bodies ), std::placeholders::_1, 100.0,
                         numerical_derivatives::order8 );

    // Create scaling object, to be used for all one-way doppler partials in current link end.
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling > transmitterProperTimePartials =
            createDopplerProperTimePartials( transmitterDopplerProperTimeInterface, oneWayDopplerLinkEnds,
                                             observation_models::transmitter );
    std::shared_ptr< OneWayDopplerProperTimeComponentScaling > receiverProperTimePartials =
            createDopplerProperTimePartials( receiverDopplerProperTimeInterface, oneWayDopplerLinkEnds,
                                             observation_models::receiver  );

    std::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaling = std::make_shared< OneWayDopplerScaling >(
                std::bind( &linear_algebra::evaluateSecondBlockInStateVector, transmitterNumericalStateDerivativeFunction, std::placeholders::_1 ),
                std::bind( &linear_algebra::evaluateSecondBlockInStateVector, receiverNumericalStateDerivativeFunction, std::placeholders::_1 ),
                transmitterProperTimePartials,
                receiverProperTimePartials );


    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    SingleLinkObservationPartialList dopplerPartials;

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        std::string acceleratedBody;
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state ||
                ( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state ) )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else
        {
            throw std::runtime_error( "Error when making one way doppler partials, could not identify parameter " +
                                      std::to_string(
                                          initialDynamicalParameters.at( i )->getParameterName( ).first ) );
        }


        // Create position one-way doppler partial for current body
        std::shared_ptr< ObservationPartial< 1 > > currentDopplerPartialWrtPosition = createOneWayDopplerPartialWrtBodyState(
                    oneWayDopplerLinkEnds, bodies, acceleratedBody, oneWayDopplerScaling, lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current doppler and current body)
        if( currentDopplerPartialWrtPosition != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            dopplerPartials[ currentPair ] = currentDopplerPartialWrtPosition;
        }

        // Increment current index by size of body initial state (6).
        currentIndex += 6;
    }

    // Iterate over all double parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator =
         doubleParametersToEstimate.begin( ); parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position one-way doppler partial for current parameter
        std::shared_ptr< ObservationPartial< 1 > > currentDopplerPartial = createOneWayDopplerPartialWrtParameter(
                    oneWayDopplerLinkEnds, bodies, parameterIterator->second,
                    oneWayDopplerScaling, lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current doppler and current parameter)
        if( currentDopplerPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            dopplerPartials[ currentPair ] = currentDopplerPartial;
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

        std::shared_ptr< ObservationPartial< 1 > > currentDopplerPartial;
        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentDopplerPartial = createOneWayDopplerPartialWrtParameter(
                        oneWayDopplerLinkEnds, bodies, parameterIterator->second, oneWayDopplerScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentDopplerPartial = createObservationPartialWrtLinkProperty< 1 >(
                        oneWayDopplerLinkEnds, observation_models::one_way_doppler, parameterIterator->second, useBiasPartials );
        }


        // Check if partial is non-nullptr (i.e. whether dependency exists between current doppler and current parameter)
        if( currentDopplerPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            dopplerPartials[ currentPair ] = currentDopplerPartial;
        }
    }

    // Return complete set of partials and scaling object.
    return std::make_pair( dopplerPartials, oneWayDopplerScaling );
}

//! Function to generate one-way doppler partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate one-way doppler partials for all parameters that are to be estimated, for all sets of link ends.
 *  The one-way doppler partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  \param observationModelList List of all observation models (must be one-way Doppler) for which partials are to be created,
 *  with map key being the link ends
 *  \param bodies List of all bodies, for creating one-way doppler partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Map of SingleLinkObservationPartialList, representing all necessary one-way doppler partials of a single link end,
 *  and OneWayDopplerScaling, object, used for scaling the position partial members of all OneWayDopplerPartials in link end.
 */
template< typename ObservationScalarType, typename TimeType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
std::shared_ptr< PositionPartialScaling > > > createOneWayDopplerPartials(
        const std::map< observation_models::LinkEnds,
        std::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > > observationModelList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate,
        const bool useBiasPartials = true )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
            std::shared_ptr< PositionPartialScaling > > > dopplerPartials;

    // Iterate over all link ends.
    for( auto linkIterator: observationModelList )
    {
        // Check if required link end types are present
        if( ( linkIterator.first.count( observation_models::receiver ) == 0 ) ||
                ( linkIterator.first.count( observation_models::transmitter ) == 0 ) )
        {
            throw std::runtime_error( "Error when making 1-way doppler partials, did not find both receiver and transmitter in link ends" );

        }


        std::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
                dopplerObservationModel =
                std::dynamic_pointer_cast< observation_models::OneWayDopplerObservationModel
                < ObservationScalarType, TimeType > >( linkIterator.second );

        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > singleLinkLightTimeCorrections;
        if( dopplerObservationModel == nullptr )
        {
            throw std::runtime_error( "Error when making one-way Doppler partials. Type is inconsistent" );
        }
        else
        {
            singleLinkLightTimeCorrections = dopplerObservationModel->getLightTimeCalculator( )->getLightTimeCorrection( );
        }

        // Create doppler partials for current link ends
        dopplerPartials[ linkIterator.first ] = createOneWayDopplerPartials< ObservationScalarType >(
                    linkIterator.first, bodies, parametersToEstimate,
                    dopplerObservationModel->getTransmitterProperTimeRateCalculator( ),
                    dopplerObservationModel->getReceiverProperTimeRateCalculator( ),
                    singleLinkLightTimeCorrections, useBiasPartials );
    }

    // Return complete set of link ends.
    return dopplerPartials;
}




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
        const observation_models::LinkEnds& twoWayDopplerLinkEnds,
        const std::shared_ptr< observation_models::TwoWayDopplerObservationModel< ParameterType, TimeType > >
        twoWayObservationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > lightTimeCorrections =
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >( ) )

{

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
        currentLightTimeCorrections.clear( );
        if( lightTimeCorrections.size( ) > 0 )
        {
            currentLightTimeCorrections = lightTimeCorrections.at( i );
        }

        // Define links for current one-way Doppler link
        currentLinkEnds.clear( );
        currentLinkEnds[ observation_models::transmitter ] = twoWayDopplerLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i, 3 ) );
        currentLinkEnds[ observation_models::receiver ] = twoWayDopplerLinkEnds.at(
                    observation_models::getNWayLinkEnumFromIndex( i + 1, 3 ) );

        std::shared_ptr< observation_models::DopplerProperTimeRateInterface > transmitterDopplerProperTimeInterface;
        std::shared_ptr< observation_models::DopplerProperTimeRateInterface > receiverDopplerProperTimeInterface;
        if( i == 0 )
        {
            transmitterDopplerProperTimeInterface =
                    twoWayObservationModel->getUplinkDopplerCalculator( )->getTransmitterProperTimeRateCalculator( );
            receiverDopplerProperTimeInterface =
                    twoWayObservationModel->getUplinkDopplerCalculator( )->getReceiverProperTimeRateCalculator( );
            oneWayDopplerModels.push_back(
                        observation_models::getSizeOneObservationFunctionAtDoublePrecisionFromObservationModel<
                        ParameterType, TimeType >( twoWayObservationModel->getUplinkDopplerCalculator( ) ) );

        }
        else if( i == 1 )
        {
            transmitterDopplerProperTimeInterface =
                    twoWayObservationModel->getDownlinkDopplerCalculator( )->getTransmitterProperTimeRateCalculator( );
            receiverDopplerProperTimeInterface =
                    twoWayObservationModel->getDownlinkDopplerCalculator( )->getReceiverProperTimeRateCalculator( );
            oneWayDopplerModels.push_back(
                        observation_models::getSizeOneObservationFunctionAtDoublePrecisionFromObservationModel<
                        ParameterType, TimeType >( twoWayObservationModel->getDownlinkDopplerCalculator( ) ) );
        }


        // Create one-way Doppler partials for current link
        constituentOneWayDopplerPartials.push_back(
                    createOneWayDopplerPartials(
                        currentLinkEnds, bodies, parametersToEstimate, transmitterDopplerProperTimeInterface,
                        receiverDopplerProperTimeInterface, currentLightTimeCorrections, false ) );
        constituentOneWayRangePartials.push_back(
                    createOneWayRangePartials( currentLinkEnds, bodies, parametersToEstimate,
                                               currentLightTimeCorrections ) );
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
                oneWayDopplerScalers, oneWayRangeScalings, oneWayDopplerModels );

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
        if( isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentTwoWayDopplerPartial = createObservationPartialWrtLinkProperty< 1 >(
                        twoWayDopplerLinkEnds, observation_models::two_way_doppler, parameterIterator->second );
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

//! Function to generate two-way Doppler partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate two-way Doppler partials for all parameters that are to be estimated, for all sets of link ends.
 *  The two-way Doppler partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  The two-way Doppler partials are built from one-way range partials of the constituent links
 *  \param observationModelList List of all two-way Doppler models (as a function of LinkEnds) for which partials are to be
 *  created.
 *  \param bodies List of all bodies, for creating two-way Doppler partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction used (empty by default). First vector entry is
 *  index of link in 2-way link ends (up and downlink), second vector is list of light-time corrections.
 *  \return Map of SingleLinkObservationPartialList, representing all necessary two-way Doppler partials of a single link end,
 *  and TwoWayDopplerScaling, object, used for scaling the position partial members of all TwoWayDopplerPartials in link end.
 */
template< typename ParameterType, typename TimeType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > >
createTwoWayDopplerPartials(
        const std::map< observation_models::LinkEnds,
        std::shared_ptr< observation_models::ObservationModel< 1, ParameterType, TimeType > > > observationModelList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    std::map< observation_models::LinkEnds,
            std::pair< SingleLinkObservationPartialList, std::shared_ptr< PositionPartialScaling > > > partialMap;
    std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;

    // Iterate over all sets of link ends, and  create associated two-way Doppler partials
    for( typename std::map< observation_models::LinkEnds,
         std::shared_ptr< observation_models::ObservationModel< 1, ParameterType, TimeType > > >::const_iterator modelIterator =
         observationModelList.begin( ); modelIterator != observationModelList.end( ); modelIterator++ )
    {
        // Retrieve light-time corrections
        if( lightTimeCorrections.count( modelIterator->first ) > 0 )
        {
            currentLightTimeCorrections = lightTimeCorrections.at( modelIterator->first );
        }
        else
        {
            currentLightTimeCorrections.clear( );
        }

        // Create two-way Doppler partials for current LinkEnds
        partialMap[ modelIterator->first ] = createTwoWayDopplerPartials< ParameterType, TimeType >(
                    modelIterator->first, std::dynamic_pointer_cast< observation_models::TwoWayDopplerObservationModel<
                    ParameterType, TimeType > >( modelIterator->second ), bodies, parametersToEstimate,
                    currentLightTimeCorrections );
    }

    return partialMap;
}




}

}

#endif // TUDAT_CREATEONEWAYDOPPLERPARTIALS_H

