/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Mathematics/Interpolators/interpolator.h"

#include "Tudat/Astrodynamics/ObservationModels/ObservableCorrections/lightTimeCorrection.h"
#include "Tudat/SimulationSetup/EstimationSetup/createCartesianStatePartials.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCalculator.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/oneWayDopplerPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/initialTranslationalState.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"
#include "Tudat/SimulationSetup/EstimationSetup/createLightTimeCorrectionPartials.h"

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
 *  \param bodyMap List of all bodies, for creating one-way doppler partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param oneWayDopplerScaler Object scale position partials to one-way doppler partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way doppler partial object wrt a single parameter (is NULL if no parameter dependency exists).
 */
template< typename ParameterType >
boost::shared_ptr< ObservationPartial< 1 > > createOneWayDopplerPartialWrtParameter(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const boost::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaler,
        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    if( boost::dynamic_pointer_cast< OneWayDopplerScaling >( oneWayDopplerScaler ) == NULL )
    {
        throw std::runtime_error( "Error, expected one-way doppler scaling when making one-way doppler partial" );
    }
    boost::shared_ptr< ObservationPartial< 1 > > oneWayDopplerPartial;

    {
        std::map< observation_models::LinkEndType, boost::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( oneWayDopplerLinkEnds, bodyMap, parameterToEstimate );

        // Create one-doppler partials if any position partials are created (i.e. if any dependency exists).
        boost::shared_ptr< OneWayDopplerPartial > testOneWayDopplerPartial  = boost::make_shared< OneWayDopplerPartial >(
                    boost::dynamic_pointer_cast< OneWayDopplerScaling >( oneWayDopplerScaler ),
                    positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );

        if( positionPartials.size( ) > 0 || testOneWayDopplerPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) > 0
                || oneWayDopplerScaler->getProperTimeParameterDependencySize( parameterToEstimate->getParameterName( ) ) > 0 )
        {
            oneWayDopplerPartial = testOneWayDopplerPartial;
        }
    }

    // Return doppler partial object (NULL if no dependency exists).
    return oneWayDopplerPartial;
}

//! Function to generate one-way doppler partial wrt a position of a body.
/*!
 *  Function to generate one-way doppler partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param oneWayDopplerLinkEnds Link ends (transmitter and receiever) for which one-way doppler partials are to be calculated
 *  (i.e. for which one-way doppler observations are to be processed).
 *  \param bodyMap List of all bodies, for creating one-way doppler partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param oneWayDopplerScaler Object scale position partials to one-way doppler partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return One-way doppler partial object wrt a current position of a body (is NULL if no parameter dependency exists).
 */
boost::shared_ptr< OneWayDopplerPartial > createOneWayDopplerPartialWrtBodyState(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::string bodyToEstimate,
        const boost::shared_ptr< PositionPartialScaling > oneWayDopplerScaler,
        const std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) );

//! Function to create an object that computes the scaling of the state partials to obtain proper time rate partials
/*!
 * Function to create an object that computes the scaling of the state partials to obtain proper time rate partials. A single
 * scaling object is used for a single link end of the one-way Doppler partials
 * \param dopplerProperTimeInterface Object that is used to computed proper-time rate in one-way Doppler modelkkl
 * \param oneWayDopplerLinkEnds Link ends of observable
 * \param linkEndAtWhichPartialIsComputed Link end for which proper-time partials are to be created
 * \return Scaling object for proper-time rate partials
 */
boost::shared_ptr< OneWayDopplerProperTimeComponentScaling > createDopplerProperTimePartials(
        const boost::shared_ptr< observation_models::DopplerProperTimeRateInterface > dopplerProperTimeInterface,
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
 *  \param bodyMap List of all bodies, for creating one-way doppler partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param transmitterDopplerProperTimeInterface Proper time rate calculator for transmitter
 *  \param receiverDopplerProperTimeInterface Proper time rate calculator for receiver
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary one-way doppler partials of a single link end, and OneWayDopplerScaling, object, used for
 *  scaling the position partial members of all OneWayDopplerPartials in link end.
 */
template< typename ParameterType >
std::pair< SingleLinkObservationPartialList, boost::shared_ptr< PositionPartialScaling > > createOneWayDopplerPartials(
        const observation_models::LinkEnds oneWayDopplerLinkEnds,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const boost::shared_ptr< observation_models::DopplerProperTimeRateInterface > transmitterDopplerProperTimeInterface,
        const boost::shared_ptr< observation_models::DopplerProperTimeRateInterface > receiverDopplerProperTimeInterface,
        const std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrections =
        std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > >( ) )
{
    std::vector< boost::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        lightTimeCorrectionPartialObjects = observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections );
    }

    boost::function< Eigen::Vector6d( const double )> transmitterNumericalStateDerivativeFunction =
            boost::bind( &numerical_derivatives::computeCentralDifference< Eigen::Vector6d, double >,
                observation_models::getLinkEndCompleteEphemerisFunction< double, double >(
                    oneWayDopplerLinkEnds.at( observation_models::transmitter ), bodyMap ), _1, 100.0,
                numerical_derivatives::order8 );
    boost::function< Eigen::Vector6d( const double )> receiverNumericalStateDerivativeFunction =
            boost::bind( numerical_derivatives::computeCentralDifference< Eigen::Vector6d, double >,
                observation_models::getLinkEndCompleteEphemerisFunction< double, double >(
                    oneWayDopplerLinkEnds.at( observation_models::receiver ), bodyMap ), _1, 100.0,
                numerical_derivatives::order8 );

    // Create scaling object, to be used for all one-way doppler partials in current link end.
    boost::shared_ptr< OneWayDopplerProperTimeComponentScaling > transmitterProperTimePartials =
           createDopplerProperTimePartials( transmitterDopplerProperTimeInterface, oneWayDopplerLinkEnds,
                                            observation_models::transmitter );
    boost::shared_ptr< OneWayDopplerProperTimeComponentScaling > receiverProperTimePartials =
            createDopplerProperTimePartials( receiverDopplerProperTimeInterface, oneWayDopplerLinkEnds,
                                             observation_models::receiver  );

    boost::shared_ptr< OneWayDopplerScaling > oneWayDopplerScaling = boost::make_shared< OneWayDopplerScaling >(
            boost::bind( &linear_algebra::evaluateSecondBlockInStateVector, transmitterNumericalStateDerivativeFunction, _1 ),
            boost::bind( &linear_algebra::evaluateSecondBlockInStateVector, receiverNumericalStateDerivativeFunction, _1 ),
                transmitterProperTimePartials,
                receiverProperTimePartials );


    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    SingleLinkObservationPartialList dopplerPartials;

    std::vector< boost::shared_ptr< estimatable_parameters::EstimatableParameter<
            Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            parametersToEstimate->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        std::string acceleratedBody;
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else
        {
            throw std::runtime_error( "Error when making one way doppler partials, could not identify parameter " +
                       boost::lexical_cast< std::string >(
                                          initialDynamicalParameters.at( i )->getParameterName( ).first ) );
        }


        // Create position one-way doppler partial for current body
        boost::shared_ptr< ObservationPartial< 1 > > currentDopplerPartialWrtPosition = createOneWayDopplerPartialWrtBodyState(
                    oneWayDopplerLinkEnds, bodyMap, acceleratedBody, oneWayDopplerScaling, lightTimeCorrectionPartialObjects );

        // Check if partial is non-null (i.e. whether dependency exists between current doppler and current body)
        if( currentDopplerPartialWrtPosition != NULL )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            dopplerPartials[ currentPair ] = currentDopplerPartialWrtPosition;
        }

        // Increment current index by size of body initial state (6).
        currentIndex += 6;
    }

    // Iterate over all double parameters that are to be estimated.
    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator =
         doubleParametersToEstimate.begin( ); parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position one-way doppler partial for current parameter
        boost::shared_ptr< ObservationPartial< 1 > > currentDopplerPartial = createOneWayDopplerPartialWrtParameter(
                    oneWayDopplerLinkEnds, bodyMap, parameterIterator->second,
                    oneWayDopplerScaling, lightTimeCorrectionPartialObjects );

        // Check if partial is non-null (i.e. whether dependency exists between current doppler and current parameter)
        if( currentDopplerPartial != NULL )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            dopplerPartials[ currentPair ] = currentDopplerPartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate =
            parametersToEstimate->getVectorParameters( );
    for( std::map< int, boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator =
         vectorParametersToEstimate.begin( ); parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        boost::shared_ptr< ObservationPartial< 1 > > currentDopplerPartial;
        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentDopplerPartial = createOneWayDopplerPartialWrtParameter(
                        oneWayDopplerLinkEnds, bodyMap, parameterIterator->second, oneWayDopplerScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentDopplerPartial = createObservationPartialWrtLinkProperty< 1 >(
                        oneWayDopplerLinkEnds, observation_models::one_way_doppler, parameterIterator->second );
        }


        // Check if partial is non-null (i.e. whether dependency exists between current doppler and current parameter)
        if( currentDopplerPartial != NULL )
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
 *  \param bodyMap List of all bodies, for creating one-way doppler partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary one-way doppler partials of a single link end,
 *  and OneWayDopplerScaling, object, used for scaling the position partial members of all OneWayDopplerPartials in link end.
 */
template< typename ObservationScalarType, typename TimeType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
boost::shared_ptr< PositionPartialScaling > > > createOneWayDopplerPartials(
        const std::map< observation_models::LinkEnds,
        boost::shared_ptr< observation_models::ObservationModel< 1, ObservationScalarType, TimeType > > > observationModelList,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< estimatable_parameters::EstimatableParameterSet< ObservationScalarType > > parametersToEstimate )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationPartialList,
            boost::shared_ptr< PositionPartialScaling > > > dopplerPartials;

    // Iterate over all link ends.
    for( auto linkIterator: observationModelList )
    {
        // Check if required link end types are present
        if( ( linkIterator.first.count( observation_models::receiver ) == 0 ) ||
                ( linkIterator.first.count( observation_models::transmitter ) == 0 ) )
        {
            throw std::runtime_error( "Error when making 1-way doppler partials, did not find both receiver and transmitter in link ends" );

        }


         boost::shared_ptr< observation_models::OneWayDopplerObservationModel< ObservationScalarType, TimeType > >
                 dopplerObservationModel =
                 boost::dynamic_pointer_cast< observation_models::OneWayDopplerObservationModel
                 < ObservationScalarType, TimeType > >( linkIterator.second );

        std::vector< boost::shared_ptr< observation_models::LightTimeCorrection > > singleLinkLightTimeCorrections;
        if( dopplerObservationModel == NULL )
        {
            throw std::runtime_error( "Error when making one-way Doppler partials. Type is inconsistent" );
        }
        else
        {
            singleLinkLightTimeCorrections = dopplerObservationModel->getLightTimeCalculator( )->getLightTimeCorrection( );
        }

        // Create doppler partials for current link ends
        dopplerPartials[ linkIterator.first ] = createOneWayDopplerPartials< ObservationScalarType >(
                    linkIterator.first, bodyMap, parametersToEstimate,
                    dopplerObservationModel->getTransmitterProperTimeRateCalculator( ),
                    dopplerObservationModel->getReceiverProperTimeRateCalculator( ),
                    singleLinkLightTimeCorrections );
    }

    // Return complete set of link ends.
    return dopplerPartials;
}

}

}

#endif // TUDAT_CREATEONEWAYDOPPLERPARTIALS_H

