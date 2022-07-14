/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEDIRECTOBSERVATIONPARTIALS_H
#define TUDAT_CREATEDIRECTOBSERVATIONPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/interpolator.h"

#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayLinkObservationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayRangePartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerPartial.h"
#include "tudat/astro/observation_models.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"
#include "tudat/simulation/estimation_setup/createPositionPartialScaling.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace observation_partials
{



//! Function to generate observation partial wrt a single  parameter.
/*!
 *  Function to generate observation partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param oneWayLinkEnds Link ends (transmitter and receiever) for which observation partials are to be calculated
 *  (i.e. for which  observation observations are to be processed).
 *  \param bodies List of all bodies, for creating observation partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param positionPartialScaler Object scale position partials to observation partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return observation partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType, int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtParameter(
        const observation_models::LinkEnds oneWayLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    std::shared_ptr< ObservationPartial< ObservationSize > > observationPartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( oneWayLinkEnds, bodies, parameterToEstimate );

        // Create observation partials if any position partials are created (i.e. if any dependency exists).
        std::shared_ptr< DirectObservationPartial< ObservationSize > > testObservationPartial  =
                std::make_shared< DirectObservationPartial< ObservationSize > >(
                    positionPartialScaler, positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );
        if( positionPartials.size( ) > 0 || testObservationPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) > 0 ||
                testObservationPartial->useLinkIndependentPartials())
        {
            observationPartial = testObservationPartial;
        }
    }

    // Return observation partial object (nullptr if no dependency exists).
    return observationPartial;
}

//! Function to generate observation partial wrt a position of a body.
/*!
 *  Function to generate observation partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param oneWayLinkEnds Link ends (transmitter and receiever) for which observation partials are to be calculated
 *  (i.e. for which observation observations are to be processed).
 *  \param bodies List of all bodies, for creating observation partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param positionPartialScaler Object scale position partials to observation partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return observation partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtBodyPosition(
        const observation_models::LinkEnds oneWayLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >( ) )
{
    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyState( oneWayLinkEnds, bodies, bodyToEstimate );

    // Create observation partials if any position partials are created (i.e. if any dependency exists).
    std::shared_ptr< DirectObservationPartial< ObservationSize > > observationPartial;
    if( positionPartials.size( ) > 0 )
    {
        observationPartial = std::make_shared< DirectObservationPartial< ObservationSize > >(
                    positionPartialScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    // Return observation partial object (nullptr if no dependency exists).
    return observationPartial;
}

//! Function to generate observation partial wrt a rotational state of a body.
/*!
 *  Function to generate observation partial wrt a rotational state of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param oneWayLinkEnds Link ends (transmitter and receiever) for which observation partials are to be calculated
 *  (i.e. for which observation observations are to be processed).
 *  \param bodies List of all bodies, for creating observation partial.
 *  \param bodyToEstimate Name of body wrt rotational state of which a partial is to be created.
 *  \param positionPartialScaler Object scale position partials to observation partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return observation partial object wrt a current rotational state of a body (is nullptr if no parameter dependency exists).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtBodyRotationalState(
        const observation_models::LinkEnds oneWayLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > >&
        lightTimeCorrectionPartialObjects  )
{
    // Create position partials of link ends for current body position
    std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
            createCartesianStatePartialsWrtBodyRotationalState( oneWayLinkEnds, bodies, bodyToEstimate );

    // Create observation partials if any position partials are created (i.e. if any dependency exists).
    std::shared_ptr< DirectObservationPartial< ObservationSize > > observationPartial;
    if( positionPartials.size( ) > 0 )
    {
        observationPartial = std::make_shared< DirectObservationPartial< ObservationSize > >(
                    positionPartialScaler, positionPartials, std::make_pair(
                        estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
                    lightTimeCorrectionPartialObjects );
    }

    // Return observation partial object (nullptr if no dependency exists).
    return observationPartial;
}


//! Function to generate observation partials and associated scaler for single link ends.
/*!
 *  Function to generate observation partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends set.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter and receiever linkEndType) that are to be used.
 *  \param oneWayLinkEnds Link ends (transmitter and receiever) for which observation partials are to be calculated
 *  (i.e. for which observation observations are to be processed).
 *  \param bodies List of all bodies, for creating observation partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary observation partials of a single link end, and ObservationPartial< ObservationSize >, object, used for
 *  scaling the position partial members of all ObservationPartials in link end.
 */
template< typename ParameterType, int ObservationSize, typename TimeType = double  >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >,
std::shared_ptr< PositionPartialScaling > >
createSingleLinkObservationPartials(
        const observation_models::LinkEnds oneWayLinkEnds,
        const observation_models::ObservableType observableType,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >& lightTimeCorrections =
        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > >( ),
        const bool useBiasPartials = true,
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > observationModel = nullptr )
{
    std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > lightTimeCorrectionPartialObjects;

    if( lightTimeCorrections.size( ) > 0 )
    {
        lightTimeCorrectionPartialObjects = observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections );
    }

    // Create scaling object, to be used for all observation partials in current link end.
    //    std::shared_ptr< ObservationPartialScalingCreator< ObservationSize > > observationPartialScalingCreator =
    //            std::make_shared< ObservationPartialScalingCreator< ObservationSize > >;
    std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionScaling =
            ObservationPartialScalingCreator< ObservationSize >::createPositionScalingObject(
                oneWayLinkEnds, observableType, bodies, observationModel );

    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > observationPartials;

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

            // Create position observation partial for current body
            std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial =
                    createObservationPartialWrtBodyPosition< ObservationSize >(
                        oneWayLinkEnds, bodies, acceleratedBody, positionScaling, lightTimeCorrectionPartialObjects );

            // Check if partial is non-null (i.e. whether dependency exists between current observable and current body)
            if( currentObservationPartial != nullptr )
            {
                // Add partial to the list.
                currentPair = std::pair< int, int >( currentIndex, 6 );
                observationPartials[ currentPair ] = currentObservationPartial;
            }

            // Increment current index by size of body initial state (6).
            currentIndex += 6;
        }
        else if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_rotational_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;

            // Create position observation partial for current body
            std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial =
                    createObservationPartialWrtBodyRotationalState< ObservationSize >(
                        oneWayLinkEnds, bodies, acceleratedBody, positionScaling, lightTimeCorrectionPartialObjects );

            // Check if partial is non-null (i.e. whether dependency exists between current observable and current body)
            if( currentObservationPartial != nullptr )
            {
                // Add partial to the list.
                currentPair = std::pair< int, int >( currentIndex, 7 );
                observationPartials[ currentPair ] = currentObservationPartial;
            }

            // Increment current index by size of body initial state (6).
            currentIndex += 7;
        }
        else
        {
            throw std::runtime_error( "Error when making observation partials, could not identify parameter " +
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
        // Create position observation partial for current parameter
        std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial =
                createObservationPartialWrtParameter< double, ObservationSize >(
                    oneWayLinkEnds, bodies, parameterIterator->second,
                    positionScaling, lightTimeCorrectionPartialObjects );


        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentObservationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            observationPartials[ currentPair ] = currentObservationPartial;
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
        // Create observation partial for current parameter
        std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentObservationPartial = createObservationPartialWrtParameter< Eigen::VectorXd, ObservationSize >(
                        oneWayLinkEnds, bodies, parameterIterator->second, positionScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentObservationPartial = createObservationPartialWrtLinkProperty< ObservationSize >(
                        oneWayLinkEnds, observableType, parameterIterator->second, useBiasPartials );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentObservationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            observationPartials[ currentPair ] = currentObservationPartial;
        }
    }

    // Return complete set of partials and scaling object.
    return std::make_pair( observationPartials, positionScaling );
}

template< typename ParameterType, int ObservationSize, typename TimeType = double >
std::map< observation_models::LinkEnds,
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > > createSingleLinkObservationPartialsList(
        const std::vector< observation_models::LinkEnds >& linkEndsList,
        const observation_models::ObservableType observableType,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        std::map< observation_models::LinkEnds, std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >
        lightTimeCorrections =
        std::map< observation_models::LinkEnds, std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ),
        const bool useBiasPartials = true,
        const std::vector< std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > > observationModelList =
        std::vector< std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > >( ) )
{
    bool useObservationModel = false;
    if( observationModelList.size( ) != 0 )
    {
        if( linkEndsList.size( ) != observationModelList.size( ) )
        {
            throw std::runtime_error( "Error when creating observaion partials, list of observation models provided, but size is incompatible." );
        }
        useObservationModel = true;
    }
    // Declare return list.
    std::map< observation_models::LinkEnds,
            std::pair<
            std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >,
            std::shared_ptr< PositionPartialScaling > >
            > observationPartials;

    // Iterate over all link ends.
    std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > currentObservationModel = nullptr;
    for( unsigned int i = 0; i < linkEndsList.size( ); i++ )
    {
        observation_models::LinkEnds linkEnds = linkEndsList.at( i );

        std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > singleLinkLightTimeCorrections;
        if( lightTimeCorrections.count( linkEnds ) > 0 )
        {
            if( lightTimeCorrections.at( linkEnds ).size( ) > 1 )
            {
                std::cerr << "Error when making observation partials, light time corrections for "
                          << lightTimeCorrections.at( linkEnds ).size( ) << " links found" << std::endl;
            }
            else if( lightTimeCorrections.at( linkEnds ).size( ) == 1 )
            {
                singleLinkLightTimeCorrections = lightTimeCorrections.at( linkEnds ).at( 0 );
            }
        }

        if( useObservationModel == true )
        {
            currentObservationModel = observationModelList.at( i );
        }

        // Create observation partials for current link ends
        observationPartials[ linkEnds ] = createSingleLinkObservationPartials< ParameterType, ObservationSize >(
                    linkEnds, observableType, bodies, parametersToEstimate,
                    singleLinkLightTimeCorrections, useBiasPartials, currentObservationModel );
    }

    // Return complete set of link ends.
    return observationPartials;
}

//! Function to generate observation partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate observation partials for all parameters that are to be estimated, for all sets of link ends.
 *  The observation partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter and receiever linkEndType)
 *  that are to be used.
 *  \param bodies List of all bodies, for creating observation partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \param useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation bias parameters
 *  \return Map of SingleLinkObservationPartialList, representing all necessary observation partials of a single link end,
 *  and ObservationPartial< ObservationSize >, object, used for scaling the position partial members of all ObservationPartials in link end.
 */
template< typename ParameterType, typename TimeType, int ObservationSize >
std::map< observation_models::LinkEnds,
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > ,
std::shared_ptr< PositionPartialScaling > > > createSingleLinkObservationPartialsList(
        const std::map< observation_models::LinkEnds,
        std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > > observationModelList,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool useBiasPartials = true )
{

    std::map< observation_models::LinkEnds, std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >
            lightTimeCorrections = getLightTimeCorrectionsList( observationModelList );

    std::vector< observation_models::LinkEnds > linkEndsList = utilities::createVectorFromMapKeys( observationModelList );

    observation_models::ObservableType observableType = observation_models::undefined_observation_model;
    // Iterate over all link ends.
    for( auto it : observationModelList )
    {
        if( observableType == observation_models::undefined_observation_model )
        {
            observableType = it.second->getObservableType( );
        }
        else if( observableType != it.second->getObservableType( ) )
        {
            throw std::runtime_error( "Error when creating single link observation partials, input models are inconsistent" );
        }
    }

    // Return complete set of link ends.
    return createSingleLinkObservationPartialsList< ParameterType, ObservationSize >(
                linkEndsList, observableType, bodies, parametersToEstimate, lightTimeCorrections, useBiasPartials,
                utilities::createVectorFromMapValues( observationModelList ) );
}

}

}


#endif // TUDAT_CREATEDIRECTOBSERVATIONPARTIALS_H
