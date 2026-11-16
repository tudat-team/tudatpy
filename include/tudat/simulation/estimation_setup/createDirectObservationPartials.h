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

#include <Eigen/Core>
#include <map>
#include <memory>
#include <vector>

#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"
#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationModel.h"
#include "tudat/astro/observation_models/oneWayDopplerObservationModel.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/orbit_determination/observation_partials/angularPositionPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/directObservationPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayDopplerPartial.h"
#include "tudat/astro/orbit_determination/observation_partials/oneWayRangePartial.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/simulation/estimation_setup/createClockPartials.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"
#include "tudat/simulation/estimation_setup/createObservationBiasPartial.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/createPositionPartialScaling.h"

namespace tudat
{

namespace observation_partials
{

//! Function to create a direct observation partial w.r.t. a parameter, if the observable depends on that parameter.
/*!
 *  Function to create a direct observation partial w.r.t. a parameter, if the observable depends on that parameter: through
 *  the link-end positions, through the light-time corrections, or directly (as defined by the position partial scaling).
 *  \param positionPartialScaler Object to scale position partials to observation partials for current link ends.
 *  \param positionPartials Partials of the link-end positions w.r.t. the parameter.
 *  \param parameterIdentifier Identifier of the parameter w.r.t. which the partial is to be created.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used.
 *  \return Observation partial object w.r.t. the parameter (nullptr if no dependency exists).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createDirectObservationPartial(
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > >& positionPartials,
        const estimatable_parameters::EstimatebleParameterIdentifier& parameterIdentifier,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
                lightTimeCorrectionPartialObjects )
{
    std::shared_ptr< DirectObservationPartial< ObservationSize > > observationPartial =
            std::make_shared< DirectObservationPartial< ObservationSize > >(
                    positionPartialScaler, positionPartials, parameterIdentifier, lightTimeCorrectionPartialObjects );
    if( positionPartials.size( ) > 0 || observationPartial->getNumberOfLighTimeCorrectionPartialsFunctions( ) > 0 ||
        observationPartial->useLinkIndependentPartials( ) )
    {
        return observationPartial;
    }
    return nullptr;
}

//! Function to generate observation partial wrt a single  parameter.
/*!
 *  Function to generate observation partial wrt a single  parameter, for a single link ends (which
 * must contain a transmitter and receiever  linkEndType). \tparam ParameterType Type of parameter
 * (double for size 1, VectorXd for larger size). \param oneWayLinkEnds Link ends (transmitter and
 * receiever) for which observation partials are to be calculated (i.e. for which  observation
 * observations are to be processed). \param bodies List of all bodies, for creating observation
 * partial. \param parameterToEstimate Object of current parameter that is to be estimated. \param
 * positionPartialScaler Object scale position partials to observation partials for current link
 * ends. \param lightTimeCorrectionPartialObjects List of light time correction partials to be used
 * (empty by default) \return observation partial object wrt a single parameter (is nullptr if no
 * parameter dependency exists).
 */
template< typename ParameterType, int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtParameter(
        const observation_models::LinkEnds oneWayLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
                lightTimeCorrectionPartialObjects = {} )
{
    return createDirectObservationPartial< ObservationSize >(
            positionPartialScaler,
            createCartesianStatePartialsWrtParameter( oneWayLinkEnds, bodies, parameterToEstimate ),
            parameterToEstimate->getParameterName( ),
            lightTimeCorrectionPartialObjects );
}

//! Function to generate observation partial wrt a position of a body.
/*!
 *  Function to generate observation partial wrt a position of a body, for a single link ends (which
 * must contain a transmitter and receiever  linkEndType). \param oneWayLinkEnds Link ends
 * (transmitter and receiever) for which observation partials are to be calculated (i.e. for which
 * observation observations are to be processed). \param bodies List of all bodies, for creating
 * observation partial. \param bodyToEstimate Name of body wrt position of which a partial is to be
 * created. \param positionPartialScaler Object scale position partials to observation partials for
 * current link ends. \param lightTimeCorrectionPartialObjects List of light time correction
 * partials to be used (empty by default) \param bodiesWithEstimatedState Names of all bodies with an estimated
 * translational state (see createCartesianStatePartialsWrtBodyState) \return observation partial object wrt a current position
 * of a body (is nullptr if no parameter dependency exists).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtBodyPosition(
        const observation_models::LinkEnds oneWayLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
                lightTimeCorrectionPartialObjects = {},
        const std::vector< std::string >& bodiesWithEstimatedState = std::vector< std::string >( ) )
{
    return createDirectObservationPartial< ObservationSize >(
            positionPartialScaler,
            createCartesianStatePartialsWrtBodyState( oneWayLinkEnds, bodies, bodyToEstimate, bodiesWithEstimatedState ),
            std::make_pair( estimatable_parameters::initial_body_state, std::make_pair( bodyToEstimate, "" ) ),
            lightTimeCorrectionPartialObjects );
}

//! Function to generate observation partial wrt a rotational state of a body.
/*!
 *  Function to generate observation partial wrt a rotational state of a body, for a single link
 * ends (which must contain a transmitter and receiever  linkEndType). \param oneWayLinkEnds Link
 * ends (transmitter and receiever) for which observation partials are to be calculated (i.e. for
 * which observation observations are to be processed). \param bodies List of all bodies, for
 * creating observation partial. \param bodyToEstimate Name of body wrt rotational state of which a
 * partial is to be created. \param positionPartialScaler Object scale position partials to
 * observation partials for current link ends. \param lightTimeCorrectionPartialObjects List of
 * light time correction partials to be used (empty by default) \return observation partial object
 * wrt a current rotational state of a body (is nullptr if no parameter dependency exists).
 */
template< int ObservationSize >
std::shared_ptr< ObservationPartial< ObservationSize > > createObservationPartialWrtBodyRotationalState(
        const observation_models::LinkEnds oneWayLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionPartialScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
                lightTimeCorrectionPartialObjects )
{
    return createDirectObservationPartial< ObservationSize >(
            positionPartialScaler,
            createCartesianStatePartialsWrtBodyRotationalState( oneWayLinkEnds, bodies, bodyToEstimate ),
            std::make_pair( estimatable_parameters::initial_rotational_body_state, std::make_pair( bodyToEstimate, "" ) ),
            lightTimeCorrectionPartialObjects );
}

//! Function to generate observation partials and associated scaler for single link ends.
/*!
 *  Function to generate observation partials and associated scaler for all parameters that are to
 * be estimated, for a single link ends set. The set of parameters and bodies that are to be
 * estimated, as well as the set of link ends (each of which must contain a transmitter and
 * receiever linkEndType) that are to be used. \param oneWayLinkEnds Link ends (transmitter and
 * receiever) for which observation partials are to be calculated (i.e. for which observation
 * observations are to be processed). \param bodies List of all bodies, for creating observation
 * partials. \param parametersToEstimate Set of parameters that are to be estimated \param
 * lightTimeCorrections List of light time correction partials to be used (empty by default) \param
 * useBiasPartials Boolean to denote whether this function should create partials w.r.t. observation
 * bias parameters \return Set of observation partials with associated indices in complete vector of
 * parameters that are estimated, representing all  necessary observation partials of a single link
 * end, and ObservationPartial< ObservationSize >, object, used for scaling the position partial
 * members of all ObservationPartials in link end.
 */
template< typename ParameterType, int ObservationSize, typename TimeType = double >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > >,
           std::shared_ptr< PositionPartialScaling > >
createSingleLinkObservationPartials(
        const std::shared_ptr< observation_models::ObservationModel< ObservationSize, ParameterType, TimeType > > observationModel,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const bool isPartialForDifferencedObservable = false,
        const bool isPartialForConcatenatedObservable = false )
{
    observation_models::LinkEnds oneWayLinkEnds = observationModel->getLinkEnds( );
    observation_models::ObservableType observableType = observationModel->getObservableType( );

    // Create light-time correction partials for each light-time leg
    std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > > lightTimeCorrectionPartialObjects;
    for( auto lightTimeCorrections : observation_models::getLightTimeCorrections( observationModel ) )
    {
        lightTimeCorrectionPartialObjects.push_back( observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections ) );
    }

    // Create scaling object, to be used for all observation partials in current link end.
    //    std::shared_ptr< ObservationPartialScalingCreator< ObservationSize > >
    //    observationPartialScalingCreator =
    //            std::make_shared< ObservationPartialScalingCreator< ObservationSize > >;
    std::shared_ptr< DirectPositionPartialScaling< ObservationSize > > positionScaling =
            ObservationPartialScalingCreator< ObservationSize >::createPositionScalingObject(
                    oneWayLinkEnds, observableType, bodies, observationModel, parametersToEstimate );

    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< ObservationSize > > > observationPartials;

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::Matrix< ParameterType, Eigen::Dynamic, 1 > > > >
            initialDynamicalParameters = parametersToEstimate->getEstimatedInitialStateParameters( );

    // Retrieve all bodies with estimated translational state
    std::vector< std::string > bodiesWithEstimatedState =
            estimatable_parameters::getListOfBodiesToEstimate( parametersToEstimate )[ propagators::translational_state ];

    // Iterate over list of bodies of which the partials of the accelerations acting on them are
    // required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        std::string acceleratedBody;
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::initial_body_state ||
            initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;

            // Create position observation partial for current body
            std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial =
                    createObservationPartialWrtBodyPosition< ObservationSize >( oneWayLinkEnds,
                                                                                bodies,
                                                                                acceleratedBody,
                                                                                positionScaling,
                                                                                lightTimeCorrectionPartialObjects,
                                                                                bodiesWithEstimatedState );

            // Check if partial is non-null (i.e. whether dependency exists between current
            // observable and current body)
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

            // Check if partial is non-null (i.e. whether dependency exists between current
            // observable and current body)
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
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator parameterIterator =
                 doubleParametersToEstimate.begin( );
         parameterIterator != doubleParametersToEstimate.end( );
         parameterIterator++ )
    {
        // Create position observation partial for current parameter
        std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial =
                createObservationPartialWrtParameter< double, ObservationSize >(
                        oneWayLinkEnds, bodies, parameterIterator->second, positionScaling, lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current
        // observable and current parameter)
        if( currentObservationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            observationPartials[ currentPair ] = currentObservationPartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > > vectorParametersToEstimate =
            parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >::iterator parameterIterator =
                 vectorParametersToEstimate.begin( );
         parameterIterator != vectorParametersToEstimate.end( );
         parameterIterator++ )
    {
        // Create observation partial for current parameter
        std::shared_ptr< ObservationPartial< ObservationSize > > currentObservationPartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first ) )
        {
            currentObservationPartial = createObservationPartialWrtParameter< Eigen::VectorXd, ObservationSize >(
                    oneWayLinkEnds, bodies, parameterIterator->second, positionScaling, lightTimeCorrectionPartialObjects );
        }
        else
        {
            std::function< std::shared_ptr< ObservationPartial< ObservationSize > >( const std::string& ) >
                    partialWrtStateCreationFunction = std::bind( &createObservationPartialWrtBodyPosition< ObservationSize >,
                                                                 oneWayLinkEnds,
                                                                 bodies,
                                                                 std::placeholders::_1,
                                                                 positionScaling,
                                                                 lightTimeCorrectionPartialObjects,
                                                                 std::vector< std::string >( ) );
            currentObservationPartial =
                    createObservationPartialWrtLinkProperty< ObservationSize >( oneWayLinkEnds,
                                                                                observableType,
                                                                                parameterIterator->second,
                                                                                bodies,
                                                                                isPartialForDifferencedObservable,
                                                                                isPartialForConcatenatedObservable,
                                                                                observationPartials,
                                                                                partialWrtStateCreationFunction,
                                                                                observationModel->getObservationBiasCalculator( ) );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current
        // observable and current parameter)
        if( currentObservationPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, parameterIterator->second->getParameterSize( ) );
            observationPartials[ currentPair ] = currentObservationPartial;
        }
    }

    // Return complete set of partials and scaling object.
    return std::make_pair( observationPartials, positionScaling );
}

}  // namespace observation_partials

}  // namespace tudat

#endif  // TUDAT_CREATEDIRECTOBSERVATIONPARTIALS_H
