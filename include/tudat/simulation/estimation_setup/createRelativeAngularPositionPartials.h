/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATERELATIVEANGULARPOSITIONPARTIALS_H
#define TUDAT_CREATERELATIVEANGULARPOSITIONPARTIALS_H

#include <vector>
#include <map>

#include <memory>

#include <Eigen/Core>

#include "tudat/math/interpolators/interpolator.h"

#include "tudat/simulation/estimation_setup/createCartesianStatePartials.h"
#include "tudat/astro/orbit_determination/observation_partials/relativeAngularPositionPartial.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrectionPartials.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/initialTranslationalState.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"


namespace tudat
{

namespace observation_partials
{

//! Function to generate relative angular position partial wrt a position of a body.
/*!
 *  Function to generate relative angular position partial wrt a position of a body, for a single link ends (which must contain a
 *  transmitter and receiever  linkEndType).
 *  \param relativeAngularPositionLinkEnds Link ends (transmitter and receiever) for which partials are to be calculated
 *  (i.e. for which relative angular position observation observations are to be processed).
 *  \param bodies List of all bodies, for creating angular position partial.
 *  \param bodyToEstimate Name of body wrt position of which a partial is to be created.
 *  \param relativeAngularPositionScaler Object scale position partials to relative angular position partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Relative angular position partial object wrt a current position of a body (is nullptr if no parameter dependency exists).
 */
std::shared_ptr< RelativeAngularPositionPartial > createRelativeAngularPositionPartialWrtBodyPosition(
        const observation_models::LinkEnds relativeAngularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::string bodyToEstimate,
        const std::shared_ptr< RelativeAngularPositionScaling > relativeAngularPositionScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) );

//! Function to generate relative angular position partial wrt a single  parameter.
/*!
 *  Function to generate relative angular position partial wrt a single  parameter, for a single link ends (which must contain a
 *  transmitter and receiver linkEndType).
 *  \tparam ParameterType Type of parameter (double for size 1, VectorXd for larger size).
 *  \param relativeAngularPositionLinkEnds Link ends (transmitter and receiever) for which relative angular position partials are to be
 *  calculated (i.e. for which  relative angular position observations are to be processed).
 *  \param bodies List of all bodies, for creating angular position partial.
 *  \param parameterToEstimate Object of current parameter that is to be estimated.
 *  \param relativeAngularPositionScaler Object scale position partials to relative angular position partials for current link ends.
 *  \param lightTimeCorrectionPartialObjects List of light time correction partials to be used (empty by default)
 *  \return Relative angular position partial object wrt a single parameter (is nullptr if no parameter dependency exists).
 */
template< typename ParameterType >
std::shared_ptr< RelativeAngularPositionPartial > createRelativeAngularPositionPartialWrtParameter(
        const observation_models::LinkEnds relativeAngularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameter< ParameterType > > parameterToEstimate,
        const std::shared_ptr< RelativeAngularPositionScaling > relativeAngularPositionScaler,
        const std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >&
        lightTimeCorrectionPartialObjects =
        std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > >( ) )
{
    std::shared_ptr< RelativeAngularPositionPartial > relativeAngularPositionPartial;

    {
        std::map< observation_models::LinkEndType, std::shared_ptr< CartesianStatePartial > > positionPartials =
                createCartesianStatePartialsWrtParameter( relativeAngularPositionLinkEnds, bodies, parameterToEstimate );

        std::shared_ptr< RelativeAngularPositionPartial > testRelativeAngularPositionPartial = std::make_shared< RelativeAngularPositionPartial >(
                    relativeAngularPositionScaler, positionPartials, parameterToEstimate->getParameterName( ),
                    lightTimeCorrectionPartialObjects );
        // Create angular position partials if any position partials are created (i.e. if any dependency exists).
        if( positionPartials.size( ) > 0 || testRelativeAngularPositionPartial->getNumberOfLightTimeCorrectionPartialsFunctions( ) )
        {
            relativeAngularPositionPartial = testRelativeAngularPositionPartial;
        }
    }
    return relativeAngularPositionPartial;
}

//! Function to generate relative angular position partials and associated scaler for single link end.
/*!
 *  Function to generate relative angular position partials and associated scaler for all parameters that are to be estimated,
 *  for a single link ends.
 *  The set of parameters and bodies that are to be estimated, as well as the set of link ends
 *  (each of which must contain a transmitter, transmitter2 and receiver linkEndType) that are to be used.
 *  \param relativeAngularPositionLinkEnds Link ends (transmitter, transmitter2 and receiver) for which relative angular position partials are to be
 *  calculated (i.e. for which relative angular position observations are to be processed).
 *  \param bodies List of all bodies, for creating relative angular position partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states of
 *  requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Set of observation partials with associated indices in complete vector of parameters that are estimated,
 *  representing all  necessary relative angular position partials of a single link end, and RelativeAngularPositionScaling, object, used for
 *  scaling the position partial members of all RelativeAngularPositionPartials in link end.
 */
template< typename ParameterType >
std::pair< std::map< std::pair< int, int >, std::shared_ptr< ObservationPartial< 2 > > >,
std::shared_ptr< PositionPartialScaling > >
createRelativeAngularPositionPartials(
        const observation_models::LinkEnds relativeAngularPositionLinkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >& lightTimeCorrections =
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > >( ) )

{
    std::vector< std::vector< std::shared_ptr< observation_partials::LightTimeCorrectionPartial > > > lightTimeCorrectionPartialObjects;
    if( lightTimeCorrections.size( ) > 0 )
    {
        if( lightTimeCorrections.size( ) != 2 )
        {
            throw std::runtime_error( "Error when making relative angular position partials, light time corrections for "
                                      + std::to_string( lightTimeCorrections.size( ) ) + " links found, instead of 2.");
        }
        lightTimeCorrectionPartialObjects.push_back(
                observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections[ 0 ] ) );
        lightTimeCorrectionPartialObjects.push_back(
                observation_partials::createLightTimeCorrectionPartials( lightTimeCorrections[ 1 ] ) );
    }

    // Create scaling object, to be used for all relative angular position partials in current link end.
    std::shared_ptr< RelativeAngularPositionScaling > relativeAngularPositionScaling =
            std::make_shared< RelativeAngularPositionScaling >( );

    SingleLinkObservationTwoPartialList relativeAngularPositionPartials;


    // Initialize vector index variables.
    int currentIndex = 0;
    std::pair< int, int > currentPair = std::pair< int, int >( currentIndex, 1 );

    std::vector< std::shared_ptr< estimatable_parameters::EstimatableParameter<
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
        else if( initialDynamicalParameters.at( i )->getParameterName( ).first == estimatable_parameters::arc_wise_initial_body_state )
        {
            acceleratedBody = initialDynamicalParameters.at( i )->getParameterName( ).second.first;
        }
        else
        {
            throw std::runtime_error( "Error when making relative angular position partials, could not identify parameter" );
        }

        // Create position relative angular position partial for current body
        std::shared_ptr< RelativeAngularPositionPartial > currentRelativeAngularPositionPartial = createRelativeAngularPositionPartialWrtBodyPosition(
                    relativeAngularPositionLinkEnds, bodies, acceleratedBody, relativeAngularPositionScaling,
                    lightTimeCorrectionPartialObjects );

        // Check if partial is non-nullptr (i.e. whether dependency exists between current relative angular position and current body)
        if( currentRelativeAngularPositionPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( currentIndex, 6 );
            relativeAngularPositionPartials[ currentPair ] = currentRelativeAngularPositionPartial;
        }

        // Increment current index by size of body initial state (6).
        currentIndex += 6;
    }

    // Iterate over all double parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > > doubleParametersToEstimate =
            parametersToEstimate->getDoubleParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > >::iterator
         parameterIterator = doubleParametersToEstimate.begin( );
         parameterIterator != doubleParametersToEstimate.end( ); parameterIterator++ )
    {
        // Create position relative angular position partial for current parameter
        std::shared_ptr< RelativeAngularPositionPartial > currentRelativeAngularPositionPartial = createRelativeAngularPositionPartialWrtParameter(
                    relativeAngularPositionLinkEnds, bodies, parameterIterator->second, relativeAngularPositionScaling,
                    lightTimeCorrectionPartialObjects );

        if( currentRelativeAngularPositionPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first, 1 );
            relativeAngularPositionPartials[ currentPair ] = currentRelativeAngularPositionPartial;
        }
    }

    // Iterate over all vector parameters that are to be estimated.
    std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > >
            vectorParametersToEstimate = parametersToEstimate->getVectorParameters( );
    for( std::map< int, std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd  > > >::iterator
         parameterIterator = vectorParametersToEstimate.begin( );
         parameterIterator != vectorParametersToEstimate.end( ); parameterIterator++ )
    {

        // Create position relative angular position partial for current parameter
        std::shared_ptr< ObservationPartial< 2 > > currentRelativeAngularPositionPartial;

        if( !isParameterObservationLinkProperty( parameterIterator->second->getParameterName( ).first )  )
        {
            currentRelativeAngularPositionPartial = createRelativeAngularPositionPartialWrtParameter(
                        relativeAngularPositionLinkEnds, bodies, parameterIterator->second, relativeAngularPositionScaling,
                        lightTimeCorrectionPartialObjects );
        }
        else
        {
            currentRelativeAngularPositionPartial = createObservationPartialWrtLinkProperty< 2 >(
                        relativeAngularPositionLinkEnds, observation_models::relative_angular_position, parameterIterator->second );
        }

        // Check if partial is non-nullptr (i.e. whether dependency exists between current observable and current parameter)
        if( currentRelativeAngularPositionPartial != nullptr )
        {
            // Add partial to the list.
            currentPair = std::pair< int, int >( parameterIterator->first,
                                                 parameterIterator->second->getParameterSize( ) );
            relativeAngularPositionPartials[ currentPair ] = currentRelativeAngularPositionPartial;
        }

    }
    return std::make_pair( relativeAngularPositionPartials, relativeAngularPositionScaling );
}

//! Function to generate relative angular position partials for all parameters that are to be estimated, for all sets of link ends.
/*!
 *  Function to generate relative angular position partials for all parameters that are to be estimated, for all sets of link ends.
 *  The relative angular position partials are generated per set of link ends. The set of parameters and bodies that are to be
 *  estimated, as well as the set of link ends (each of which must contain a transmitter, transmitter2 and receiver linkEndType)
 *  that are to be used.
 *  \param linkEnds Vector of all link ends for which relative angular position partials are to be calculated (i.e. for which one-way
 *  relative angular position observations are  to be processed).
 *  \param bodies List of all bodies, for creating relative angular position partials.
 *  \param parametersToEstimate Set of parameters that are to be estimated (in addition to initial states
 *  of requested bodies)
 *  \param lightTimeCorrections List of light time correction partials to be used (empty by default)
 *  \return Map of SingleLinkObservationPartialList, representing all necessary relative angular position partials of a single link
 *  end, and RelativeAngularPositionScaling, object, used for scaling the position partial members of all RelativeAngularPositionPartials in
 *  link end.
 */
template< typename ParameterType >
std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationTwoPartialList,
std::shared_ptr< PositionPartialScaling > > >
createRelativeAngularPositionPartials(
        const std::vector< observation_models::LinkEnds > linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const std::shared_ptr< estimatable_parameters::EstimatableParameterSet< ParameterType > > parametersToEstimate,
        const std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >& lightTimeCorrections =
        std::map< observation_models::LinkEnds,
        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > >( ) )
{
    // Declare return list.
    std::map< observation_models::LinkEnds, std::pair< SingleLinkObservationTwoPartialList,
            std::shared_ptr< PositionPartialScaling > > > relativeAngularPositionPartials;

    // Iterate over all link ends.
    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        // Check if required link end types are present
        if( ( linkEnds[ i ].count( observation_models::receiver ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter ) == 0 ) ||
                ( linkEnds[ i ].count( observation_models::transmitter2 ) == 0 ) )
        {
            throw std::runtime_error( "Error when making relative angular position partials, did not find receiver, transmitter and transmitter2 in link ends" );

        }

        std::vector< std::vector< std::shared_ptr< observation_models::LightTimeCorrection > > > currentLightTimeCorrections;
        if( lightTimeCorrections.count( linkEnds.at( i ) ) > 0 )
        {
            if( lightTimeCorrections.at( linkEnds.at( i ) ).size( ) != 2 )
            {
                std::cerr << "Error when making relative angular position partials, light time corrections for " <<
                           lightTimeCorrections.at( linkEnds.at( i ) ).size( ) << " links found" << std::endl;
            }
            currentLightTimeCorrections = lightTimeCorrections.at( linkEnds.at( i ) ); //.at( 0 )
        }

        // Create relative angular position partials for current link ends
        relativeAngularPositionPartials[ linkEnds[ i ] ] = createRelativeAngularPositionPartials(
                    linkEnds[ i ], bodies, parametersToEstimate, currentLightTimeCorrections );
    }
    return relativeAngularPositionPartials;
}

}

}

#endif // TUDAT_CREATERELATIVEANGULARPOSITIONPARTIALS_H
