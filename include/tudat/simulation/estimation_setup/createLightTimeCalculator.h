/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATELIGHTTIMECALCULATOR_H
#define TUDAT_CREATELIGHTTIMECALCULATOR_H

#include "tudat/astro/ephemerides/compositeEphemeris.h"
#include "tudat/astro/observation_models/lightTimeSolution.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/createLightTimeCorrection.h"
#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"

namespace tudat
{
namespace observation_models
{

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
    const std::shared_ptr< ephemerides::Ephemeris >& transmitterCompleteEphemeris,
    const std::shared_ptr< ephemerides::Ephemeris >& receiverCompleteEphemeris,
    const simulation_setup::SystemOfBodies& bodies,
    const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections,
    const LinkEnds& linkEnds,
    const LinkEndType& transmittingLinkEndType,
    const LinkEndType& receivingLinkEndType,
    const ObservableType observableType,
    const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
    = std::make_shared< LightTimeConvergenceCriteria >( )  )
{
    std::vector< std::shared_ptr< LightTimeCorrection > > lightTimeCorrectionFunctions;

    // Create lighttime correction functions from lightTimeCorrections
    for( unsigned int i = 0; i < lightTimeCorrections.size( ); i++ )
    {
        std::shared_ptr< LightTimeCorrection > lightTimeCorrectionFunction = createLightTimeCorrections(
            lightTimeCorrections[ i ], bodies, linkEnds, transmittingLinkEndType,
            receivingLinkEndType, observableType );

        if ( lightTimeCorrectionFunction != nullptr )
        {
            lightTimeCorrectionFunctions.push_back( lightTimeCorrectionFunction );
        }
    }

    // Create light time calculator.
    return std::make_shared< LightTimeCalculator< ObservationScalarType, TimeType > >
        ( transmitterCompleteEphemeris, receiverCompleteEphemeris, lightTimeCorrectionFunctions,
          lightTimeConvergenceCriteria );
}

//! Function to create a light-time calculation object
/*!
 *  Function to create a light-time calculation object from light time correction settings environment and link end
 *  state functions.
 *  \param transmitterCompleteEphemeris Function returning the transmitter Cartesian state as a function of time.
 *  \param receiverCompleteEphemeris Function returning the receiver Cartesian state as a function of time.
 *  \param bodies List of body objects that comprises the environment
 *  \param lightTimeCorrections List of light time corrections (w.r.t. Euclidean distance) that are applied when computing
 *  light time.
 *  \param transmittingLinkEnd Identifier for transmitting link end.
 *  \param receivingLinkEnd Identifier for receiving link end.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
        const std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) >& transmitterCompleteEphemeris,
        const std::function< Eigen::Matrix< ObservationScalarType, 6, 1 >( const TimeType ) >& receiverCompleteEphemeris,
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections,
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType,
        const ObservableType observableType,
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( )  )
{
    return createLightTimeCalculator< ObservationScalarType, TimeType >(
        std::make_shared< ephemerides::CustomEphemeris< TimeType, ObservationScalarType > >( transmitterCompleteEphemeris, bodies.getFrameOrigin( ), bodies.getFrameOrientation( ) ),
        std::make_shared< ephemerides::CustomEphemeris< TimeType, ObservationScalarType > >( transmitterCompleteEphemeris, bodies.getFrameOrigin( ), bodies.getFrameOrientation( ) ),
            bodies, lightTimeCorrections, linkEnds, transmittingLinkEndType, receivingLinkEndType, observableType, lightTimeConvergenceCriteria );
}

//! Function to create a light-time calculation object
/*!
 *  Function to create a light-time calculation object from light time correction settings environment and link end
 *  identifiers.
 *  \param transmittingLinkEnd Identifier for transmitting link end.
 *  \param receivingLinkEnd Identifier for receiving link end.
 *  \param bodies List of body objects that comprises the environment
 *  \param lightTimeCorrections List of light time corrections (w.r.t. Euclidean distance) that are applied when computing
 *  light time.
 */
template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< observation_models::LightTimeCalculator< ObservationScalarType, TimeType > >
createLightTimeCalculator(
        const LinkEnds& linkEnds,
        const LinkEndType& transmittingLinkEndType,
        const LinkEndType& receivingLinkEndType,
        const simulation_setup::SystemOfBodies& bodies,
        const ObservableType observableType = undefined_observation_model,
        const std::vector< std::shared_ptr< LightTimeCorrectionSettings > >& lightTimeCorrections =
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ),
        const std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ) )
{

    // Get link end state functions and create light time calculator.
    return createLightTimeCalculator< ObservationScalarType, TimeType >(
                simulation_setup::getLinkEndCompleteEphemeris< TimeType, ObservationScalarType >(
                    linkEnds.at( transmittingLinkEndType ), bodies ),
                simulation_setup::getLinkEndCompleteEphemeris< TimeType, ObservationScalarType >(
                    linkEnds.at( receivingLinkEndType ), bodies ),
                bodies, lightTimeCorrections, linkEnds, transmittingLinkEndType, receivingLinkEndType,
                observableType, lightTimeConvergenceCriteria );
}

template< typename ObservationScalarType = double, typename TimeType = double >
std::shared_ptr< MultiLegLightTimeCalculator< ObservationScalarType, TimeType > > createMultiLegLightTimeCalculator(
        const LinkEnds& linkEnds,
        const simulation_setup::SystemOfBodies& bodies,
        const ObservableType observableType = undefined_observation_model,
        const std::vector< std::vector< std::shared_ptr< LightTimeCorrectionSettings > > >& lightTimeCorrections
            = std::vector< std::vector< std::shared_ptr< LightTimeCorrectionSettings > > >( ),
        const std::vector< std::shared_ptr< LightTimeConvergenceCriteria > >& singleLegsLightTimeConvergenceCriteria
            = std::vector< std::shared_ptr< LightTimeConvergenceCriteria > >( ),
        const std::shared_ptr< LightTimeConvergenceCriteria > multiLegLightTimeConvergenceCriteria
            = std::make_shared< LightTimeConvergenceCriteria >( ) )
{

    // Check if link ends contain a receiver and a transmitter
    if( linkEnds.count( receiver ) == 0 )
    {
        throw std::runtime_error( "Error when making multi-leg light time calculator, no receiver found" );
    }
    if( linkEnds.count( transmitter ) == 0 )
    {
        throw std::runtime_error( "Error when making multi-leg light time calculator, no transmitter found" );
    }

    // Check link end consistency.
    for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( );
         linkEndIterator++ )
    {
        if( ( linkEndIterator->first != transmitter ) && ( linkEndIterator->first != receiver ) )
        {
            int linkEndIndex = static_cast< int >( linkEndIterator->first );
            LinkEndType previousLinkEndType = static_cast< LinkEndType >( linkEndIndex - 1 );

            if( linkEnds.count( previousLinkEndType ) == 0 )
            {
                throw std::runtime_error( "Error when making multi-leg light time calculator, did not find link end type " +
                                          std::to_string( previousLinkEndType ) );
            }
        }
    }

    // Check consistency of convergence criteria size, corrections size, and number of link ends
    unsigned int numberOfLinks = linkEnds.size( ) - 1;
    if ( !( lightTimeCorrections.size( ) == numberOfLinks || lightTimeCorrections.empty( ) ) )
    {
        throw std::runtime_error(
                "Error when making multi-leg light time calculator: size of single-leg light time corrections (" +
                std::to_string( lightTimeCorrections.size( ) ) +
                ") are inconsistent with number of links (" + std::to_string( numberOfLinks ) + ")." );
    }
    else if ( !( singleLegsLightTimeConvergenceCriteria.size( ) == numberOfLinks || singleLegsLightTimeConvergenceCriteria.empty( ) ) )
    {
        throw std::runtime_error(
                "Error when making multi-leg light time calculator: size of single-leg convergence criteria (" +
                std::to_string( singleLegsLightTimeConvergenceCriteria.size( ) ) +
                ") are inconsistent with number of links (" + std::to_string( numberOfLinks ) + ")." );
    }

    // Define light-time calculator list
    std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > lightTimeCalculators;
    // Boolean indicating whether any of the used light time corrections require multi-leg iterations
    bool multiLegIterationsRequired = false;
    // Iterate over all link ends and create light-time calculators
    LinkEnds::const_iterator transmitterIterator = linkEnds.begin( );
    LinkEnds::const_iterator receiverIterator = linkEnds.begin( );
    receiverIterator++;
    for( unsigned int i = 0; i < linkEnds.size( ) - 1; i++ )
    {
        // Get convergence criteria of current leg
        std::shared_ptr< LightTimeConvergenceCriteria > currentConvergenceCriteria;
        if ( singleLegsLightTimeConvergenceCriteria.empty( ) )
        {
            currentConvergenceCriteria = std::make_shared< LightTimeConvergenceCriteria >( );
        }
        else
        {
            currentConvergenceCriteria = singleLegsLightTimeConvergenceCriteria.at( i );
        }

        // Get light time corrections of current leg
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > currentLightTimeCorrections;
        if ( lightTimeCorrections.empty( ) )
        { }
        else
        {
            currentLightTimeCorrections = lightTimeCorrections.at( i );
        }

        lightTimeCalculators.push_back(
                createLightTimeCalculator< ObservationScalarType, TimeType >(
                        linkEnds, transmitterIterator->first, receiverIterator->first,
                        bodies, observableType, currentLightTimeCorrections,
                        currentConvergenceCriteria ) );

        // Check whether any of the corrections requires multi-leg iterations
        if ( !multiLegIterationsRequired )
        {
            for ( unsigned j = 0; j < currentLightTimeCorrections.size( ); ++j )
            {
                if ( requiresMultiLegIterations( currentLightTimeCorrections.at( j )->getCorrectionType( ) ) )
                {
                    multiLegIterationsRequired = true;
                    break;
                }
            }
        }

        transmitterIterator++;
        receiverIterator++;
    }

    // Create multi-leg light time calculator
    std::shared_ptr< observation_models::MultiLegLightTimeCalculator< ObservationScalarType, TimeType > >
            multiLegLightTimeCalculator = std::make_shared< observation_models::MultiLegLightTimeCalculator<
                    ObservationScalarType, TimeType > >(
                            lightTimeCalculators, multiLegLightTimeConvergenceCriteria, multiLegIterationsRequired );

    return multiLegLightTimeCalculator;
}


} // namespace observation_models

} // namespace tudat

#endif // TUDAT_CREATELIGHTTIMECALCULATOR_H
