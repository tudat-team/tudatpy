/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>

#include <functional>
#include <boost/make_shared.hpp>


#include "tudat/simulation/estimation_setup/createObservationViability.h"

namespace tudat
{

namespace observation_models
{


//! Function to filter list of observationViabilitySettings, so that only those relevant for single set of link ends are retained
ObservationViabilitySettingsList filterObservationViabilitySettings(
        const ObservationViabilitySettingsList& observationViabilitySettings,
        const LinkEnds& linkEnds )
{
    ObservationViabilitySettingsList filteredViabilitySettings;

    // Iterate over all viability settings
    for( unsigned int i = 0; i < observationViabilitySettings.size( ); i++ )
    {
        // Iterate over all link ends
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            // Check if present viabilitytt setting is relevant
            if( linkEndIterator->second == observationViabilitySettings.at( i )->getAssociatedLinkEnd( ) ||
                    ( ( observationViabilitySettings.at( i )->getAssociatedLinkEnd( ).second == "" ) &&
                      ( observationViabilitySettings.at( i )->getAssociatedLinkEnd( ).first == linkEndIterator->second.first ) ) )
            {
                filteredViabilitySettings.push_back( observationViabilitySettings.at( i ) );
                break;
            }
        }
    }

    return filteredViabilitySettings;
}


//! Function to create an object to check if a minimum elevation angle condition is met for an observation
std::shared_ptr< MinimumElevationAngleCalculator > createMinimumElevationAngleCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings,
        const std::string& stationName )
{
    if( observationViabilitySettings->observationViabilityType_ != minimum_elevation_angle )
    {
        throw std::runtime_error( "Error when making minimum elevation angle calculator, inconsistent input" );
    }

    // If specific link end is specified
    std::string groundStationNameToUse;
    if( observationViabilitySettings->getAssociatedLinkEnd( ).second != "" )
    {
        groundStationNameToUse = observationViabilitySettings->getAssociatedLinkEnd( ).second;
        if( groundStationNameToUse != stationName )
        {
            throw std::runtime_error( "Error when making minimum elevation angle calculator, inconsistent station input" );
        }
    }
    else
    {
        groundStationNameToUse = stationName;
    }

    if( bodies.count( observationViabilitySettings->getAssociatedLinkEnd( ).first ) == 0 )
    {
        throw std::runtime_error( "Error when making minimum elevation angle calculator, body " +
                                  observationViabilitySettings->getAssociatedLinkEnd( ).first + " not found." );
    }

    // Retrieve pointing angles calculator
    std::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAngleCalculator =
            bodies.at( observationViabilitySettings->getAssociatedLinkEnd( ).first )->
            getGroundStation( groundStationNameToUse )->getPointingAnglesCalculator( );

    // Create check object
    double minimumElevationAngle = observationViabilitySettings->getDoubleParameter( );
    return std::make_shared< MinimumElevationAngleCalculator >(
                getLinkStateAndTimeIndicesForLinkEnd(
                    linkEnds,observationType, observationViabilitySettings->getAssociatedLinkEnd( ) ),
                minimumElevationAngle, pointingAngleCalculator );
}

//! Function to create an object to check if a body avoidance angle condition is met for an observation
std::shared_ptr< BodyAvoidanceAngleCalculator > createBodyAvoidanceAngleCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings )
{
    if( observationViabilitySettings->observationViabilityType_ != body_avoidance_angle )
    {
        throw std::runtime_error( "Error when making body avoidance angle calculator, inconsistent input" );
    }

    if( bodies.count( observationViabilitySettings->getStringParameter( ) ) == 0 )
    {
        throw std::runtime_error( "Error when making body avoidance angle calculator, body " +
                                  observationViabilitySettings->getStringParameter( ) + " not found." );
    }

    // Create state function of body to be avoided.
    std::function< Eigen::Vector6d( const double ) > stateFunctionOfBodyToAvoid =
            std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                         bodies.at( observationViabilitySettings->getStringParameter( ) ), std::placeholders::_1 );

    // Create check object
    double bodyAvoidanceAngle = observationViabilitySettings->getDoubleParameter( );
    return std::make_shared< BodyAvoidanceAngleCalculator >(
                getLinkStateAndTimeIndicesForLinkEnd(
                    linkEnds,observationType, observationViabilitySettings->getAssociatedLinkEnd( ) ),
                bodyAvoidanceAngle, stateFunctionOfBodyToAvoid, observationViabilitySettings->getStringParameter( ) );
}

//! Function to create an object to check if a body occultation condition is met for an observation
std::shared_ptr< OccultationCalculator > createOccultationCalculator(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::shared_ptr< ObservationViabilitySettings > observationViabilitySettings )
{
    if( observationViabilitySettings->observationViabilityType_ != body_occultation )
    {
        throw std::runtime_error( "Error when making occultation calculator, inconsistent input" );
    }

    if( bodies.count( observationViabilitySettings->getStringParameter( ) ) == 0 )
    {
        throw std::runtime_error( "Error when making occultation calculator, body " +
                                  observationViabilitySettings->getStringParameter( ) + " not found." );
    }

    // Create state function of occulting body.
    std::function< Eigen::Vector6d( const double ) > stateOfOccultingBody =
            std::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                         bodies.at( observationViabilitySettings->getStringParameter( ) ), std::placeholders::_1 );

    // Create check object
    if( bodies.at( observationViabilitySettings->getStringParameter( ) )->getShapeModel( ) == nullptr )
    {
        throw std::runtime_error( "Error when makig occultation calculator, no shape model found for " +
                                  observationViabilitySettings->getStringParameter( ) );
    }
    double occultingBodyRadius =
            bodies.at( observationViabilitySettings->getStringParameter( ) )->getShapeModel( )->getAverageRadius( );
    return std::make_shared< OccultationCalculator >(
                getLinkStateAndTimeIndicesForLinkEnd(
                    linkEnds, observationType, observationViabilitySettings->getAssociatedLinkEnd( ) ),
                stateOfOccultingBody, occultingBodyRadius );
}

//! Function to create an list of obervation viability conditions for a single set of link ends
std::vector< std::shared_ptr< ObservationViabilityCalculator > > createObservationViabilityCalculators(
        const simulation_setup::SystemOfBodies& bodies,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings )
{
    std::vector< std::shared_ptr< ObservationViabilityCalculator > > linkViabilityCalculators;

    std::vector< std::shared_ptr< ObservationViabilitySettings > > relevantObservationViabilitySettings =
            filterObservationViabilitySettings( observationViabilitySettings, linkEnds );

    for( unsigned int i = 0; i < relevantObservationViabilitySettings.size( ); i++ )
    {

        switch( relevantObservationViabilitySettings.at( i )->observationViabilityType_ )
        {
        case minimum_elevation_angle:
        {
            // Create list of ground stations for which elevation angle check is to be made.
            std::vector< std::string > listOfGroundStations;
            for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( );
                 linkEndIterator != linkEnds.end( ); linkEndIterator++ )
            {
                if( linkEndIterator->second.first == relevantObservationViabilitySettings.at( i )->getAssociatedLinkEnd( ).first )
                {
                    if( std::find( listOfGroundStations.begin( ), listOfGroundStations.end( ), linkEndIterator->second.second ) ==
                            listOfGroundStations.end( ) )
                    {
                        listOfGroundStations.push_back( linkEndIterator->second.second );
                    }
                }
            }

            // Create elevation angle check separately for eah ground station: check requires different pointing angles calculator
            for( unsigned int j = 0; j < listOfGroundStations.size( ); j++ )
            {
                linkViabilityCalculators.push_back(
                            createMinimumElevationAngleCalculator(
                                bodies, linkEnds, observationType, relevantObservationViabilitySettings.at( i ),
                                listOfGroundStations.at( j ) ) );
            }
            break;
        }
        case body_avoidance_angle:

            linkViabilityCalculators.push_back(
                        createBodyAvoidanceAngleCalculator(
                            bodies, linkEnds, observationType, relevantObservationViabilitySettings.at( i ) ) );
            break;
        case body_occultation:

            linkViabilityCalculators.push_back(
                        createOccultationCalculator(
                            bodies, linkEnds, observationType, relevantObservationViabilitySettings.at( i ) ) );
            break;
        default:
            throw std::runtime_error(
                        "Error when making observation viability calculator, type not recognized " +
                        std::to_string(
                            relevantObservationViabilitySettings.at( i )->observationViabilityType_ ) );
        }

    }

    return linkViabilityCalculators;
}

//! Function to create an list of obervation viability conditions for a number of sets of link ends, for a single observable type
std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > >
createObservationViabilityCalculators(
        const simulation_setup::SystemOfBodies& bodies,
        const std::vector< LinkEnds > linkEnds,
        const ObservableType observationType,
        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings )
{
    std::map< LinkEnds, std::vector< std::shared_ptr< ObservationViabilityCalculator > > > viabilityCalculators;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        viabilityCalculators[ linkEnds.at( i ) ] =
                createObservationViabilityCalculators(
                    bodies, linkEnds.at( i ), observationType, observationViabilitySettings );
    }

    return viabilityCalculators;
}


} // namespace observation_models

} // namespace tudat

