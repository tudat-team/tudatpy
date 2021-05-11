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


#include "tudat/simulation/estimation_setup/createObservationModel.h"

namespace tudat
{

namespace observation_models
{

std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortObservationModelSettingsByType(
        const std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings )
{
    std::map< ObservableType, std::vector< std::shared_ptr< ObservationModelSettings > > > sortedObservationModelSettings;
    for( unsigned int i = 0; i < observationModelSettings.size( ); i++ )
    {
        sortedObservationModelSettings[ observationModelSettings.at( i )->observableType_ ].push_back(
                observationModelSettings.at( i ) );
    }
    return sortedObservationModelSettings;
}

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

//! Function to retrieve the link end indices in link end states/times that are to be used in viability calculation
std::vector< std::pair< int, int > > getLinkEndIndicesForObservationViability(
        const LinkEnds& linkEnds, const ObservableType observableType,  const LinkEndId linkEndToCheck )
{
    std::vector< std::pair< int, int > > linkEndIndices;

    switch( observableType )
    {
    case one_way_range:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) ||
                ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) && ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way range viability indices" );
        }
        break;
    case one_way_doppler:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                                  ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way doppler viability indices" );
        }
        break;
    case two_way_doppler:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                                  ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }

        if( linkEnds.at( reflector1 ) == linkEndToCheck || ( ( linkEnds.at( reflector1 ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }

        if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }

        if( linkEndIndices.size( ) == 0 )
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way doppler viability indices" );
        }
        break;
    case one_way_differenced_range:
        if( linkEnds.at( transmitter ) == linkEndToCheck || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                              linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
            linkEndIndices.push_back( std::make_pair( 2, 3 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
            linkEndIndices.push_back( std::make_pair( 3, 2 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant 1-way differenced range viability indices" );
        }
        break;
    case n_way_range:
    {
        std::vector< int > matchingLinkEndIndices = getNWayLinkEndIndicesFromLinkEndId( linkEndToCheck, linkEnds );
        if( matchingLinkEndIndices.size( ) > 0 )
        {
            for( unsigned int i = 0; i < matchingLinkEndIndices.size( ); i++ )
            {
                if( matchingLinkEndIndices.at( i ) == 0 )
                {
                    linkEndIndices.push_back( std::make_pair( 0, 1 ) );
                }
                else if( matchingLinkEndIndices.at( i ) == static_cast< int >( linkEnds.size( ) )  - 1 )
                {
                    linkEndIndices.push_back( std::make_pair( 2 * ( linkEnds.size( ) - 1 ) - 1,
                                                              2 * ( linkEnds.size( ) - 1 ) - 2 ) );
                }
                else
                {
                    linkEndIndices.push_back(
                                std::make_pair( 2 * matchingLinkEndIndices.at( i ),
                                                2 * matchingLinkEndIndices.at( i ) + 1 ) );
                    linkEndIndices.push_back(
                                std::make_pair( 2 * matchingLinkEndIndices.at( i ) - 1,
                                                2 * matchingLinkEndIndices.at( i ) - 2 ) );
                }
            }
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant n-way range viability indices" );
        }
        break;
    }
    case angular_position:
        if( ( linkEnds.at( transmitter ) == linkEndToCheck ) || ( ( linkEnds.at( transmitter ).first == linkEndToCheck.first ) &&
                                                                  ( linkEndToCheck.second == "" ) ) )
        {
            linkEndIndices.push_back( std::make_pair( 0, 1 ) );
        }
        else if( linkEnds.at( receiver ) == linkEndToCheck || ( ( linkEnds.at( receiver ).first == linkEndToCheck.first ) &&
                                                                linkEndToCheck.second == "" ) )
        {
            linkEndIndices.push_back( std::make_pair( 1, 0 ) );
        }
        else
        {
            throw std::runtime_error( "Error, parsed irrelevant angular position viability indices" );
        }
        break;
    case position_observable:

        throw std::runtime_error( "Error, parsed irrelevant position observable viability indices" );
        break;
    case euler_angle_313_observable:

        throw std::runtime_error( "Error, parsed irrelevant euler angle observable viability indices" );
        break;
    case velocity_observable:

        throw std::runtime_error( "Error, parsed irrelevant position observable viability indices" );
        break;
    default:
        throw std::runtime_error( "Error, observable type " + std::to_string(
                                      observableType ) + " not recognized when making viability link ends" );

    }

    return linkEndIndices;
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
        if( groundStationNameToUse != stationName )
        {
            throw std::runtime_error( "Error when making minimum elevation angle calculator, inconsistent station input" );
        }
        groundStationNameToUse = observationViabilitySettings->getAssociatedLinkEnd( ).second;
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
                getLinkEndIndicesForObservationViability(
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
                getLinkEndIndicesForObservationViability(
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
                getLinkEndIndicesForObservationViability(
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

//! Function to create a list of obervation viability conditions for any number of sets of link ends and observable types
PerObservableObservationViabilityCalculatorList
createObservationViabilityCalculators(
        const simulation_setup::SystemOfBodies& bodies,
        const std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable,
        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings )
{
    PerObservableObservationViabilityCalculatorList viabilityCalculators;

    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = linkEndsPerObservable.begin( );
         observableIterator != linkEndsPerObservable.end( ); observableIterator++ )
    {
        for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
        {
            viabilityCalculators[ observableIterator->first ][ observableIterator->second.at( i ) ] =
                    createObservationViabilityCalculators(
                        bodies, observableIterator->second.at( i ), observableIterator->first, observationViabilitySettings );
        }
    }

    return viabilityCalculators;
}


//void addObservationViabilitySettings(
//        std::shared_ptr< ObservationModelSettings > observationModelSettings,
//        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& viabilitySettingsList )
//{
//    std::vector< std::shared_ptr< ObservationViabilitySettings > > currentSettings =
//            std::move( observationModelSettings->viabilitySettingsList_ );
//    currentSettings.insert( std::end( currentSettings ), std::begin( viabilitySettingsList ), std::end( viabilitySettingsList ) );
//    observationModelSettings->viabilitySettingsList_ = std::move( currentSettings );
//}

//void addObservationViabilitySettings(
//        std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings,
//        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& viabilitySettingsList )
//{
//    for( unsigned int i = 0; i < observationModelSettings.size( ); i++ )
//    {
//        addObservationViabilitySettings( observationModelSettings.at( i ), viabilitySettingsList );
//    }
//}

//void addObservationViabilitySettings(
//        std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings,
//        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& viabilitySettingsList,
//        const ObservableType observationType )
//{
//    for( unsigned int i = 0; i < observationModelSettings.size( ); i++ )
//    {
//        if( observationModelSettings.at( i )->observableType_ == observationType )
//        {
//            addObservationViabilitySettings( observationModelSettings.at( i ), viabilitySettingsList );
//        }
//    }
//}

//void addObservationViabilitySettings(
//        std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings,
//        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& viabilitySettingsList,
//        const ObservableType observationType,
//        const std::vector< LinkEnds > linkEnds )
//{
//    for( unsigned int i = 0; i < observationModelSettings.size( ); i++ )
//    {
//        if( observationModelSettings.at( i )->observableType_ == observationType )
//        {
//            if( std::find( linkEnds.begin( ), linkEnds.end( ), observationModelSettings.at( i )->linkEnds_ ) !=
//                    linkEnds.end( ) )
//            {
//                addObservationViabilitySettings( observationModelSettings.at( i ), viabilitySettingsList );
//            }
//        }
//    }
//}

//void addObservationViabilitySettings(
//        std::vector< std::shared_ptr< ObservationModelSettings > >& observationModelSettings,
//        const std::vector< std::shared_ptr< ObservationViabilitySettings > >& viabilitySettingsList,
//        const std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerType )
//{
//    for( auto it : linkEndsPerType )
//    {
//        addObservationViabilitySettings(
//                    observationModelSettings, viabilitySettingsList, it.first, it.second );
//    }
//}

} // namespace observation_models

} // namespace tudat

