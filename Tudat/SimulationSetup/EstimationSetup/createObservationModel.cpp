/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <map>

#include <boost/function.hpp>
#include <boost/make_shared.hpp>


#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"

namespace tudat
{

namespace observation_models
{


//! Function to create list of observation models sorted by observable type and link ends from list only sorted in link ends.
SortedObservationSettingsMap convertUnsortedToSortedObservationSettingsMap(
        const ObservationSettingsMap& unsortedObservationSettingsMap )
{
    SortedObservationSettingsMap sortedObservationSettingsMap;

    for( ObservationSettingsMap::const_iterator iterator = unsortedObservationSettingsMap.begin( );
         iterator != unsortedObservationSettingsMap.end( ); iterator++ )
    {
        sortedObservationSettingsMap[ iterator->second->observableType_ ][ iterator->first ] =
                iterator->second;
    }
    return sortedObservationSettingsMap;
}


ObservationViabilitySettingsList filterObservationViabilitySettings(
        const ObservationViabilitySettingsList observationViabilitySettings,
        const LinkEnds linkEnds )
{
    ObservationViabilitySettingsList filteredViabilitySettings;

    for( unsigned int i = 0; i < observationViabilitySettings.size( ); i++ )
    {
        for( LinkEnds::const_iterator linkEndIterator = linkEnds.begin( ); linkEndIterator != linkEnds.end( ); linkEndIterator++ )
        {
            if( linkEndIterator->second == observationViabilitySettings.at( i )->associatedLinkEnd_ ||
                    ( ( observationViabilitySettings.at( i )->associatedLinkEnd_.second == "" ) &&
                      ( observationViabilitySettings.at( i )->associatedLinkEnd_.first == linkEndIterator->second.first ) ) )
            {
                filteredViabilitySettings.push_back( observationViabilitySettings.at( i ) );
                break;
            }
        }
    }

    return filteredViabilitySettings;
}

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
            std::cerr<<"Warning, parsed irrelevant 1-way range viability indices"<<std::endl;
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
            std::cerr<<"Warning, parsed irrelevant 1-way range viability indices"<<std::endl;
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
            std::cerr<<"Warning, parsed irrelevant 1-way range rate viability indices"<<std::endl;
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
                    linkEndIndices.push_back( std::make_pair( 2 * ( linkEnds.size( ) - 1 ) - 1, 2 * ( linkEnds.size( ) - 1 ) - 2 ) );
                }
                else
                {
                    linkEndIndices.push_back( std::make_pair( 2 * matchingLinkEndIndices.at( i ), 2 * matchingLinkEndIndices.at( i ) + 1 ) );
                    linkEndIndices.push_back( std::make_pair( 2 * matchingLinkEndIndices.at( i ) - 1, 2 * matchingLinkEndIndices.at( i ) - 2 ) );
                }
            }
        }
        else
        {
            std::cerr<<"Warning, parsed irrelevant n-way range viability indices"<<std::endl;
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
            std::cerr<<"Warning, parsed irrelevant angular position viability indices"<<std::endl;
        }
        break;
    case position_observable:

        std::cerr<<"Warning, parsed irrelevant position observable viability indices"<<std::endl;
        break;
    default:
        std::cerr<<"Error, observable type "<<observableType<<" not recognized when making viability link ends"<<std::endl;

    }

    return linkEndIndices;
}

boost::shared_ptr< MinimumElevationAngleCalculator > createMinimumElevationAngleCalculator(
        const simulation_setup::NamedBodyMap& bodyMap,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const boost::shared_ptr< ObservationViabilitySettings > observationViabilitySettings )
{
    if( observationViabilitySettings->observationViabilityType_ != minimum_elevation_angle )
    {
        std::cerr<<"Error when making minimum elevation angle calculator, inconsistent input"<<std::endl;
    }

    std::string groundStationName;
    if( observationViabilitySettings->associatedLinkEnd_.second != "" )
    {
        groundStationName = observationViabilitySettings->associatedLinkEnd_.second;
    }
    else
    {
        for( LinkEnds::const_iterator it = linkEnds.begin( ); it != linkEnds.end( ); it++ )
        {
            if( it->second.first == observationViabilitySettings->associatedLinkEnd_.first )
            {
                groundStationName = it->second.second;
                break;
            }
        }
    }

    boost::shared_ptr< ground_stations::PointingAnglesCalculator > pointingAngleCalculator =
            bodyMap.at( observationViabilitySettings->associatedLinkEnd_.first )->
            getGroundStation( groundStationName )->getPointingAnglesCalculator( );

    double minimumElevationAngle = observationViabilitySettings->doubleParameter_;

    return boost::make_shared< MinimumElevationAngleCalculator >(
                getLinkEndIndicesForObservationViability( linkEnds,observationType, observationViabilitySettings->associatedLinkEnd_ ),
                minimumElevationAngle, pointingAngleCalculator );
}

boost::shared_ptr< BodyAvoidanceAngleCalculator > createBodyAvoidanceAngleCalculator(
        const simulation_setup::NamedBodyMap& bodyMap,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const boost::shared_ptr< ObservationViabilitySettings > observationViabilitySettings )
{
    if( observationViabilitySettings->observationViabilityType_ != body_avoidance_angle )
    {
        std::cerr<<"Error when making body avoidance angle calculator, inconsistent input"<<std::endl;
    }

    boost::function< Eigen::Vector6d( const double ) > stateFunctionOfBodyToAvoid =
            boost::bind( &simulation_setup::Body::getStateInBaseFrameFromEphemeris< double, double >,
                         bodyMap.at( observationViabilitySettings->stringParameter_ ), _1 );

    double bodyAvoidanceAngle = observationViabilitySettings->doubleParameter_;

    return boost::make_shared< BodyAvoidanceAngleCalculator >(
                getLinkEndIndicesForObservationViability( linkEnds,observationType, observationViabilitySettings->associatedLinkEnd_ ),
                bodyAvoidanceAngle, stateFunctionOfBodyToAvoid, observationViabilitySettings->stringParameter_ );
}

std::vector< boost::shared_ptr< ObservationViabilityCalculator > > createObservationViabilityCalculators(
        const simulation_setup::NamedBodyMap& bodyMap,
        const LinkEnds linkEnds,
        const ObservableType observationType,
        const std::vector< boost::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings )
{
    std::vector< boost::shared_ptr< ObservationViabilityCalculator > > linkViabilityCalculators;

    std::vector< boost::shared_ptr< ObservationViabilitySettings > > relevantObservationViabilitySettings =
            filterObservationViabilitySettings( observationViabilitySettings, linkEnds );

    for( unsigned int i = 0; i < relevantObservationViabilitySettings.size( ); i++ )
    {
        switch( relevantObservationViabilitySettings.at( i )->observationViabilityType_ )
        {
        case minimum_elevation_angle:

            linkViabilityCalculators.push_back( createMinimumElevationAngleCalculator(
                                                    bodyMap, linkEnds, observationType, relevantObservationViabilitySettings.at( i ) ) );
            break;
        case body_avoidance_angle:

            linkViabilityCalculators.push_back( createBodyAvoidanceAngleCalculator(
                                                    bodyMap, linkEnds, observationType, relevantObservationViabilitySettings.at( i ) ) );
            break;
        default:
            std::cerr<<"Error when making observation viability calculator, type not recognized "<<
                       relevantObservationViabilitySettings.at( i )->observationViabilityType_<<std::endl;
        }

    }

    return linkViabilityCalculators;
}

std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationViabilityCalculator > > >
createObservationViabilityCalculators(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< LinkEnds > linkEnds,
        const ObservableType observationType,
        const std::vector< boost::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings )
{
    std::map< LinkEnds, std::vector< boost::shared_ptr< ObservationViabilityCalculator > > > viabilityCalculators;

    for( unsigned int i = 0; i < linkEnds.size( ); i++ )
    {
        viabilityCalculators[ linkEnds.at( i ) ] =
                createObservationViabilityCalculators(
                    bodyMap, linkEnds.at( i ), observationType, observationViabilitySettings );
    }

    return viabilityCalculators;
}

PerObservableObservationViabilityCalculatorList
createObservationViabilityCalculators(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable,
        const std::vector< boost::shared_ptr< ObservationViabilitySettings > >& observationViabilitySettings )
{
    PerObservableObservationViabilityCalculatorList viabilityCalculators;

    for( std::map< ObservableType, std::vector< LinkEnds > >::const_iterator observableIterator = linkEndsPerObservable.begin( );
         observableIterator != linkEndsPerObservable.end( ); observableIterator++ )
    {
        for( unsigned int i = 0; i < observableIterator->second.size( ); i++ )
        {
            viabilityCalculators[ observableIterator->first ][ observableIterator->second.at( i ) ] =
                    createObservationViabilityCalculators(
                        bodyMap, observableIterator->second.at( i ), observableIterator->first, observationViabilitySettings );
        }
    }

    return viabilityCalculators;
}

} // namespace observation_models

} // namespace tudat

