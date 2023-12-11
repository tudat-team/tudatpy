/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
//
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/utilities.h"

#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::orbit_determination;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::numerical_integrators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::statistics;
using namespace tudat::ground_stations;



BOOST_AUTO_TEST_SUITE( test_observation_dependent_variables )


void compareAgainstReference(
    const std::shared_ptr<ObservationCollection< > > simulatedObservations,
    const std::vector< std::shared_ptr< ObservationDependentVariableSettings > >& dependentVariableSettingsList,
    const std::vector< std::map< double, Eigen::VectorXd > >& referenceReceiverDependentVariableResults,
    const ObservableType observableType,
    const double expectedTimeOffset )
{
    BOOST_CHECK_EQUAL( dependentVariableSettingsList.size( ), referenceReceiverDependentVariableResults.size( ) );

    for( unsigned int i = 0; i < dependentVariableSettingsList.size( ); i++ )
    {
        std::map<double, Eigen::VectorXd> computedDependentVariables = observation_models::getDependentVariableResultList(
            simulatedObservations, dependentVariableSettingsList.at( i ), observableType );
        std::map< double, Eigen::VectorXd > referenceDependentVariables = referenceReceiverDependentVariableResults.at( i );

        BOOST_CHECK_EQUAL( computedDependentVariables.size( ), referenceDependentVariables.size( ) );

        if( referenceDependentVariables.size( ) > 0 )
        {
            int variableSize = referenceDependentVariables.begin( )->second.rows( );
            auto referenceIterator = referenceDependentVariables.begin( );
            auto computedIterator = computedDependentVariables.begin( );
            for ( auto it: referenceDependentVariables )
            {

                BOOST_CHECK_CLOSE_FRACTION( computedIterator->first - referenceIterator->first, expectedTimeOffset, 4.0 * std::numeric_limits< double >::epsilon( ) );

//                double currentTime = it.first;

//                std::cout<<i<<" "<<computedIterator->second<<" "<<referenceIterator->second<<std::endl;

                for ( int j = 0; j < variableSize; j++ )
                {
                    BOOST_CHECK_SMALL( std::fabs( computedIterator->second( j ) - referenceIterator->second( j ) ),
                                       std::numeric_limits<double>::epsilon( ) * referenceIterator->second.norm( ) );
                }
                referenceIterator++;
                computedIterator++;
            }
        }
    }
}

double computeLineSegmentToCenterOfMassDistance(
    const Eigen::Vector3d lineSegmentStart,
    const Eigen::Vector3d lineSegmentEnd,
    const Eigen::Vector3d pointLocation )
{
    Eigen::Vector3d lineDirection = lineSegmentEnd - lineSegmentStart;
    Eigen::Vector3d startToPoint = pointLocation - lineSegmentStart;
    Eigen::Vector3d endToPoint = pointLocation - lineSegmentEnd;

    double startInnerProduct = startToPoint.dot( lineDirection );
    double endInnerProduct = endToPoint.dot( lineDirection );

    double distance = TUDAT_NAN;
    if( startInnerProduct * endInnerProduct > 0 )
    {
        if( startToPoint.norm( ) < endToPoint.norm( ) )
        {
            distance = startToPoint.norm( );
        }
        else
        {
            distance = endToPoint.norm( );
        }
    }
    else
    {
        double angle = linear_algebra::computeAngleBetweenVectors( startToPoint, lineDirection );
        distance = std::sin( angle ) * startToPoint.norm( );
    }
    return distance;
}

//! Test whether the observation dependent variables are computed correctly
/*
 *  In this test, the calculation of observation dependent variables is checked against theoretical expectiations,
 *  for a link between an Earth ground station and a moon orbiter. The check is done for a one-way range observable,
 *  with transmitter/receiver as station/spacecraft (and the other way around). It is then checked whether the
 *  corresponding link in other observables yields identical results
 */
BOOST_AUTO_TEST_CASE( testObservationDependentVariables )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, "Earth" );

    // Add spacecraft orbiting Moon in Keplerian orbit
    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings = keplerEphemerisSettings(
        keplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    

    // Creatre ground stations
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 1.0 ).finished( ), geodetic_position );

    // Add relevant systems for DSN ovservable (X-band link; 3GHz transmission frequency)
    bodies.at( "Earth" )->getGroundStation( "Station1" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    bodies.at( "Earth" )->getGroundStation( "Station1" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );
    bodies.at( "Earth" )->getGroundStation( "Station1" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 3.0E9 ) );

    bodies.at( "Earth" )->getGroundStation( "Station2" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    bodies.at( "Earth" )->getGroundStation( "Station2" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );
    bodies.at( "Earth" )->getGroundStation( "Station2" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 3.0E9 ) );

    bodies.at( "MoonOrbiter" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    bodies.at( "MoonOrbiter" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );

    // Define relevant sets of link ends.

    // Station to spacecraft (1-way)
    std::vector< LinkEnds > stationTransmitterOneWayLinkEnds;

    // Spacecraft to station (1-way)
    std::vector< LinkEnds > stationReceiverOneWayLinkEnds;

    // Station->spacecraft->station (2-way)
    std::vector< LinkEnds > stationReceiverTwoWayLinkEnds;

    // Spacecraft->station->spacecraft (2-way)
    std::vector< LinkEnds > stationRetransmitterTwoWayLinkEnds;

    // Relative observation link ends: orbiter and Mars to station
    std::vector< LinkEnds > stationReceiverRelativeLinkEnds;

    // Relative observation link ends: Mars and orbiter to station
     std::vector< LinkEnds > stationReceiverOppositeRelativeLinkEnds;

    // Define link ends to/from each of the two ground stations for above list
    LinkEnds linkEnds;
    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        stationTransmitterOneWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        stationReceiverOneWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        stationReceiverTwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        stationRetransmitterTwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        linkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "Mars", "" ) );
        stationReceiverRelativeLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Mars", "" ) );
        stationReceiverOppositeRelativeLinkEnds.push_back( linkEnds );

    }


    // Station 2->spacecraft->station 1 (2-way)
    std::vector< LinkEnds > stationReceiverThreeWayLinkEnds;

    // Station 1->spacecraft->station 2 (2-way)
    std::vector< LinkEnds > stationTransmitterThreeWayLinkEnds;

    linkEnds.clear( );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );
    stationReceiverThreeWayLinkEnds.push_back( linkEnds );

    linkEnds.clear( );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );
    linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    stationTransmitterThreeWayLinkEnds.push_back( linkEnds );

    // List of dependent variables, computed for 1st observation model, and check for each subsequent one
    std::vector< std::map< double, Eigen::VectorXd > > referenceReceiverDependentVariableResults;
    std::vector< std::map< double, Eigen::VectorXd > > referenceTransmitterDependentVariableResults;

    // Run analysis for each observable (only compare against theory for 1st one)
    for( unsigned int currentObservableTestCase = 0; currentObservableTestCase < 9; currentObservableTestCase++ )
    {
        // Check if observable is differenced
        bool isDifferencedObservable = false;

        // Check geometry type: 0: 1-way; 1: 2-/3-way; 2: relative observation
        int geometryType = -1;

        // Set current observable
        ObservableType currentObservableType;
        if( currentObservableTestCase == 0 )
        {
            currentObservableType = one_way_range;
            geometryType = 0;
        }
        else if( currentObservableTestCase == 1 )
        {
            currentObservableType = angular_position;
            geometryType = 0;
        }
        else if( currentObservableTestCase == 2 )
        {
            currentObservableType = one_way_doppler;
            geometryType = 0;
        }
        else if( currentObservableTestCase == 3 )
        {
            currentObservableType = n_way_range;
            geometryType = 1;
        }
        else if( currentObservableTestCase == 4 )
        {
            currentObservableType = two_way_doppler;
            geometryType = 1;
        }
        else if( currentObservableTestCase == 5 )
        {
            currentObservableType = one_way_differenced_range;
            geometryType = 0;
            isDifferencedObservable = true;
        }
        else if( currentObservableTestCase == 6 )
        {
            currentObservableType = n_way_differenced_range;
            geometryType = 1;
            isDifferencedObservable = true;
        }
        else if( currentObservableTestCase == 7 )
        {
            currentObservableType = dsn_n_way_averaged_doppler;
            geometryType = 1;
            isDifferencedObservable = true;
        }
        else if( currentObservableTestCase == 8 )
        {
            currentObservableType = relative_angular_position;
            geometryType = 2;
        }

        // For geometry type, set number of link end cases (each using different link ends and/or reference link end)
        int numberOfLinkEndCases = -1;
        if( geometryType == 0 )
        {
            numberOfLinkEndCases = 2;
        }
        else if( geometryType == 1 )
        {
            numberOfLinkEndCases = 6;
        }
        else if( geometryType == 2 )
        {
            numberOfLinkEndCases = 2;
        }

        // If differenced observable, check for both start and end
        if( isDifferencedObservable )
        {
            numberOfLinkEndCases *= 2;
        }

        // Iterate over all link end settings
        for( int currentLinkEndCase = 0; currentLinkEndCase < numberOfLinkEndCases; currentLinkEndCase++ ) //numberOfLinkEndCases
        {
            // Define to check against which one-way range the results should be compared
            bool compareAgainstReceiver = -1;

            // Define properties of link to check
            LinkEndType referenceLinkEnd = unidentified_link_end;
            IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined;
            LinkEndType originatingLinkEndRole = unidentified_link_end;

            // If observable is differenced, check for both interval start and end
            if( isDifferencedObservable )
            {
                if( currentLinkEndCase < numberOfLinkEndCases / 2 )
                {
                    integratedObservableHandling = interval_start;
                }
                else
                {
                    integratedObservableHandling = interval_end;
                }
            }

            // Set current link ends.
            LinkEnds currentLinkEnds;
            switch( geometryType )
            {
            case 0:
            {
                switch ( currentLinkEndCase % 2 )
                {
                case 0:
                    currentLinkEnds = stationReceiverOneWayLinkEnds.at( 0 );
                    referenceLinkEnd = receiver;
                    originatingLinkEndRole = transmitter;
                    compareAgainstReceiver = true;
                    break;
                case 1:
                    currentLinkEnds = stationTransmitterOneWayLinkEnds.at( 0 );
                    referenceLinkEnd = transmitter;
                    originatingLinkEndRole = receiver;
                    compareAgainstReceiver = false;
                    break;
                default:
                    throw std::runtime_error( "Error in observation dependent variable unit test A " );
                }
                break;
            }
            case 1:
            {
                switch ( currentLinkEndCase % 2 )
                {
                case 0:
                    currentLinkEnds = stationReceiverTwoWayLinkEnds.at( 0 );
                    referenceLinkEnd = receiver;
                    originatingLinkEndRole = retransmitter;
                    compareAgainstReceiver = true;
                    break;
                case 1:
                    currentLinkEnds = stationReceiverTwoWayLinkEnds.at( 0 );
                    referenceLinkEnd = transmitter;
                    originatingLinkEndRole = retransmitter;
                    compareAgainstReceiver = false;
                    break;
                case 2:
                    currentLinkEnds = stationTransmitterThreeWayLinkEnds.at( 0 );
                    referenceLinkEnd = transmitter;
                    originatingLinkEndRole = retransmitter;
                    compareAgainstReceiver = false;
                    break;
                case 3:
                    currentLinkEnds = stationReceiverThreeWayLinkEnds.at( 0 );
                    referenceLinkEnd = receiver;
                    originatingLinkEndRole = retransmitter;
                    compareAgainstReceiver = true;
                    break;
                case 4:
                    currentLinkEnds = stationRetransmitterTwoWayLinkEnds.at( 0 );
                    referenceLinkEnd = retransmitter;
                    originatingLinkEndRole = receiver;
                    compareAgainstReceiver = false;
                    break;
                case 5:
                    currentLinkEnds = stationRetransmitterTwoWayLinkEnds.at( 0 );
                    referenceLinkEnd = retransmitter;
                    originatingLinkEndRole = transmitter;
                    compareAgainstReceiver = true;
                    break;
                default:
                    throw std::runtime_error( "Error in observation dependent variable unit test B " );
                }
                break;
            }
            case 2:
            {
                switch ( currentLinkEndCase )
                {
                case 0:
                    currentLinkEnds = stationReceiverRelativeLinkEnds.at( 0 );
                    referenceLinkEnd = receiver;
                    originatingLinkEndRole = transmitter;
                    compareAgainstReceiver = true;
                    break;
                case 1:
                    currentLinkEnds = stationReceiverOppositeRelativeLinkEnds.at( 0 );
                    referenceLinkEnd = receiver;
                    originatingLinkEndRole = transmitter2;
                    compareAgainstReceiver = true;
                    break;
                default:
                    throw std::runtime_error( "Error in observation dependent variable unit test C " );
                }
                break;
            }
            }

            // Skip DSN n-way differenced observable if reference link end is not receiver
            if( !( currentObservableType == dsn_n_way_averaged_doppler && referenceLinkEnd != receiver ) )
            {

                std::cout<<currentObservableTestCase<<" "<<currentLinkEnds.size( )<<" "<<geometryType<<" "<<currentLinkEndCase<<std::endl;

                // Define link ends for current observable
                std::map<ObservableType, std::vector<LinkEnds> > linkEndsPerObservable;
                linkEndsPerObservable[ currentObservableType ].push_back( currentLinkEnds );

                // Define observation settings for each observable/link ends combination
                std::vector<std::shared_ptr<ObservationModelSettings> > observationSettingsList;
                for ( std::map<ObservableType, std::vector<LinkEnds> >::iterator linkEndIterator = linkEndsPerObservable.begin( );
                      linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
                {
                    ObservableType currentObservable = linkEndIterator->first;
                    std::vector<LinkEnds> currentLinkEndsList = linkEndIterator->second;

                    for ( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
                    {
                        // Create observation settings
                        if( currentObservableType == n_way_differenced_range )
                        {
                            observationSettingsList.push_back( std::make_shared<NWayDifferencedRangeObservationSettings>(
                                currentLinkEndsList.at( i )));
                        }
                        else if( currentObservableType == dsn_n_way_averaged_doppler )
                        {
                            observationSettingsList.push_back( std::make_shared<DsnNWayAveragedDopplerObservationSettings>(
                                currentLinkEndsList.at( i )));
                        }
                        else
                        {
                            observationSettingsList.push_back( std::make_shared<ObservationModelSettings>(
                                currentObservable, currentLinkEndsList.at( i )));
                        }
                    }
                }

                // Create observation simulators
                std::vector<std::shared_ptr<ObservationSimulatorBase<double, double> > > observationSimulators =
                    createObservationSimulators( observationSettingsList, bodies );

                // Define ancilliary settings
                std::shared_ptr< ObservationAncilliarySimulationSettings > ancilliarySettings = nullptr;
                double integrationTime = 60.0;
                double referenceTimeShift = 0.0;
                if( currentObservableType == dsn_n_way_averaged_doppler )
                {
                    ancilliarySettings =
                    getDsnNWayAveragedDopplerAncillarySettings(
                        std::vector< FrequencyBands >{ x_band, x_band }, x_band, 7.0e9, integrationTime );
                }
                else if( isDifferencedObservable )
                {
                    ancilliarySettings = std::make_shared< ObservationAncilliarySimulationSettings >( );
                    ancilliarySettings->setAncilliaryDoubleData( doppler_integration_time, integrationTime );
                }

                // For differenced observables; shift reference time by half the integration time
                if( isDifferencedObservable )
                {
                    if( integratedObservableHandling == interval_start )
                    {
                        referenceTimeShift = integrationTime / 2.0;
                    }
                    else if( integratedObservableHandling == interval_end )
                    {
                        referenceTimeShift = -integrationTime / 2.0;
                    }
                }

                // Define osbervation times.
                std::vector<double> baseTimeList;
                double observationTimeStart = initialEphemerisTime + 1000.0;
                double observationInterval = 100.0;
                for ( unsigned int i = 0; i < 3; i++ )
                {
                    for ( unsigned int j = 0; j < 432; j++ )
                    {
                        baseTimeList.push_back( observationTimeStart + referenceTimeShift + static_cast< double >( i ) * 86400.0 +
                                                static_cast< double >( j ) * observationInterval );
                    }
                }

                // Define observation simulation settings (observation type, link end, times and reference link end)
                std::vector<std::shared_ptr<ObservationSimulationSettings<double> > > measurementSimulationInput;

                for ( std::map<ObservableType, std::vector<LinkEnds> >::iterator linkEndIterator = linkEndsPerObservable.begin( );
                      linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
                {
                    ObservableType currentObservable = linkEndIterator->first;
                    std::vector<LinkEnds> currentLinkEndsList = linkEndIterator->second;
                    for ( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
                    {
                        measurementSimulationInput.push_back(
                            std::make_shared<TabulatedObservationSimulationSettings<> >(
                                currentObservable, currentLinkEndsList.at( i ), baseTimeList, referenceLinkEnd,
                                std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ), nullptr,
                                ancilliarySettings ) );
                    }
                }

                // Define settings for dependent variables
                std::vector<std::shared_ptr<ObservationDependentVariableSettings> > dependentVariableList;

                std::shared_ptr<ObservationDependentVariableSettings> elevationAngleSettings =
                    std::make_shared<StationAngleObservationDependentVariableSettings>(
                        station_elevation_angle, LinkEndId( std::make_pair( "Earth", "Station1" )),
                        referenceLinkEnd, integratedObservableHandling, originatingLinkEndRole );
                std::shared_ptr<ObservationDependentVariableSettings> azimuthAngleSettings =
                    std::make_shared<StationAngleObservationDependentVariableSettings>(
                        station_azimuth_angle, LinkEndId( std::make_pair( "Earth", "Station1" )),
                        referenceLinkEnd, integratedObservableHandling, originatingLinkEndRole );
                std::shared_ptr<ObservationDependentVariableSettings> targetRangeSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        target_range, referenceLinkEnd, originatingLinkEndRole, integratedObservableHandling );
                std::shared_ptr<ObservationDependentVariableSettings> targetInverseRangeSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        target_range, originatingLinkEndRole, referenceLinkEnd, integratedObservableHandling );                
                std::shared_ptr<ObservationDependentVariableSettings> linkBodyCenterDistanceSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        link_body_center_distance, referenceLinkEnd, originatingLinkEndRole, integratedObservableHandling, "Moon" );
                std::shared_ptr<ObservationDependentVariableSettings> linkBodyCenterDistanceInverseSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        link_body_center_distance, originatingLinkEndRole, referenceLinkEnd, integratedObservableHandling, "Moon" );
                std::shared_ptr<ObservationDependentVariableSettings> linkLimbDistanceSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        link_limb_distance, referenceLinkEnd, originatingLinkEndRole, integratedObservableHandling, "Moon" );
                std::shared_ptr<ObservationDependentVariableSettings> linkLimbDistanceInverseSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        link_limb_distance, originatingLinkEndRole, referenceLinkEnd, integratedObservableHandling, "Moon" );
                std::shared_ptr<ObservationDependentVariableSettings> moonAvoidanceAngleSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        body_avoidance_angle_variable, referenceLinkEnd, originatingLinkEndRole, integratedObservableHandling, "Moon" );
                std::shared_ptr<ObservationDependentVariableSettings> orbitalPlaneAngleSettings =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        link_angle_with_orbital_plane, referenceLinkEnd, originatingLinkEndRole, integratedObservableHandling, "Moon" );

                dependentVariableList.push_back( elevationAngleSettings );
                dependentVariableList.push_back( azimuthAngleSettings );
                dependentVariableList.push_back( targetRangeSettings );
                dependentVariableList.push_back( targetInverseRangeSettings );
                dependentVariableList.push_back( linkBodyCenterDistanceSettings );
                dependentVariableList.push_back( linkBodyCenterDistanceInverseSettings );
                dependentVariableList.push_back( linkLimbDistanceSettings );
                dependentVariableList.push_back( linkLimbDistanceInverseSettings );
                dependentVariableList.push_back( moonAvoidanceAngleSettings );
                dependentVariableList.push_back( orbitalPlaneAngleSettings );

                addDependentVariablesToObservationSimulationSettings(
                    measurementSimulationInput, dependentVariableList, bodies );

                // Simulate noise-free observations
                std::shared_ptr<ObservationCollection<> > idealObservationsAndTimes = simulateObservations<double, double>(
                    measurementSimulationInput, observationSimulators, bodies );

                // If first case (one-way range) compare against theoretical expectations
                if ( currentObservableTestCase == 0 )
                {
                    std::map<double, Eigen::VectorXd> elevationAngles1 = getDependentVariableResultList(
                        idealObservationsAndTimes, elevationAngleSettings, currentObservableType );
                    std::map<double, Eigen::VectorXd> azimuthAngles1 = getDependentVariableResultList(
                        idealObservationsAndTimes, azimuthAngleSettings, currentObservableType );

                    std::map<double, Eigen::VectorXd> targetRanges1 = getDependentVariableResultList(
                        idealObservationsAndTimes, targetRangeSettings, currentObservableType );
                    std::map<double, Eigen::VectorXd> targetInverseRanges1 = getDependentVariableResultList(
                        idealObservationsAndTimes, targetInverseRangeSettings, currentObservableType );

                    std::map<double, Eigen::VectorXd> linkBodyDistances = getDependentVariableResultList(
                        idealObservationsAndTimes, linkBodyCenterDistanceSettings, currentObservableType );
                    std::map<double, Eigen::VectorXd> linkBodyInverseDistances = getDependentVariableResultList(
                        idealObservationsAndTimes, linkBodyCenterDistanceInverseSettings, currentObservableType );

                    std::map<double, Eigen::VectorXd> linkLimbDistances = getDependentVariableResultList(
                        idealObservationsAndTimes, linkLimbDistanceSettings, currentObservableType );
                    std::map<double, Eigen::VectorXd> linkLimbInverseDistances = getDependentVariableResultList(
                        idealObservationsAndTimes, linkLimbDistanceInverseSettings, currentObservableType );
                    std::map<double, Eigen::VectorXd> moonAvoidanceAngles = getDependentVariableResultList(
                        idealObservationsAndTimes, moonAvoidanceAngleSettings, currentObservableType );
                    std::map<double, Eigen::VectorXd> orbitalPlaneAngles = getDependentVariableResultList(
                        idealObservationsAndTimes, orbitalPlaneAngleSettings, currentObservableType );

                    // Add data to reference cases against which subsequent observables will be compared
                    if( currentLinkEndCase == 0 )
                    {
                        referenceReceiverDependentVariableResults.push_back( elevationAngles1 );
                        referenceReceiverDependentVariableResults.push_back( azimuthAngles1 );
                        referenceReceiverDependentVariableResults.push_back( targetRanges1 );
                        referenceReceiverDependentVariableResults.push_back( targetInverseRanges1 );
                        referenceReceiverDependentVariableResults.push_back( linkBodyDistances );
                        referenceReceiverDependentVariableResults.push_back( linkBodyInverseDistances );
                        referenceReceiverDependentVariableResults.push_back( linkLimbDistances );
                        referenceReceiverDependentVariableResults.push_back( linkLimbInverseDistances );
                        referenceReceiverDependentVariableResults.push_back( moonAvoidanceAngles );
                        referenceReceiverDependentVariableResults.push_back( orbitalPlaneAngles );

                    }
                    else if( currentLinkEndCase == 1 )
                    {
                        referenceTransmitterDependentVariableResults.push_back( elevationAngles1 );
                        referenceTransmitterDependentVariableResults.push_back( azimuthAngles1 );
                        referenceTransmitterDependentVariableResults.push_back( targetRanges1 );
                        referenceTransmitterDependentVariableResults.push_back( targetInverseRanges1 );
                        referenceTransmitterDependentVariableResults.push_back( linkBodyDistances );
                        referenceTransmitterDependentVariableResults.push_back( linkBodyInverseDistances );
                        referenceTransmitterDependentVariableResults.push_back( linkLimbDistances );
                        referenceTransmitterDependentVariableResults.push_back( linkLimbInverseDistances );
                        referenceTransmitterDependentVariableResults.push_back( moonAvoidanceAngles );
                        referenceTransmitterDependentVariableResults.push_back( orbitalPlaneAngles );
                    }

                    // Check size of dependent variable results vectors
                    std::map< LinkEnds, int >  linkEndIdentifiers = idealObservationsAndTimes->getLinkEndIdentifierMap( );
                    std::vector< int > linkEndIds = idealObservationsAndTimes->getConcatenatedLinkEndIds( );

                    int numberOfLinkEnds1Observations = utilities::countNumberOfOccurencesInVector<int>(
                        linkEndIds, linkEndIdentifiers.at( currentLinkEnds ) );

                    BOOST_CHECK_EQUAL( elevationAngles1.size( ), numberOfLinkEnds1Observations );
                    BOOST_CHECK_EQUAL( azimuthAngles1.size( ), numberOfLinkEnds1Observations );
                    BOOST_CHECK_EQUAL( targetRanges1.size( ), numberOfLinkEnds1Observations );

                    // Retrieve pointing angles calculators
                    std::shared_ptr< PointingAnglesCalculator> pointingAnglesCalculator1 =
                        bodies.at( "Earth" )->getGroundStation( "Station1" )->getPointingAnglesCalculator( );
                    std::shared_ptr< PointingAnglesCalculator> pointingAnglesCalculator2 =
                        bodies.at( "Earth" )->getGroundStation( "Station2" )->getPointingAnglesCalculator( );

                    // Retrieve observation model
                    std::shared_ptr<ObservationModel<1, double, double> > observationModel1 =
                        std::dynamic_pointer_cast<ObservationSimulator<1, double, double> >(
                            observationSimulators.at( 0 ))->getObservationModel( currentLinkEnds );

                    // Iterate over all times
                    std::vector<double> linkEndTimes;
                    std::vector<Eigen::Matrix<double, 6, 1> > linkEndStates;
                    for ( auto it: elevationAngles1 )
                    {
                        double currentTime = it.first;
                        double currentElevation = elevationAngles1.at( currentTime )( 0 );
                        double currentAzimuth = azimuthAngles1.at( currentTime )( 0 );
                        double targetRange = targetRanges1.at( currentTime )( 0 );
                        double targetInverseRange = targetInverseRanges1.at( currentTime )( 0 );
                        double linkBodyDistance = linkBodyDistances.at( currentTime )( 0 );
                        double linkBodyInverseDistance = linkBodyInverseDistances.at( currentTime )( 0 );
                        double linkLimbDistance = linkLimbDistances.at( currentTime )( 0 );
                        double linkLimbInverseDistance = linkLimbInverseDistances.at( currentTime )( 0 );
                        double moonAvoidanceAngle = moonAvoidanceAngles.at( currentTime )( 0 );
                        double orbitalPlaneAngle = orbitalPlaneAngles.at( currentTime )( 0 );

                        observationModel1->computeIdealObservationsWithLinkEndData(
                            currentTime, referenceLinkEnd, linkEndTimes, linkEndStates );

                        Eigen::Vector3d vectorToTarget = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 )).segment( 0, 3 );
                        Eigen::Vector6d moonState = spice_interface::getBodyCartesianStateAtEpoch(
                            "Moon", "Earth", "ECLIPJ2000", "None", ( linkEndTimes.at( 0 ) + linkEndTimes.at( 1 ) )/ 2.0 );

                        Eigen::Vector3d stationToMoon = moonState.segment< 3 >( 0 );
                        Eigen::Vector6d moonToSpacecraft = -moonState;

                        if( referenceLinkEnd == transmitter )
                        {
                            vectorToTarget *= -1.0;
                            stationToMoon -= linkEndStates.at( 0 ).segment< 3 >( 0 );
                            moonToSpacecraft = linkEndStates.at( 1 ) - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "Earth", "ECLIPJ2000", "None", linkEndTimes.at( 1 ) );
                        }
                        else
                        {
                            stationToMoon -= linkEndStates.at( 1 ).segment< 3 >( 0 );
                            moonToSpacecraft = linkEndStates.at( 0 ) - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "Earth", "ECLIPJ2000", "None", linkEndTimes.at( 0 ) );
                        }

                        Eigen::Vector3d orbitalAngularMomentum = moonToSpacecraft.segment< 3 >( 0 ).cross( moonToSpacecraft.segment< 3 >( 3 ) );

                        double referenceTime = ( referenceLinkEnd == transmitter ) ? linkEndTimes.at( 0 ) : linkEndTimes.at( 1 );

                        double elevationAngle = pointingAnglesCalculator1->calculateElevationAngleFromInertialVector(
                            vectorToTarget, referenceTime );
                        BOOST_CHECK_SMALL(( elevationAngle - currentElevation ), std::numeric_limits<double>::epsilon( ));

                        double azimuthAngle = pointingAnglesCalculator1->calculateAzimuthAngleFromInertialVector(
                            vectorToTarget, referenceTime );

                        double linkDistanceToMoon = computeLineSegmentToCenterOfMassDistance(
                            linkEndStates.at( 0 ).segment< 3 >( 0 ), linkEndStates.at( 1 ).segment< 3 >( 0 ),
                            moonState.segment< 3 >( 0 ) );

                        double manualMoonAvoidanceAngle = linear_algebra::computeAngleBetweenVectors( stationToMoon, vectorToTarget );
                        double manualOrbitalPlaneAngle = linear_algebra::computeAngleBetweenVectors( orbitalAngularMomentum, vectorToTarget ) - mathematical_constants::PI / 2.0;

                        BOOST_CHECK_SMALL( std::fabs( azimuthAngle - currentAzimuth ), std::numeric_limits<double>::epsilon( ));

                        BOOST_CHECK_SMALL( std::fabs( targetRange - vectorToTarget.norm( )),
                                          std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));
                        BOOST_CHECK_SMALL( std::fabs( targetInverseRange - targetRange ),
                                          std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));
                        BOOST_CHECK_SMALL( std::fabs( targetInverseRange - vectorToTarget.norm( )),
                                          std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));

                        BOOST_CHECK_SMALL( std::fabs( linkDistanceToMoon - linkBodyDistance ), 1.0E-4 );

                        BOOST_CHECK_SMALL( std::fabs( linkBodyDistance - linkBodyInverseDistance ),
                                          std::numeric_limits<double>::epsilon( ) * 1.0E7 );
                        BOOST_CHECK_SMALL( std::fabs( linkLimbDistance - linkLimbInverseDistance ),
                                          std::numeric_limits<double>::epsilon( ) * 1.0E7 );
                        BOOST_CHECK_SMALL( std::fabs( linkBodyDistance - linkLimbDistance - spice_interface::getAverageRadius( "Moon" ) ),
                                          std::numeric_limits<double>::epsilon( ) * 1.0E7 );
                        BOOST_CHECK_SMALL( std::fabs( moonAvoidanceAngle - manualMoonAvoidanceAngle ), 1.0E-12 );
                        BOOST_CHECK_SMALL( std::fabs( orbitalPlaneAngle - manualOrbitalPlaneAngle ), 1.0E-10 );


                    }
                }
                if( compareAgainstReceiver )
                {
                    compareAgainstReference(
                        idealObservationsAndTimes, dependentVariableList, referenceReceiverDependentVariableResults, currentObservableType, referenceTimeShift );
                }
                else
                {
                    compareAgainstReference(
                        idealObservationsAndTimes, dependentVariableList, referenceTransmitterDependentVariableResults, currentObservableType, referenceTimeShift );
                }
            }

        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

