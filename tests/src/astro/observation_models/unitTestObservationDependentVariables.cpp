/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/utilities.h"
#include "tudat/basics/testMacros.h"

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
        std::map<double, Eigen::VectorXd> computedDependentVariables = simulatedObservations->getDependentVariableHistory( dependentVariableSettingsList.at( i ) );
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
    for( unsigned int currentObservableTestCase = 0; currentObservableTestCase < 8; currentObservableTestCase++ )
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
//        else if( currentObservableTestCase == 7 )
//        {
//            currentObservableType = dsn_n_way_averaged_doppler;
//            geometryType = 1;
//            isDifferencedObservable = true;
//        }
        else if( currentObservableTestCase == 7 )
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
                std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariableList;

                std::shared_ptr< ObservationDependentVariableSettings > elevationAngleSettings = elevationAngleDependentVariable(
                        referenceLinkEnd, LinkEndId( std::make_pair( "Earth", "Station1" ) ), originatingLinkEndRole, LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > linkEndTypeElevationAngleSettings = elevationAngleDependentVariable(
                        referenceLinkEnd, LinkEndId( "", "" ), originatingLinkEndRole, LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > azimuthAngleSettings = azimuthAngleDependentVariable(
                        referenceLinkEnd, LinkEndId( std::make_pair( "Earth", "Station1" ) ), originatingLinkEndRole, LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > linkEndTypeAzimuthAngleSettings = azimuthAngleDependentVariable(
                        referenceLinkEnd, LinkEndId( std::make_pair( "", "" ) ), originatingLinkEndRole, LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > targetRangeSettings = targetRangeBetweenLinkEndsDependentVariable(
                        originatingLinkEndRole, referenceLinkEnd, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > targetInverseRangeSettings = targetRangeBetweenLinkEndsDependentVariable(
                        referenceLinkEnd, originatingLinkEndRole, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > linkBodyCenterDistanceSettings = linkBodyCenterDistanceDependentVariable(
                        "Moon", originatingLinkEndRole, referenceLinkEnd, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > linkBodyCenterDistanceInverseSettings = linkBodyCenterDistanceDependentVariable(
                        "Moon", referenceLinkEnd, originatingLinkEndRole, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > linkLimbDistanceSettings = linkLimbDistanceDependentVariable(
                        "Moon", originatingLinkEndRole, referenceLinkEnd, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > linkLimbDistanceInverseSettings = linkLimbDistanceDependentVariable(
                        "Moon", referenceLinkEnd, originatingLinkEndRole, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > moonAvoidanceAngleSettings = bodyAvoidanceAngleDependentVariable(
                        "Moon", originatingLinkEndRole, referenceLinkEnd, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > moonAvoidanceAngleSettings2 = bodyAvoidanceAngleDependentVariable(
                        "Moon", referenceLinkEnd, originatingLinkEndRole, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );
                std::shared_ptr< ObservationDependentVariableSettings > orbitalPlaneAngleSettings = linkAngleWrtOrbitalPlaneDependentVariable(
                        "Moon", originatingLinkEndRole, referenceLinkEnd, LinkEndId( "", "" ), LinkEndId( "", "" ), integratedObservableHandling );

                dependentVariableList.push_back( elevationAngleSettings );
                dependentVariableList.push_back( linkEndTypeElevationAngleSettings );
                dependentVariableList.push_back( azimuthAngleSettings );
                dependentVariableList.push_back( linkEndTypeAzimuthAngleSettings );
                dependentVariableList.push_back( targetRangeSettings );
                dependentVariableList.push_back( targetInverseRangeSettings );
                dependentVariableList.push_back( linkBodyCenterDistanceSettings );
                dependentVariableList.push_back( linkBodyCenterDistanceInverseSettings );
                dependentVariableList.push_back( linkLimbDistanceSettings );
                dependentVariableList.push_back( linkLimbDistanceInverseSettings );
                dependentVariableList.push_back( moonAvoidanceAngleSettings );
                dependentVariableList.push_back( moonAvoidanceAngleSettings2 );
                dependentVariableList.push_back( orbitalPlaneAngleSettings );


                addDependentVariablesToObservationSimulationSettings(
                    measurementSimulationInput, dependentVariableList, bodies );

                // Simulate noise-free observations
                std::shared_ptr<ObservationCollection<> > idealObservationsAndTimes = simulateObservations<double, double>(
                    measurementSimulationInput, observationSimulators, bodies );

                // If first case (one-way range) compare against theoretical expectations
                if ( currentObservableTestCase == 0 )
                {
                    std::map< double, Eigen::VectorXd > elevationAngles1 = idealObservationsAndTimes->getDependentVariableHistory( elevationAngleSettings );
                    std::map< double, Eigen::VectorXd > elevationAngles2 = idealObservationsAndTimes->getDependentVariableHistory( linkEndTypeElevationAngleSettings );

                    std::map< double, Eigen::VectorXd > azimuthAngles1 = idealObservationsAndTimes->getDependentVariableHistory( azimuthAngleSettings );
                    std::map< double, Eigen::VectorXd > azimuthAngles2 = idealObservationsAndTimes->getDependentVariableHistory( linkEndTypeAzimuthAngleSettings );

                    std::map< double, Eigen::VectorXd > targetRanges1 = idealObservationsAndTimes->getDependentVariableHistory( targetRangeSettings );
                    std::map< double, Eigen::VectorXd > targetInverseRanges1 = idealObservationsAndTimes->getDependentVariableHistory( targetInverseRangeSettings );

                    std::map< double, Eigen::VectorXd > linkBodyDistances = idealObservationsAndTimes->getDependentVariableHistory( linkBodyCenterDistanceSettings );
                    std::map< double, Eigen::VectorXd > linkBodyInverseDistances = idealObservationsAndTimes->getDependentVariableHistory( linkBodyCenterDistanceInverseSettings );

                    std::map< double, Eigen::VectorXd > linkLimbDistances = idealObservationsAndTimes->getDependentVariableHistory( linkLimbDistanceSettings );
                    std::map< double, Eigen::VectorXd > linkLimbInverseDistances = idealObservationsAndTimes->getDependentVariableHistory( linkLimbDistanceInverseSettings );

                    std::map< double, Eigen::VectorXd > moonAvoidanceAngles = idealObservationsAndTimes->getDependentVariableHistory( moonAvoidanceAngleSettings );
                    std::map< double, Eigen::VectorXd > moonAvoidanceAngles2 = idealObservationsAndTimes->getDependentVariableHistory( moonAvoidanceAngleSettings2 );

                    std::map< double, Eigen::VectorXd > orbitalPlaneAngles = idealObservationsAndTimes->getDependentVariableHistory( orbitalPlaneAngleSettings );

                    // Add data to reference cases against which subsequent observables will be compared
                    if( currentLinkEndCase == 0 )
                    {
                        referenceReceiverDependentVariableResults.push_back( elevationAngles1 );
                        referenceReceiverDependentVariableResults.push_back( elevationAngles2 );
                        referenceReceiverDependentVariableResults.push_back( azimuthAngles1 );
                        referenceReceiverDependentVariableResults.push_back( azimuthAngles2 );
                        referenceReceiverDependentVariableResults.push_back( targetRanges1 );
                        referenceReceiverDependentVariableResults.push_back( targetInverseRanges1 );
                        referenceReceiverDependentVariableResults.push_back( linkBodyDistances );
                        referenceReceiverDependentVariableResults.push_back( linkBodyInverseDistances );
                        referenceReceiverDependentVariableResults.push_back( linkLimbDistances );
                        referenceReceiverDependentVariableResults.push_back( linkLimbInverseDistances );
                        referenceReceiverDependentVariableResults.push_back( moonAvoidanceAngles );
                        referenceReceiverDependentVariableResults.push_back( moonAvoidanceAngles2 );
                        referenceReceiverDependentVariableResults.push_back( orbitalPlaneAngles );

                    }
                    else if( currentLinkEndCase == 1 )
                    {
                        referenceTransmitterDependentVariableResults.push_back( elevationAngles1 );
                        referenceTransmitterDependentVariableResults.push_back( elevationAngles2 );
                        referenceTransmitterDependentVariableResults.push_back( azimuthAngles1 );
                        referenceTransmitterDependentVariableResults.push_back( azimuthAngles2 );
                        referenceTransmitterDependentVariableResults.push_back( targetRanges1 );
                        referenceTransmitterDependentVariableResults.push_back( targetInverseRanges1 );
                        referenceTransmitterDependentVariableResults.push_back( linkBodyDistances );
                        referenceTransmitterDependentVariableResults.push_back( linkBodyInverseDistances );
                        referenceTransmitterDependentVariableResults.push_back( linkLimbDistances );
                        referenceTransmitterDependentVariableResults.push_back( linkLimbInverseDistances );
                        referenceTransmitterDependentVariableResults.push_back( moonAvoidanceAngles );
                        referenceTransmitterDependentVariableResults.push_back( moonAvoidanceAngles2 );
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
                        double currentElevation2 = elevationAngles2.at( currentTime )( 0 );
                        double currentAzimuth = azimuthAngles1.at( currentTime )( 0 );
                        double currentAzimuth2 = azimuthAngles2.at( currentTime )( 0 );
                        double targetRange = targetRanges1.at( currentTime )( 0 );
                        double targetInverseRange = targetInverseRanges1.at( currentTime )( 0 );
                        double linkBodyDistance = linkBodyDistances.at( currentTime )( 0 );
                        double linkBodyInverseDistance = linkBodyInverseDistances.at( currentTime )( 0 );
                        double linkLimbDistance = linkLimbDistances.at( currentTime )( 0 );
                        double linkLimbInverseDistance = linkLimbInverseDistances.at( currentTime )( 0 );
                        double moonAvoidanceAngle = moonAvoidanceAngles.at( currentTime )( 0 );
                        double moonAvoidanceAngle2 = moonAvoidanceAngles2.at( currentTime )( 0 );
                        double orbitalPlaneAngle = orbitalPlaneAngles.at( currentTime )( 0 );

                        observationModel1->computeIdealObservationsWithLinkEndData(
                            currentTime, referenceLinkEnd, linkEndTimes, linkEndStates );

                        Eigen::Vector3d vectorToTarget = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 )).segment( 0, 3 );
                        Eigen::Vector6d moonState = spice_interface::getBodyCartesianStateAtEpoch(
                            "Moon", "Earth", "ECLIPJ2000", "None", ( linkEndTimes.at( 0 ) + linkEndTimes.at( 1 ) )/ 2.0 );

                        Eigen::Vector3d stationToMoon = moonState.segment< 3 >( 0 );
                        Eigen::Vector3d spacecraftToMoon = moonState.segment< 3 >( 0 );
                        Eigen::Vector6d moonToSpacecraft = -moonState;

                        if( referenceLinkEnd == transmitter )
                        {
                            vectorToTarget *= -1.0;
                            stationToMoon -= linkEndStates.at( 0 ).segment< 3 >( 0 );
                            spacecraftToMoon -= linkEndStates.at( 1 ).segment< 3 >( 0 );
                            moonToSpacecraft = linkEndStates.at( 1 ) - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "Earth", "ECLIPJ2000", "None", linkEndTimes.at( 1 ) );
                        }
                        else
                        {
                            stationToMoon -= linkEndStates.at( 1 ).segment< 3 >( 0 );
                            spacecraftToMoon -= linkEndStates.at( 0 ).segment< 3 >( 0 );
                            moonToSpacecraft = linkEndStates.at( 0 ) - spice_interface::getBodyCartesianStateAtEpoch(
                                "Moon", "Earth", "ECLIPJ2000", "None", linkEndTimes.at( 0 ) );
                        }

                        Eigen::Vector3d orbitalAngularMomentum = moonToSpacecraft.segment< 3 >( 0 ).cross( moonToSpacecraft.segment< 3 >( 3 ) );

                        double referenceTime = ( referenceLinkEnd == transmitter ) ? linkEndTimes.at( 0 ) : linkEndTimes.at( 1 );

                        double elevationAngle = pointingAnglesCalculator1->calculateElevationAngleFromInertialVector(
                            vectorToTarget, referenceTime );
                        BOOST_CHECK_SMALL(( elevationAngle - currentElevation ), std::numeric_limits<double>::epsilon( ));
                        BOOST_CHECK_SMALL(( elevationAngle - currentElevation2 ), std::numeric_limits<double>::epsilon( ));

                        double azimuthAngle = pointingAnglesCalculator1->calculateAzimuthAngleFromInertialVector(
                            vectorToTarget, referenceTime );

                        double linkDistanceToMoon = computeLineSegmentToCenterOfMassDistance(
                            linkEndStates.at( 0 ).segment< 3 >( 0 ), linkEndStates.at( 1 ).segment< 3 >( 0 ),
                            moonState.segment< 3 >( 0 ) );

                        double manualMoonAvoidanceAngle = linear_algebra::computeAngleBetweenVectors( stationToMoon, vectorToTarget );
                        double manualMoonAvoidanceAngle2 = linear_algebra::computeAngleBetweenVectors( spacecraftToMoon, - vectorToTarget );
                        double manualOrbitalPlaneAngle = linear_algebra::computeAngleBetweenVectors( orbitalAngularMomentum, vectorToTarget ) - mathematical_constants::PI / 2.0;

                        BOOST_CHECK_SMALL( std::fabs( azimuthAngle - currentAzimuth ), std::numeric_limits<double>::epsilon( ));
                        BOOST_CHECK_SMALL( std::fabs( azimuthAngle - currentAzimuth2 ), std::numeric_limits<double>::epsilon( ));

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
                        BOOST_CHECK_SMALL( std::fabs( moonAvoidanceAngle2 - manualMoonAvoidanceAngle2 ), 1.0E-10 );
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

//! Test whether the interfaces to create and get observation dependent variables work properly
/*
 *  This test does not check the calculation of observation dependent variables, but rather verifies that the interfaces used to easily create observation dependent
 *  variables and retrieve the associated compute values work as expected. This check is done for various observable types and link ends.
 */
BOOST_AUTO_TEST_CASE( testObservationDependentVariablesInterface )
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
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth" );

    // Add spacecraft orbiting Moon in Keplerian orbit
    bodySettings.addSettings( "MoonOrbiter" );
    Eigen::Vector6d keplerElements = Eigen::Vector6d::Zero( );
    keplerElements( 0 ) = 2.0E6;
    keplerElements( 1 ) = 0.1;
    keplerElements( 2 ) = 1.0;
    bodySettings.at( "MoonOrbiter" )->ephemerisSettings = keplerEphemerisSettings(
            keplerElements, 0.0, spice_interface::getBodyGravitationalParameter( "Moon" ), "Moon" );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );


    // Create ground stations
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 1.0 ).finished( ), geodetic_position );

    // Add relevant systems for DSN observable (X-band link; 3GHz transmission frequency)
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
    LinkEnds oneWayLinkEnds;
    oneWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    oneWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );

    // Station->spacecraft->station (2-way)
    LinkEnds station1TwoWayLinkEnds;
    station1TwoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    station1TwoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    station1TwoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );

    LinkEnds station2TwoWayLinkEnds;
    station2TwoWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );
    station2TwoWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    station2TwoWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );

    // Relative observation link ends: orbiter and Mars to station
    LinkEnds relativeLinkEnds;
    relativeLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    relativeLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    relativeLinkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "Mars", "" ) );

    // Station 2->spacecraft->station 1 (3-way)
    LinkEnds threeWayLinkEnds;
    threeWayLinkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    threeWayLinkEnds[ retransmitter ] = LinkEndId( std::make_pair( "MoonOrbiter", "" ) );
    threeWayLinkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );

    // Define link ends per observable type
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ n_way_differenced_range ] = { threeWayLinkEnds };
    linkEndsPerObservable[ dsn_n_way_averaged_doppler ] = { station1TwoWayLinkEnds, station2TwoWayLinkEnds };
    linkEndsPerObservable[ one_way_range ] = { oneWayLinkEnds };
    linkEndsPerObservable[ relative_angular_position ] = { relativeLinkEnds };


    // Define observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;

    // 3-way range
    observationSettingsList.push_back( std::make_shared< NWayDifferencedRangeObservationSettings >( threeWayLinkEnds ) );

    // 2-way DSN Doppler (for both ground stations)
    observationSettingsList.push_back( std::make_shared< DsnNWayAveragedDopplerObservationSettings >( station1TwoWayLinkEnds ) );
    observationSettingsList.push_back( std::make_shared< DsnNWayAveragedDopplerObservationSettings >( station2TwoWayLinkEnds ) );

    // 1-way range (for both ground stations)
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, oneWayLinkEnds ) );

    // relative angular position
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( relative_angular_position, relativeLinkEnds ) );

    // Create observation simulators
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( observationSettingsList, bodies );


    // Define observation simulation times and ancillary settings
    std::map< ObservableType, std::vector< double > > baseTimeListPerObservable;
    std::map< ObservableType, std::shared_ptr< ObservationAncilliarySimulationSettings > > ancillarySettingsPerObservable;
    for ( auto linkEndIterator : linkEndsPerObservable )
    {
        ObservableType currentObservable = linkEndIterator.first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator.second;

        // Define ancilliary settings
        std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings = nullptr;
        double integrationTime = 60.0;
        double referenceTimeShift = 0.0;
        if( currentObservable == dsn_n_way_averaged_doppler )
        {
            std::vector< double > delays = std::vector< double >( { 1.0e-3 } );
            ancillarySettings = getDsnNWayAveragedDopplerAncillarySettings(
                    std::vector< FrequencyBands >{ x_band, x_band }, x_band, 7.0e9, integrationTime, delays );
        }
        else if( currentObservable == n_way_differenced_range )
        {
            ancillarySettings = std::make_shared< ObservationAncilliarySimulationSettings >( );
            ancillarySettings->setAncilliaryDoubleData( doppler_integration_time, integrationTime );

            std::vector< double > delays = std::vector< double >( { 1.0e-3 } );
            ancillarySettings->setAncilliaryDoubleVectorData( link_ends_delays, delays );
        }
        ancillarySettingsPerObservable[ currentObservable ] = ancillarySettings;

        // For differenced observables; shift reference time by half the integration time
        if( currentObservable == dsn_n_way_averaged_doppler || currentObservable == n_way_differenced_range )
        {
            referenceTimeShift = integrationTime / 2.0;
        }

        // Define observation times.
        std::vector< double > baseTimeList;
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
        baseTimeListPerObservable[ currentObservable ] = baseTimeList;
    }


    // Station elevation settings
    std::shared_ptr< ObservationDependentVariableSettings > elevationAngleSettings = elevationAngleDependentVariable( );

    // Station azimuth settings
    std::shared_ptr< ObservationDependentVariableSettings > azimuthStationSettings1 = azimuthAngleDependentVariable( unidentified_link_end, LinkEndId( "Earth", "Station1" ) );
    std::shared_ptr< ObservationDependentVariableSettings > azimuthStationSettings2 = azimuthAngleDependentVariable( );

    // Limb distance settings
    std::shared_ptr<ObservationDependentVariableSettings> limbDistanceSettings = linkLimbDistanceDependentVariable( "Moon" );

    // Avoidance angle settings
    std::shared_ptr< ObservationDependentVariableSettings > moonAvoidanceAngleSettings = bodyAvoidanceAngleDependentVariable( "Moon", unidentified_link_end, receiver );

    // Integration time settings
    std::shared_ptr< ObservationDependentVariableSettings > integrationTimeSettings = integrationTimeDependentVariable( );

    // Retransmission delays settings
    std::shared_ptr< ObservationDependentVariableSettings > retransmissionDelaysSettings = retransmissionDelaysDependentVariable( );

    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariablesList;
    dependentVariablesList.push_back( elevationAngleSettings );
    dependentVariablesList.push_back( azimuthStationSettings1 );
    dependentVariablesList.push_back( limbDistanceSettings );
    dependentVariablesList.push_back( moonAvoidanceAngleSettings );
    dependentVariablesList.push_back( integrationTimeSettings );
    dependentVariablesList.push_back( retransmissionDelaysSettings );

    std::map< unsigned int, std::map< ObservationDependentVariables, std::map< ObservableType, std::vector< unsigned int > > > > numberOfSettingsToBeCreated;

    // number of settings for test case 0: dependent variables set in simulation settings
    std::map< ObservationDependentVariables, std::map< ObservableType, std::vector< unsigned int > > > numberOfSettingsTestCase0 = {
            {
                station_elevation_angle, { { n_way_differenced_range, { 2 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 2, 2 } ) },
                                         { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
            },
            {
                station_azimuth_angle, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 2, 0 } ) },
                                         { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
            },
            {
                link_limb_distance, { { n_way_differenced_range, { 2 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 2, 2 } ) },
                                      { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
            },
            {
                body_avoidance_angle_variable, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 1, 1 } ) },
                                                 { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
            },
            {
                integration_time_dependent_variable, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 1, 1 } ) },
                                                               { one_way_range, { 0 } }, { relative_angular_position, { 0 } } }
            },
            {
                retransmission_delays_dependent_variable, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 1, 1 } ) },
                                                            { one_way_range, { 0 } }, { relative_angular_position, { 0 } } }
            }
    };


    // number of settings for test case 1: dependent variables defined after observation collection is created
    std::map< ObservationDependentVariables, std::map< ObservableType, std::vector< unsigned int > > > numberOfSettingsTestCase1 = {
        {
                station_elevation_angle, { { n_way_differenced_range, { 2 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 2, 2 } ) },
                                           { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
        },
        {
                station_azimuth_angle, { { n_way_differenced_range, { 2 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 2, 0 } ) },
                                         { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
        },
        {
                link_limb_distance, { { n_way_differenced_range, { 2 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 2, 2 } ) },
                                      { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
        },
        {
                body_avoidance_angle_variable, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 1, 0 } ) },
                                                 { one_way_range, { 1 } }, { relative_angular_position, { 2 } } }
        },
        {
                integration_time_dependent_variable, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 1, 1 } ) },
                                                             { one_way_range, { 0 } }, { relative_angular_position, { 0 } } }
        },
        {
                retransmission_delays_dependent_variable, { { n_way_differenced_range, { 1 } }, { dsn_n_way_averaged_doppler, std::vector< unsigned int >( { 1, 1 } ) },
                                                          { one_way_range, { 0 } }, { relative_angular_position, { 0 } } }
        }
    };

    numberOfSettingsToBeCreated[ 0 ] = numberOfSettingsTestCase0;
    numberOfSettingsToBeCreated[ 1 ] = numberOfSettingsTestCase1;


    std::map< ObservationDependentVariables, std::vector< std::vector< Eigen::MatrixXd > > > dependentVariablesReferenceValues;
    for ( unsigned int testCase = 0 ; testCase < 2 ; testCase++ )
    {
        // Define observation simulation settings
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        for ( auto linkEndIt : linkEndsPerObservable )
        {
            ObservableType observableType = linkEndIt.first;
            std::vector< LinkEnds > linkEndsList = linkEndIt.second;
            for ( unsigned int i = 0; i < linkEndsList.size( ); i++ )
            {
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< > >(
                        linkEndIt.first, linkEndsList.at( i ), baseTimeListPerObservable.at( observableType ), receiver,
                        std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > >( ), nullptr, ancillarySettingsPerObservable.at( observableType ) ) );
            }
        }

        if ( testCase == 0 )
        {
            // Add dependent variables to simulation settings
            addDependentVariablesToObservationSimulationSettings( measurementSimulationInput, dependentVariablesList, bodies );
        }

        // Simulate noise-free observations
        std::shared_ptr< ObservationCollection< > > idealObservationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, observationSimulators, bodies );


        if ( testCase == 1 )
        {
            // Add dependent variables after observation collection is created
            std::shared_ptr< ObservationCollectionParser > elevationAngleParser =
                    idealObservationsAndTimes->addDependentVariable( elevationAngleSettings, bodies );
            std::shared_ptr< ObservationCollectionParser > azimuthAngleParser1 =
                    idealObservationsAndTimes->addDependentVariable( azimuthStationSettings1, bodies );
            std::shared_ptr< ObservationCollectionParser > azimuthAngleParser2 =
                    idealObservationsAndTimes->addDependentVariable( azimuthStationSettings2, bodies, observationParser( n_way_differenced_range ) );
            std::shared_ptr< ObservationCollectionParser > limbDistanceParser =
                    idealObservationsAndTimes->addDependentVariable( limbDistanceSettings, bodies );
            std::shared_ptr< ObservationCollectionParser > moonAngleParser =
                    idealObservationsAndTimes->addDependentVariable( moonAvoidanceAngleSettings, bodies, observationParser( std::make_pair( "Earth", "Station1" ) ) );
            std::shared_ptr< ObservationCollectionParser > integrationTimeParser =
                    idealObservationsAndTimes->addDependentVariable( integrationTimeSettings, bodies );
            std::shared_ptr< ObservationCollectionParser > retransmissionDelaysParser =
                    idealObservationsAndTimes->addDependentVariable( retransmissionDelaysSettings, bodies );

            // Compute dependent variables
            computeResidualsAndDependentVariables< double, double >( idealObservationsAndTimes, observationSimulators, bodies );
        }

        // Define number of dependent variables settings that should be created
        std::map< ObservationDependentVariables, std::map< ObservableType, std::vector< unsigned int > > > expectedNumberOfSettings;
        if ( testCase == 0 )
        {
            expectedNumberOfSettings = numberOfSettingsTestCase0;
        }
        else if ( testCase == 1 )
        {
            expectedNumberOfSettings = numberOfSettingsTestCase1;
        }

        // Parse all dependent variable types
        for ( auto variableIt : expectedNumberOfSettings )
        {
            // Parse all observable types
            for ( auto observableIt : variableIt.second )
            {
                // Retrieve single observation sets for given observable
                std::vector< std::shared_ptr< SingleObservationSet< double, double > > > observationSets =
                        idealObservationsAndTimes->getSingleObservationSets( observationParser( observableIt.first ) );

                if ( observationSets.size( ) != observableIt.second.size( ) )
                {
                    throw std::runtime_error( "Error when comparing number of dependent variable settings effectively created, number of reference values inconsistent"
                                              "with number of observation sets for observable " + std::to_string( observableIt.first ) + "." );
                }

                for ( unsigned int i = 0 ; i < observationSets.size( ) ; i++ )
                {
                    // Retrieve relevant dependent variable settings and check numbers of settings are consistent
                    std::vector< Eigen::MatrixXd > dependentVariables = observationSets.at( i )->getAllCompatibleDependentVariables(
                            std::make_shared< ObservationDependentVariableSettings >( variableIt.first ) );
                    BOOST_CHECK( dependentVariables.size( ) == observableIt.second.at( i ) );
                }
            }
        }

        std::map< ObservableType, std::map< LinkEnds, std::vector< std::shared_ptr< SingleObservationSet< double, double > > > > > sortedObservationSets =
                idealObservationsAndTimes->getObservationsSets( );

        // Save dependent variables values from first test case
        if ( testCase == 0 )
        {
            for ( auto currentSettings : dependentVariablesList )
            {
                ObservationDependentVariables variableType = currentSettings->variableType_;

                std::vector< std::vector< Eigen::MatrixXd > > currentDependentVariablesSortedPerSet;

                for ( auto set : idealObservationsAndTimes->getSingleObservationSets( ) )
                {
                    std::vector< std::pair< int, int > > compatibleIndicesAndSizes;
                    for ( auto it : set->getDependentVariableCalculator( )->getSettingsIndicesAndSizes( ) )
                    {
                        if ( it.second->areSettingsCompatible( currentSettings ) )
                        {
                            compatibleIndicesAndSizes.push_back( it.first );
                        }
                    }

                    std::vector< Eigen::MatrixXd > currentSetDependentVariablesPerSettings;
                    Eigen::MatrixXd currentSetFullDependentVariables = set->getObservationsDependentVariablesMatrix( );
                    for ( auto indicesIt : compatibleIndicesAndSizes )
                    {
                        Eigen::MatrixXd singleDependentVariableValues = Eigen::MatrixXd::Zero( currentSetFullDependentVariables.rows( ), indicesIt.second );
                        for ( unsigned int i = 0 ; i < currentSetFullDependentVariables.rows( ) ; i++ )
                        {
                            singleDependentVariableValues.block( i, 0, 1, indicesIt.second ) =
                                    currentSetFullDependentVariables.block( i, indicesIt.first, 1, indicesIt.second );
                        }
                        currentSetDependentVariablesPerSettings.push_back( singleDependentVariableValues );
                    }

                    if ( currentSetDependentVariablesPerSettings.size( ) > 0 )
                    {
                        currentDependentVariablesSortedPerSet.push_back( currentSetDependentVariablesPerSettings );
                    }
                }

                dependentVariablesReferenceValues[ variableType ] = currentDependentVariablesSortedPerSet;

            }
        }

        // Compare dependent variable values w.r.t. first test case
        if ( testCase == 1 )
        {
            for ( auto currentSettings : dependentVariablesList )
            {
                std::vector< std::vector< std::shared_ptr< ObservationDependentVariableSettings > > > compatibleSettingsList =
                        idealObservationsAndTimes->getCompatibleDependentVariablesSettingsList( currentSettings ).first;
                std::vector< std::vector< Eigen::MatrixXd > > dependentVariableValues =
                        idealObservationsAndTimes->getAllCompatibleDependentVariables( currentSettings ).first;

                // Retrieve reference values from first test case
                std::vector< std::vector< Eigen::MatrixXd > > referenceValues =
                        dependentVariablesReferenceValues.at( currentSettings->variableType_ );


                // Check that the number of single observation sets for which the given settings are relevant is consistent
                if ( currentSettings != moonAvoidanceAngleSettings )
                {
                    BOOST_CHECK( dependentVariableValues.size( ) == referenceValues.size( ) );
                }

                std::shared_ptr< ObservationCollectionParser > currentSettingsParser =
                        idealObservationsAndTimes->getCompatibleDependentVariablesSettingsList( currentSettings ).second;
                std::vector< std::shared_ptr< SingleObservationSet< double, double > > > currentSettingsSets =
                        idealObservationsAndTimes->getSingleObservationSets( currentSettingsParser );
                BOOST_CHECK( currentSettingsSets.size( ) == referenceValues.size( ) );

                // Parse dependent variable values per single observation set
                for ( unsigned int k = 0 ; k < dependentVariableValues.size( ) ; k++ )
                {
                    // Check that the number of settings per single observation sets is consistent
                    BOOST_CHECK( dependentVariableValues.at( k ).size( ) == referenceValues.at( k ).size( ) );

                    // Current single observation set
                    std::shared_ptr< SingleObservationSet< double, double > > currentSet = currentSettingsSets.at( k );

                    for ( unsigned int j = 0 ; j < dependentVariableValues.at( k ).size( ) ; j++ )
                    {
                        // Check that the dependent variable sizes and values are consistent
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariableValues.at( k ).at( j ), referenceValues.at( k ).at( j ), 1.0e-12 );

                        // Retrieve current complete dependent variables list
                        std::shared_ptr< ObservationDependentVariableSettings > currentCompleteSettings = compatibleSettingsList.at( k ).at( j );
                        Eigen::MatrixXd dependentVariablesFromCompleteSettings =
                                currentSet->getSingleDependentVariable( currentCompleteSettings );

                        // Check that the dependent variable sizes and values are consistent
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( dependentVariableValues.at( k ).at( j ), dependentVariablesFromCompleteSettings, 1.0e-12 );
                    }
                }
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

