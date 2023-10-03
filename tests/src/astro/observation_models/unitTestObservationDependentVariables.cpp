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

//                std::cout<<i<<" "<<computedDependentVariables.at( currentTime )<<" "<<
//                         referenceDependentVariables.at( currentTime )<<std::endl;

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

//! Test whether observation noise is correctly added when simulating noisy observations
BOOST_AUTO_TEST_CASE( testObservationDependentVariables )
//int main( )
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
    double finalEphemerisTime = double( 1.0E7 + 3.0 * physical_constants::JULIAN_DAY );

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0, finalEphemerisTime + 3600.0 );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Earth",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    

    // Creatre ground stations
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );


    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 1.0 ).finished( ), geodetic_position );

    bodies.at( "Earth" )->getGroundStation( "Station1" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    bodies.at( "Earth" )->getGroundStation( "Station1" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );
    bodies.at( "Earth" )->getGroundStation( "Station1" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 3.0E9 ) );

    bodies.at( "Earth" )->getGroundStation( "Station2" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    bodies.at( "Earth" )->getGroundStation( "Station2" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );
    bodies.at( "Earth" )->getGroundStation( "Station2" )->setTransmittingFrequencyCalculator( std::make_shared< ConstantFrequencyInterpolator >( 3.0E9 ) );

    bodies.at( "Moon" )->setVehicleSystems( std::make_shared< system_models::VehicleSystems >( ) );
    bodies.at( "Moon" )->getVehicleSystems( )->setTransponderTurnaroundRatio( );

    // Define parameters.
    std::vector< LinkEnds > stationTransmitterOneWayLinkEnds;
    std::vector< LinkEnds > stationReceiverOneWayLinkEnds;

    std::vector< LinkEnds > stationReceiverTwoWayLinkEnds;
    std::vector< LinkEnds > stationRetransmitterTwoWayLinkEnds;
    std::vector< LinkEnds > stationReceiverThreeWayLinkEnds;
    std::vector< LinkEnds > stationTransmitterThreeWayLinkEnds;

    std::vector< LinkEnds > stationReceiverRelativeLinkEnds;
    std::vector< LinkEnds > stationReceiverOppositeRelativeLinkEnds;


//    std::vector< LinkEnds > stationTransmitterRelativeLinkEnds;
//    std::vector< LinkEnds > stationTransmitter2RelativeLinkEnds;


    // Define link ends to/from ground stations to Me de    oon
    LinkEnds linkEnds;
    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {

        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Moon", "" ) );
        stationTransmitterOneWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Moon", "" ) );
        stationReceiverOneWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "Moon", "" ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        stationReceiverTwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Moon", "" ) );
        linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Moon", "" ) );
        stationRetransmitterTwoWayLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Moon", "" ) );
        linkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "Mars", "" ) );
        stationReceiverRelativeLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "Moon", "" ) );
        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Mars", "" ) );
        stationReceiverOppositeRelativeLinkEnds.push_back( linkEnds );


//        linkEnds.clear( );
//        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
//        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Moon", "" ) );
//        linkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "Mars", "" ) );
//        stationTransmitterRelativeLinkEnds.push_back( linkEnds );
//
//        linkEnds.clear( );
//        linkEnds[ transmitter2 ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( i ) ) );
//        linkEnds[ receiver ] = LinkEndId( std::make_pair( "Moon", "" ) );
//        linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Mars", "" ) );
//        stationTransmitter2RelativeLinkEnds.push_back( linkEnds );
    }

    linkEnds.clear( );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "Moon", "" ) );
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );
    stationReceiverThreeWayLinkEnds.push_back( linkEnds );

    linkEnds.clear( );
    linkEnds[ receiver ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 1 ) ) );
    linkEnds[ retransmitter ] = LinkEndId( std::make_pair( "Moon", "" ) );
    linkEnds[ transmitter ] = LinkEndId( std::make_pair( "Earth", groundStationNames.at( 0 ) ) );
    stationTransmitterThreeWayLinkEnds.push_back( linkEnds );





    std::map<double, Eigen::VectorXd> referenceElevationAngles;
    std::map<double, Eigen::VectorXd> referenceAzimuthAngles;
    std::map<double, Eigen::VectorXd> referenceTargetRanges;

    std::vector< std::map< double, Eigen::VectorXd > > referenceReceiverDependentVariableResults;
    std::vector< std::map< double, Eigen::VectorXd > > referenceTransmitterDependentVariableResults;

    for( unsigned int currentObservableTestCase = 0; currentObservableTestCase < 9; currentObservableTestCase++ )
    {
        bool isDifferencedObservable = false;
        int geometryType = -1;
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

//            relative_angular_position = 9,

//            dsn_one_way_averaged_doppler = 12,


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

        if( isDifferencedObservable )
        {
            numberOfLinkEndCases *= 2;
        }

        for( int currentLinkEndCase = 0; currentLinkEndCase < numberOfLinkEndCases; currentLinkEndCase++ )
        {
            bool compareAgainstReceiver = -1;
            LinkEndType referenceLinkEnd = unidentified_link_end;
            IntegratedObservationPropertyHandling integratedObservableHandling = interval_undefined;
            LinkEndType originatingLinkEndRole = unidentified_link_end;

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
//                case 2:
//                    currentLinkEnds = stationTransmitterRelativeLinkEnds.at( 0 );
//                    referenceLinkEnd = transmitter;
//                    originatingLinkEndRole = receiver;
//                    compareAgainstReceiver = false;
//                    break;
//                case 3:
//                    currentLinkEnds = stationTransmitter2RelativeLinkEnds.at( 0 );
//                    referenceLinkEnd = transmitter2;
//                    originatingLinkEndRole = receiver;
//                    compareAgainstReceiver = false;
//                    break;
                default:
                    throw std::runtime_error( "Error in observation dependent variable unit test C " );
                }
                break;
            }
            }

            if( !( currentObservableType == dsn_n_way_averaged_doppler && referenceLinkEnd != receiver ) )
            {

                std::cout<<currentLinkEnds.size( )<<" "<<geometryType<<" "<<currentLinkEndCase<<std::endl;
                // Define (arbitrary) link ends for each observable
                std::map<ObservableType, std::vector<LinkEnds> > linkEndsPerObservable;
                linkEndsPerObservable[ currentObservableType ].push_back( currentLinkEnds );

                // Define observation settings for each observable/link ends combination
                std::vector<std::shared_ptr<ObservationModelSettings> > observationSettingsList;
                for ( std::map<ObservableType, std::vector<LinkEnds> >::iterator
                          linkEndIterator = linkEndsPerObservable.begin( );
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

    //            std::vector<std::shared_ptr<ObservationViabilitySettings> > viabilitySettingsList;
    //            viabilitySettingsList.push_back( elevationAngleViabilitySettings(
    //                std::make_pair( "Earth", "Station1" ), 25.0 * mathematical_constants::PI / 180.0 ));
    //            viabilitySettingsList.push_back( elevationAngleViabilitySettings(
    //                std::make_pair( "Earth", "Station2" ), 25.0 * mathematical_constants::PI / 180.0 ));
    //
    //            addViabilityToObservationSimulationSettings(
    //                measurementSimulationInput, viabilitySettingsList );

                // Define settings for dependent variables
                std::vector<std::shared_ptr<ObservationDependentVariableSettings> > dependentVariableList;

                std::shared_ptr<ObservationDependentVariableSettings> elevationAngleSettings1 =
                    std::make_shared<StationAngleObservationDependentVariableSettings>(
                        station_elevation_angle, LinkEndId( std::make_pair( "Earth", "Station1" )),
                        referenceLinkEnd, integratedObservableHandling, originatingLinkEndRole );
                std::shared_ptr<ObservationDependentVariableSettings> azimuthAngleSettings1 =
                    std::make_shared<StationAngleObservationDependentVariableSettings>(
                        station_azimuth_angle, LinkEndId( std::make_pair( "Earth", "Station1" )),
                        referenceLinkEnd, integratedObservableHandling, originatingLinkEndRole );
                std::shared_ptr<ObservationDependentVariableSettings> targetRangeSettings1 =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        target_range, referenceLinkEnd, originatingLinkEndRole, integratedObservableHandling );
                std::shared_ptr<ObservationDependentVariableSettings> targetInverseRangeSettings1 =
                    std::make_shared<InterlinkObservationDependentVariableSettings>(
                        target_range, originatingLinkEndRole, referenceLinkEnd, integratedObservableHandling );

                dependentVariableList.push_back( elevationAngleSettings1 );
                dependentVariableList.push_back( azimuthAngleSettings1 );
                dependentVariableList.push_back( targetRangeSettings1 );
                dependentVariableList.push_back( targetInverseRangeSettings1 );



                addDependentVariablesToObservationSimulationSettings(
                    measurementSimulationInput, dependentVariableList, bodies );

                // Simulate noise-free observations
                std::shared_ptr<ObservationCollection<> > idealObservationsAndTimes = simulateObservations<double, double>(
                    measurementSimulationInput, observationSimulators, bodies );

                if ( currentObservableTestCase == 0 )
                {

                    std::map<double, Eigen::VectorXd> elevationAngles1 = getDependentVariableResultList(
                        idealObservationsAndTimes, elevationAngleSettings1, currentObservableType );
                    std::map<double, Eigen::VectorXd> azimuthAngles1 = getDependentVariableResultList(
                        idealObservationsAndTimes, azimuthAngleSettings1, currentObservableType );
                    std::map<double, Eigen::VectorXd> targetRanges1 = getDependentVariableResultList(
                        idealObservationsAndTimes, targetRangeSettings1, currentObservableType );
                    std::map<double, Eigen::VectorXd> targetInverseRanges1 = getDependentVariableResultList(
                        idealObservationsAndTimes, targetInverseRangeSettings1, currentObservableType );

                    if( currentLinkEndCase == 0 )
                    {
                        referenceReceiverDependentVariableResults.push_back( elevationAngles1 );
                        referenceReceiverDependentVariableResults.push_back( azimuthAngles1 );
                        referenceReceiverDependentVariableResults.push_back( targetRanges1 );
                        referenceReceiverDependentVariableResults.push_back( targetInverseRanges1 );
                    }
                    else if( currentLinkEndCase == 1 )
                    {
                        referenceTransmitterDependentVariableResults.push_back( elevationAngles1 );
                        referenceTransmitterDependentVariableResults.push_back( azimuthAngles1 );
                        referenceTransmitterDependentVariableResults.push_back( targetRanges1 );
                        referenceTransmitterDependentVariableResults.push_back( targetInverseRanges1 );
                    }

                    std::map<LinkEnds, int>
                        linkEndIdentifiers = idealObservationsAndTimes->getLinkEndIdentifierMap( );
                    std::vector<int> linkEndIds = idealObservationsAndTimes->getConcatenatedLinkEndIds( );

                    int numberOfLinkEnds1Observations = utilities::countNumberOfOccurencesInVector<int>(
                        linkEndIds, linkEndIdentifiers.at( currentLinkEnds ));

                    BOOST_CHECK_EQUAL( elevationAngles1.size( ), numberOfLinkEnds1Observations );
                    BOOST_CHECK_EQUAL( azimuthAngles1.size( ), numberOfLinkEnds1Observations );
                    BOOST_CHECK_EQUAL( targetRanges1.size( ), numberOfLinkEnds1Observations );
                    BOOST_CHECK_EQUAL( targetInverseRanges1.size( ), numberOfLinkEnds1Observations );

                    std::shared_ptr< PointingAnglesCalculator> pointingAnglesCalculator1 =
                        bodies.at( "Earth" )->getGroundStation( "Station1" )->getPointingAnglesCalculator( );
                    std::shared_ptr< PointingAnglesCalculator> pointingAnglesCalculator2 =
                        bodies.at( "Earth" )->getGroundStation( "Station2" )->getPointingAnglesCalculator( );

                    std::shared_ptr<ObservationModel<1, double, double> > observationModel1 =
                        std::dynamic_pointer_cast<ObservationSimulator<1, double, double> >(
                            observationSimulators.at( 0 ))->getObservationModel( currentLinkEnds );

                    std::vector<double> linkEndTimes;
                    std::vector<Eigen::Matrix<double, 6, 1> > linkEndStates;

                    for ( auto it: elevationAngles1 )
                    {
                        double currentTime = it.first;
                        double currentElevation = elevationAngles1.at( currentTime )( 0 );
                        double currentAzimuth = azimuthAngles1.at( currentTime )( 0 );
                        double targetRange = targetRanges1.at( currentTime )( 0 );
                        double targetInverseRange = targetInverseRanges1.at( currentTime )( 0 );

                        observationModel1->computeIdealObservationsWithLinkEndData(
                            currentTime, referenceLinkEnd, linkEndTimes, linkEndStates );
                        Eigen::Vector3d vectorToTarget = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 )).segment( 0, 3 );
                        if( referenceLinkEnd == transmitter )
                        {
                            vectorToTarget *= -1.0;
                        }
                        double referenceTime = ( referenceLinkEnd == transmitter ) ? linkEndTimes.at( 0 ) : linkEndTimes.at( 1 );

                        double elevationAngle = pointingAnglesCalculator1->calculateElevationAngleFromInertialVector(
                            vectorToTarget, referenceTime );
                        BOOST_CHECK_SMALL(( elevationAngle - currentElevation ), std::numeric_limits<double>::epsilon( ));

                        double azimuthAngle = pointingAnglesCalculator1->calculateAzimuthAngleFromInertialVector(
                            vectorToTarget, referenceTime );
                        BOOST_CHECK_SMALL(( azimuthAngle - currentAzimuth ), std::numeric_limits<double>::epsilon( ));

                        BOOST_CHECK_SMALL(( targetRange - vectorToTarget.norm( )),
                                          std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));
                        BOOST_CHECK_SMALL(( targetInverseRange - targetRange ),
                                          std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));
                        BOOST_CHECK_SMALL(( targetInverseRange - vectorToTarget.norm( )),
                                          std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));

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

