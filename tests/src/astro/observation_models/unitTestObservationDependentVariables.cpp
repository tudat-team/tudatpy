/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "tudat/simulation/estimation.h"

//namespace tudat
//{
//namespace unit_tests
//{

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

//BOOST_AUTO_TEST_SUITE( test_observation_noise_models )

//// Function to conver
//double ignoreInputeVariable( std::function< double( ) > inputFreeFunction, const double dummyInput )
//{
//    return inputFreeFunction( );
//}

//! Test whether observation noise is correctly added when simulating noisy observations
int main( )
{
    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );

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
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::vector< LinkEnds > twoWayLinkEnds;

    // Define link ends to/from ground stations to Moon
    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Moon", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Moon", "" );
        stationReceiverLinkEnds.push_back( linkEnds );

        twoWayLinkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ retransmitter ] = std::make_pair( "Moon", "" );
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        twoWayLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrary) link ends for each observable
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );

    // Define observation settings for each observable/link ends combination
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Create observation settings
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                   currentObservable, currentLinkEndsList.at( i ) ) );
        }
    }

    // Create observation simulators
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > >  observationSimulators =
            createObservationSimulators( observationSettingsList, bodies );

    // Define osbervation times.
    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double observationInterval = 10.0;
    for( unsigned int i = 0; i < 14; i++ )
    {
        for( unsigned int j = 0; j < 4320; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    // Define observation simulation settings (observation type, link end, times and reference link end)
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                        std::make_shared< TabulatedObservationSimulationSettings< > >(
                            currentObservable, currentLinkEndsList.at( i ), baseTimeList, receiver ) );
        }
    }

    std::vector< std::shared_ptr< observation_models::ObservationViabilitySettings > > viabilitySettingsList;
    viabilitySettingsList.push_back( elevationAngleViabilitySettings(
                                         std::make_pair( "Earth", "Station1" ), 25.0 * mathematical_constants::PI / 180.0 ) );
    viabilitySettingsList.push_back( elevationAngleViabilitySettings(
                                         std::make_pair( "Earth", "Station2" ), 25.0 * mathematical_constants::PI / 180.0 ) );

    addViabilityToObservationSimulationSettings(
            measurementSimulationInput, viabilitySettingsList );

    // Define settings for dependent variables
    std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariableList;

    std::shared_ptr< ObservationDependentVariableSettings > elevationAngleSettings1 =
            std::make_shared< StationAngleObservationDependentVariableSettings >(
                station_elevation_angle, std::make_pair( "Earth", "Station1" ) );
    std::shared_ptr< ObservationDependentVariableSettings > azimuthAngleSettings1 =
            std::make_shared< StationAngleObservationDependentVariableSettings >(
                station_azimuth_angle, std::make_pair( "Earth", "Station1" ) );

//    std::shared_ptr< ObservationDependentVariableSettings > elevationAngleSettings2 =
//            std::make_shared< StationAngleObservationDependentVariableSettings >(
//                station_elevation_angle, std::make_pair( "Earth", "Station2" ) );

    dependentVariableList.push_back( elevationAngleSettings1 );
    dependentVariableList.push_back( azimuthAngleSettings1 );

    addDependentVariablesToObservationSimulationSettings(
                measurementSimulationInput, dependentVariableList, bodies );



    // Simulate noise-free observations
    std::shared_ptr< ObservationCollection< > > idealObservationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, observationSimulators, bodies );

    std::cout<<"Simulated observations "<<std::endl;
    std::map< double, Eigen::VectorXd > elevationAngles;
    std::map< double, Eigen::VectorXd > azimuthAngles;

    std::cout<<"Getting elevation angle"<<std::endl;
    elevationAngles = getDependentVariableResultList(
                idealObservationsAndTimes, elevationAngleSettings1,
                one_way_range );
    input_output::writeDataMapToTextFile( elevationAngles,
                                          "elevationAngles1_range.dat",
                                          "/home/dominic/Software/Tudat30Bundle/test-output/" );

    std::cout<<"Getting azimuth angle"<<std::endl;
    azimuthAngles = getDependentVariableResultList(
                idealObservationsAndTimes, azimuthAngleSettings1,
                one_way_range );
    input_output::writeDataMapToTextFile( azimuthAngles,
                                          "azimuthAngles1_range.dat",
                                          "/home/dominic/Software/Tudat30Bundle/test-output/" );

//    elevationAngles = getDependentVariableResultList(
//                idealObservationsAndTimes, elevationAngleSettings1,
//                one_way_range, stationReceiverLinkEnds[ 0 ] );
//    input_output::writeDataMapToTextFile( elevationAngles,
//                                          "elevationAngles2_range.dat",
//                                          "/home/dominic/Software/Tudat30Bundle/test-output/" );

//    elevationAngles = getDependentVariableResultList(
//                idealObservationsAndTimes, elevationAngleSettings2,
//                one_way_doppler );
//    input_output::writeDataMapToTextFile( elevationAngles,
//                                          "elevationAngles1_doppler.dat",
//                                          "/home/dominic/Software/Tudat30Bundle/test-output/" );

//    elevationAngles = getDependentVariableResultList(
//                idealObservationsAndTimes, elevationAngleSettings2,
//                one_way_doppler, stationReceiverLinkEnds[ 1 ] );
//    input_output::writeDataMapToTextFile( elevationAngles,
//                                          "elevationAngles2_doppler.dat",
//                                          "/home/dominic/Software/Tudat30Bundle/test-output/" );


}

//class ObservationDependentVariableWrapper
//{
//  ObservableType observableType_;

//  LinkEnds linkEnds;

//  std::map< double, Eigen::VectorXd > dependentVariables_;

//  std::vector< std::shared_ptr< ObservationDependentVariableSettings > > dependentVariableList;

//  std::vector< int > sizeIndices_;

//  std::map< double, Eigen::VectorXd > getSingleDependentVariables( std::shared_ptr< ObservationDependentVariableSettings > );



//};

//BOOST_AUTO_TEST_SUITE_END( )

//}

//}

