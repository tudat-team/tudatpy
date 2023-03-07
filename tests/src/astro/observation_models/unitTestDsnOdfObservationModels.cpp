/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
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

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_dsn_odf_observation_models )

BOOST_AUTO_TEST_CASE( testDsnNWayAveragedDopplerModel )
{

    spice_interface::loadStandardSpiceKernels(
            { "/Users/pipas/Documents/messenger-spice/msgr_040803_150430_150430_od431sc_2.bsp" } );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_040803_080216_120401.bsp" );
//    spice_interface::loadSpiceKernelInTudat( "/Users/pipas/Documents/messenger-spice/msgr_040803_150430_150430_od431sc_2.bsp" );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mercury" );
    bodiesToCreate.push_back( "Venus" );

    // Define light-time perturbing bodies
    std::vector< std::string > lightTimePerturbingBodies = { "Earth", "Sun", "Mercury", "Venus" };

    // Specify initial time
    double initialEphemerisTime = 234224663.2 - 10.0 * 86400.0;
    double finalEphemerisTime = 234349306.2 + 10.0 * 86400.0;
    double ephemerisTimeStep = 300.0;
    double buffer = 10.0 * ephemerisTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings bodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    bodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create spacecraft
    std::string spacecraftName = "Messenger";
    bodySettings.addSettings( spacecraftName );
    bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared< InterpolatedSpiceEphemerisSettings >(
                    initialEphemerisTime - buffer, finalEphemerisTime + buffer, ephemerisTimeStep, "SSB", "ECLIPJ2000",
                    std::make_shared< interpolators::LagrangeInterpolatorSettings >( 6 ), spacecraftName );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Read and process ODF file data
    std::shared_ptr< ProcessedOdfFileContents > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents >(
                    input_output::readOdfFile( "/Users/pipas/Documents/dsn_trk-2-18/odf07155.dat" ), // "/Users/pipas/Documents/mro-rawdata-odf/mromagr2009_332_1945xmmmv1.odf"
                    bodies.getBody( "Earth" ), true, spacecraftName );

    std::cout << std::endl << "Spacecraft name: " << processedOdfFileContents->getSpacecraftName( ) << std::endl;
    std::vector< std::string > groundStationsNames = processedOdfFileContents->getGroundStationsNames( );
    std::cout << "Ground stations: ";
    for ( unsigned int i = 0; i < groundStationsNames.size( ); ++i )
    {
        std::cout << groundStationsNames.at( i ) << " ";
    }
    std::cout << std::endl;

    std::cout << std::setprecision( 15 );
    std::cout << "Start time: " << processedOdfFileContents->getStartAndEndTime( bodies ).first << std::endl;
    std::cout << "End time: " << processedOdfFileContents->getStartAndEndTime( bodies ).second << std::endl;

    std::vector< observation_models::ObservableType > observableTypes = processedOdfFileContents->getProcessedObservableTypes( );
    std::cout << "Observable types: ";
    for ( unsigned int i = 0; i < observableTypes.size( ); ++i )
    {
        std::cout << observableTypes.at( i ) << " ";
    }
    std::cout << std::endl;

    // Create observed observation collection
    std::shared_ptr< observation_models::ObservationCollection< > > observedObservationCollection =
            observation_models::createOdfObservedObservationCollection(
                    processedOdfFileContents, bodies );

    std::cout << std::endl << "Observation type start and size:" << std::endl;
    std::map< observation_models::ObservableType, std::pair< int, int > > observationTypeStartAndSize =
            observedObservationCollection->getObservationTypeStartAndSize( );
    for ( auto it = observationTypeStartAndSize.begin(); it != observationTypeStartAndSize.end(); ++it )
    {
        std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << std::endl;
    }

    std::cout << "Observed observables: "<< observedObservationCollection->getObservationVector( ) << std::endl;

    // Create computed observation collection
    std::vector< std::shared_ptr< observation_models::ObservationModelSettings > > observationModelSettingsList;

    std::vector< std::shared_ptr< observation_models::LightTimeCorrectionSettings > > lightTimeCorrectionSettings =
            { std::make_shared< observation_models::FirstOrderRelativisticLightTimeCorrectionSettings >(
                    lightTimePerturbingBodies ) };

    std::map < observation_models::ObservableType, std::vector< observation_models::LinkEnds > > linkEndsPerObservable =
            observedObservationCollection->getLinkEndsPerObservableType( );
    for ( auto it = linkEndsPerObservable.begin(); it != linkEndsPerObservable.end(); ++it )
    {
        for ( unsigned int i = 0; i < it->second.size(); ++i )
        {
//            observationModelSettingsList.push_back(
//                    std::make_shared< observation_models::ObservationModelSettings >(
//                            it->first, it->second.at( i ), lightTimeCorrectionSettings, nullptr, nullptr ) );
            if ( it->first == observation_models::dsn_n_way_averaged_doppler )
            {
                observationModelSettingsList.push_back(
                    std::make_shared< observation_models::DsnNWayAveragedDopplerObservationSettings >(
                            it->second.at( i ), lightTimeCorrectionSettings, nullptr, nullptr ) );
            }
        }
    }

    std::vector< std::shared_ptr< observation_models::ObservationSimulatorBase< double, double > > >
            observationSimulators = observation_models::createObservationSimulators(
                    observationModelSettingsList, bodies );


    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > observationSimulationSettings =
            observation_models::createOdfObservationSimulationSettingsList< double >(
                    observedObservationCollection );


    std::shared_ptr< observation_models::ObservationCollection< double, double > >
            simulatedObservationCollection = simulation_setup::simulateObservations< double, double >(
                    observationSimulationSettings, observationSimulators, bodies );

    std::cout << std::endl << "Observation type start and size:" << std::endl;
    observationTypeStartAndSize = simulatedObservationCollection->getObservationTypeStartAndSize( );
    for ( auto it = observationTypeStartAndSize.begin(); it != observationTypeStartAndSize.end(); ++it )
    {
        std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << std::endl;
    }

    std::cout << simulatedObservationCollection->getObservationVector( ) << std::endl;


//    Eigen::Matrix< double, Eigen::Dynamic, 1 > simulatedObservations =
//            simulatedObservationCollection->getObservationVector( );
//
//    Eigen::MatrixXd observations ( simulatedObservationCollection->getObservationVector( ).size( ), 4 );
//    observations.col( 1 ) = simulatedObservationCollection->getObservationVector( );
//    observations.col( 3 ) = observedObservationCollection->getObservationVector( );
//
//    for ( unsigned int i = 0; i < simulatedObservationCollection->getObservationVector( ).size( ); ++i )
//    {
//        observations( i, 0 ) = simulatedObservationCollection->getConcatenatedTimeVector( ).at( i );
//        observations( i, 2 ) = observedObservationCollection->getConcatenatedTimeVector( ).at( i );
//    }
//
//    std::ofstream file("/Users/pipas/tudatpy-testing/observations.txt");
//    file << std::setprecision( 15 ) << observations;
//    file.close();

}

BOOST_AUTO_TEST_SUITE_END( )

}

}