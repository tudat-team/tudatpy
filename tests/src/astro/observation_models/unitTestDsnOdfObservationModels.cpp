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
#include "tudat/astro/orbit_determination/processOdfFile.h"

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

    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );

    // Define light-time perturbing bodies
    std::vector< std::string > lightTimePerturbingBodies = { "Earth", "Sun" };

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create observed observation collection
    std::shared_ptr< observation_models::ObservationCollection< > > observedObservationCollection =
            orbit_determination::createOdfObservedObservationCollection(
                    orbit_determination::processOdfFileContents( input_output::readOdfFile(
//                            "/Users/pipas/Documents/mro-rawdata-odf/mromagr2009_332_1945xmmmv1.odf" ) ),
                            "/Users/pipas/Documents/dsn_trk-2-18/odf07155.dat" ) ),
                    bodies );

    std::cout << std::endl << "Observation type start and size:" << std::endl;
    std::map< observation_models::ObservableType, std::pair< int, int > > observationTypeStartAndSize =
            observedObservationCollection->getObservationTypeStartAndSize( );
    for ( auto it = observationTypeStartAndSize.begin(); it != observationTypeStartAndSize.end(); ++it )
    {
        std::cout << it->first << " " << std::get<0>(it->second) << " " << std::get<1>(it->second) << std::endl;
    }

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
            orbit_determination::createOdfObservationSimulationSettingsList< double >(
                    observedObservationCollection );


    std::shared_ptr< observation_models::ObservationCollection< double, double > >
            simulatedObservationCollection = simulation_setup::simulateObservations< double, double >(
                    observationSimulationSettings, observationSimulators, bodies );


}

BOOST_AUTO_TEST_SUITE_END( )

}

}