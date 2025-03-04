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

#include <boost/test/unit_test.hpp>
#include <limits>
#include <string>

#include "tudat/basics/testMacros.h"
#include "tudat/io/readOdfFile.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;

using namespace tudat;

BOOST_AUTO_TEST_SUITE( test_process_odf_data )

BOOST_AUTO_TEST_CASE( testProcessOdfData )
{
    // Define ODF data paths
    std::string odFile = tudat::paths::getTudatTestDataPath( ) + "/odf07155.odf";

    // Laod raw ODF data
    std::vector< std::shared_ptr< input_output::OdfRawFileContents > > rawOdfDataVector;
    std::shared_ptr< input_output::OdfRawFileContents > rawOdfContents = std::make_shared< input_output::OdfRawFileContents >( odFile );

    // Process ODF file data
    std::string spacecraftName = "MESSENGER";
    std::shared_ptr< ProcessedOdfFileContents< Time > > processedOdfFileContents =
            std::make_shared< ProcessedOdfFileContents< Time > >( rawOdfContents, spacecraftName, true );

    std::vector< observation_models::ObservableType > obsType = processedOdfFileContents->getProcessedObservableTypes( );

    // Create data structure that handles Observed Data in Tudat
    std::shared_ptr< observation_models::ObservationCollection< long double, Time > > observedObservationCollection =
            observation_models::createOdfObservedObservationCollection< long double, Time >( processedOdfFileContents,
                                                                                             { dsn_n_way_averaged_doppler, n_way_range } );

    auto observationSets = observedObservationCollection->getObservationsSets( );

    // Check the observations for NWayRange and DsnNWayAveragedDoppler
    for( const auto &observableTypeEntry: observationSets )
    {
        observation_models::ObservableType observableType = observableTypeEntry.first;
        const auto &linkEndsMap = observableTypeEntry.second;

        for( const auto &linkEndsEntry: linkEndsMap )
        {
            const observation_models::LinkEnds &linkEnds = linkEndsEntry.first;
            const auto &observationSetVector = linkEndsEntry.second;

            for( const auto &observationSet: observationSetVector )
            {
                // Get the observations and times
                auto observations = observationSet->getObservations( );
                auto observationTimes = observationSet->getObservationTimes( );
                auto ancillarySettings = observationSet->getAncilliarySettings( );

                if( !observations.empty( ) && !observationTimes.empty( ) )
                {
                    // Check NWayRange observable
                    if( observableType == observation_models::n_way_range )
                    {
                        BOOST_CHECK_EQUAL( observation_models::getObservableName( observableType ), "NWayRange" );
                        BOOST_CHECK_EQUAL( linkEnds.at( transmitter ).stationName_, "DSS-14" );
                        BOOST_CHECK_EQUAL( linkEnds.at( retransmitter ).bodyName_, "MESSENGER" );
                        BOOST_CHECK_EQUAL( linkEnds.at( receiver ).stationName_, "DSS-14" );
                        BOOST_CHECK_CLOSE_FRACTION( double( observationTimes.at( 0 ) ), 234262616.184812933, 1e-9 );
                        BOOST_CHECK_CLOSE_FRACTION( observations[ 0 ].transpose( )( 0 ), 333589.366953747, 1e-9 );

                        // Check ancillary settings
                        if( ancillarySettings != nullptr )
                        {
                            BOOST_CHECK_EQUAL(
                                    ancillarySettings->getDoubleData( ).at( observation_models::reception_reference_frequency_band ),
                                    1.000000000 );
                            BOOST_CHECK_CLOSE_FRACTION(
                                    ancillarySettings->getDoubleData( ).at( observation_models::sequential_range_lowest_ranging_component ),
                                    14.000000000,
                                    1e-9 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::link_ends_delays ).size( ),
                                               3 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::frequency_bands ).size( ),
                                               2 );
                        }
                    }

                    // Check DsnNWayAveragedDoppler observable
                    if( observableType == observation_models::dsn_n_way_averaged_doppler )
                    {
                        BOOST_CHECK_EQUAL( observation_models::getObservableName( observableType ), "DsnNWayAveragedDoppler" );
                        BOOST_CHECK_EQUAL( linkEnds.at( transmitter ).stationName_, "DSS-14" );
                        BOOST_CHECK_EQUAL( linkEnds.at( retransmitter ).bodyName_, "MESSENGER" );
                        BOOST_CHECK_EQUAL( linkEnds.at( receiver ).stationName_, "DSS-14" );
                        BOOST_CHECK_CLOSE_FRACTION( double( observationTimes.at( 0 ) ), 234262457.184812993, 1e-9 );
                        BOOST_CHECK_CLOSE_FRACTION( observations[ 0 ].transpose( )( 0 ), 1.563486099, 1e-9 );

                        // Check ancillary settings
                        if( ancillarySettings != nullptr )
                        {
                            BOOST_CHECK_EQUAL(
                                    ancillarySettings->getDoubleData( ).at( observation_models::reception_reference_frequency_band ),
                                    1.000000000 );
                            BOOST_CHECK_CLOSE_FRACTION(
                                    ancillarySettings->getDoubleData( ).at( observation_models::doppler_integration_time ),
                                    60.000000000,
                                    1e-9 );
                            BOOST_CHECK_CLOSE_FRACTION(
                                    ancillarySettings->getDoubleData( ).at( observation_models::doppler_reference_frequency ),
                                    7177641534.000000000,
                                    1e-9 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::link_ends_delays ).size( ),
                                               3 );
                            BOOST_CHECK_EQUAL( ancillarySettings->getDoubleVectorData( ).at( observation_models::frequency_bands ).size( ),
                                               2 );
                        }
                    }
                }
                break;
            }
            break;
        }
    }

    BOOST_AUTO_TEST_SUITE_END( )
}
