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

#include "tudat/interface/json/estimation/observation.h"
#include "tudat/interface/json/estimation/parameter.h"
#include "tudat/interface/json/propagation/acceleration.h"
#include "tudat/interface/json/support/valueConversions.h"
#include "tudat/interface/json/tests/unitTestSupport.h"
#include "tudat/simulation/estimation_setup/createEstimatableParameters.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/propagation_setup/createAccelerationModels.h"

namespace tudat
{

namespace unit_tests
{

#define INPUT( filename )                                                                          \
    ( json_interface::inputDirectory( ) / boost::filesystem::path( __FILE__ ).stem( ) / filename ) \
            .string( )

BOOST_AUTO_TEST_SUITE( test_json_observation )

// Test 1: sphericalHarmonicGravity
BOOST_AUTO_TEST_CASE( test_json_acceleration_sphericalHarmonicGravity )
{
    using namespace tudat;
    using namespace tudat::simulation_setup;
    // using namespace tudat::unit_tests;

    using namespace tudat::simulation_setup;
    using namespace tudat::observation_models;
    using namespace tudat::json_interface;
    using namespace tudat::estimatable_parameters;

    // Create ground stations from geodetic positions.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    // Create list of link ends in which station is receiver and in which station is transmitter
    // (with spacecraft other link end).
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;
    std::vector< LinkEnds > twoWayLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );
        stationReceiverLinkEnds.push_back( linkEnds );

        linkEnds.clear( );
        linkEnds[ receiver ] = std::make_pair( "Vehicle", "" );
        linkEnds[ reflector ] = std::make_pair( "Earth", groundStationNames.at( i ) );
        linkEnds[ transmitter ] = std::make_pair( "Vehicle", "" );

        twoWayLinkEnds.push_back( linkEnds );
    }

    // Define (arbitrarily) link ends to be used for 1-way range, 1-way doppler and angular position
    // observables
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    linkEndsPerObservable[ position_observable ].push_back( stationTransmitterLinkEnds[ 0 ] );

    linkEndsPerObservable[ one_way_differenced_range ].push_back( stationTransmitterLinkEnds[ 0 ] );

    linkEndsPerObservable[ n_way_range ].push_back( twoWayLinkEnds[ 0 ] );
    linkEndsPerObservable[ n_way_range ].push_back( twoWayLinkEnds[ 1 ] );

    //    n_way_range = 5,
    //    two_way_doppler = 6

    // Iterate over all observable types and associated link ends, and creatin settings for
    // observation
    observation_models::ObservationSettingsListPerLinkEnd observationSettingsMap;
    // std::map< LinkEnds, std::vector< std::shared_ptr < ObservationSettings > > >
    // observationSettingsMap;

    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator =
                 linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Define bias and light-time correction settings
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            std::shared_ptr< LightTimeCorrectionSettings > lightTimeCorrections;
            if( currentObservable == one_way_range )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                        Eigen::Vector1d::Constant( 1.0 ), true );
                std::vector< std::string > perturbingBodies = { "Mars", "Moon" };
                lightTimeCorrections =
                        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                perturbingBodies );
            }

            if( currentObservable == angular_position && i == 1 )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                        Eigen::Vector2d::Constant( 1.0E-6 ), false );
                std::vector< std::string > perturbingBodies = { "Mars", "Moon" };
                lightTimeCorrections =
                        std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                perturbingBodies );
            }

            if( currentObservable == angular_position && i == 1 )
            {
                std::vector< double > arcStartTimes;
                arcStartTimes.push_back( 0.0 );
                arcStartTimes.push_back( 1.0E6 );
                arcStartTimes.push_back( 3.0e6 );

                std::vector< Eigen::VectorXd > observationBiases;
                observationBiases.push_back( Eigen::Vector2d::Constant( 1.0E-6 ) );
                observationBiases.push_back( Eigen::Vector2d::Constant( 1.0E-7 ) );
                observationBiases.push_back( Eigen::Vector2d::Constant( 2.0E-6 ) );
                observationBiases[ 2 ]( 1 ) = 1.0E5;

                biasSettings = std::make_shared< ArcWiseConstantObservationBiasSettings >(
                        arcStartTimes, observationBiases, transmitter, true );
            }

            if( currentObservable == one_way_doppler && i == 0 )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                        Eigen::Vector1d::Constant( 1.0E-12 ), true );
            }

            if( currentObservable == one_way_doppler && i == 1 )
            {
                std::vector< double > arcStartTimes;
                arcStartTimes.push_back( 0.0 );
                arcStartTimes.push_back( 1.0E6 );
                arcStartTimes.push_back( 3.0e6 );

                std::vector< Eigen::VectorXd > observationBiases;
                observationBiases.push_back( Eigen::Vector1d::Constant( 1.0E-12 ) );
                observationBiases.push_back( Eigen::Vector1d::Constant( 1.0E-13 ) );
                observationBiases.push_back( Eigen::Vector1d::Constant( 2.0E-12 ) );

                biasSettings = std::make_shared< ArcWiseConstantObservationBiasSettings >(
                        arcStartTimes, observationBiases, transmitter, true );
            }

            if( currentObservable == one_way_range && i == 0 )
            {
                std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >(
                        Eigen::Vector1d::Constant( 1.0E-6 ), false ) );
                biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >(
                        Eigen::Vector1d::Constant( 2.0 ), true ) );
                biasSettings =
                        std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
            }

            std::shared_ptr< ObservationSettings > currentObservationSettings;
            if( currentObservable == one_way_doppler && i == 0 )
            {
                currentObservationSettings = std::make_shared< OneWayDopplerObservationSettings >(
                        lightTimeCorrections,
                        std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >(
                                "Earth" ),
                        std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Sun" ),
                        biasSettings );
            }
            else if( currentObservable == one_way_differenced_range )
            {
                currentObservationSettings =
                        std::make_shared< OneWayDifferencedRangeRateObservationSettings >(
                                [ = ]( const double ) { return 60.0; } );
            }
            else if( currentObservable == n_way_range )
            {
                std::vector< double > retransmissionTimes;
                retransmissionTimes.push_back( 1.0E-3 );
                if( i == 0 )
                {
                    currentObservationSettings = std::make_shared< NWayRangeObservationSettings >(
                            lightTimeCorrections,
                            3,
                            [ = ]( const double ) { return retransmissionTimes; },
                            biasSettings );
                }
                else if( i == 1 )
                {
                    std::vector< std::string > perturbingBodies = { "Mars", "Moon" };

                    std::vector< std::shared_ptr< ObservationSettings > >
                            oneWayRangeObservationSettings;
                    oneWayRangeObservationSettings.push_back( std::make_shared<
                                                              ObservationSettings >(
                            one_way_range,
                            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                    perturbingBodies ) ) );
                    oneWayRangeObservationSettings.push_back(
                            std::make_shared< ObservationSettings >(
                                    one_way_range,
                                    nullptr,
                                    std::make_shared< ConstantObservationBiasSettings >(
                                            Eigen::Vector1d::Constant( 1.0 ), true ) ) );

                    currentObservationSettings = std::make_shared< NWayRangeObservationSettings >(
                            oneWayRangeObservationSettings,
                            [ = ]( const double ) { return retransmissionTimes; },
                            biasSettings );
                }
            }
            else
            {
                currentObservationSettings = std::make_shared< ObservationSettings >(
                        currentObservable, lightTimeCorrections, biasSettings );
            }

            // Define settings for observable, no light-time corrections, and biases for selected
            // links
            observationSettingsMap[ currentLinkEndsList.at( i ) ].push_back(
                    currentObservationSettings );
        }
    }

    nlohmann::json jsonObject;
    to_json( jsonObject, observationSettingsMap );

    std::string fileName = INPUT( "observationOutput.json" );

    std::ofstream outputFile( fileName );
    outputFile << jsonObject.dump( 2 );
    outputFile.close( );

    observation_models::ObservationSettingsListPerLinkEnd observationSettingsMapFromJson;

    from_json( jsonObject, observationSettingsMapFromJson );

    BOOST_CHECK_EQUAL_JSON( observationSettingsMap, observationSettingsMapFromJson );

    nlohmann::json jsonObjectFromFile = parseJSONFile( fileName );

    observation_models::ObservationSettingsListPerLinkEnd observationSettingsMapFromFile =
            parseJSONFile< observation_models::ObservationSettingsListPerLinkEnd >( fileName );

    BOOST_CHECK_EQUAL_JSON( observationSettingsMap, observationSettingsMapFromFile );
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
