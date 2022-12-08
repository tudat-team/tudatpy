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
#include <boost/make_shared.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat;


BOOST_AUTO_TEST_SUITE( test_n_way_range_model )

std::vector< double > getRetransmissionDelays( const double evaluationTime, const int numberOfRetransmitters )
{
    std::vector< double > retransmissionDelays;
    if( evaluationTime < 1.5E5 )
    {
        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( 0.0 );
        }
    }
    else
    {
        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( evaluationTime * 5.0E-17 * static_cast< double >( i + 1 ) );
        }
    }
    return retransmissionDelays;
}

BOOST_AUTO_TEST_CASE( testNWayRateRangeModel )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Moon" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    

    // Create ground stations
    std::pair< std::string, std::string > earthStationStation = std::pair< std::string, std::string >( "Earth", "EarthStation" );
    std::pair< std::string, std::string > earthStationStation2 = std::pair< std::string, std::string >( "Earth", "EarthStation2" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );
    createGroundStation( bodies.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "EarthStation", ( Eigen::Vector3d( ) << 1.0, 0.1, -1.4 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "EarthStation2", ( Eigen::Vector3d( ) << -30.0, 1.2, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( earthStationStation );
    groundStations.push_back( earthStationStation2 );
    groundStations.push_back( mslStation );

    // Define list of observation times for which to check model
    std::vector< double > observationTimes;
    observationTimes.push_back( 0.0 );
    observationTimes.push_back( 2.0E5 );
    observationTimes.push_back( 3.0E5 );

    // Define link ends for observations.

    // Test 2-way model
    {
        // Define link ends for 2-way model and constituent one-way models
        LinkEnds twoWayLinkEnds;
        twoWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation"  );
        twoWayLinkEnds[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );
        twoWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2"  );
        int numberOfLinkEnds = twoWayLinkEnds.size( );

        LinkEnds uplinkLinkEnds;
        uplinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation" );
        uplinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );

        LinkEnds downlinkLinkEnds;
        downlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );
        downlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2" );

        // Create light-time correction settings.
        std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
        lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                   lightTimePerturbingBodies ) );


        // Create observation settings for 2-way model and constituent one-way models
        std::shared_ptr< ObservationModelSettings > uplinkObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_range, uplinkLinkEnds, lightTimeCorrectionSettings );
        std::shared_ptr< ObservationModelSettings > downlinkObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_range, downlinkLinkEnds, lightTimeCorrectionSettings );
        std::vector< std::shared_ptr< ObservationModelSettings > > twoWayLinkSettings;
        twoWayLinkSettings.push_back( uplinkObservableSettings );
        twoWayLinkSettings.push_back( downlinkObservableSettings );

        std::shared_ptr< NWayRangeObservationSettings > twoWayRangeObservableSettings = std::make_shared< NWayRangeObservationSettings >
                ( twoWayLinkSettings );

        std::shared_ptr< NWayDifferencedRangeObservationSettings > twoWayRangeRateObservableSettings =
                std::make_shared< NWayDifferencedRangeObservationSettings >(
                    twoWayLinkEnds, lightTimeCorrectionSettings, nullptr );

        // Create observation models
        std::shared_ptr< ObservationModel< 1, double, double > > twoWayRangeObservationModel =
                ObservationModelCreator< 1, double, double >::createObservationModel(
                    twoWayRangeObservableSettings, bodies );

        std::shared_ptr< ObservationModel< 1, double, double > > twoWayRangeRateObservationModel =
                ObservationModelCreator< 1, double, double >::createObservationModel(
                    twoWayRangeRateObservableSettings, bodies );


        // Define link ends time and state vectors for  2-way model and constituent one-way models
        std::vector< double > rangeStartTimes;
        std::vector< Eigen::Vector6d > rangeStartStates;

        std::vector< double > rangeEndTimes;
        std::vector< Eigen::Vector6d > rangeEndStates;

        std::vector< double > rangeRateTimes;
        std::vector< Eigen::Vector6d > rangeRateStates;

        // Define variables for range measurements
        Eigen::VectorXd twoWayRangeStart;
        Eigen::VectorXd twoWayRangeEnd;
        Eigen::VectorXd twoWayRangeRate;

        // Iterate over each 2-way link end as reference link end

        std::vector< LinkEndType > referenceLinkEnds = { transmitter, receiver };
        std::vector< LinkEndType > fullLinkEnds = { transmitter, retransmitter, receiver };

        for( unsigned int observationTimeNumber = 0; observationTimeNumber < observationTimes.size( ); observationTimeNumber++ )
        {
            double observationTime = observationTimes.at( observationTimeNumber );
            for( unsigned int i = 0; i < referenceLinkEnds.size( ); i++ )
            {
                // Compute 2-way range
                twoWayRangeStart = twoWayRangeObservationModel->computeObservationsWithLinkEndData(
                            observationTime - 30.0, referenceLinkEnds.at( i ), rangeStartTimes, rangeStartStates );
                twoWayRangeEnd = twoWayRangeObservationModel->computeObservationsWithLinkEndData(
                            observationTime + 30.0, referenceLinkEnds.at( i ), rangeEndTimes, rangeEndStates );
                twoWayRangeRate = twoWayRangeRateObservationModel->computeObservationsWithLinkEndData(
                            observationTime, referenceLinkEnds.at( i ), rangeRateTimes, rangeRateStates,
                            getAveragedDopplerAncilliarySettings( 60.0 ) );

                double rangeRateError =  ( twoWayRangeRate - ( - twoWayRangeStart + twoWayRangeEnd ) / 60.0 )( 0 );
                BOOST_CHECK_SMALL( rangeRateError, 1.0E-9 );

                for( unsigned int j = 0; j < fullLinkEnds.size( ); j++ )
                {
                    std::vector< int > linkEndIndices =
                            getLinkEndIndicesForLinkEndTypeAtObservable(
                                n_way_differenced_range, fullLinkEnds.at( j ), twoWayLinkEnds.size( ) );
                    int numberOfRangeTimes = rangeStartTimes.size( );

                    for( unsigned int k = 0; k < linkEndIndices.size( ); k++ )
                    {
                        std::cout<<i<<" "<<j<<" "<<k<<std::endl;
                        BOOST_CHECK_SMALL( rangeStartTimes.at( k ) - rangeRateTimes.at( k ), 1.0E-9 );
                        BOOST_CHECK_SMALL( rangeEndTimes.at( k ) - rangeRateTimes.at( k + numberOfRangeTimes ), 1.0E-9 );

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeStartStates.at( k ), rangeRateStates.at( k ), ( 100.0 * std::numeric_limits< double >::epsilon( ) ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeEndStates.at( k ), rangeRateStates.at( k + numberOfRangeTimes ), ( 100.0 * std::numeric_limits< double >::epsilon( ) ) );
                    }
                }
            }
        }
    }

    // Test 4-way model
    {
        // Define link ends for 4-way model and constituent one-way models
        LinkEnds fourWayLinkEnds;
        fourWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation"  );
        fourWayLinkEnds[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );
        fourWayLinkEnds[ reflector2 ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2"  );
        fourWayLinkEnds[ reflector3 ] = std::make_pair< std::string, std::string >( "Moon" , ""  );
        fourWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );

        LinkEnds firstlinkLinkEnds;
        firstlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation" );
        firstlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );

        LinkEnds secondlinkLinkEnds;
        secondlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );
        secondlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2" );

        LinkEnds thirdlinkLinkEnds;
        thirdlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2" );
        thirdlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Moon" , "" );

        LinkEnds fourthlinkLinkEnds;
        fourthlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Moon" , "" );
        fourthlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );


        // Create light-time correction settings.
        std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
        lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                   lightTimePerturbingBodies ) );


        // Create observation settings for 4-way model and constituent one-way models
        std::shared_ptr< ObservationModelSettings > firstlinkObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_range, firstlinkLinkEnds, lightTimeCorrectionSettings );
        std::shared_ptr< ObservationModelSettings > secondlinkObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_range, secondlinkLinkEnds, lightTimeCorrectionSettings );
        std::shared_ptr< ObservationModelSettings > thirdlinkObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_range, thirdlinkLinkEnds, lightTimeCorrectionSettings );
        std::shared_ptr< ObservationModelSettings > fourthlinkObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_range, fourthlinkLinkEnds, lightTimeCorrectionSettings );

        std::vector< std::shared_ptr< ObservationModelSettings > > fourWayLinkSettings;
        fourWayLinkSettings.push_back( firstlinkObservableSettings );
        fourWayLinkSettings.push_back( secondlinkObservableSettings );
        fourWayLinkSettings.push_back( thirdlinkObservableSettings );
        fourWayLinkSettings.push_back( fourthlinkObservableSettings );
        std::shared_ptr< NWayRangeObservationSettings > fourWayRangeObservableSettings =
                std::make_shared< NWayRangeObservationSettings >
                ( fourWayLinkSettings );

        std::shared_ptr< NWayDifferencedRangeObservationSettings > fourWayRangeRateObservableSettings =
                std::make_shared< NWayDifferencedRangeObservationSettings >(
                    fourWayLinkEnds, lightTimeCorrectionSettings, nullptr );


        // Create observation models
        std::shared_ptr< ObservationModel< 1, double, double > > fourWayRangeObservationModel =
                ObservationModelCreator< 1, double, double >::createObservationModel(
                    fourWayRangeObservableSettings, bodies );

        std::shared_ptr< ObservationModel< 1, double, double > > fourWayRangeRateObservationModel =
                ObservationModelCreator< 1, double, double >::createObservationModel(
                    fourWayRangeRateObservableSettings, bodies );


        // Define link ends time and state vectors for  4-way model and constituent one-way models
        std::vector< double > rangeStartTimes;
        std::vector< Eigen::Vector6d > rangeStartStates;

        std::vector< double > rangeEndTimes;
        std::vector< Eigen::Vector6d > rangeEndStates;

        std::vector< double > rangeRateTimes;
        std::vector< Eigen::Vector6d > rangeRateStates;

        // Define variables for range measurements
        Eigen::VectorXd fourWayRangeStart;
        Eigen::VectorXd fourWayRangeEnd;
        Eigen::VectorXd fourWayRangeRate;

        // Iterate over each 4-way link end as reference link end

        std::vector< LinkEndType > referenceLinkEnds = { transmitter, receiver };
        std::vector< LinkEndType > fullLinkEnds = { transmitter, retransmitter, retransmitter2, retransmitter3, receiver };

        for( unsigned int observationTimeNumber = 0; observationTimeNumber < observationTimes.size( ); observationTimeNumber++ )
        {
            double observationTime = observationTimes.at( observationTimeNumber );
            for( unsigned int i = 0; i < referenceLinkEnds.size( ); i++ )
            {
                // Compute 4-way range
                fourWayRangeStart = fourWayRangeObservationModel->computeObservationsWithLinkEndData(
                            observationTime - 30.0, referenceLinkEnds.at( i ), rangeStartTimes, rangeStartStates );
                fourWayRangeEnd = fourWayRangeObservationModel->computeObservationsWithLinkEndData(
                            observationTime + 30.0, referenceLinkEnds.at( i ), rangeEndTimes, rangeEndStates );
                fourWayRangeRate = fourWayRangeRateObservationModel->computeObservationsWithLinkEndData(
                            observationTime, referenceLinkEnds.at( i ), rangeRateTimes, rangeRateStates,
                            getAveragedDopplerAncilliarySettings( 60.0 ) );

                double rangeRateError =  ( fourWayRangeRate - ( - fourWayRangeStart + fourWayRangeEnd ) / 60.0 )( 0 );
                BOOST_CHECK_SMALL( rangeRateError, 1.0E-9 );

                for( unsigned int j = 0; j < fullLinkEnds.size( ); j++ )
                {
                    std::vector< int > linkEndIndices =
                            getLinkEndIndicesForLinkEndTypeAtObservable(
                                n_way_differenced_range, fullLinkEnds.at( j ), fourWayLinkEnds.size( ) );
                    int numberOfRangeTimes = rangeStartTimes.size( );

                    for( unsigned int k = 0; k < linkEndIndices.size( ); k++ )
                    {
                        BOOST_CHECK_SMALL( rangeStartTimes.at( k ) - rangeRateTimes.at( k ), 1.0E-9 );
                        BOOST_CHECK_SMALL( rangeEndTimes.at( k ) - rangeRateTimes.at( k + numberOfRangeTimes ), 1.0E-9 );

                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeStartStates.at( k ), rangeRateStates.at( k ), ( 100.0 * std::numeric_limits< double >::epsilon( ) ) );
                        TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeEndStates.at( k ), rangeRateStates.at( k + numberOfRangeTimes ), ( 100.0 * std::numeric_limits< double >::epsilon( ) ) );
                    }
                }
            }
        }
    }
    //        // Test 4-way model
    //        {
    //            // Define link ends for 4-way model and constituent one-way models
    //            LinkEnds fourWayLinkEnds;
    //            fourWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation"  );
    //            fourWayLinkEnds[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );
    //            fourWayLinkEnds[ reflector2 ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2"  );
    //            fourWayLinkEnds[ reflector3 ] = std::make_pair< std::string, std::string >( "Moon" , ""  );
    //            fourWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );

    //            LinkEnds firstlinkLinkEnds;
    //            firstlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation" );
    //            firstlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );

    //            LinkEnds secondlinkLinkEnds;
    //            secondlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );
    //            secondlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2" );

    //            LinkEnds thirdlinkLinkEnds;
    //            thirdlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2" );
    //            thirdlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Moon" , "" );

    //            LinkEnds fourthlinkLinkEnds;
    //            fourthlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Moon" , "" );
    //            fourthlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation" );

    //            // Create light-time correction settings.
    //            std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    //            std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    //            lightTimeCorrectionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
    //                                                       lightTimePerturbingBodies ) );


    //            // Create observation settings for 4-way model and constituent one-way models
    //            std::shared_ptr< ObservationModelSettings > firstlinkObservableSettings = std::make_shared< ObservationModelSettings >
    //                    ( one_way_range, firstlinkLinkEnds, lightTimeCorrectionSettings );
    //            std::shared_ptr< ObservationModelSettings > secondlinkObservableSettings = std::make_shared< ObservationModelSettings >
    //                    ( one_way_range, secondlinkLinkEnds, lightTimeCorrectionSettings );
    //            std::shared_ptr< ObservationModelSettings > thirdlinkObservableSettings = std::make_shared< ObservationModelSettings >
    //                    ( one_way_range, thirdlinkLinkEnds, lightTimeCorrectionSettings );
    //            std::shared_ptr< ObservationModelSettings > fourthlinkObservableSettings = std::make_shared< ObservationModelSettings >
    //                    ( one_way_range, fourthlinkLinkEnds, lightTimeCorrectionSettings );

    //            std::vector< std::shared_ptr< ObservationModelSettings > > fourWayLinkSettings;
    //            fourWayLinkSettings.push_back( firstlinkObservableSettings );
    //            fourWayLinkSettings.push_back( secondlinkObservableSettings );
    //            fourWayLinkSettings.push_back( thirdlinkObservableSettings );
    //            fourWayLinkSettings.push_back( fourthlinkObservableSettings );
    //            std::shared_ptr< NWayRangeObservationSettings > fourWayObservableSettings =
    //                    std::make_shared< NWayRangeObservationSettings >(
    //                        fourWayLinkSettings, std::bind( &getRetransmissionDelays, std::placeholders::_1, 3 ) );

    //            // Create observation models
    //            std::shared_ptr< ObservationModel< 1, double, double > > firstlinkObservationModel =
    //                    ObservationModelCreator< 1, double, double >::createObservationModel(
    //                        firstlinkObservableSettings, bodies );
    //            std::shared_ptr< ObservationModel< 1, double, double > > secondlinkObservationModel =
    //                    ObservationModelCreator< 1, double, double >::createObservationModel(
    //                        secondlinkObservableSettings, bodies );
    //            std::shared_ptr< ObservationModel< 1, double, double > > thirdlinkObservationModel =
    //                    ObservationModelCreator< 1, double, double >::createObservationModel(
    //                        thirdlinkObservableSettings, bodies );
    //            std::shared_ptr< ObservationModel< 1, double, double > > fourthlinkObservationModel =
    //                    ObservationModelCreator< 1, double, double >::createObservationModel(
    //                        fourthlinkObservableSettings, bodies );
    //            std::shared_ptr< ObservationModel< 1, double, double > > fourWayObservationModel =
    //                    ObservationModelCreator< 1, double, double >::createObservationModel(
    //                        fourWayObservableSettings, bodies );

    //            // Define link ends time and state vectors for  2-way model and constituent one-way models
    //            std::vector< double > firstlinkLinkEndTimes;
    //            std::vector< Eigen::Vector6d > firstlinkLinkEndStates;

    //            std::vector< double > secondlinkLinkEndTimes;
    //            std::vector< Eigen::Vector6d > secondlinkLinkEndStates;

    //            std::vector< double > thirdlinkLinkEndTimes;
    //            std::vector< Eigen::Vector6d > thirdlinkLinkEndStates;

    //            std::vector< double > fourthlinkLinkEndTimes;
    //            std::vector< Eigen::Vector6d > fourthlinkLinkEndStates;

    //            std::vector< double > fourWayLinkEndTimes;
    //            std::vector< Eigen::Vector6d > fourWayLinkEndStates;

    //            // Define variables for range measurements
    //            Eigen::VectorXd fourWayRange;
    //            Eigen::VectorXd firstlinkRange;
    //            Eigen::VectorXd secondlinkRange;
    //            Eigen::VectorXd thirdlinkRange;
    //            Eigen::VectorXd fourthlinkRange;


    //            // Iterate over each 4-way link end as reference link end
    //            double observationTime = observationTimes.at( observationTimeNumber );
    //            double firstlinkObservationTime, secondlinkObservationTime, thirdlinkObservationTime, fourthlinkObservationTime;
    //            LinkEndType sublinkReferenceLinkEnd = transmitter;
    //            for( LinkEnds::const_iterator linkEndIterator = fourWayLinkEnds.begin( ); linkEndIterator != fourWayLinkEnds.end( );
    //                 linkEndIterator++ )
    //            {
    //                // Compute 4-way range
    //                fourWayRange = fourWayObservationModel->computeObservationsWithLinkEndData(
    //                            observationTime, linkEndIterator->first, fourWayLinkEndTimes, fourWayLinkEndStates );

    //                // Set observation times of constituent one-way ranges
    //                firstlinkObservationTime = fourWayLinkEndTimes.at( 0 );
    //                secondlinkObservationTime = fourWayLinkEndTimes.at( 2 );
    //                thirdlinkObservationTime = fourWayLinkEndTimes.at( 4 );
    //                fourthlinkObservationTime = fourWayLinkEndTimes.at( 6 );

    //                // Compute constituent one-way ranges
    //                firstlinkRange = firstlinkObservationModel->computeObservationsWithLinkEndData(
    //                            firstlinkObservationTime, sublinkReferenceLinkEnd, firstlinkLinkEndTimes, firstlinkLinkEndStates );
    //                secondlinkRange = secondlinkObservationModel->computeObservationsWithLinkEndData(
    //                            secondlinkObservationTime, sublinkReferenceLinkEnd, secondlinkLinkEndTimes, secondlinkLinkEndStates );
    //                thirdlinkRange = thirdlinkObservationModel->computeObservationsWithLinkEndData(
    //                            thirdlinkObservationTime, sublinkReferenceLinkEnd, thirdlinkLinkEndTimes, thirdlinkLinkEndStates );
    //                fourthlinkRange = fourthlinkObservationModel->computeObservationsWithLinkEndData(
    //                            fourthlinkObservationTime, sublinkReferenceLinkEnd, fourthlinkLinkEndTimes, fourthlinkLinkEndStates );

    //                // Check validity of retransmission delay
    //                std::vector< double > retransmissionDelays = getRetransmissionDelays(
    //                            observationTime, 3 );
    //                BOOST_CHECK_SMALL(
    //                            std::fabs( fourWayLinkEndTimes.at( 2 ) - fourWayLinkEndTimes.at( 1 ) - retransmissionDelays.at( 0 ) ),
    //                            observationTime  * std::numeric_limits< double >::epsilon( )  );
    //                BOOST_CHECK_SMALL(
    //                            std::fabs( fourWayLinkEndTimes.at( 4 ) - fourWayLinkEndTimes.at( 3 ) - retransmissionDelays.at( 1 ) ),
    //                            observationTime  * std::numeric_limits< double >::epsilon( )  );
    //                BOOST_CHECK_SMALL(
    //                            std::fabs( fourWayLinkEndTimes.at( 6 ) - fourWayLinkEndTimes.at( 5 ) - retransmissionDelays.at( 2 ) ),
    //                            observationTime  * std::numeric_limits< double >::epsilon( )  );

    //                // Check if link end times are consistent
    //                BOOST_CHECK_CLOSE_FRACTION( firstlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 0 ),
    //                                            std::numeric_limits< double >::epsilon( ) );
    //                BOOST_CHECK_CLOSE_FRACTION( firstlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 1 ),
    //                                            std::numeric_limits< double >::epsilon( ) );

    //                BOOST_CHECK_CLOSE_FRACTION( secondlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 2 ),
    //                                            std::numeric_limits< double >::epsilon( ) );
    //                BOOST_CHECK_CLOSE_FRACTION( secondlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 3 ),
    //                                            std::numeric_limits< double >::epsilon( ) );


    //                BOOST_CHECK_CLOSE_FRACTION( thirdlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 4 ),
    //                                            std::numeric_limits< double >::epsilon( ) );
    //                BOOST_CHECK_CLOSE_FRACTION( thirdlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 5 ),
    //                                            std::numeric_limits< double >::epsilon( ) );

    //                BOOST_CHECK_CLOSE_FRACTION( fourthlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 6 ),
    //                                            std::numeric_limits< double >::epsilon( ) );
    //                BOOST_CHECK_CLOSE_FRACTION( fourthlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 7 ),
    //                                            std::numeric_limits< double >::epsilon( ) );

    //                // Check if link end states are consistent
    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( firstlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 0 ),
    //                                                   std::numeric_limits< double >::epsilon( ) );
    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( firstlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 1 ),
    //                                                   std::numeric_limits< double >::epsilon( ) );

    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( secondlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 2 ),
    //                                                   std::numeric_limits< double >::epsilon( ) );
    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( secondlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 3 ),
    //                                                   std::numeric_limits< double >::epsilon( ) );

    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( thirdlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 4 ),
    //                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( thirdlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 5 ),
    //                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fourthlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 6 ),
    //                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
    //                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fourthlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 7 ),
    //                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

    //                // Check if range observations from 4-way and consituent one-ay are equal
    //                BOOST_CHECK_CLOSE_FRACTION(
    //                            ( firstlinkRange + secondlinkRange + thirdlinkRange + fourthlinkRange )( 0 ) +
    //                              ( retransmissionDelays.at( 0 ) + retransmissionDelays.at( 1 ) + retransmissionDelays.at( 2 ) ) *
    //                              physical_constants::SPEED_OF_LIGHT, fourWayRange( 0 ),
    //                            2.0 * std::numeric_limits< double >::epsilon( ) );
    //            }
    //        }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}

