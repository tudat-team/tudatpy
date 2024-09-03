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


BOOST_AUTO_TEST_SUITE( test_n_way_range_model )

std::vector< double > getRetransmissionDelays( const double evaluationTime, const int numberOfRetransmitters )
{
    std::vector< double > retransmissionDelays;
    if( evaluationTime < 1.5E5 )
    {
        // Transmission delay
        retransmissionDelays.push_back( 0.0 );
        // Retransmission delay
        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( 0.0 );
        }
        // Reception delay
        retransmissionDelays.push_back( 0.0 );
    }
    else
    {
        // Transmission delay
        retransmissionDelays.push_back( 1.0e-5 );
        // Retransmission delay
        for( int i = 0; i < numberOfRetransmitters; i++ )
        {
            retransmissionDelays.push_back( evaluationTime * 5.0E-12 * static_cast< double >( i + 1 ) );
        }
        // Reception delay
        retransmissionDelays.push_back( 2.0e-5 );
    }
    return retransmissionDelays;
}

BOOST_AUTO_TEST_CASE( testNWayRangeModel )
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
    observationTimes.push_back( 1.0E5 );
    observationTimes.push_back( 2.0E5 );
    observationTimes.push_back( 3.0E5 );

    // Define link ends for observations.
    for( unsigned int observationTimeNumber = 0; observationTimeNumber < observationTimes.size( ); observationTimeNumber++ )
    {
        // Test 2-way model
        {
            // Define link ends for 2-way model and constituent one-way models
            LinkEnds twoWayLinkEnds;
            twoWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation"  );
            twoWayLinkEnds[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );
            twoWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2"  );

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

            // Create light time convergence criteria
            std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( true );

            // Create observation settings for 2-way model and constituent one-way models
            std::shared_ptr< ObservationModelSettings > uplinkObservableSettings = std::make_shared< ObservationModelSettings >
                    ( one_way_range, uplinkLinkEnds, lightTimeCorrectionSettings );
            std::shared_ptr< ObservationModelSettings > downlinkObservableSettings = std::make_shared< ObservationModelSettings >
                    ( one_way_range, downlinkLinkEnds, lightTimeCorrectionSettings );
            std::vector< std::shared_ptr< ObservationModelSettings > > twoWayLinkSettings;
            twoWayLinkSettings.push_back( uplinkObservableSettings );
            twoWayLinkSettings.push_back( downlinkObservableSettings );
            std::shared_ptr< NWayRangeObservationSettings > twoWayObservableSettings = std::make_shared< NWayRangeObservationSettings >
                    ( twoWayLinkSettings, nullptr, lightTimeConvergenceCriteria );

            // Create observation models
            std::shared_ptr< ObservationModel< 1, double, double > > uplinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        uplinkObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > downlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        downlinkObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > twoWayObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        twoWayObservableSettings, bodies );

            // Define link ends time and state vectors for  2-way model and constituent one-way models
            std::vector< double > uplinkLinkEndTimes;
            std::vector< Eigen::Vector6d > uplinkLinkEndStates;

            std::vector< double > downlinkLinkEndTimes;
            std::vector< Eigen::Vector6d > downlinkLinkEndStates;

            std::vector< double > twoWayLinkEndTimes;
            std::vector< Eigen::Vector6d > twoWayLinkEndStates;

            // Define variables for range measurements
            Eigen::VectorXd twoWayRange;
            Eigen::VectorXd uplinkRange;
            Eigen::VectorXd downlinkRange;

            // Iterate over each 2-way link end as reference link end
            double observationTime = observationTimes.at( observationTimeNumber );
            double uplinkObservationTime = TUDAT_NAN, downlinkObservationTime = TUDAT_NAN;
            LinkEndType uplinkReferenceLinkEnd = unidentified_link_end , downlinkReferenceLinkEnd = unidentified_link_end;
            for( LinkEnds::const_iterator linkEndIterator = twoWayLinkEnds.begin( ); linkEndIterator != twoWayLinkEnds.end( );
                 linkEndIterator++ )
            {
                std::vector< double > retransmissionDelays = getRetransmissionDelays( observationTime, 1 );
                // Only check reference link ends other than transmitter/ receiver if the associated delay is not zero.
                // Light time calculation is not implemented for case with non-zero delay
                if ( linkEndIterator->first != transmitter && linkEndIterator->first != receiver )
                {
                    if ( retransmissionDelays.at( 1 ) != 0 )
                    {
                        continue;
                    }
                }

                // Compute 2-way range
                twoWayRange = twoWayObservationModel->computeObservationsWithLinkEndData(
                            observationTime, linkEndIterator->first, twoWayLinkEndTimes, twoWayLinkEndStates,
                            getNWayRangeAncilliarySettings( retransmissionDelays ) );

                // Set observation times/reference link ends of constituent one-way ranges
                if( linkEndIterator->first == transmitter )
                {
                    uplinkObservationTime = twoWayLinkEndTimes.at( 0 );
                    downlinkObservationTime = twoWayLinkEndTimes.at( 2 );
                    uplinkReferenceLinkEnd = transmitter;
                    downlinkReferenceLinkEnd = transmitter;
                }
                else if( linkEndIterator->first == reflector1 )
                {
                    uplinkObservationTime = twoWayLinkEndTimes.at( 1 );
                    downlinkObservationTime = twoWayLinkEndTimes.at( 2 );
                    uplinkReferenceLinkEnd = receiver;
                    downlinkReferenceLinkEnd = transmitter;
                }
                else if( linkEndIterator->first == receiver )
                {
                    uplinkObservationTime = twoWayLinkEndTimes.at( 1 );
                    downlinkObservationTime = twoWayLinkEndTimes.at( 3 );
                    uplinkReferenceLinkEnd = receiver;
                    downlinkReferenceLinkEnd = receiver;
                }

                // Compute constituent one-way ranges
                uplinkRange = uplinkObservationModel->computeObservationsWithLinkEndData(
                            uplinkObservationTime, uplinkReferenceLinkEnd, uplinkLinkEndTimes, uplinkLinkEndStates );
                downlinkRange = downlinkObservationModel->computeObservationsWithLinkEndData(
                            downlinkObservationTime, downlinkReferenceLinkEnd, downlinkLinkEndTimes, downlinkLinkEndStates );

                // Check validity of retransmission delay
                BOOST_CHECK_SMALL(
                            std::fabs( twoWayLinkEndTimes.at( 2 ) - twoWayLinkEndTimes.at( 1 ) - retransmissionDelays.at( 1 ) ),
                            observationTime * std::numeric_limits< double >::epsilon( ) );
                // Check validity of transmission delay
                if ( linkEndIterator->first == transmitter )
                {
                    BOOST_CHECK_SMALL(
                            std::fabs( twoWayLinkEndTimes.at( 0 ) - observationTime - retransmissionDelays.at( 0 ) ),
                            observationTime * std::numeric_limits< double >::epsilon( ) );
                }
                // Check validity of reception delay
                else if ( linkEndIterator->first == receiver )
                {
                    BOOST_CHECK_SMALL(
                            std::fabs( observationTime - twoWayLinkEndTimes.at( 3 ) - retransmissionDelays.at( 2 ) ),
                            observationTime * std::numeric_limits< double >::epsilon( ) );
                }

                // Check number of multi-leg iterations
                int numIter = std::dynamic_pointer_cast< NWayRangeObservationModel< double, double > >(
                            twoWayObservationModel )->getMultiLegLightTimeCalculator( )->getNumberOfMultiLegIterations( );
                bool iterateMultipleLegs = std::dynamic_pointer_cast< NWayRangeObservationModel< double, double > >(
                            twoWayObservationModel )->getMultiLegLightTimeCalculator( )->getIterateMultiLegLightTime( );
                BOOST_CHECK_EQUAL( numIter, 0 );
                BOOST_CHECK_EQUAL( iterateMultipleLegs, false );

                // Check if link end times are consistent
                BOOST_CHECK_CLOSE_FRACTION( uplinkLinkEndTimes.at( 0 ), twoWayLinkEndTimes.at( 0 ),
                                            std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( uplinkLinkEndTimes.at( 1 ), twoWayLinkEndTimes.at( 1 ),
                                            std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_CLOSE_FRACTION( downlinkLinkEndTimes.at( 0 ), twoWayLinkEndTimes.at( 2 ),
                                            std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( downlinkLinkEndTimes.at( 1 ), twoWayLinkEndTimes.at( 3 ),
                                            std::numeric_limits< double >::epsilon( ) );

                // Check if link end states are consistent
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( uplinkLinkEndStates.at( 0 ), twoWayLinkEndStates.at( 0 ),
                                                   std::numeric_limits< double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( uplinkLinkEndStates.at( 1 ), twoWayLinkEndStates.at( 1 ),
                                                   std::numeric_limits< double >::epsilon( ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( downlinkLinkEndStates.at( 0 ), twoWayLinkEndStates.at( 2 ),
                                                   std::numeric_limits< double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( downlinkLinkEndStates.at( 1 ), twoWayLinkEndStates.at( 3 ),
                                                   std::numeric_limits< double >::epsilon( ) );

                // Check if range observations from 2-way and constituent one-way are equal
                double delaySum = retransmissionDelays.at( 0 ) + retransmissionDelays.at( 1 ) + retransmissionDelays.at( 2 );
                BOOST_CHECK_CLOSE_FRACTION(
                            ( uplinkRange + downlinkRange )( 0 ) + delaySum *
                            physical_constants::SPEED_OF_LIGHT,
                            twoWayRange( 0 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
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

            // Create light time convergence criteria
            bool iterateMultipleLegs = observationTimeNumber % 2;
            std::shared_ptr< LightTimeConvergenceCriteria > lightTimeConvergenceCriteria =
                std::make_shared< LightTimeConvergenceCriteria >( true, iterateMultipleLegs );


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
            std::shared_ptr< NWayRangeObservationSettings > fourWayObservableSettings =
                    std::make_shared< NWayRangeObservationSettings >(
                        fourWayLinkSettings, nullptr, lightTimeConvergenceCriteria );

            // Create observation models
            std::shared_ptr< ObservationModel< 1, double, double > > firstlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        firstlinkObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > secondlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        secondlinkObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > thirdlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        thirdlinkObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > fourthlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        fourthlinkObservableSettings, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > fourWayObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        fourWayObservableSettings, bodies );

            // Define link ends time and state vectors for  2-way model and constituent one-way models
            std::vector< double > firstlinkLinkEndTimes;
            std::vector< Eigen::Vector6d > firstlinkLinkEndStates;

            std::vector< double > secondlinkLinkEndTimes;
            std::vector< Eigen::Vector6d > secondlinkLinkEndStates;

            std::vector< double > thirdlinkLinkEndTimes;
            std::vector< Eigen::Vector6d > thirdlinkLinkEndStates;

            std::vector< double > fourthlinkLinkEndTimes;
            std::vector< Eigen::Vector6d > fourthlinkLinkEndStates;

            std::vector< double > fourWayLinkEndTimes;
            std::vector< Eigen::Vector6d > fourWayLinkEndStates;

            // Define variables for range measurements
            Eigen::VectorXd fourWayRange;
            Eigen::VectorXd firstlinkRange;
            Eigen::VectorXd secondlinkRange;
            Eigen::VectorXd thirdlinkRange;
            Eigen::VectorXd fourthlinkRange;


            // Iterate over each 4-way link end as reference link end
            double observationTime = observationTimes.at( observationTimeNumber );
            double firstlinkObservationTime, secondlinkObservationTime, thirdlinkObservationTime, fourthlinkObservationTime;
            LinkEndType sublinkReferenceLinkEnd = transmitter;
            for( LinkEnds::const_iterator linkEndIterator = fourWayLinkEnds.begin( ); linkEndIterator != fourWayLinkEnds.end( );
                 linkEndIterator++ )
            {
                std::vector< double > retransmissionDelays = getRetransmissionDelays( observationTime, 3 );
                // Only check reference link ends other than transmitter/ receiver if the associated delay is not zero.
                // Light time calculation is not implemented for case with non-zero delay
                if ( linkEndIterator->first != transmitter && linkEndIterator->first != receiver )
                {
                    if ( retransmissionDelays.at( 1 ) != 0 || retransmissionDelays.at( 2 ) != 0 || retransmissionDelays.at( 3 ) != 0 )
                    {
                        continue;
                    }
                }

                // Compute 4-way range
                fourWayRange = fourWayObservationModel->computeObservationsWithLinkEndData(
                            observationTime, linkEndIterator->first, fourWayLinkEndTimes, fourWayLinkEndStates,
                            getNWayRangeAncilliarySettings( retransmissionDelays ) );

                // Set observation times of constituent one-way ranges
                firstlinkObservationTime = fourWayLinkEndTimes.at( 0 );
                secondlinkObservationTime = fourWayLinkEndTimes.at( 2 );
                thirdlinkObservationTime = fourWayLinkEndTimes.at( 4 );
                fourthlinkObservationTime = fourWayLinkEndTimes.at( 6 );

                // Compute constituent one-way ranges
                firstlinkRange = firstlinkObservationModel->computeObservationsWithLinkEndData(
                            firstlinkObservationTime, sublinkReferenceLinkEnd, firstlinkLinkEndTimes, firstlinkLinkEndStates );
                secondlinkRange = secondlinkObservationModel->computeObservationsWithLinkEndData(
                            secondlinkObservationTime, sublinkReferenceLinkEnd, secondlinkLinkEndTimes, secondlinkLinkEndStates );
                thirdlinkRange = thirdlinkObservationModel->computeObservationsWithLinkEndData(
                            thirdlinkObservationTime, sublinkReferenceLinkEnd, thirdlinkLinkEndTimes, thirdlinkLinkEndStates );
                fourthlinkRange = fourthlinkObservationModel->computeObservationsWithLinkEndData(
                            fourthlinkObservationTime, sublinkReferenceLinkEnd, fourthlinkLinkEndTimes, fourthlinkLinkEndStates );

                // Check validity of retransmission delay
                BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 2 ) - fourWayLinkEndTimes.at( 1 ) - retransmissionDelays.at( 1 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );
                BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 4 ) - fourWayLinkEndTimes.at( 3 ) - retransmissionDelays.at( 2 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );
                BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 6 ) - fourWayLinkEndTimes.at( 5 ) - retransmissionDelays.at( 3 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );
                // Check validity of transmission delay
                if ( linkEndIterator->first == transmitter )
                {
                    BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 0 ) - observationTime - retransmissionDelays.at( 0 ) ),
                            observationTime * std::numeric_limits< double >::epsilon( ) );
                }
                // Check validity of reception delay
                else if ( linkEndIterator->first == receiver )
                {
                    BOOST_CHECK_SMALL(
                            std::fabs( observationTime - fourWayLinkEndTimes.at( 7 ) - retransmissionDelays.at( 4 ) ),
                            observationTime * std::numeric_limits< double >::epsilon( ) );
                }

                 // Check number of multi-leg iterations
                int numIter = std::dynamic_pointer_cast< NWayRangeObservationModel< double, double > >(
                            fourWayObservationModel )->getMultiLegLightTimeCalculator( )->getNumberOfMultiLegIterations( );
                bool iterateMultipleLegs = std::dynamic_pointer_cast< NWayRangeObservationModel< double, double > >(
                            fourWayObservationModel )->getMultiLegLightTimeCalculator( )->getIterateMultiLegLightTime( );
                BOOST_CHECK_EQUAL( numIter, 0 );
                BOOST_CHECK_EQUAL( iterateMultipleLegs, false );

                // Check if link end times are consistent
                BOOST_CHECK_CLOSE_FRACTION( firstlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 0 ),
                                            std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( firstlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 1 ),
                                            std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_CLOSE_FRACTION( secondlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 2 ),
                                            std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( secondlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 3 ),
                                            std::numeric_limits< double >::epsilon( ) );


                BOOST_CHECK_CLOSE_FRACTION( thirdlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 4 ),
                                            std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( thirdlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 5 ),
                                            std::numeric_limits< double >::epsilon( ) );

                BOOST_CHECK_CLOSE_FRACTION( fourthlinkLinkEndTimes.at( 0 ), fourWayLinkEndTimes.at( 6 ),
                                            std::numeric_limits< double >::epsilon( ) );
                BOOST_CHECK_CLOSE_FRACTION( fourthlinkLinkEndTimes.at( 1 ), fourWayLinkEndTimes.at( 7 ),
                                            std::numeric_limits< double >::epsilon( ) );

                // Check if link end states are consistent
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( firstlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 0 ),
                                                   std::numeric_limits< double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( firstlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 1 ),
                                                   std::numeric_limits< double >::epsilon( ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( secondlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 2 ),
                                                   std::numeric_limits< double >::epsilon( ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( secondlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 3 ),
                                                   std::numeric_limits< double >::epsilon( ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( thirdlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 4 ),
                                                   ( 25.0 * std::numeric_limits< double >::epsilon( ) ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( thirdlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 5 ),
                                                   ( 25.0 * std::numeric_limits< double >::epsilon( ) ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fourthlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 6 ),
                                                   ( 25.0 * std::numeric_limits< double >::epsilon( ) ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fourthlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 7 ),
                                                   ( 25.0 * std::numeric_limits< double >::epsilon( ) ) );

                // Check if range observations from 4-way and constituent one-way are equal
                double delaySum = retransmissionDelays.at( 0 ) + retransmissionDelays.at( 1 ) + retransmissionDelays.at( 2 ) +
                        retransmissionDelays.at( 3 ) + retransmissionDelays.at( 4 );
                BOOST_CHECK_CLOSE_FRACTION(
                            ( firstlinkRange + secondlinkRange + thirdlinkRange + fourthlinkRange )( 0 ) +
                              delaySum * physical_constants::SPEED_OF_LIGHT, fourWayRange( 0 ),
                            2.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }
    }
}


BOOST_AUTO_TEST_CASE( testTwoWayRangeModelTimeScaleBias )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
        getDefaultBodySettings(
            bodiesToCreate );

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
    observationTimes.push_back( 1.0E8 );
    observationTimes.push_back( 2.0E8 );
    observationTimes.push_back( 3.0E8 );

    // Define link ends for observations.
    for( unsigned int observationTimeNumber = 0; observationTimeNumber < observationTimes.size( ); observationTimeNumber++ )
    {
        // Define link ends for 2-way model and constituent one-way models
        LinkEnds twoWayLinkEnds;
        twoWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation"  );
        twoWayLinkEnds[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars" , "MarsStation"  );
        twoWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , "EarthStation2"  );

        // Create observation models
        std::shared_ptr< NWayRangeObservationSettings > twoWayObservableSettings = std::make_shared< NWayRangeObservationSettings >
            ( twoWayLinkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( )  );
        std::shared_ptr< ObservationModel< 1, double, double > > twoWayObservationModel =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                twoWayObservableSettings, bodies );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;
        double unbiasedTwoWayRange = twoWayObservationModel->computeObservationsWithLinkEndData(
            observationTimes.at( observationTimeNumber ), receiver,
            linkEndTimes, linkEndStates )( 0 );

        std::shared_ptr< ObservationBiasSettings > biasSettings = twoWayTimeScaleRangeBias( );
        std::shared_ptr< NWayRangeObservationSettings > twoWayObservableSettingsWithBias = std::make_shared< NWayRangeObservationSettings >
            ( twoWayLinkEnds, std::vector< std::shared_ptr< LightTimeCorrectionSettings > >( ), biasSettings );
        std::shared_ptr< ObservationModel< 1, double, double > > twoWayObservationModelWithBias =
            ObservationModelCreator< 1, double, double >::createObservationModel(
                twoWayObservableSettingsWithBias, bodies );

        double biasedTwoWayRange = twoWayObservationModelWithBias->computeObservationsWithLinkEndData(
            observationTimes.at( observationTimeNumber ), receiver,
            linkEndTimes, linkEndStates )( 0 );

        std::cout<<biasedTwoWayRange - unbiasedTwoWayRange<<std::endl;
        std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > timeScaleConverter = earth_orientation::defaultTimeConverter;
        double transmissionTimeDifference = timeScaleConverter->getCurrentTime< Time >( basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, linkEndTimes.at( 0 ),
                                                    bodies.at( "Earth" )->getGroundStationMap( ).at( "EarthStation" )->getNominalStationState( )->getNominalCartesianPosition( ) ) - Time( linkEndTimes.at( 0 ) );
        double receptionTimeDifference = timeScaleConverter->getCurrentTime< Time >( basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, linkEndTimes.at( 3 ),
                                                    bodies.at( "Earth" )->getGroundStationMap( ).at( "EarthStation2" )->getNominalStationState( )->getNominalCartesianPosition( ) ) - Time( linkEndTimes.at( 3 ) );
        BOOST_CHECK_SMALL( ( receptionTimeDifference -transmissionTimeDifference ) * physical_constants::SPEED_OF_LIGHT - (biasedTwoWayRange - unbiasedTwoWayRange), 1.0E-4 );
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}

