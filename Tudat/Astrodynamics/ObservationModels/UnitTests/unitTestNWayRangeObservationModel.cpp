/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>
#include <boost/make_shared.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"

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
    std::map< std::string, boost::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > earthStationStation = std::pair< std::string, std::string >( "Earth", "EarthStation" );
    std::pair< std::string, std::string > earthStationStation2 = std::pair< std::string, std::string >( "Earth", "EarthStation2" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );
    createGroundStation( bodyMap.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "EarthStation", ( Eigen::Vector3d( ) << 1.0, 0.1, -1.4 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "EarthStation2", ( Eigen::Vector3d( ) << -30.0, 1.2, 2.1 ).finished( ),
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
            twoWayLinkEnds[ transmitter ] = std::make_pair( "Earth" , "EarthStation"  );
            twoWayLinkEnds[ reflector1 ] = std::make_pair( "Mars" , "MarsStation"  );
            twoWayLinkEnds[ receiver ] = std::make_pair( "Earth" , "EarthStation2"  );

            LinkEnds uplinkLinkEnds;
            uplinkLinkEnds[ transmitter ] = std::make_pair( "Earth" , "EarthStation" );
            uplinkLinkEnds[ receiver ] = std::make_pair( "Mars" , "MarsStation" );

            LinkEnds downlinkLinkEnds;
            downlinkLinkEnds[ transmitter ] = std::make_pair( "Mars" , "MarsStation" );
            downlinkLinkEnds[ receiver ] = std::make_pair( "Earth" , "EarthStation2" );

            // Create light-time correction settings.
            std::vector< std::string > lightTimePerturbingBodies = boost::assign::list_of( "Sun" );
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
            lightTimeCorrectionSettings.push_back( boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                       lightTimePerturbingBodies ) );


            // Create observation settings for 2-way model and constituent one-way models
            boost::shared_ptr< ObservationSettings > uplinkObservableSettings = boost::make_shared< ObservationSettings >
                    ( one_way_range, lightTimeCorrectionSettings );
            boost::shared_ptr< ObservationSettings > downlinkObservableSettings = boost::make_shared< ObservationSettings >
                    ( one_way_range, lightTimeCorrectionSettings );
            std::vector< boost::shared_ptr< ObservationSettings > > twoWayLinkSettings;
            twoWayLinkSettings.push_back( uplinkObservableSettings );
            twoWayLinkSettings.push_back( downlinkObservableSettings );
            boost::shared_ptr< NWayRangeObservationSettings > twoWayObservableSettings = boost::make_shared< NWayRangeObservationSettings >
                    ( twoWayLinkSettings, boost::bind( &getRetransmissionDelays, _1, 1 ) );

            // Create observation models
            boost::shared_ptr< ObservationModel< 1, double, double > > uplinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        uplinkLinkEnds, uplinkObservableSettings, bodyMap );
            boost::shared_ptr< ObservationModel< 1, double, double > > downlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        downlinkLinkEnds, downlinkObservableSettings, bodyMap );
            boost::shared_ptr< ObservationModel< 1, double, double > > twoWayObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        twoWayLinkEnds, twoWayObservableSettings, bodyMap );

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
                // Compute 2-way range
                twoWayRange = twoWayObservationModel->computeObservationsWithLinkEndData(
                            observationTime, linkEndIterator->first, twoWayLinkEndTimes, twoWayLinkEndStates );

                // Set observation times/reference link ends of constituent one-way ranges
                if( linkEndIterator->first == transmitter )
                {
                    uplinkObservationTime = observationTime;
                    downlinkObservationTime = twoWayLinkEndTimes.at( 2 );
                    uplinkReferenceLinkEnd = transmitter;
                    downlinkReferenceLinkEnd = transmitter;
                }
                else if( linkEndIterator->first == reflector1 )
                {
                    uplinkObservationTime = observationTime;
                    downlinkObservationTime = twoWayLinkEndTimes.at( 2 );
                    uplinkReferenceLinkEnd = receiver;
                    downlinkReferenceLinkEnd = transmitter;
                }
                else if( linkEndIterator->first == receiver )
                {
                    uplinkObservationTime = twoWayLinkEndTimes.at( 1 );
                    downlinkObservationTime = observationTime;
                    uplinkReferenceLinkEnd = receiver;
                    downlinkReferenceLinkEnd = receiver;
                }

                // Compute constituent one-way ranges
                uplinkRange = uplinkObservationModel->computeObservationsWithLinkEndData(
                            uplinkObservationTime, uplinkReferenceLinkEnd, uplinkLinkEndTimes, uplinkLinkEndStates );
                downlinkRange = downlinkObservationModel->computeObservationsWithLinkEndData(
                            downlinkObservationTime, downlinkReferenceLinkEnd, downlinkLinkEndTimes, downlinkLinkEndStates );

                // Check validity of retransmission delay
                std::vector< double > retransmissionDelays = getRetransmissionDelays(
                            observationTime, 1 );
                BOOST_CHECK_SMALL(
                            std::fabs( twoWayLinkEndTimes.at( 2 ) - twoWayLinkEndTimes.at( 1 ) - retransmissionDelays.at( 0 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );

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

                // Check if range observations from 2-way and consituent one-ay are equal
                BOOST_CHECK_CLOSE_FRACTION(
                            ( uplinkRange + downlinkRange )( 0 ) + retransmissionDelays.at( 0 ) *
                            physical_constants::SPEED_OF_LIGHT,
                            twoWayRange( 0 ), 2.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }

        // Test 4-way model
        {
            // Define link ends for 4-way model and constituent one-way models
            LinkEnds fourWayLinkEnds;
            fourWayLinkEnds[ transmitter ] = std::make_pair( "Earth" , "EarthStation"  );
            fourWayLinkEnds[ reflector1 ] = std::make_pair( "Mars" , "MarsStation"  );
            fourWayLinkEnds[ reflector2 ] = std::make_pair( "Earth" , "EarthStation2"  );
            fourWayLinkEnds[ reflector3 ] = std::make_pair( "Moon" , ""  );
            fourWayLinkEnds[ receiver ] = std::make_pair( "Mars" , "MarsStation"  );

            LinkEnds firstlinkLinkEnds;
            firstlinkLinkEnds[ transmitter ] = std::make_pair( "Earth" , "EarthStation" );
            firstlinkLinkEnds[ receiver ] = std::make_pair( "Mars" , "MarsStation" );

            LinkEnds secondlinkLinkEnds;
            secondlinkLinkEnds[ transmitter ] = std::make_pair( "Mars" , "MarsStation" );
            secondlinkLinkEnds[ receiver ] = std::make_pair( "Earth" , "EarthStation2" );

            LinkEnds thirdlinkLinkEnds;
            thirdlinkLinkEnds[ transmitter ] = std::make_pair( "Earth" , "EarthStation2" );
            thirdlinkLinkEnds[ receiver ] = std::make_pair( "Moon" , "" );

            LinkEnds fourthlinkLinkEnds;
            fourthlinkLinkEnds[ transmitter ] = std::make_pair( "Moon" , "" );
            fourthlinkLinkEnds[ receiver ] = std::make_pair( "Mars" , "MarsStation" );

            // Create light-time correction settings.
            std::vector< std::string > lightTimePerturbingBodies = boost::assign::list_of( "Sun" );
            std::vector< boost::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
            lightTimeCorrectionSettings.push_back( boost::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >(
                                                       lightTimePerturbingBodies ) );


            // Create observation settings for 4-way model and constituent one-way models
            boost::shared_ptr< ObservationSettings > firstlinkObservableSettings = boost::make_shared< ObservationSettings >
                    ( one_way_range, lightTimeCorrectionSettings );
            boost::shared_ptr< ObservationSettings > secondlinkObservableSettings = boost::make_shared< ObservationSettings >
                    ( one_way_range, lightTimeCorrectionSettings );
            boost::shared_ptr< ObservationSettings > thirdlinkObservableSettings = boost::make_shared< ObservationSettings >
                    ( one_way_range, lightTimeCorrectionSettings );
            boost::shared_ptr< ObservationSettings > fourthlinkObservableSettings = boost::make_shared< ObservationSettings >
                    ( one_way_range, lightTimeCorrectionSettings );

            std::vector< boost::shared_ptr< ObservationSettings > > fourWayLinkSettings;
            fourWayLinkSettings.push_back( firstlinkObservableSettings );
            fourWayLinkSettings.push_back( secondlinkObservableSettings );
            fourWayLinkSettings.push_back( thirdlinkObservableSettings );
            fourWayLinkSettings.push_back( fourthlinkObservableSettings );
            boost::shared_ptr< NWayRangeObservationSettings > fourWayObservableSettings =
                    boost::make_shared< NWayRangeObservationSettings >(
                        fourWayLinkSettings, boost::bind( &getRetransmissionDelays, _1, 3 ) );

            // Create observation models
            boost::shared_ptr< ObservationModel< 1, double, double > > firstlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        firstlinkLinkEnds, firstlinkObservableSettings, bodyMap );
            boost::shared_ptr< ObservationModel< 1, double, double > > secondlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        secondlinkLinkEnds, secondlinkObservableSettings, bodyMap );
            boost::shared_ptr< ObservationModel< 1, double, double > > thirdlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        thirdlinkLinkEnds, thirdlinkObservableSettings, bodyMap );
            boost::shared_ptr< ObservationModel< 1, double, double > > fourthlinkObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        fourthlinkLinkEnds, fourthlinkObservableSettings, bodyMap );
            boost::shared_ptr< ObservationModel< 1, double, double > > fourWayObservationModel =
                    ObservationModelCreator< 1, double, double >::createObservationModel(
                        fourWayLinkEnds, fourWayObservableSettings, bodyMap );

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
                // Compute 4-way range
                fourWayRange = fourWayObservationModel->computeObservationsWithLinkEndData(
                            observationTime, linkEndIterator->first, fourWayLinkEndTimes, fourWayLinkEndStates );

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
                std::vector< double > retransmissionDelays = getRetransmissionDelays(
                            observationTime, 3 );
                BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 2 ) - fourWayLinkEndTimes.at( 1 ) - retransmissionDelays.at( 0 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );
                BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 4 ) - fourWayLinkEndTimes.at( 3 ) - retransmissionDelays.at( 1 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );
                BOOST_CHECK_SMALL(
                            std::fabs( fourWayLinkEndTimes.at( 6 ) - fourWayLinkEndTimes.at( 5 ) - retransmissionDelays.at( 2 ) ),
                            observationTime  * std::numeric_limits< double >::epsilon( )  );

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
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( thirdlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 5 ),
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fourthlinkLinkEndStates.at( 0 ), fourWayLinkEndStates.at( 6 ),
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( fourthlinkLinkEndStates.at( 1 ), fourWayLinkEndStates.at( 7 ),
                                                   ( 10.0 * std::numeric_limits< double >::epsilon( ) ) );

                // Check if range observations from 4-way and consituent one-ay are equal
                BOOST_CHECK_CLOSE_FRACTION(
                            ( firstlinkRange + secondlinkRange + thirdlinkRange + fourthlinkRange )( 0 ) +
                              ( retransmissionDelays.at( 0 ) + retransmissionDelays.at( 1 ) + retransmissionDelays.at( 2 ) ) *
                              physical_constants::SPEED_OF_LIGHT, fourWayRange( 0 ),
                            2.0 * std::numeric_limits< double >::epsilon( ) );
            }
        }
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}

