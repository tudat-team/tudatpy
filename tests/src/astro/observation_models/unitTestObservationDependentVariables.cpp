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
//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MAIN

#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/utilities.h"

#include "tudat/simulation/estimation.h"
//
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
//
//BOOST_AUTO_TEST_SUITE( test_observation_dependent_variables )
//
////! Test whether observation noise is correctly added when simulating noisy observations
//BOOST_AUTO_TEST_CASE( testObservationDependentVariables )
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
    groundStationNames.push_back( "Station3" );//            relative_angular_position = 9,


    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 1.0 ).finished( ), geodetic_position );

    // Define parameters.
    std::vector< LinkEnds > stationTransmitterOneWayLinkEnds;
    std::vector< LinkEnds > stationReceiverOneWayLinkEnds;

    std::vector< LinkEnds > stationReceiverTwoWayLinkEnds;
    std::vector< LinkEnds > stationRetransmitterTwoWayLinkEnds;
    std::vector< LinkEnds > stationReceiverThreeWayLinkEnds;
    std::vector< LinkEnds > stationTransmitterThreeWayLinkEnds;

    // Define link ends to/from ground stations to Moon
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

    for( unsigned int observableTestCase = 0; observableTestCase < 1; observableTestCase++ )
    {
        int geometryType = -1;
        ObservableType currentObservableType;
        if( observableTestCase == 0 )
        {
            currentObservableType = one_way_range;
            geometryType = 0;
        }
        else if( observableTestCase == 1 )
        {
            currentObservableType = angular_position;
            geometryType = 0;
        }
        else if( observableTestCase == 2 )
        {
            currentObservableType = one_way_doppler;
            geometryType = 0;
        }
        else if( observableTestCase == 3 )
        {
            currentObservableType = n_way_range;
            geometryType = 1;
        }
        else if( observableTestCase == 4 )
        {
            currentObservableType = two_way_doppler;
            geometryType = 1;
        }
        
//            n_way_range = 5,
//            two_way_doppler = 6,

//            relative_angular_position = 9,

//            dsn_one_way_averaged_doppler = 12,
//            one_way_differenced_range = 4,

//            dsn_n_way_averaged_doppler = 13
//            n_way_differenced_range = 10,

//            position_observable = 2,
//            velocity_observable = 8,
//            relative_position_observable = 11,
//            euler_angle_313_observable = 7,

        int numberOfLinkEndCases = -1;
        if( geometryType == 0 )
        {
            numberOfLinkEndCases = 1;
        }
        else if( geometryType == 1 )
        {
            numberOfLinkEndCases = 4;
        }

        for( int currentLinkEndCase = 0; currentLinkEndCase < numberOfLinkEndCases; currentLinkEndCase++ )
        {

            LinkEnds currentLinkEnds;
            switch( geometryType )
            {
            case 0:
            {
                switch ( currentLinkEndCase )
                {
                case 0:
                    currentLinkEnds = stationReceiverOneWayLinkEnds.at( 0 );
                    break;
                case 1:
                    currentLinkEnds = stationTransmitterOneWayLinkEnds.at( 0 );
                    break;
                default:
                    throw std::runtime_error( "Error in observation dependent variable unit test A " );
                }
                break;
            }
            case 1:
            {
                switch ( currentLinkEndCase )
                {
                case 0:
                    currentLinkEnds = stationReceiverTwoWayLinkEnds.at( 0 );
                    break;
                case 1:
                    currentLinkEnds = stationReceiverThreeWayLinkEnds.at( 0 );
                    break;
                case 2:
                    currentLinkEnds = stationTransmitterThreeWayLinkEnds.at( 0 );
                    break;
                case 3:
                    currentLinkEnds = stationRetransmitterTwoWayLinkEnds.at( 0 );
                    break;
                default:
                    throw std::runtime_error( "Error in observation dependent variable unit test A " );
                }
                break;
            }
            }
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
                    observationSettingsList.push_back( std::make_shared<ObservationModelSettings>(
                        currentObservable, currentLinkEndsList.at( i )));
                }
            }

            // Create observation simulators
            std::vector<std::shared_ptr<ObservationSimulatorBase<double, double> > > observationSimulators =
                createObservationSimulators( observationSettingsList, bodies );

            // Define osbervation times.
            std::vector<double> baseTimeList;
            double observationTimeStart = initialEphemerisTime + 1000.0;
            double observationInterval = 10.0;
            for ( unsigned int i = 0; i < 14; i++ )
            {
                for ( unsigned int j = 0; j < 4320; j++ )
                {
                    baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                            static_cast< double >( j ) * observationInterval );
                }
            }

            // Define observation simulation settings (observation type, link end, times and reference link end)
            std::vector<std::shared_ptr<ObservationSimulationSettings<double> > > measurementSimulationInput;
            for (
                std::map<ObservableType, std::vector<LinkEnds> >::iterator
                    linkEndIterator = linkEndsPerObservable.begin( );
                linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
            {
                ObservableType currentObservable = linkEndIterator->first;
                std::vector<LinkEnds> currentLinkEndsList = linkEndIterator->second;
                for ( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
                {
                    measurementSimulationInput.push_back(
                        std::make_shared<TabulatedObservationSimulationSettings<> >(
                            currentObservable, currentLinkEndsList.at( i ), baseTimeList, receiver ));
                }
            }

            std::vector<std::shared_ptr<observation_models::ObservationViabilitySettings> > viabilitySettingsList;
            viabilitySettingsList.push_back( elevationAngleViabilitySettings(
                std::make_pair( "Earth", "Station1" ), 25.0 * mathematical_constants::PI / 180.0 ));
            viabilitySettingsList.push_back( elevationAngleViabilitySettings(
                std::make_pair( "Earth", "Station2" ), 25.0 * mathematical_constants::PI / 180.0 ));

            addViabilityToObservationSimulationSettings(
                measurementSimulationInput, viabilitySettingsList );

            // Define settings for dependent variables
            std::vector<std::shared_ptr<ObservationDependentVariableSettings> > dependentVariableList;

            std::shared_ptr<ObservationDependentVariableSettings> elevationAngleSettings1 =
                std::make_shared<StationAngleObservationDependentVariableSettings>(
                    station_elevation_angle, LinkEndId( std::make_pair( "Earth", "Station1" )));
            std::shared_ptr<ObservationDependentVariableSettings> azimuthAngleSettings1 =
                std::make_shared<StationAngleObservationDependentVariableSettings>(
                    station_azimuth_angle, LinkEndId( std::make_pair( "Earth", "Station1" )));
            std::shared_ptr<ObservationDependentVariableSettings> targetRangeSettings1 =
                std::make_shared<InterlinkObservationDependentVariableSettings>(
                    target_range, transmitter, receiver );
            std::shared_ptr<ObservationDependentVariableSettings> targetInverseRangeSettings1 =
                std::make_shared<InterlinkObservationDependentVariableSettings>(
                    target_range, receiver, transmitter );

            std::shared_ptr<ObservationDependentVariableSettings> elevationAngleSettings2 =
                std::make_shared<StationAngleObservationDependentVariableSettings>(
                    station_elevation_angle, LinkEndId( std::make_pair( "Earth", "Station2" )));

            dependentVariableList.push_back( elevationAngleSettings1 );
            dependentVariableList.push_back( azimuthAngleSettings1 );
            dependentVariableList.push_back( elevationAngleSettings2 );

            addDependentVariablesToObservationSimulationSettings(
                measurementSimulationInput, dependentVariableList, bodies );

            addDependentVariablesToObservationSimulationSettings(
                measurementSimulationInput, { targetRangeSettings1 }, bodies, currentObservableType, currentLinkEnds );
            addDependentVariablesToObservationSimulationSettings(
                measurementSimulationInput, { targetInverseRangeSettings1 }, bodies, currentObservableType, currentLinkEnds );


            // Simulate noise-free observations
            std::shared_ptr<ObservationCollection<> > idealObservationsAndTimes = simulateObservations<double, double>(
                measurementSimulationInput, observationSimulators, bodies );

            if ( observableTestCase == 0 )
            {

                std::map<double, Eigen::VectorXd> elevationAngles1 = getDependentVariableResultList(
                    idealObservationsAndTimes, elevationAngleSettings1, currentObservableType );
                std::map<double, Eigen::VectorXd> elevationAngles2 = getDependentVariableResultList(
                    idealObservationsAndTimes, elevationAngleSettings2, currentObservableType );
                std::map<double, Eigen::VectorXd> azimuthAngles1 = getDependentVariableResultList(
                    idealObservationsAndTimes, azimuthAngleSettings1, currentObservableType );
                std::map<double, Eigen::VectorXd> targetRanges1 = getDependentVariableResultList(
                    idealObservationsAndTimes, targetRangeSettings1, currentObservableType );
                std::map<double, Eigen::VectorXd> targetInverseRanges1 = getDependentVariableResultList(
                    idealObservationsAndTimes, targetInverseRangeSettings1, currentObservableType );

                referenceElevationAngles = elevationAngles1;
                referenceAzimuthAngles = azimuthAngles1;
                referenceTargetRanges = targetRanges1;

                std::map<observation_models::LinkEnds, int>
                    linkEndIdentifiers = idealObservationsAndTimes->getLinkEndIdentifierMap( );
                std::vector<int> linkEndIds = idealObservationsAndTimes->getConcatenatedLinkEndIds( );

                int numberOfLinkEnds1Observations = utilities::countNumberOfOccurencesInVector<int>(
                    linkEndIds, linkEndIdentifiers.at( stationReceiverOneWayLinkEnds[ 0 ] ));
                int numberOfLinkEnds2Observations = utilities::countNumberOfOccurencesInVector<int>(
                    linkEndIds, linkEndIdentifiers.at( stationReceiverOneWayLinkEnds[ 1 ] ));

                BOOST_CHECK_EQUAL( elevationAngles1.size( ), numberOfLinkEnds1Observations );
                BOOST_CHECK_EQUAL( azimuthAngles1.size( ), numberOfLinkEnds1Observations );
                BOOST_CHECK_EQUAL( targetRanges1.size( ), numberOfLinkEnds1Observations );
                BOOST_CHECK_EQUAL( targetInverseRanges1.size( ), numberOfLinkEnds1Observations );
                BOOST_CHECK_EQUAL( elevationAngles2.size( ), numberOfLinkEnds2Observations );

                std::shared_ptr<ground_stations::PointingAnglesCalculator> pointingAnglesCalculator1 =
                    bodies.at( "Earth" )->getGroundStation( "Station1" )->getPointingAnglesCalculator( );
                std::shared_ptr<ground_stations::PointingAnglesCalculator> pointingAnglesCalculator2 =
                    bodies.at( "Earth" )->getGroundStation( "Station2" )->getPointingAnglesCalculator( );

                std::shared_ptr<ObservationModel<1, double, double> > observationModel1 =
                    std::dynamic_pointer_cast<ObservationSimulator<1, double, double> >(
                        observationSimulators.at( 0 ))->getObservationModel( stationReceiverOneWayLinkEnds[ 0 ] );
                std::shared_ptr<ObservationModel<1, double, double> > observationModel2 =
                    std::dynamic_pointer_cast<ObservationSimulator<1, double, double> >(
                        observationSimulators.at( 0 ))->getObservationModel( stationReceiverOneWayLinkEnds[ 1 ] );

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
                        currentTime, receiver, linkEndTimes, linkEndStates );
                    Eigen::Vector3d vectorToTarget = ( linkEndStates.at( 0 ) - linkEndStates.at( 1 )).segment( 0, 3 );

                    double elevationAngle = pointingAnglesCalculator1->calculateElevationAngleFromInertialVector(
                        vectorToTarget, linkEndTimes.at( 1 ));
                    BOOST_CHECK_SMALL(( elevationAngle - currentElevation ), std::numeric_limits<double>::epsilon( ));

                    double azimuthAngle = pointingAnglesCalculator1->calculateAzimuthAngleFromInertialVector(
                        vectorToTarget, linkEndTimes.at( 1 ));
                    BOOST_CHECK_SMALL(( azimuthAngle - currentAzimuth ), std::numeric_limits<double>::epsilon( ));

                    BOOST_CHECK_SMALL(( targetRange - vectorToTarget.norm( )),
                                      std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));
                    BOOST_CHECK_SMALL(( targetInverseRange - targetRange ),
                                      std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));
                    BOOST_CHECK_SMALL(( targetInverseRange - vectorToTarget.norm( )),
                                      std::numeric_limits<double>::epsilon( ) * vectorToTarget.norm( ));

                }
            }
            else if ( observableTestCase < 3 )
            {
                std::map<double, Eigen::VectorXd> elevationAngles = getDependentVariableResultList(
                    idealObservationsAndTimes, elevationAngleSettings1, currentObservableType );
                std::map<double, Eigen::VectorXd> azimuthAngles = getDependentVariableResultList(
                    idealObservationsAndTimes, azimuthAngleSettings1, currentObservableType );
                std::map<double, Eigen::VectorXd> targetRanges = getDependentVariableResultList(
                    idealObservationsAndTimes, targetRangeSettings1, currentObservableType );

                for ( auto it: referenceElevationAngles )
                {
                    BOOST_CHECK(( elevationAngles.count( it.first ) > 0 ));
                    BOOST_CHECK(( azimuthAngles.count( it.first ) > 0 ));
                    BOOST_CHECK(( targetRanges.count( it.first ) > 0 ));

                    double currentTime = it.first;

                    BOOST_CHECK_SMALL( std::fabs(
                        elevationAngles.at( currentTime )( 0 ) - referenceElevationAngles.at( currentTime )( 0 )),
                                       std::numeric_limits<double>::epsilon( ));
                    BOOST_CHECK_SMALL( std::fabs(
                        azimuthAngles.at( currentTime )( 0 ) - referenceAzimuthAngles.at( currentTime )( 0 )),
                                       std::numeric_limits<double>::epsilon( ));
                    BOOST_CHECK_SMALL(
                        std::fabs( targetRanges.at( currentTime )( 0 ) - referenceTargetRanges.at( currentTime )( 0 )),
                        std::numeric_limits<double>::epsilon( ) * referenceTargetRanges.at( currentTime )( 0 ));

                }
            }
        }
    }
}
//
//BOOST_AUTO_TEST_SUITE_END( )
//
//}
//
//}
//
