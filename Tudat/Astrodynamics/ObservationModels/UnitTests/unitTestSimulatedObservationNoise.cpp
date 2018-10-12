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

#include "Tudat/SimulationSetup/tudatEstimationHeader.h"

namespace tudat
{
namespace unit_tests
{

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

BOOST_AUTO_TEST_SUITE( test_observation_noise_models )

// Function to conver
double ignoreInputeVariable( std::function< double( ) > inputFreeFunction, const double dummyInput )
{
    return inputFreeFunction( );
}

//! Test whether observation noise is correctly added when simulating noisy observations
BOOST_AUTO_TEST_CASE( testObservationNoiseModels )
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
    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0, finalEphemerisTime + 3600.0 );
    bodySettings[ "Earth" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Earth",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY ) );

    NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Creatre ground stations: same position, but different representation
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodyMap.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ), geodetic_position );
    createGroundStation( bodyMap.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ), geodetic_position );

    // Define parameters.
    std::vector< LinkEnds > stationReceiverLinkEnds;
    std::vector< LinkEnds > stationTransmitterLinkEnds;

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
    }

    // Define (arbitrary) link ends for each observable
    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    linkEndsPerObservable[ one_way_doppler ].push_back( stationReceiverLinkEnds[ 1 ] );
    linkEndsPerObservable[ one_way_doppler ].push_back( stationTransmitterLinkEnds[ 2 ] );

    linkEndsPerObservable[ angular_position ].push_back( stationReceiverLinkEnds[ 2 ] );
    linkEndsPerObservable[ angular_position ].push_back( stationTransmitterLinkEnds[ 1 ] );

    // Set range biases for two range links
    double rangeBias1 = 3.0;
    double rangeBias2 = -7.2;

    // Define observation settings for each observable/link ends combination
    observation_models::ObservationSettingsMap observationSettingsMap;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Add range bias for first two 1-way range observations
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if( ( currentObservable == one_way_range ) && ( i == 0 ) )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                            Eigen::Vector1d::Constant( rangeBias1 ), false );
            }
            else if( ( currentObservable == one_way_range ) && ( i == 1 ) )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >(
                            Eigen::Vector1d::Constant( rangeBias2 ), false );
            }

            // Create observation settings
            observationSettingsMap.insert(
                        std::make_pair( currentLinkEndsList.at( i ),
                                        std::make_shared< ObservationSettings >(
                                            currentObservable, std::shared_ptr< LightTimeCorrectionSettings >( ),
                                            biasSettings ) ) );
        }
    }

    // Create observation simulators
    std::map< ObservableType,
            std::shared_ptr< ObservationSimulatorBase< double, double > > >  observationSimulators =
            createObservationSimulators( observationSettingsMap, bodyMap );

    // Define osbervation times. NOTE: These times are not checked w.r.t. visibility and are used for testing purposes only.
    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double observationInterval = 5.0;
    for( unsigned int i = 0; i < 3; i++ )
    {
        for( unsigned int j = 0; j < 10000; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }

    // Define observation simulation settings (observation type, link end, times and reference link end)
    std::map< ObservableType, std::map< LinkEnds, std::shared_ptr< ObservationSimulationTimeSettings< double > > > >
            measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput[ currentObservable ][ currentLinkEndsList.at( i ) ] =
                    std::make_shared< TabulatedObservationSimulationTimeSettings< double > >(
                        receiver, baseTimeList );
        }
    }


    typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > >
            SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    // Simulate noise-free observations
    PodInputDataType idealObservationsAndTimes = simulateObservations< double, double >(
                measurementSimulationInput, observationSimulators );

    std::map< ObservableType, std::map< LinkEnds, std::vector< double > > > observationDifference;
    std::map< ObservableType, std::map< LinkEnds, double > > meanObservationDifference;
    std::map< ObservableType, std::map< LinkEnds, double > > standardDeviationObservationDifference;

    // Test noise simulation, with single, constant, distribution for all observables and link ends
    {
        // Define (arbitrary) noise properties
        double constantOffset = 12.0;
        double constantStandardDeviation = 0.5;

        // Create noise function
        std::function< double( ) > inputFreeNoiseFunction = createBoostContinuousRandomVariableGeneratorFunction(
                    normal_boost_distribution, { constantOffset, constantStandardDeviation }, 0.0 );
        std::function< double( const double ) > noiseFunction =
                std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                           inputFreeNoiseFunction, std::placeholders::_1 );

        // Simulate noisy observables
        PodInputDataType constantNoiseObservationsAndTimes = simulateObservationsWithNoise< double, double >(
                    measurementSimulationInput, observationSimulators, noiseFunction );

        // Compare ideal and noise observations for each combination of observable/link ends
        for( PodInputDataType::const_iterator dataIterator = constantNoiseObservationsAndTimes.begin( );
             dataIterator != constantNoiseObservationsAndTimes.end( ); dataIterator++ )
        {
            for( SingleObservablePodInputType::const_iterator innerDataIterator = dataIterator->second.begin( );
                 innerDataIterator != dataIterator->second.end( ); innerDataIterator++ )
            {
                // Compute mean and standard deviation of difference bewteen noisy and ideal observations.
                Eigen::VectorXd dataDifference = innerDataIterator->second.first -
                        idealObservationsAndTimes.at( dataIterator->first ).at( innerDataIterator->first ).first;
                meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                        computeAverageOfVectorComponents( dataDifference );
                standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                        computeStandardDeviationOfVectorComponents( dataDifference );

                // Compare with imposed mean and standard deviation of noise.
                BOOST_CHECK_CLOSE_FRACTION(
                            meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ],
                        constantOffset, 1.0E-2 );
                BOOST_CHECK_CLOSE_FRACTION(
                            standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ],
                        constantStandardDeviation, 1.0E-2 );
            }
        }
    }

    // Test noise simulation, with difference distribution for each observable.
    {

        // Define (arbitrary) noise properties for observables
        std::map< ObservableType, double > constantOffsets;
        constantOffsets[ one_way_range ] = -200.0;
        constantOffsets[ one_way_doppler ] = -2.8E-5;
        constantOffsets[ angular_position ] = 3.0E-4;

        std::map< ObservableType, double > constantStandardDeviations;
        constantStandardDeviations[ one_way_range ] = 2.4;
        constantStandardDeviations[ one_way_doppler ] = 7.5E-8;
        constantStandardDeviations[ angular_position ] = 6.3E-6;

        // Create noise function for each observable
        std::map< ObservableType, std::function< double( const double ) > > noiseFunctionPerObservable;
        for( std::map< ObservableType, double >::const_iterator typeIterator = constantOffsets.begin( );
             typeIterator != constantOffsets.end( ); typeIterator++ )
        {
            noiseFunctionPerObservable[ typeIterator->first ] =
                    std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                               createBoostContinuousRandomVariableGeneratorFunction(
                                   normal_boost_distribution,
            { constantOffsets.at( typeIterator->first ), constantStandardDeviations.at( typeIterator->first ) },
                                   0.0 ), std::placeholders::_1 );
        }

        // Simulate noisy observables
        PodInputDataType noisyPerObservableObservationsAndTimes = simulateObservationsWithNoise< double, double >(
                    measurementSimulationInput, observationSimulators, noiseFunctionPerObservable );

        // Compare ideal and noise observations for each combination of observable/link ends
        for( PodInputDataType::const_iterator dataIterator = noisyPerObservableObservationsAndTimes.begin( );
             dataIterator != noisyPerObservableObservationsAndTimes.end( ); dataIterator++ )
        {
            for( SingleObservablePodInputType::const_iterator innerDataIterator = dataIterator->second.begin( );
                 innerDataIterator != dataIterator->second.end( ); innerDataIterator++ )
            {
                // Compute mean and standard deviation of difference bewteen noisy and ideal observations.
                Eigen::VectorXd dataDifference = innerDataIterator->second.first -
                        idealObservationsAndTimes.at( dataIterator->first ).at( innerDataIterator->first ).first;
                meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                        computeAverageOfVectorComponents( dataDifference );
                standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                        computeStandardDeviationOfVectorComponents( dataDifference );

                // Compare with imposed mean and standard deviation of noise.
                BOOST_CHECK_CLOSE_FRACTION(
                            meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ],
                        constantOffsets[ dataIterator->first ], 1.0E-2 );
                BOOST_CHECK_CLOSE_FRACTION(
                            standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ],
                        constantStandardDeviations[ dataIterator->first ], 1.0E-2 );
            }
        }

    }

    // Test noise simulation, with difference distribution for each observable and set of link ends.
    {

        // Define (arbitrary) noise properties for observable, per link ends
        std::map< ObservableType, std::map< LinkEnds, double > > constantOffsetsPerStation;
        constantOffsetsPerStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 0 ) ] = 2.4;
        constantOffsetsPerStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 1 ) ] = -65.3;
        constantOffsetsPerStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 2 ) ] = 54.1;

        constantOffsetsPerStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 0 ) ] = 4.3E-6;
        constantOffsetsPerStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 1 ) ] = -3.4E-5;

        constantOffsetsPerStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 0 ) ] = 5.3E-7;
        constantOffsetsPerStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 1 ) ] = 3.33E-6;

        std::map< ObservableType, std::map< LinkEnds, double > > constantStandardDeviationsStation;
        constantStandardDeviationsStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 0 ) ] = 0.65;
        constantStandardDeviationsStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 1 ) ] = 1.34;
        constantStandardDeviationsStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 2 ) ] = 4.33;

        constantStandardDeviationsStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 0 ) ] = 2.6E-8;
        constantStandardDeviationsStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 1 ) ] = 2.2E-8;;

        constantStandardDeviationsStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 0 ) ] = 1.2E-12;
        constantStandardDeviationsStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 1 ) ] = 4.3E-10;

        // Create noise function for each observable and link ends combination
        std::map< ObservableType, std::function< double( const double ) > > noiseFunctionPerObservable;
        std::map< ObservableType, std::map< LinkEnds, std::function< double( const double ) > > > noiseFunctionPerLinkEnd;
        for( std::map< ObservableType, std::map< LinkEnds, double > >::const_iterator
             typeIterator = constantOffsetsPerStation.begin( );
             typeIterator != constantOffsetsPerStation.end( );
             typeIterator++ )
        {
            for( std::map< LinkEnds, double >::const_iterator linkEndIterator = typeIterator->second.begin( );
                 linkEndIterator != typeIterator->second.end( ); linkEndIterator++ )
            {
                noiseFunctionPerLinkEnd[ typeIterator->first ][ linkEndIterator->first ] =
                        std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                                   createBoostContinuousRandomVariableGeneratorFunction(
                                       normal_boost_distribution,
                { constantOffsetsPerStation.at( typeIterator->first ).at( linkEndIterator->first ), constantStandardDeviationsStation.at( typeIterator->first ).at( linkEndIterator->first ) },
                                       0.0 ), std::placeholders::_1 );
            }
        }

        // Simulate noisy observables
        PodInputDataType noisyPerObservableAndLinkEndsObservationsAndTimes = simulateObservationsWithNoise< double, double >(
                    measurementSimulationInput, observationSimulators, noiseFunctionPerLinkEnd );

        // Compare ideal and noise observations for each combination of observable/link ends
        for( PodInputDataType::const_iterator dataIterator = noisyPerObservableAndLinkEndsObservationsAndTimes.begin( );
             dataIterator != noisyPerObservableAndLinkEndsObservationsAndTimes.end( ); dataIterator++ )
        {
            for( SingleObservablePodInputType::const_iterator innerDataIterator = dataIterator->second.begin( );
                 innerDataIterator != dataIterator->second.end( ); innerDataIterator++ )
            {
                // Compute mean and standard deviation of difference bewteen noisy and ideal observations.
                Eigen::VectorXd dataDifference = innerDataIterator->second.first -
                        idealObservationsAndTimes.at( dataIterator->first ).at( innerDataIterator->first ).first;
                meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                        computeAverageOfVectorComponents( dataDifference );
                standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ] =
                        computeStandardDeviationOfVectorComponents( dataDifference );

                // Compare with imposed mean and standard deviation of noise.
                BOOST_CHECK_CLOSE_FRACTION(
                            meanObservationDifference[ dataIterator->first ][ innerDataIterator->first ],
                        constantOffsetsPerStation[ dataIterator->first ][ innerDataIterator->first ], 1.0E-2 );
                BOOST_CHECK_CLOSE_FRACTION(
                            standardDeviationObservationDifference[ dataIterator->first ][ innerDataIterator->first ],
                        constantStandardDeviationsStation[ dataIterator->first ][ innerDataIterator->first ], 1.0E-2 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}

