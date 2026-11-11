/*    Copyright (c) 2010-2019, Delft University of Technology
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
#include "tudat/simulation/environment_setup/createBodiesFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include <string>

#include <boost/test/included/unit_test.hpp>

#include "tudat/astro/ephemerides/constantEphemeris.h"
#include "tudat/astro/observation_models/angularPositionObservationModel.h"
#include "tudat/astro/observation_models/transmissionFrequencyInterface.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/math/statistics/basicStatistics.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::estimatable_parameters;
using namespace tudat::interpolators;
using namespace tudat::spice_interface;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::ephemerides;
using namespace tudat::propagators;
using namespace tudat::basic_astrodynamics;
using namespace tudat::coordinate_conversions;
using namespace tudat::statistics;

BOOST_AUTO_TEST_SUITE( test_observation_noise_models )

//! Check the conversion from Gaussian RA*cos(DEC) noise to RA noise at each observation time.
BOOST_AUTO_TEST_CASE( testAngularPositionNoiseScaling )
{
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Target" );
    linkEnds[ receiver ] = LinkEndId( "Observer" );
    const std::vector< double > declinations = { 0.0, mathematical_constants::PI / 3.0, -mathematical_constants::PI / 3.0, 1.4 };
    const std::shared_ptr< ObservationSimulationSettings<> > settings =
            tabulatedObservationSimulationSettings( angular_position, linkEnds, declinations );

    // An independent generator with the same seed provides the unscaled Gaussian samples.
    const double noiseAmplitude = 0.01;
    const int seed = noiseSeed;
    const std::function< Eigen::VectorXd( const double ) > referenceNoiseFunction =
            getIndependentGaussianNoiseFunction( noiseAmplitude, 0.0, seed, 2 );
    addGaussianNoiseToAngularPositionObservationSimulationSettings<>( { settings }, noiseAmplitude );

    for( unsigned int i = 0; i < declinations.size( ); i++ )
    {
        const double observationTime = static_cast< double >( i );
        const double declination = declinations.at( i );
        const Eigen::VectorXd unscaledNoise = referenceNoiseFunction( observationTime );
        Eigen::Vector2d observation( 0.5, declination );
        Eigen::VectorXd dependentVariables;
        addNoiseAndDependentVariableToObservation< 2 >( observation,
                                                        observationTime,
                                                        dependentVariables,
                                                        {},
                                                        {},
                                                        nullptr,
                                                        angular_position,
                                                        settings->getObservationNoiseFunction( ),
                                                        nullptr,
                                                        settings->getScaleAngularPositionNoise( ) );

        // DEC must be taken before adding its noise, and only the RA sample is scaled.
        BOOST_CHECK_SMALL( ( observation( 0 ) - 0.5 ) * std::cos( declination ) - unscaledNoise( 0 ), 1.0e-15 );
        BOOST_CHECK_SMALL( observation( 1 ) - declination - unscaledNoise( 1 ), 1.0e-15 );
    }

    const std::shared_ptr< ObservationSimulationSettings<> > rangeSettings =
            tabulatedObservationSimulationSettings( one_way_range, linkEnds, std::vector< double >{ 0.0 } );
    BOOST_CHECK_THROW( rangeSettings->setScaleAngularPositionNoise( true ), std::runtime_error );

    // Applying the RA scaling to an already normalized angular position model must be rejected.
    const std::shared_ptr< LightTimeCalculator<> > lightTimeCalculator = std::make_shared< LightTimeCalculator<> >(
            std::make_shared< ConstantEphemeris >( ( Eigen::Vector6d( ) << 1.0e8, 1.0e8, 1.0e8, 0.0, 0.0, 0.0 ).finished( ) ),
            std::make_shared< ConstantEphemeris >( Eigen::Vector6d::Zero( ) ) );
    const std::shared_ptr< ObservationModel< 2 > > normalizedModel =
            std::make_shared< AngularPositionObservationModel<> >( linkEnds, lightTimeCalculator, nullptr, true );
    BOOST_CHECK_THROW( simulateObservationWithCheck< 2 >(
                               0.0, normalizedModel, receiver, {}, settings->getObservationNoiseFunction( ), nullptr, nullptr, true ),
                       std::runtime_error );
}

//! Check the RA noise RMS through the normal observation simulation interfaces at five declinations.
BOOST_AUTO_TEST_CASE( testSimulatedAngularPositionNoiseScaling )
{
    const int numberOfObservations = 1000;
    const double noiseAmplitude = 1.0e-5;
    const double earthRadius = 6371.0e3;
    const Eigen::Vector3d stationPosition( earthRadius, 0.0, 0.0 );
    std::vector< double > observationTimes;
    for( int i = 0; i < numberOfObservations; i++ )
    {
        observationTimes.push_back( 1.0e6 + static_cast< double >( i ) );
    }

    LinkEnds linkEnds;
    linkEnds[ transmitter ] = LinkEndId( "Target" );
    linkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
    const int originalNoiseSeed = noiseSeed;

    for( const double declinationInDegrees : { 0.0, 20.0, 40.0, 60.0, 80.0 } )
    {
        BOOST_TEST_CONTEXT( "Declination: " << declinationInDegrees << " degrees" )
        {
            const double declination = declinationInDegrees * mathematical_constants::PI / 180.0;
            const double rightAscension = 0.5;
            const Eigen::Vector3d relativePosition = 1.0e11 *
                    Eigen::Vector3d( std::cos( declination ) * std::cos( rightAscension ),
                                     std::cos( declination ) * std::sin( rightAscension ),
                                     std::sin( declination ) );
            Eigen::Vector6d targetState = Eigen::Vector6d::Zero( );
            targetState.segment< 3 >( 0 ) = stationPosition + relativePosition;

            // Compute the true angles directly from the static station-to-target geometry.
            const Eigen::Vector3d lineOfSight = targetState.segment< 3 >( 0 ) - stationPosition;
            const double trueRightAscension = std::atan2( lineOfSight.y( ), lineOfSight.x( ) );
            const double trueDeclination = std::atan2( lineOfSight.z( ), lineOfSight.head< 2 >( ).norm( ) );

            BodyListSettings bodySettings;
            bodySettings.addSettings( "Earth" );
            bodySettings.addSettings( "Target" );
            bodySettings.at( "Earth" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( Eigen::Vector6d::Zero( ) );
            bodySettings.at( "Earth" )->rotationModelSettings =
                    std::make_shared< SimpleRotationModelSettings >( "ECLIPJ2000", "IAU_Earth", Eigen::Quaterniond::Identity( ), 0.0, 0.0 );
            bodySettings.at( "Earth" )->shapeModelSettings = std::make_shared< SphericalBodyShapeSettings >( earthRadius );
            bodySettings.at( "Target" )->ephemerisSettings = std::make_shared< ConstantEphemerisSettings >( targetState );
            SystemOfBodies bodies = createSystemOfBodies( bodySettings );
            createGroundStation( bodies.at( "Earth" ), "Station", stationPosition, cartesian_position );

            const std::vector< std::shared_ptr< ObservationModelSettings > > observationModelSettings = { angularPositionSettings(
                    linkEnds ) };
            const std::vector< std::shared_ptr< ObservationSimulatorBase<> > > observationSimulators =
                    createObservationSimulators( observationModelSettings, bodies );
            const std::vector< std::shared_ptr< ObservationSimulationSettings<> > > simulationSettings = {
                tabulatedObservationSimulationSettings( angular_position, linkEnds, observationTimes, receiver )
            };

            // Use the same Gaussian samples at each declination to isolate the effect of cos(DEC).
            noiseSeed = 12345;
            addGaussianNoiseToAngularPositionObservationSimulationSettings( simulationSettings, noiseAmplitude );
            const std::shared_ptr< ObservationCollection<> > simulatedObservations =
                    simulateObservations( simulationSettings, observationSimulators, bodies );
            const Eigen::VectorXd angularPositions = simulatedObservations->getSingleLinkObservations( angular_position, linkEnds );
            BOOST_REQUIRE_EQUAL( angularPositions.size( ), 2 * numberOfObservations );

            double sumSquaredRaNoise = 0.0;
            double sumSquaredDecNoise = 0.0;
            for( int i = 0; i < numberOfObservations; i++ )
            {
                const double raNoise = trueRightAscension - angularPositions( 2 * i );
                const double decNoise = trueDeclination - angularPositions( 2 * i + 1 );
                sumSquaredRaNoise += raNoise * raNoise;
                sumSquaredDecNoise += decNoise * decNoise;
            }
            const double raNoiseRms = std::sqrt( sumSquaredRaNoise / numberOfObservations );
            const double decNoiseRms = std::sqrt( sumSquaredDecNoise / numberOfObservations );
            const double expectedRaNoiseRms = noiseAmplitude / std::cos( trueDeclination );
            BOOST_TEST_MESSAGE( "DEC = " << declinationInDegrees << " deg; RA noise RMS = " << raNoiseRms
                                         << "; expected = " << expectedRaNoiseRms << "; DEC noise RMS = " << decNoiseRms );

            // Allow sampling variation in the RMS of 1000 independent Gaussian observations.
            BOOST_CHECK_CLOSE_FRACTION( raNoiseRms, expectedRaNoiseRms, 0.1 );
            BOOST_CHECK_CLOSE_FRACTION( decNoiseRms, noiseAmplitude, 0.1 );
        }
    }
    noiseSeed = originalNoiseSeed;
}

// Function to conver
double ignoreInputeVariable( std::function< double( ) > inputFreeFunction, const double dummyInput )
{
    return inputFreeFunction( );
}

//! Test whether observation noise is correctly added when simulating noisy observations
BOOST_AUTO_TEST_CASE( testObservationNoiseModels )
{
    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = double( 1.0E7 );
    double finalEphemerisTime = double( 1.0E7 + 3.0 * physical_constants::JULIAN_DAY );

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, initialEphemerisTime - 3600.0, finalEphemerisTime + 3600.0 );
    bodySettings.at( "Earth" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
            "ECLIPJ2000",
            "IAU_Earth",
            spice_interface::computeRotationQuaternionBetweenFrames( "ECLIPJ2000", "IAU_Earth", initialEphemerisTime ),
            initialEphemerisTime,
            2.0 * mathematical_constants::PI / ( physical_constants::JULIAN_DAY ) );

    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Creatre ground stations: same position, but different representation
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

    // Define link ends to/from ground stations to Moon
    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Moon", "" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Moon", "" );
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
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;

        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            // Add range bias for first two 1-way range observations
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if( ( currentObservable == one_way_range ) && ( i == 0 ) )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Constant( rangeBias1 ), false );
            }
            else if( ( currentObservable == one_way_range ) && ( i == 1 ) )
            {
                biasSettings = std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Constant( rangeBias2 ), false );
            }

            // Create observation settings
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                    currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ) );
        }
    }

    // Create observation simulators
    std::vector< std::shared_ptr< ObservationSimulatorBase< double, double > > > observationSimulators =
            createObservationSimulators( observationSettingsList, bodies );

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

    std::map< double, Eigen::VectorXd > targetAngles = getTargetAnglesAndRange(
            bodies, std::make_pair< std::string, std::string >( "Earth", "Station1" ), "Moon", baseTimeList, true );

    // Define observation simulation settings (observation type, link end, times and reference link end)
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkEnds > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( );
         linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkEnds > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings<> >(
                    currentObservable, currentLinkEndsList.at( i ), baseTimeList, receiver ) );
        }
    }

    // Simulate noise-free observations
    std::shared_ptr< ObservationCollection<> > idealObservationsAndTimes =
            simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );

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
                           inputFreeNoiseFunction,
                           std::placeholders::_1 );

        // Simulate noisy observables
        addNoiseFunctionToObservationSimulationSettings( measurementSimulationInput, noiseFunction );
        std::shared_ptr< ObservationCollection<> > constantNoiseObservationsAndTimes =
                simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );

        // Compare ideal and noise observations for each combination of observable/link ends
        for( auto observableIterator : linkEndsPerObservable )
        {
            ObservableType currentObservable = observableIterator.first;
            std::vector< LinkEnds > linkEndsList = observableIterator.second;
            for( unsigned int k = 0; k < linkEndsList.size( ); k++ )
            {
                LinkEnds currentLinkEnds = linkEndsList.at( k );

                // Compute mean and standard deviation of difference bewteen noisy and ideal observations.
                Eigen::VectorXd dataDifference =
                        constantNoiseObservationsAndTimes->getSingleLinkObservations( currentObservable, currentLinkEnds ) -
                        idealObservationsAndTimes->getSingleLinkObservations( currentObservable, currentLinkEnds );

                meanObservationDifference[ currentObservable ][ currentLinkEnds ] = computeAverageOfVectorComponents( dataDifference );
                standardDeviationObservationDifference[ currentObservable ][ currentLinkEnds ] =
                        computeStandardDeviationOfVectorComponents( dataDifference );

                // Compare with imposed mean and standard deviation of noise.
                BOOST_CHECK_CLOSE_FRACTION( meanObservationDifference[ currentObservable ][ currentLinkEnds ], constantOffset, 1.0E-2 );
                BOOST_CHECK_CLOSE_FRACTION(
                        standardDeviationObservationDifference[ currentObservable ][ currentLinkEnds ], constantStandardDeviation, 1.0E-2 );
            }
        }
    }

    // Test noise simulation, with difference distribution for each observable.
    {
        // Define (arbitrary) noise properties for observables
        std::map< ObservableType, double > constantOffsets;
        constantOffsets[ one_way_range ] = -200.0;
        constantOffsets[ one_way_doppler ] = -2.8E3;
        constantOffsets[ angular_position ] = 3.0E-4;

        std::map< ObservableType, double > constantStandardDeviations;
        constantStandardDeviations[ one_way_range ] = 2.4;
        constantStandardDeviations[ one_way_doppler ] = 7.5;
        constantStandardDeviations[ angular_position ] = 6.3E-6;

        clearNoiseFunctionFromObservationSimulationSettings( measurementSimulationInput );
        // Create noise function for each observable
        std::map< ObservableType, std::function< double( const double ) > > noiseFunctionPerObservable;
        for( std::map< ObservableType, double >::const_iterator typeIterator = constantOffsets.begin( );
             typeIterator != constantOffsets.end( );
             typeIterator++ )
        {
            std::function< double( const double ) > noiseFunction =
                    std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                               createBoostContinuousRandomVariableGeneratorFunction(
                                       normal_boost_distribution,
                                       { constantOffsets.at( typeIterator->first ), constantStandardDeviations.at( typeIterator->first ) },
                                       0.0 ),
                               std::placeholders::_1 );

            addNoiseFunctionToObservationSimulationSettings( measurementSimulationInput, noiseFunction, typeIterator->first );
        }

        // Simulate noisy observables
        std::shared_ptr< ObservationCollection<> > noisyPerObservableObservationsAndTimes =
                simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );

        // Compare ideal and noise observations for each combination of observable/link ends
        for( auto observableIterator : linkEndsPerObservable )
        {
            ObservableType currentObservable = observableIterator.first;
            std::vector< LinkEnds > linkEndsList = observableIterator.second;
            for( unsigned int k = 0; k < linkEndsList.size( ); k++ )
            {
                LinkEnds currentLinkEnds = linkEndsList.at( k );

                // Compute mean and standard deviation of difference bewteen noisy and ideal observations.
                Eigen::VectorXd dataDifference =
                        noisyPerObservableObservationsAndTimes->getSingleLinkObservations( currentObservable, currentLinkEnds ) -
                        idealObservationsAndTimes->getSingleLinkObservations( currentObservable, currentLinkEnds );

                meanObservationDifference[ currentObservable ][ currentLinkEnds ] = computeAverageOfVectorComponents( dataDifference );
                standardDeviationObservationDifference[ currentObservable ][ currentLinkEnds ] =
                        computeStandardDeviationOfVectorComponents( dataDifference );

                // Compare with imposed mean and standard deviation of noise.
                BOOST_CHECK_CLOSE_FRACTION(
                        meanObservationDifference[ currentObservable ][ currentLinkEnds ], constantOffsets[ currentObservable ], 1.0E-2 );
                BOOST_CHECK_CLOSE_FRACTION( standardDeviationObservationDifference[ currentObservable ][ currentLinkEnds ],
                                            constantStandardDeviations[ currentObservable ],
                                            1.0E-2 );
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

        constantOffsetsPerStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 0 ) ] = 4.3E2;
        constantOffsetsPerStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 1 ) ] = -3.4E3;

        constantOffsetsPerStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 0 ) ] = 5.3E-7;
        constantOffsetsPerStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 1 ) ] = 3.33E-6;

        std::map< ObservableType, std::map< LinkEnds, double > > constantStandardDeviationsStation;
        constantStandardDeviationsStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 0 ) ] = 0.65;
        constantStandardDeviationsStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 1 ) ] = 1.34;
        constantStandardDeviationsStation[ one_way_range ][ linkEndsPerObservable[ one_way_range ].at( 2 ) ] = 4.33;

        constantStandardDeviationsStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 0 ) ] = 2.6;
        constantStandardDeviationsStation[ one_way_doppler ][ linkEndsPerObservable[ one_way_doppler ].at( 1 ) ] = 2.2;

        constantStandardDeviationsStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 0 ) ] = 1.2E-12;
        constantStandardDeviationsStation[ angular_position ][ linkEndsPerObservable[ angular_position ].at( 1 ) ] = 4.3E-10;

        clearNoiseFunctionFromObservationSimulationSettings( measurementSimulationInput );

        // Create noise function for each observable and link ends combination
        std::map< ObservableType, std::function< double( const double ) > > noiseFunctionPerObservable;
        std::map< ObservableType, std::map< LinkEnds, std::function< double( const double ) > > > noiseFunctionPerLinkEnd;
        for( std::map< ObservableType, std::map< LinkEnds, double > >::const_iterator typeIterator = constantOffsetsPerStation.begin( );
             typeIterator != constantOffsetsPerStation.end( );
             typeIterator++ )
        {
            for( std::map< LinkEnds, double >::const_iterator linkEndIterator = typeIterator->second.begin( );
                 linkEndIterator != typeIterator->second.end( );
                 linkEndIterator++ )
            {
                std::function< double( const double ) > noiseFunction =
                        std::bind( &utilities::evaluateFunctionWithoutInputArgumentDependency< double, const double >,
                                   createBoostContinuousRandomVariableGeneratorFunction(
                                           normal_boost_distribution,
                                           { constantOffsetsPerStation.at( typeIterator->first ).at( linkEndIterator->first ),
                                             constantStandardDeviationsStation.at( typeIterator->first ).at( linkEndIterator->first ) },
                                           0.0 ),
                                   std::placeholders::_1 );

                addNoiseFunctionToObservationSimulationSettings(
                        measurementSimulationInput, noiseFunction, typeIterator->first, linkEndIterator->first );
            }
        }

        // Simulate noisy observables
        std::shared_ptr< ObservationCollection<> > noisyPerObservableAndLinkEndsObservationsAndTimes =
                simulateObservations< double, double >( measurementSimulationInput, observationSimulators, bodies );

        // Compare ideal and noise observations for each combination of observable/link ends
        for( auto observableIterator : linkEndsPerObservable )
        {
            ObservableType currentObservable = observableIterator.first;
            std::vector< LinkEnds > linkEndsList = observableIterator.second;
            for( unsigned int k = 0; k < linkEndsList.size( ); k++ )
            {
                LinkEnds currentLinkEnds = linkEndsList.at( k );

                // Compute mean and standard deviation of difference bewteen noisy and ideal observations.
                Eigen::VectorXd dataDifference =
                        noisyPerObservableAndLinkEndsObservationsAndTimes->getSingleLinkObservations( currentObservable, currentLinkEnds ) -
                        idealObservationsAndTimes->getSingleLinkObservations( currentObservable, currentLinkEnds );
                meanObservationDifference[ currentObservable ][ currentLinkEnds ] = computeAverageOfVectorComponents( dataDifference );
                standardDeviationObservationDifference[ currentObservable ][ currentLinkEnds ] =
                        computeStandardDeviationOfVectorComponents( dataDifference );

                // Compare with imposed mean and standard deviation of noise.
                BOOST_CHECK_CLOSE_FRACTION( meanObservationDifference[ currentObservable ][ currentLinkEnds ],
                                            constantOffsetsPerStation[ currentObservable ][ currentLinkEnds ],
                                            1.0E-2 );
                BOOST_CHECK_CLOSE_FRACTION( standardDeviationObservationDifference[ currentObservable ][ currentLinkEnds ],
                                            constantStandardDeviationsStation[ currentObservable ][ currentLinkEnds ],
                                            1.0E-2 );
            }
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
