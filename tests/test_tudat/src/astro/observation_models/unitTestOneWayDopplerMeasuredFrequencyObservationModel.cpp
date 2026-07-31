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

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/simulation/estimation_setup/createObservationModelFactory.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::earth_orientation;

BOOST_AUTO_TEST_SUITE( test_one_way_doppler_measured_frequency )

/*
 * Test 1: Core observable verification
 *
 * Verifies that the one-way Doppler measured frequency observable satisfies
 *     f_rx = f_tx * (1 + D)
 * where D is the fractional one-way Doppler computed independently.
 * Also verifies that link-end times and states from the measured frequency model
 * are consistent with those from the plain one-way Doppler model.
 */
BOOST_AUTO_TEST_CASE( testOneWayDopplerMeasuredFrequencyCoreObservable )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Sun" );
    bodiesToCreate.push_back( "Mars" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Set frequency of transmitter (Mars) to X-Band: 8.4 GHz
    double transmitterFrequency = 8.4e9;
    bodies.getBody( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmitterFrequency ) );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    // Create light-time correction settings
    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create one-way Doppler measured frequency observation settings via the convenience function
    std::shared_ptr< ObservationModelSettings > measuredFreqSettings =
            oneWayDopplerMeasuredFrequencySettings( linkEnds, lightTimeCorrectionSettings, basic_astrodynamics::tdb_scale );

    // Create a plain one-way Doppler observation model for cross-validation.
    // normalizeWithSpeedOfLight = true gives fractional Doppler D = (dt_rx/dt_tx - 1).
    std::shared_ptr< ObservationModelSettings > plainDopplerSettings = std::make_shared< OneWayDopplerObservationModelSettings >(
            linkEnds, lightTimeCorrectionSettings, nullptr, nullptr, nullptr, std::make_shared< LightTimeConvergenceCriteria >( ), true );

    // Create observation models.
    std::shared_ptr< ObservationModel< 1, double, Time > > measuredFreqModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( measuredFreqSettings, bodies );
    std::shared_ptr< ObservationModel< 1, double, Time > > plainDopplerModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( plainDopplerSettings, bodies );

    // Compute observations at mid-epoch
    Time observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

    std::vector< double > linkEndTimesMeasFreq, linkEndTimesDoppler;
    std::vector< Eigen::Vector6d > linkEndStatesMeasFreq, linkEndStatesDoppler;

    double measuredFrequency = measuredFreqModel->computeObservationsWithLinkEndData(
            observationTime, receiver, linkEndTimesMeasFreq, linkEndStatesMeasFreq )( 0 );

    double dopplerObservable = plainDopplerModel->computeObservationsWithLinkEndData(
            observationTime, receiver, linkEndTimesDoppler, linkEndStatesDoppler )( 0 );

    // Verify f_rx = f_tx * (1 + D)
    double expectedMeasuredFrequency = transmitterFrequency * ( 1.0 + dopplerObservable );
    // std::cout << std::setprecision( 16 )
    //           << "Measured frequency:  " << measuredFrequency << std::endl
    //           << "Expected frequency:  " << expectedMeasuredFrequency << std::endl
    //           << "Doppler (fractional): " << dopplerObservable << std::endl
    //           << "Difference: " << measuredFrequency - expectedMeasuredFrequency << std::endl;

    // The measured frequency and expected value should agree to numerical precision
    BOOST_CHECK_CLOSE_FRACTION( measuredFrequency, expectedMeasuredFrequency, 4.0 * std::numeric_limits< double >::epsilon( ) );

    // Verify that link-end times agree between the two models
    BOOST_CHECK_CLOSE_FRACTION(
            linkEndTimesMeasFreq.at( 0 ), linkEndTimesDoppler.at( 0 ), 4.0 * std::numeric_limits< double >::epsilon( ) );
    BOOST_CHECK_CLOSE_FRACTION(
            linkEndTimesMeasFreq.at( 1 ), linkEndTimesDoppler.at( 1 ), 4.0 * std::numeric_limits< double >::epsilon( ) );

    // Verify that link-end states agree between the two models
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            linkEndStatesMeasFreq.at( 0 ), linkEndStatesDoppler.at( 0 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
    TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
            linkEndStatesMeasFreq.at( 1 ), linkEndStatesDoppler.at( 1 ), ( 4.0 * std::numeric_limits< double >::epsilon( ) ) );
}

/*
 * Test 2: Reference link-end validation
 *
 * Verifies that the model throws a std::runtime_error when invoked with
 * 'transmitter' as the reference link end (only 'receiver' is valid).
 */
BOOST_AUTO_TEST_CASE( testOneWayDopplerMeasuredFrequencyInvalidLinkEnd )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mars" };

    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    double transmitterFrequency = 8.4e9;
    bodies.getBody( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmitterFrequency ) );

    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    std::shared_ptr< ObservationModelSettings > measuredFreqSettings = oneWayDopplerMeasuredFrequencySettings( linkEnds );

    std::shared_ptr< ObservationModel< 1, double, Time > > measuredFreqModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( measuredFreqSettings, bodies );

    Time observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
    std::vector< double > linkEndTimes;
    std::vector< Eigen::Vector6d > linkEndStates;

    // Using 'transmitter' as reference link end must throw
    BOOST_CHECK_THROW( measuredFreqModel->computeObservationsWithLinkEndData( observationTime, transmitter, linkEndTimes, linkEndStates ),
                       std::runtime_error );
}

/*
 * Test 3: Observation biases
 *
 * Verifies that absolute and relative observation biases are applied correctly:
 *     biased = absoluteBias + (1 + relativeBias) * ideal
 */
BOOST_AUTO_TEST_CASE( testOneWayDopplerMeasuredFrequencyBiases )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mars" };

    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    double transmitterFrequency = 8.4e9;
    bodies.getBody( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmitterFrequency ) );

    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Define bias values
    double absoluteBias = 1.0e3;   // 1 kHz absolute bias
    double relativeBias = 1.0e-6;  // 1 ppm relative bias

    // Create bias settings
    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( absoluteBias ), true ) );
    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( relativeBias ), false ) );
    std::shared_ptr< ObservationBiasSettings > biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );

    // Create biased observation model via the convenience function
    std::shared_ptr< ObservationModelSettings > biasedSettings =
            oneWayDopplerMeasuredFrequencySettings( linkEnds, lightTimeCorrectionSettings, basic_astrodynamics::tdb_scale, biasSettings );

    std::shared_ptr< ObservationModel< 1, double, Time > > biasedModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( biasedSettings, bodies );

    Time observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

    // Compute ideal (unbiased) and biased observations
    double idealObservation = biasedModel->computeIdealObservations( observationTime, receiver )( 0 );
    double biasedObservation = biasedModel->computeObservations( observationTime, receiver )( 0 );

    // Verify: biased = absoluteBias + (1 + relativeBias) * ideal
    double expectedBiased = absoluteBias + ( 1.0 + relativeBias ) * idealObservation;

    std::cout << std::setprecision( 16 ) << "Ideal observation:  " << idealObservation << std::endl
              << "Biased observation: " << biasedObservation << std::endl
              << "Expected biased:    " << expectedBiased << std::endl
              << "Difference: " << biasedObservation - expectedBiased << std::endl;

    BOOST_CHECK_CLOSE_FRACTION( biasedObservation, expectedBiased, 1.0E-15 );
}

/*
 * Test 4: Consistency across multiple observation times
 *
 * Verifies that the relationship f_rx = f_tx * (1 + D) holds at several
 * observation epochs spanning the simulation window.
 */
BOOST_AUTO_TEST_CASE( testOneWayDopplerMeasuredFrequencyMultipleEpochs )
{
    spice_interface::loadStandardSpiceKernels( );

    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Mars" };

    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings defaultBodySettings =
            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    defaultBodySettings.at( "Earth" )->groundStationSettings = getDsnStationSettings( );

    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    double transmitterFrequency = 8.4e9;
    bodies.getBody( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmitterFrequency ) );

    LinkEnds linkEnds;
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "DSS-13" );
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    std::vector< std::string > lightTimePerturbingBodies = { "Sun" };
    std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
    lightTimeCorrectionSettings.push_back(
            std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( lightTimePerturbingBodies ) );

    // Create both models
    std::shared_ptr< ObservationModelSettings > measuredFreqSettings =
            oneWayDopplerMeasuredFrequencySettings( linkEnds, lightTimeCorrectionSettings, basic_astrodynamics::tdb_scale );

    std::shared_ptr< ObservationModelSettings > plainDopplerSettings = std::make_shared< OneWayDopplerObservationModelSettings >(
            linkEnds, lightTimeCorrectionSettings, nullptr, nullptr, nullptr, std::make_shared< LightTimeConvergenceCriteria >( ), true );

    std::shared_ptr< ObservationModel< 1, double, Time > > measuredFreqModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( measuredFreqSettings, bodies );
    std::shared_ptr< ObservationModel< 1, double, Time > > plainDopplerModel =
            ObservationModelCreator< 1, double, Time >::createObservationModel( plainDopplerSettings, bodies );

    // Test at 5 different epochs spread across the simulation window
    int numberOfTestEpochs = 5;
    for( int i = 0; i < numberOfTestEpochs; i++ )
    {
        Time observationTime =
                initialEphemerisTime + ( i + 1 ) * ( finalEphemerisTime - initialEphemerisTime ) / ( numberOfTestEpochs + 1 );

        std::vector< double > linkEndTimesMF, linkEndTimesD;
        std::vector< Eigen::Vector6d > linkEndStatesMF, linkEndStatesD;

        double measuredFrequency =
                measuredFreqModel->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimesMF, linkEndStatesMF )( 0 );

        double dopplerObservable =
                plainDopplerModel->computeObservationsWithLinkEndData( observationTime, receiver, linkEndTimesD, linkEndStatesD )( 0 );

        double expectedMeasuredFrequency = transmitterFrequency * ( 1.0 + dopplerObservable );

        std::cout << std::setprecision( 16 ) << "Epoch " << i << ": f_rx = " << measuredFrequency
                  << ", f_tx*(1+D) = " << expectedMeasuredFrequency << ", diff = " << measuredFrequency - expectedMeasuredFrequency
                  << std::endl;

        BOOST_CHECK_CLOSE_FRACTION( measuredFrequency, expectedMeasuredFrequency, 4.0 * std::numeric_limits< double >::epsilon( ) );
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
