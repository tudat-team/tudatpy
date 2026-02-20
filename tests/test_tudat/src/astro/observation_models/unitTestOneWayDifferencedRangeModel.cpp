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

#include "tudat/io/basicInputOutput.h"

#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/observation_models/oneWayRangeObservationModel.h"
#include "tudat/simulation/estimation_setup/createObservationModel.h"
#include "tudat/simulation/environment_setup/defaultBodies.h"
#include "tudat/simulation/environment_setup/createBodies.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat;
using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;

BOOST_AUTO_TEST_SUITE( test_one_way_doppler_model )

double integrationTimeFunction( const double currentObservationTime )
{
    return 60.0 + 30.0 * ( currentObservationTime - 3.0 * 86400.0 ) / ( 7.0 * 86400.0 );
}
//
// BOOST_AUTO_TEST_CASE( testOneWayDoppplerModel )
//{
//    // Load Spice kernels
//    spice_interface::loadStandardSpiceKernels( );
//
//    // Define bodies to use.
//    std::vector< std::string > bodiesToCreate;
//    bodiesToCreate.push_back( "Earth" );
//    bodiesToCreate.push_back( "Sun" );
//    bodiesToCreate.push_back( "Mars" );
//
//    // Specify initial time
//    double initialEphemerisTime = 0.0;
//    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
//    double maximumTimeStep = 3600.0;
//    double buffer = 10.0 * maximumTimeStep;
//
//    // Create bodies settings needed in simulation
//    BodyListSettings defaultBodySettings =
//            getDefaultBodySettings( bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
//
//    // Create bodies
//    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );
//
//    // Define link ends for observations.
//    LinkEnds linkEnds;
//    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "" );
//    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars", "" );
//
//    // Create range rate observation settings and model
//    std::shared_ptr< ObservationModelSettings > rangeRateObservableSettings =
//            std::make_shared< ObservationModelSettings >( one_way_differenced_range, linkEnds );
//    std::shared_ptr< ObservationModel< 1, double, double > > rangeRateObservationModel =
//            ObservationModelCreator< 1, double, double >::createObservationModel( rangeRateObservableSettings, bodies );
//
//    // Create range rate observation settings and model
//    std::shared_ptr< ObservationModelSettings > rangeObservableSettings =
//            std::make_shared< ObservationModelSettings >( one_way_range, linkEnds );
//    std::shared_ptr< ObservationModel< 1, double, double > > rangeObservationModel =
//            ObservationModelCreator< 1, double, double >::createObservationModel( rangeObservableSettings, bodies );
//
//    // Test observable for both fixed link ends
//    for( unsigned testCase = 0; testCase < 2; testCase++ )
//    {
//        for( double observationTime = 86400.0; observationTime <= 86400.0; observationTime += observationTime )
//        {
//            std::cout << "TEST: ************************************* " << testCase << " " << observationTime << std::endl;
//            double dopplerCountInterval = integrationTimeFunction( observationTime );
//            double arcStartObservationTime = observationTime - dopplerCountInterval / 2.0;
//            double arcEndObservationTime = observationTime + dopplerCountInterval / 2.0;
//            std::vector< double > rangeRateLinkEndTimes;
//            std::vector< Eigen::Vector6d > rangeRateLinkEndStates;
//
//            std::vector< double > rangeStartLinkEndTimes;
//            std::vector< Eigen::Vector6d > rangeStartLinkEndStates;
//
//            std::vector< double > rangeEndLinkEndTimes;
//            std::vector< Eigen::Vector6d > rangeEndLinkEndStates;
//
//            // Define link end
//            LinkEndType referenceLinkEnd;
//            if( testCase == 0 )
//            {
//                referenceLinkEnd = transmitter;
//            }
//            else
//            {
//                referenceLinkEnd = receiver;
//            }
//
//            // Compute observable
//            double rangeRateObservable = rangeRateObservationModel->computeObservationsWithLinkEndData(
//                    observationTime,
//                    referenceLinkEnd,
//                    rangeRateLinkEndTimes,
//                    rangeRateLinkEndStates,
//                    getAveragedDopplerAncillarySettings( dopplerCountInterval ) )( 0 );
//
//            double arcEndRange = rangeObservationModel->computeObservationsWithLinkEndData(
//                    arcEndObservationTime, referenceLinkEnd, rangeEndLinkEndTimes, rangeEndLinkEndStates )( 0 );
//            double arcStartRange = rangeObservationModel->computeObservationsWithLinkEndData(
//                    arcStartObservationTime, referenceLinkEnd, rangeStartLinkEndTimes, rangeStartLinkEndStates )( 0 );
//
//            for( unsigned linkCompare = 0; linkCompare < 2; linkCompare++ )
//            {
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        rangeStartLinkEndStates.at( linkCompare ), rangeRateLinkEndStates.at( linkCompare ), 1.0E-15 );
//                BOOST_CHECK_CLOSE_FRACTION( rangeStartLinkEndTimes.at( linkCompare ), rangeRateLinkEndTimes.at( linkCompare ), 1.0E-15 );
//
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION(
//                        rangeEndLinkEndStates.at( linkCompare ), rangeRateLinkEndStates.at( linkCompare + 2 ), 1.0E-15 );
//                BOOST_CHECK_CLOSE_FRACTION( rangeEndLinkEndTimes.at( linkCompare ), rangeRateLinkEndTimes.at( linkCompare + 2 ), 1.0E-15
//                );
//            }
//
//            double manualDifferencedRange = ( arcEndRange - arcStartRange ) / dopplerCountInterval;
//
//            // Test numerical derivative against Doppler observable
//            BOOST_CHECK_SMALL( std::fabs( manualDifferencedRange - rangeRateObservable ), 1.0E-4 * dopplerCountInterval );
//        }
//    }
//}

BOOST_AUTO_TEST_CASE( testTwoWayRangeWithFrequencyCorrections )
{
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies to use.
    std::vector< std::string > bodiesToCreate;
    bodiesToCreate.push_back( "Earth" );
    bodiesToCreate.push_back( "Mars" );
    bodiesToCreate.push_back( "Sun" );

    // Specify initial time
    double initialEphemerisTime = 0.0;
    double finalEphemerisTime = initialEphemerisTime + 7.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 10.0 * maximumTimeStep;

    // Create bodies settings needed in simulation
    BodyListSettings defaultBodySettings = getDefaultBodySettings( bodiesToCreate );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );
    bodies.at( "Mars" )->getVehicleSystems( )->setDefaultTransponderTurnaroundRatio( );

    // Create ground stations
    std::pair< std::string, std::string > earthStationStation = std::pair< std::string, std::string >( "Earth", "EarthStation" );

    createGroundStation( bodies.at( "Earth" ),
                         "EarthStation",
                         ( Eigen::Vector3d( ) << 1.0, 0.1, -1.4 ).finished( ),
                         coordinate_conversions::geodetic_position );

    double transmittingFrequencyEarth = 8.0E6;
    bodies.at( "Earth" )
            ->getGroundStation( "EarthStation" )
            ->setTransmittingFrequencyCalculator(
                    std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmittingFrequencyEarth ) );

    double transmittingFrequencyMars = 2.0E6;
    bodies.at( "Mars" )->getVehicleSystems( )->setTransmittedFrequencyCalculator(
            std::make_shared< ground_stations::ConstantFrequencyInterpolator >( transmittingFrequencyMars ) );

    // Define list of observation times for which to check model
    std::vector< double > observationTimes;
    observationTimes.push_back( 1.0E8 );
    observationTimes.push_back( 2.0E8 );
    observationTimes.push_back( 3.0E8 );

    // Define link ends for observations.
    for( unsigned int observationTimeNumber = 0; observationTimeNumber < observationTimes.size( ); observationTimeNumber++ )
    {
        for( unsigned int referenceNumber = 0; referenceNumber < 2; referenceNumber++ )
        {
            // Define link ends for 2-way model and constituent one-way models
            LinkEnds oneWayLinkEnds;
            if( referenceNumber == 0 )
            {
                oneWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation" );
                oneWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars", "" );
            }
            else
            {
                oneWayLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth", "EarthStation" );
                oneWayLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );
            }
            // Create light-time correction settings.
            std::vector< std::shared_ptr< LightTimeCorrectionSettings > > lightTimeCorrectionSettings;
            lightTimeCorrectionSettings.push_back( std::make_shared< JakowskiIonosphericCorrectionSettings >( ) );

            // Create observation settings
            std::shared_ptr< ObservationModelSettings > observableSettingsUncorrected =
                    std::make_shared< ObservationModelSettings >( one_way_differenced_range, oneWayLinkEnds );
            std::shared_ptr< ObservationModelSettings > observableSettingsCorrected =
                    std::make_shared< ObservationModelSettings >( one_way_differenced_range, oneWayLinkEnds, lightTimeCorrectionSettings );

            // Create observation models.
            std::shared_ptr< ObservationModel< 1, double, double > > observationModelCorrected =
                    ObservationModelCreator< 1, double, double >::createObservationModel( observableSettingsCorrected, bodies );
            std::shared_ptr< ObservationModel< 1, double, double > > observationModelUncorrected =
                    ObservationModelCreator< 1, double, double >::createObservationModel( observableSettingsUncorrected, bodies );

            // Compute correction
            std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings =
                    std::make_shared< ObservationAncillarySimulationSettings >( );

            double integrationTime = 7200.0;
            ancillarySettings->setAncillaryDoubleVectorData( frequency_bands, { x_band, x_band } );
            ancillarySettings->setAncillaryDoubleData( doppler_integration_time, integrationTime );

            std::vector< double > linkEndTimes;
            std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;

            double testTime = observationTimes.at( observationTimeNumber );
            double uncorrectedObservation = observationModelUncorrected->computeIdealObservationsWithLinkEndData(
                    testTime, receiver, linkEndTimes, linkEndStates, ancillarySettings )( 0 );

            // Compute corrected range
            double correctedObservation = observationModelCorrected->computeIdealObservations( testTime, receiver, ancillarySettings )( 0 );

            if( referenceNumber == 0 )
            {
                ancillarySettings->setIntermediateDoubleData( transmitter_frequency_intermediate, transmittingFrequencyEarth );
            }
            else
            {
                ancillarySettings->setIntermediateDoubleData( transmitter_frequency_intermediate, transmittingFrequencyMars );
            }
            std::shared_ptr< LightTimeCorrection > ionosphereCorrectionModel =
                    createLightTimeCorrections( std::make_shared< JakowskiIonosphericCorrectionSettings >( ),
                                                bodies,
                                                oneWayLinkEnds,
                                                transmitter,
                                                receiver,
                                                one_way_range );

            std::vector< double > subLinkEndTimes = { linkEndTimes.at( 0 ), linkEndTimes.at( 1 ) };
            std::vector< Eigen::Matrix< double, 6, 1 > > subLinkEndStates = { linkEndStates.at( 0 ), linkEndStates.at( 1 ) };
            double ionosphereCorrectionStart = ionosphereCorrectionModel->calculateLightTimeCorrectionWithMultiLegLinkEndStates(
                    subLinkEndStates, subLinkEndTimes, 0, ancillarySettings );

            subLinkEndTimes = { linkEndTimes.at( 2 ), linkEndTimes.at( 3 ) };
            subLinkEndStates = { linkEndStates.at( 2 ), linkEndStates.at( 3 ) };
            double ionosphereCorrectionEnd = ionosphereCorrectionModel->calculateLightTimeCorrectionWithMultiLegLinkEndStates(
                    subLinkEndStates, subLinkEndTimes, 0, ancillarySettings );

            //            std::cout<<( ionosphereCorrectionStart - ionosphereCorrectionEnd ) * physical_constants::SPEED_OF_LIGHT /
            //            integrationTime<<std::endl; std::cout<<correctedObservation - uncorrectedObservation<<std::endl;

            // Compare observation difference with direct correction (limited by double precision over several AU)
            BOOST_CHECK_CLOSE_FRACTION(
                    ( ionosphereCorrectionStart - ionosphereCorrectionEnd ) * physical_constants::SPEED_OF_LIGHT / integrationTime,
                    ( correctedObservation - uncorrectedObservation ),
                    1.0E-4 );
        }
    }
}
BOOST_AUTO_TEST_SUITE_END( )

}  // namespace unit_tests

}  // namespace tudat
