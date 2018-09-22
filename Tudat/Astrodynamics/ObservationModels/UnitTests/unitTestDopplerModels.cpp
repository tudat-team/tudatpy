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

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ObservationModels/oneWayRangeObservationModel.h"
#include "Tudat/SimulationSetup/EstimationSetup/createObservationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createBodies.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

namespace tudat
{
namespace unit_tests
{

using namespace tudat::observation_models;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::simulation_setup;
using namespace tudat::orbital_element_conversions;
using namespace tudat::coordinate_conversions;
using namespace tudat::unit_conversions;


BOOST_AUTO_TEST_SUITE( test_doppler_models )


BOOST_AUTO_TEST_CASE( testOneWayDoppplerModel )
{
    // Load Spice kernels
    std::string kernelsPath = input_output::getSpiceKernelPath( );
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
    std::map< std::string, std::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    // Create ground station
    const Eigen::Vector3d stationCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    createGroundStation( bodyMap.at( "Earth" ), "Station1", stationCartesianPosition, cartesian_position );

    // Create Spacecraft
    Eigen::Vector6d spacecraftOrbitalElements;
    spacecraftOrbitalElements( semiMajorAxisIndex ) = 10000.0E3;
    spacecraftOrbitalElements( eccentricityIndex ) = 0.33;
    spacecraftOrbitalElements( inclinationIndex ) = convertDegreesToRadians( 65.3 );
    spacecraftOrbitalElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    spacecraftOrbitalElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    spacecraftOrbitalElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodyMap[ "Spacecraft" ] = std::make_shared< Body >( );
    bodyMap[ "Spacecraft" ]->setEphemeris(
                createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                         spacecraftOrbitalElements, 0.0, earthGravitationalParameter, "Earth" ), "Spacecraft" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
    linkEnds[ receiver ] = std::make_pair( "Mars" , ""  );

    // Create observation settings
    std::shared_ptr< ObservationSettings > observableSettings = std::make_shared< ObservationSettings >
            ( one_way_doppler );

    // Create observation model.
    std::shared_ptr< ObservationModel< 1, double, double> > observationModel =
            ObservationModelCreator< 1, double, double>::createObservationModel(
                linkEnds, observableSettings, bodyMap );

    std::shared_ptr< OneWayDopplerObservationModel< double, double> > dopplerObservationModel =
            std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >( observationModel );

    // Test observable for both fixed link ends
    for( unsigned testCase = 0; testCase < 2; testCase++ )
    {

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Vector6d > linkEndStates;

        // Define link end
        LinkEndType referenceLinkEnd;
        if( testCase == 0 )
        {
            referenceLinkEnd = transmitter;
        }
        else
        {
            referenceLinkEnd = receiver;
        }

        // Compute observable
        double dopplerObservable = observationModel->computeObservationsWithLinkEndData(
                    observationTime, referenceLinkEnd, linkEndTimes, linkEndStates )( 0 );

        // Creare independent light time calculator object
        std::shared_ptr< LightTimeCalculator< double, double > > lightTimeCalculator =
                createLightTimeCalculator( linkEnds[ transmitter ], linkEnds[ receiver ], bodyMap );
        Eigen::Vector6d transmitterState, receiverState;
        // Compute light time
        double lightTime = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, observationTime, testCase );

        // Compare light time calculator link end conditions with observation model
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( receiverState, linkEndStates.at( 1 ), 1.0E-15 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transmitterState, linkEndStates.at( 0 ), 1.0E-15 );

            if( testCase == 0 )
            {
                BOOST_CHECK_SMALL( std::fabs( observationTime  - linkEndTimes.at( 0 ) ), 1.0E-12 );
                BOOST_CHECK_SMALL( std::fabs( observationTime + lightTime - linkEndTimes.at( 1 ) ), 1.0E-10 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( observationTime - linkEndTimes.at( 1 ) ), 1.0E-12 );
                BOOST_CHECK_SMALL( std::fabs( observationTime - lightTime - linkEndTimes.at( 0 ) ), 1.0E-10 );
            }
        }

        // Compute numerical partial derivative of light time.
        double timePerturbation = 100.0;
        double upPerturbedLightTime = lightTimeCalculator->calculateLightTime( linkEndTimes.at( 1 ) + timePerturbation, true );
        double downPerturbedLightTime = lightTimeCalculator->calculateLightTime( linkEndTimes.at( 1 ) - timePerturbation, true );

        double lightTimeSensitivity = -( upPerturbedLightTime - downPerturbedLightTime ) / ( 2.0 * timePerturbation );

        // Test numerical derivative against Doppler observable
        BOOST_CHECK_SMALL( std::fabs( lightTimeSensitivity - dopplerObservable ), 1.0E-14 );
    }

    // Test observation biases
    {
        // Create observation and bias settings
        std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
        biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( 1.0E-6 ), true ) );
        biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( 2.5E-4 ), false ) );
        std::shared_ptr< ObservationBiasSettings > biasSettings = std::make_shared< MultipleObservationBiasSettings >(
                    biasSettingsList );

        std::shared_ptr< ObservationSettings > biasedObservableSettings = std::make_shared< ObservationSettings >
                ( one_way_doppler, std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings );

        // Create observation model
        std::shared_ptr< ObservationModel< 1, double, double> > biasedObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEnds, biasedObservableSettings, bodyMap );

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double unbiasedObservation = biasedObservationModel->computeIdealObservations(
                    observationTime, receiver )( 0 );
        double biasedObservation = biasedObservationModel->computeObservations(
                    observationTime, receiver )( 0 );
        BOOST_CHECK_CLOSE_FRACTION( biasedObservation, 1.0E-6 + ( 1.0 + 2.5E-4 ) * unbiasedObservation, 1.0E-15 );

    }

    // Test proper time rates
    {
        // Define link ends for observations.
        LinkEnds linkEndsStationSpacecraft;
        linkEndsStationSpacecraft[ transmitter ] = std::make_pair( "Earth" , "Station1"  );
        linkEndsStationSpacecraft[ receiver ] = std::make_pair( "Spacecraft" , ""  );

        // Create observation settings
        std::shared_ptr< ObservationSettings > observableSettingsWithoutCorrections = std::make_shared< ObservationSettings >
                ( one_way_doppler );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithoutCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEndsStationSpacecraft, observableSettingsWithoutCorrections, bodyMap );

        // Create observation settings
        std::shared_ptr< ObservationSettings > observableSettingsWithCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  std::shared_ptr< LightTimeCorrectionSettings >( ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEndsStationSpacecraft, observableSettingsWithCorrections, bodyMap );

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double observationWithoutCorrections = observationModelWithoutCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );
        double observationWithCorrections = observationModelWithCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );

        std::shared_ptr< RotationalEphemeris > earthRotationModel =
                bodyMap.at( "Earth" )->getRotationalEphemeris( );
        Eigen::Vector3d groundStationVelocityVector =
                earthRotationModel->getDerivativeOfRotationToTargetFrame( observationTime ) *
                ( earthRotationModel->getRotationToBaseFrame( observationTime ) * stationCartesianPosition );

        std::shared_ptr< Ephemeris > spacecraftEphemeris =
                bodyMap.at( "Spacecraft" )->getEphemeris( );
        Eigen::Vector6d spacecraftState =
                spacecraftEphemeris->getCartesianState( observationTime );

        long double groundStationProperTimeRate = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( groundStationVelocityVector.norm( ), 2 ) +
                      earthGravitationalParameter / stationCartesianPosition.norm( ) ) );
        long double spacecraftProperTimeRate = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( spacecraftState.segment( 3, 3 ).norm( ), 2 ) +
                      earthGravitationalParameter / spacecraftState.segment( 0, 3 ).norm( ) ) );

        long double manualDopplerValue =
                groundStationProperTimeRate *
                ( 1.0L + static_cast< long double >( observationWithoutCorrections ) ) /
                spacecraftProperTimeRate - 1.0L;

        BOOST_CHECK_SMALL( std::fabs( static_cast< double >( manualDopplerValue ) - observationWithCorrections ),
                           static_cast< double >( std::numeric_limits< long double >::epsilon( ) ) );
    }
}


BOOST_AUTO_TEST_CASE( testTwoWayDoppplerModel )
{
    // Load Spice kernels
    std::string kernelsPath = input_output::getSpiceKernelPath( );
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
    std::map< std::string, std::shared_ptr< BodySettings > > defaultBodySettings =
            getDefaultBodySettings(
                bodiesToCreate, initialEphemerisTime - buffer, finalEphemerisTime + buffer );

    // Create bodies
    NamedBodyMap bodyMap = createBodies( defaultBodySettings );

    // Create ground stations
    const Eigen::Vector3d stationCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    createGroundStation( bodyMap.at( "Earth" ), "Station1", stationCartesianPosition, cartesian_position );

    // Set station with unrealistic position to force stronger proper time effect
    const Eigen::Vector3d stationCartesianPosition2( 4324532.0, 157372.0, -9292843.0 );
    createGroundStation( bodyMap.at( "Earth" ), "Station2", stationCartesianPosition2, cartesian_position );

    // Create Spacecraft
    Eigen::Vector6d spacecraftOrbitalElements;
    spacecraftOrbitalElements( semiMajorAxisIndex ) = 10000.0E3;
    spacecraftOrbitalElements( eccentricityIndex ) = 0.33;
    spacecraftOrbitalElements( inclinationIndex ) = convertDegreesToRadians( 65.3 );
    spacecraftOrbitalElements( argumentOfPeriapsisIndex )
            = convertDegreesToRadians( 235.7 );
    spacecraftOrbitalElements( longitudeOfAscendingNodeIndex )
            = convertDegreesToRadians( 23.4 );
    spacecraftOrbitalElements( trueAnomalyIndex ) = convertDegreesToRadians( 0.0 );
    double earthGravitationalParameter = bodyMap.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodyMap[ "Spacecraft" ] = std::make_shared< Body >( );
    bodyMap[ "Spacecraft" ]->setEphemeris(
                createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                         spacecraftOrbitalElements, 0.0, earthGravitationalParameter, "Earth" ), "Spacecraft" ) );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    {
        // Define link ends for observations.
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
        linkEnds[ reflector1 ] = std::make_pair( "Mars" , ""  );
        linkEnds[ receiver ] = std::make_pair( "Earth" , ""  );


        LinkEnds uplinkLinkEnds;
        uplinkLinkEnds[ transmitter ] = std::make_pair( "Earth" , ""  );
        uplinkLinkEnds[ receiver ] = std::make_pair( "Mars" , ""  );

        LinkEnds downlinkLinkEnds;
        downlinkLinkEnds[ transmitter ] = std::make_pair( "Mars" , ""  );
        downlinkLinkEnds[ receiver ] = std::make_pair( "Earth" , ""  );

        // Create observation settings
        std::shared_ptr< TwoWayDopplerObservationModel< double, double> > twoWayDopplerObservationModel =
                std::dynamic_pointer_cast< TwoWayDopplerObservationModel< double, double> >(
                    ObservationModelCreator< 1, double, double>::createObservationModel(
                        linkEnds, std::make_shared< ObservationSettings >( two_way_doppler ), bodyMap ) );
        std::shared_ptr< ObservationModel< 1, double, double > > twoWayRangeObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEnds, std::make_shared< ObservationSettings >( n_way_range ), bodyMap );

        std::shared_ptr< ObservationModel< 1, double, double > > uplinkDopplerObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    uplinkLinkEnds, std::make_shared< ObservationSettings >( one_way_doppler ), bodyMap );
        std::shared_ptr< ObservationModel< 1, double, double > > downlinkDopplerObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    downlinkLinkEnds, std::make_shared< ObservationSettings >( one_way_doppler ), bodyMap );


        // Creare independent light time calculator objects
        std::shared_ptr< LightTimeCalculator< double, double > > uplinkLightTimeCalculator =
                createLightTimeCalculator( linkEnds[ transmitter ], linkEnds[ reflector1 ], bodyMap );
        std::shared_ptr< LightTimeCalculator< double, double > > downlinkLightTimeCalculator =
                createLightTimeCalculator( linkEnds[ reflector1 ], linkEnds[ receiver ], bodyMap );

        // Test observable for both fixed link ends
        for( unsigned testCase = 0; testCase < 3; testCase++ )
        {

            double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
            std::vector< double > linkEndTimes;
            std::vector< Eigen::Vector6d > linkEndStates;

            std::vector< double > rangeLinkEndTimes;
            std::vector< Eigen::Vector6d > rangeLinkEndStates;


            // Define link end
            LinkEndType referenceLinkEnd, uplinkReferenceLinkEnd, downlinkReferenceLinkEnd;
            int transmitterReferenceTimeIndex, receiverReferenceTimeIndex;
            if( testCase == 0 )
            {
                referenceLinkEnd = transmitter;
                uplinkReferenceLinkEnd =  transmitter;
                downlinkReferenceLinkEnd = transmitter;
                transmitterReferenceTimeIndex = 0;
                receiverReferenceTimeIndex = 2;
            }
            else if( testCase == 1 )
            {
                referenceLinkEnd = reflector1;
                uplinkReferenceLinkEnd =  receiver;
                downlinkReferenceLinkEnd = transmitter;
                transmitterReferenceTimeIndex = 1;
                receiverReferenceTimeIndex = 2;
            }
            else
            {
                referenceLinkEnd = receiver;
                uplinkReferenceLinkEnd =  receiver;
                downlinkReferenceLinkEnd = receiver;
                transmitterReferenceTimeIndex = 1;
                receiverReferenceTimeIndex = 3;
            }

            // Compute observables
            double dopplerObservable = twoWayDopplerObservationModel->computeObservationsWithLinkEndData(
                        observationTime, referenceLinkEnd, linkEndTimes, linkEndStates )( 0 );
            double uplinkDopplerObservable = uplinkDopplerObservationModel->computeObservations(
                        linkEndTimes.at( transmitterReferenceTimeIndex ), uplinkReferenceLinkEnd )( 0 );
            double downlinkDopplerObservable = downlinkDopplerObservationModel->computeObservations(
                        linkEndTimes.at( receiverReferenceTimeIndex ), downlinkReferenceLinkEnd  )( 0 );
            twoWayRangeObservationModel->computeObservationsWithLinkEndData(
                        observationTime, referenceLinkEnd, rangeLinkEndTimes, rangeLinkEndStates )( 0 );


            // Compare light time calculator link end conditions with observation model
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeLinkEndStates.at( 3 ), linkEndStates.at( 3 ), 1.0E-15 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeLinkEndStates.at( 2 ), linkEndStates.at( 2 ), 1.0E-15 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeLinkEndStates.at( 1 ), linkEndStates.at( 1 ), 1.0E-15 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( rangeLinkEndStates.at( 0 ), linkEndStates.at( 0 ), 1.0E-15 );

                BOOST_CHECK_SMALL( std::fabs( rangeLinkEndTimes.at( 3 ) - linkEndTimes.at( 3 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL( std::fabs( rangeLinkEndTimes.at( 2 ) - linkEndTimes.at( 2 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL( std::fabs( rangeLinkEndTimes.at( 1 ) - linkEndTimes.at( 1 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL( std::fabs( rangeLinkEndTimes.at( 0 ) - linkEndTimes.at( 0 ) ), 1.0E-15 );
            }

            // Compute numerical partial derivative of light time.
            double timePerturbation = 100.0;
            double upPerturbedLightTime =
                    uplinkLightTimeCalculator->calculateLightTime( linkEndTimes.at( 1 ) + timePerturbation, true );
            double downPerturbedLightTime =
                    uplinkLightTimeCalculator->calculateLightTime( linkEndTimes.at( 1 ) - timePerturbation, true );

            double uplinkLightTimeSensitivity = -( upPerturbedLightTime - downPerturbedLightTime ) / ( 2.0 * timePerturbation );

            upPerturbedLightTime = downlinkLightTimeCalculator->calculateLightTime( linkEndTimes.at( 3 ) + timePerturbation, true );
            downPerturbedLightTime = downlinkLightTimeCalculator->calculateLightTime( linkEndTimes.at( 3 ) - timePerturbation, true );

            double downlinkLightTimeSensitivity = -( upPerturbedLightTime - downPerturbedLightTime ) / ( 2.0 * timePerturbation );

            // Test numerical derivative against Doppler observable
            BOOST_CHECK_SMALL( std::fabs( uplinkLightTimeSensitivity + downlinkLightTimeSensitivity +
                                          downlinkLightTimeSensitivity * uplinkLightTimeSensitivity - dopplerObservable ), 1.0E-14 );
            BOOST_CHECK_SMALL( std::fabs( ( uplinkDopplerObservable + 1 ) * ( downlinkDopplerObservable + 1 ) -
                                          ( dopplerObservable + 1 ) ), std::numeric_limits< double >::epsilon( ) );

        }
    }

    // Test proper time rates in two-way link where effects should cancel (no retransmission delays; transmitter and receiver are
    // same station)
    for( unsigned test = 0; test < 2; test++ )
    {
        std::string receivingStation;
        Eigen::Vector3d receivingStationPosition;

        if( test == 0 )
        {
            receivingStation = "Station1";
            receivingStationPosition = stationCartesianPosition;
        }
        else
        {
            receivingStation = "Station2";
            receivingStationPosition = stationCartesianPosition2;
        }

        // Define link ends for observations.
        LinkEnds linkEndsStationSpacecraft;
        linkEndsStationSpacecraft[ transmitter ] = std::make_pair( "Earth" , "Station1"  );
        linkEndsStationSpacecraft[ reflector1 ] = std::make_pair( "Spacecraft" , ""  );
        linkEndsStationSpacecraft[ receiver ] = std::make_pair( "Earth" , receivingStation  );

        LinkEnds uplinkLinkEndsStationSpacecraft;
        uplinkLinkEndsStationSpacecraft[ transmitter ] = std::make_pair( "Earth" , "Station1"  );
        uplinkLinkEndsStationSpacecraft[ receiver ] = std::make_pair( "Spacecraft" , ""  );

        LinkEnds downlinkLinkEndsStationSpacecraft;
        downlinkLinkEndsStationSpacecraft[ receiver ] = std::make_pair( "Earth" , receivingStation  );
        downlinkLinkEndsStationSpacecraft[ transmitter ] = std::make_pair( "Spacecraft" , ""  );

        // Create observation settings
        std::shared_ptr< ObservationSettings > observableSettingsWithoutCorrections = std::make_shared< ObservationSettings >
                ( two_way_doppler );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithoutCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEndsStationSpacecraft, observableSettingsWithoutCorrections, bodyMap );

        // Create observation settings
        std::shared_ptr< OneWayDopplerObservationSettings > oneWayObservableSettingsWithCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  std::shared_ptr< LightTimeCorrectionSettings >( ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );

        std::shared_ptr< ObservationSettings > twoWayObservableSettingsWithCorrections =
                std::make_shared< TwoWayDopplerObservationSettings >
                ( oneWayObservableSettingsWithCorrections, oneWayObservableSettingsWithCorrections );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEndsStationSpacecraft, twoWayObservableSettingsWithCorrections, bodyMap );

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double observationWithoutCorrections = observationModelWithoutCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );

        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;

        double observationWithCorrections = observationModelWithCorrections->computeIdealObservationsWithLinkEndData(
                    observationTime, receiver, linkEndTimes, linkEndStates ).x( );

        std::shared_ptr< RotationalEphemeris > earthRotationModel =
                bodyMap.at( "Earth" )->getRotationalEphemeris( );

        Eigen::Vector3d groundStationVelocityVectorAtTransmission =
                earthRotationModel->getDerivativeOfRotationToTargetFrame( linkEndTimes.at( 0 ) ) *
                ( earthRotationModel->getRotationToBaseFrame( linkEndTimes.at( 0 ) ) * stationCartesianPosition );

        Eigen::Vector3d groundStationVelocityVectorAtReception =
                earthRotationModel->getDerivativeOfRotationToTargetFrame( linkEndTimes.at( 2 ) ) *
                ( earthRotationModel->getRotationToBaseFrame( linkEndTimes.at( 2 ) ) * receivingStationPosition );



        long double groundStationProperTimeRateAtTransmission = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( groundStationVelocityVectorAtTransmission.norm( ), 2 ) +
                      earthGravitationalParameter / stationCartesianPosition.norm( ) ) );

        long double groundStationProperTimeRateAtReception = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( groundStationVelocityVectorAtReception.norm( ), 2 ) +
                      earthGravitationalParameter / receivingStationPosition.norm( ) ) );

        if( test == 0 )
        {
            BOOST_CHECK_SMALL( std::fabs( observationWithCorrections - observationWithoutCorrections ),
                               static_cast< double >( std::numeric_limits< long double >::epsilon( ) ) );
            BOOST_CHECK_SMALL( std::fabs( static_cast< double >(
                                              groundStationProperTimeRateAtTransmission - groundStationProperTimeRateAtReception ) ),
                               static_cast< double >( std::numeric_limits< long double >::epsilon( ) ) );
        }
        else

        {
            long double properTimeRatioDeviation =
                    groundStationProperTimeRateAtTransmission / groundStationProperTimeRateAtReception - 1.0L;
            long double observableDifference =
                    observationWithCorrections - ( observationWithoutCorrections + properTimeRatioDeviation +
                                                   observationWithoutCorrections * properTimeRatioDeviation );
            BOOST_CHECK_SMALL( std::fabs( static_cast< double >( observableDifference ) ),
                               static_cast< double >( std::numeric_limits< long double >::epsilon( ) ) );
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}


