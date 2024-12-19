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
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/astro/basic_astro/unitConversions.h"

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
    std::string kernelsPath = paths::getSpiceKernelPath( );
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
            getDefaultBodySettings(
                bodiesToCreate, "SSB", "ECLIPJ2000" );
    defaultBodySettings.get( "Sun" )->gravityFieldSettings = centralGravitySettings( spice_interface::getBodyGravitationalParameter( "Sun") * 10000.0 );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create ground station
    const Eigen::Vector3d stationCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    createGroundStation( bodies.at( "Earth" ), "Station1", stationCartesianPosition, cartesian_position );

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
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodies.createEmptyBody( "Spacecraft" );;
    bodies.at( "Spacecraft" )->setEphemeris(
                createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                         spacecraftOrbitalElements, 0.0, earthGravitationalParameter, "Earth", "ECLIPJ2000" ), "Spacecraft" ) );
    bodies.processBodyFrameDefinitions( );


    // Define link ends for observations.
    LinkEnds linkEnds;
    linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , ""  );
    linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , ""  );

    // Create observation settings
    for( int useCorrections = 0; useCorrections < 2; useCorrections++ )
    {
        double toleranceScaling = 1.0;
        std::vector< std::shared_ptr< LightTimeCorrectionSettings > > correctionSettings;
        if( useCorrections == 1 )
        {
            correctionSettings.push_back( std::make_shared< FirstOrderRelativisticLightTimeCorrectionSettings >( std::vector< std::string >( { "Sun" } ) ) );
            toleranceScaling *= 100.0;
        }
        std::shared_ptr<ObservationModelSettings> observableSettings = std::make_shared<ObservationModelSettings>
            ( one_way_doppler, linkEnds, correctionSettings );

        // Create observation model.
        std::shared_ptr<ObservationModel<1, double, double> > observationModel =
            ObservationModelCreator<1, double, double>::createObservationModel(
                observableSettings, bodies );

        std::shared_ptr<OneWayDopplerObservationModel<double, double> > dopplerObservationModel =
            std::dynamic_pointer_cast<OneWayDopplerObservationModel<double, double> >( observationModel );

        // Test observable for both fixed link ends
        for ( unsigned testCase = 2; testCase < 4; testCase++ )
        {
            std::cout<<testCase<<" "<<useCorrections<<std::endl;
            double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;
            std::vector<double> linkEndTimes;
            std::vector<Eigen::Vector6d> linkEndStates;

            // Define link end
            LinkEndType referenceLinkEnd;
            if ( testCase == 0 || testCase == 2 )
            {
                referenceLinkEnd = transmitter;
            }
            else
            {
                referenceLinkEnd = receiver;
            }

            double scalingTerm = physical_constants::SPEED_OF_LIGHT;
            bool useNormalization = false;
            if ( testCase > 1 )
            {
                scalingTerm = 1.0;
                useNormalization = true;
                dopplerObservationModel->setNormalizeWithSpeedOfLight( useNormalization );
            }

            // Compute observable
            double dopplerObservable = observationModel->computeObservationsWithLinkEndData(
                observationTime, referenceLinkEnd, linkEndTimes, linkEndStates )( 0 );

            // Creare independent light time calculator object
            std::shared_ptr<LightTimeCalculator<double, double> > lightTimeCalculator =
                createLightTimeCalculator( linkEnds, transmitter, receiver, bodies, undefined_observation_model, correctionSettings );
            Eigen::Vector6d transmitterState, receiverState;
            // Compute light time
            double lightTime = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                receiverState, transmitterState, observationTime, ( !( testCase == 0 || testCase == 2 )));

            // Compare light time calculator link end conditions with observation model
            {
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( receiverState, linkEndStates.at( 1 ), 1.0E-15 );
                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transmitterState, linkEndStates.at( 0 ), 1.0E-15 );

                if ( testCase == 0 || testCase == 2 )
                {
                    BOOST_CHECK_SMALL( std::fabs( observationTime - linkEndTimes.at( 0 )), 1.0E-12 );
                    BOOST_CHECK_SMALL( std::fabs( observationTime + lightTime - linkEndTimes.at( 1 )), 1.0E-10 );
                }
                else
                {
                    BOOST_CHECK_SMALL( std::fabs( observationTime - linkEndTimes.at( 1 )), 1.0E-12 );
                    BOOST_CHECK_SMALL( std::fabs( observationTime - lightTime - linkEndTimes.at( 0 )), 1.0E-10 );
                }
            }

            // Compute numerical partial derivative of light time.
            double timePerturbation = 100.0;
            double upPerturbedLightTime =
                lightTimeCalculator->calculateLightTime( linkEndTimes.at( 1 ) + timePerturbation, true );
            double downPerturbedLightTime =
                lightTimeCalculator->calculateLightTime( linkEndTimes.at( 1 ) - timePerturbation, true );

            double lightTimeSensitivity = -( upPerturbedLightTime - downPerturbedLightTime ) / ( 2.0 * timePerturbation );

            // Test numerical derivative against Doppler observable
            BOOST_CHECK_CLOSE_FRACTION( scalingTerm * lightTimeSensitivity, dopplerObservable, 1.0E-8 * toleranceScaling );

            if( useCorrections && testCase == 3 )
            {
                std::shared_ptr< LightTimeCorrection > correction = lightTimeCalculator->getLightTimeCorrection( ).at( 0 );
                double nominalLightTimeCorrection = correction->calculateLightTimeCorrection(
                    transmitterState, receiverState, linkEndTimes.at( 0 ), linkEndTimes.at( 1 ) );

                Eigen::Vector6d transmitterStateUp, receiverStateUp;
                // Compute light time
                double lightTimeUp = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                    receiverStateUp, transmitterStateUp, observationTime + timePerturbation, true );
                double lightTimeCorrectionUp = correction->calculateLightTimeCorrection(
                    transmitterStateUp, receiverStateUp, observationTime + timePerturbation - lightTimeUp, observationTime + timePerturbation );

                Eigen::Vector6d transmitterStateDown, receiverStateDown;
                // Compute light time
                double lightTimeDown = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                    receiverStateDown, transmitterStateDown, observationTime - timePerturbation, true );
                double lightTimeCorrectionDown = correction->calculateLightTimeCorrection(
                    transmitterStateDown, receiverStateDown, observationTime - timePerturbation - lightTimeUp, observationTime - timePerturbation );

                Eigen::Matrix< double, 3, 1 > lightTimeCorrectionWrtReceiver = correction->calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                    transmitterState, receiverState, linkEndTimes.at( 0 ), linkEndTimes.at( 1 ), receiver );
                Eigen::Matrix< double, 3, 1 > lightTimeCorrectionWrtTransmitter = correction->calculateLightTimeCorrectionPartialDerivativeWrtLinkEndPosition(
                    transmitterState, receiverState, linkEndTimes.at( 0 ), linkEndTimes.at( 1 ), transmitter );

                BOOST_CHECK_CLOSE_FRACTION( lightTimeCorrectionWrtReceiver.dot( receiverState.segment( 3, 3 ) ) + lightTimeCorrectionWrtTransmitter.dot( transmitterState.segment( 3, 3 ) ),
                                            ( lightTimeCorrectionUp  - lightTimeCorrectionDown ) / ( 2.0 * timePerturbation ), 1.0E-3 );

            }
        }
    }

    // Test observation biases
    {
        // Create observation and bias settings
        std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
        biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( 1.0E2 ), true ) );
        biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( 2.5E-4 ), false ) );
        std::shared_ptr< ObservationBiasSettings > biasSettings = std::make_shared< MultipleObservationBiasSettings >(
                    biasSettingsList );

        std::shared_ptr< ObservationModelSettings > biasedObservableSettings = std::make_shared< ObservationModelSettings >
                ( one_way_doppler, linkEnds, std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings );

        // Create observation model
        std::shared_ptr< ObservationModel< 1, double, double> > biasedObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    biasedObservableSettings, bodies );

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double unbiasedObservation = biasedObservationModel->computeIdealObservations(
                    observationTime, receiver )( 0 );
        double biasedObservation = biasedObservationModel->computeObservations(
                    observationTime, receiver )( 0 );
        BOOST_CHECK_CLOSE_FRACTION( biasedObservation, 1.0E2 + ( 1.0 + 2.5E-4 ) * unbiasedObservation, 1.0E-15 );

    }

    // Test proper time rates
    {
        // Define link ends for observations.
        LinkEnds linkEndsStationSpacecraft;
        linkEndsStationSpacecraft[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "Station1"  );
        linkEndsStationSpacecraft[ receiver ] = std::make_pair< std::string, std::string >( "Spacecraft" , ""  );

        // Create observation settings
        std::shared_ptr< ObservationModelSettings > observableSettingsWithoutCorrections = std::make_shared< ObservationModelSettings >
                ( one_way_doppler, linkEndsStationSpacecraft );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithoutCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    observableSettingsWithoutCorrections, bodies );

        // Create observation settings
        std::shared_ptr< ObservationModelSettings > observableSettingsWithCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  linkEndsStationSpacecraft,
                   std::shared_ptr< LightTimeCorrectionSettings >( ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    observableSettingsWithCorrections, bodies );

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double observationWithoutCorrections = observationModelWithoutCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );
        double observationWithCorrections = observationModelWithCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );

        std::shared_ptr< RotationalEphemeris > earthRotationModel =
                bodies.at( "Earth" )->getRotationalEphemeris( );
        Eigen::Vector6d groundStationEarthFixedState = Eigen::Vector6d::Zero( );
        groundStationEarthFixedState.segment( 0, 3 ) = stationCartesianPosition;
        Eigen::Vector6d groundStationGeocentricState = ephemerides::transformStateToInertialOrientation(
            groundStationEarthFixedState, observationTime, earthRotationModel );

        Eigen::Vector6d groundStationState = groundStationGeocentricState;
        groundStationState += bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( observationTime );

        Eigen::Vector6d spacecraftState = bodies.at( "Spacecraft" )->getStateInBaseFrameFromEphemeris( observationTime );
        Eigen::Vector6d spacecraftGeocentricState = bodies.at( "Spacecraft" )->getStateInBaseFrameFromEphemeris( observationTime ) -
            bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( observationTime );

        long double groundStationProperTimeRate = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( groundStationState.segment( 3, 3 ).norm( ), 2 ) +
                      earthGravitationalParameter / groundStationGeocentricState.segment( 0, 3 ).norm( ) ) );
        long double spacecraftProperTimeRate = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( spacecraftState.segment( 3, 3 ).norm( ), 2 ) +
                      earthGravitationalParameter / spacecraftGeocentricState.segment( 0, 3 ).norm( ) ) );

        long double manualDopplerValue =
                ( groundStationProperTimeRate *
                  ( 1.0L + static_cast< long double >( observationWithoutCorrections ) / physical_constants::SPEED_OF_LIGHT ) /
                  spacecraftProperTimeRate - 1.0L ) * physical_constants::SPEED_OF_LIGHT;

        BOOST_CHECK_SMALL( std::fabs( static_cast< double >( manualDopplerValue ) - observationWithCorrections ), 1.0E-6 );
    }
}


BOOST_AUTO_TEST_CASE( testTwoWayDoppplerModel )
{
    // Load Spice kernels
    std::string kernelsPath = paths::getSpiceKernelPath( );
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
            getDefaultBodySettings(
                bodiesToCreate, "SSB" );

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies( defaultBodySettings );

    // Create ground stations
    const Eigen::Vector3d stationCartesianPosition( 1917032.190, 6029782.349, -801376.113 );
    createGroundStation( bodies.at( "Earth" ), "Station1", stationCartesianPosition, cartesian_position );

    // Set station with unrealistic position to force stronger proper time effect
    const Eigen::Vector3d stationCartesianPosition2( 4324532.0, 157372.0, -9292843.0 );
    createGroundStation( bodies.at( "Earth" ), "Station2", stationCartesianPosition2, cartesian_position );

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
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    bodies.createEmptyBody( "Spacecraft" );;
    bodies.at( "Spacecraft" )->setEphemeris(
                createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                         spacecraftOrbitalElements, 0.0, earthGravitationalParameter, "Earth" ), "Spacecraft" ) );
    bodies.processBodyFrameDefinitions( );



    {
        // Define link ends for observations.
        LinkEnds linkEnds;
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , ""  );
        linkEnds[ reflector1 ] = std::make_pair< std::string, std::string >( "Mars" , ""  );
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , ""  );


        LinkEnds uplinkLinkEnds;
        uplinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , ""  );
        uplinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Mars" , ""  );

        LinkEnds downlinkLinkEnds;
        downlinkLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Mars" , ""  );
        downlinkLinkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Earth" , ""  );

        // Create observation settings
        std::shared_ptr< TwoWayDopplerObservationModel< double, double> > twoWayDopplerObservationModel =
                std::dynamic_pointer_cast< TwoWayDopplerObservationModel< double, double> >(
                    ObservationModelCreator< 1, double, double>::createObservationModel(
                        std::make_shared< ObservationModelSettings >( two_way_doppler, linkEnds ), bodies ) );
        std::shared_ptr< ObservationModel< 1, double, double > > twoWayRangeObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    twoWayRangeSimple( linkEnds ), bodies );

        std::shared_ptr< ObservationModel< 1, double, double > > uplinkDopplerObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    std::make_shared< ObservationModelSettings >( one_way_doppler, uplinkLinkEnds ), bodies );
        std::shared_ptr< ObservationModel< 1, double, double > > downlinkDopplerObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    std::make_shared< ObservationModelSettings >( one_way_doppler, downlinkLinkEnds ), bodies );


        // Creare independent light time calculator objects
        std::shared_ptr< LightTimeCalculator< double, double > > uplinkLightTimeCalculator =
                createLightTimeCalculator( linkEnds, transmitter, reflector1, bodies );
        std::shared_ptr< LightTimeCalculator< double, double > > downlinkLightTimeCalculator =
                createLightTimeCalculator( linkEnds, reflector1, receiver, bodies );

        // Test observable for both fixed link ends
        for( unsigned testCase = 0; testCase < 3; testCase++ )
        {

            for( unsigned int normalizeObservable = 0; normalizeObservable < 2; normalizeObservable++ )
            {
                if( normalizeObservable == 0 )
                {
                    twoWayDopplerObservationModel->setNormalizeWithSpeedOfLight( 0 );
                    std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >(
                                uplinkDopplerObservationModel )->setNormalizeWithSpeedOfLight( 0 );
                    std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >(
                                downlinkDopplerObservationModel )->setNormalizeWithSpeedOfLight( 0 );
                }
                else
                {
                    twoWayDopplerObservationModel->setNormalizeWithSpeedOfLight( 1 );
                    std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >(
                                uplinkDopplerObservationModel )->setNormalizeWithSpeedOfLight( 1 );
                    std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >(
                                downlinkDopplerObservationModel )->setNormalizeWithSpeedOfLight( 1 );
                }
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

                double scalingTerm = normalizeObservable ? 1.0 : physical_constants::SPEED_OF_LIGHT;

                // Test numerical derivative against Doppler observable
                BOOST_CHECK_SMALL( std::fabs( ( uplinkLightTimeSensitivity + downlinkLightTimeSensitivity +
                                                downlinkLightTimeSensitivity * uplinkLightTimeSensitivity ) * scalingTerm - dopplerObservable ),
                                   scalingTerm * 5.0E-14 );
                BOOST_CHECK_SMALL( std::fabs( ( uplinkDopplerObservable / scalingTerm + 1 ) * ( downlinkDopplerObservable / scalingTerm + 1 ) -
                                              ( dopplerObservable  / scalingTerm + 1 ) ), std::numeric_limits< double >::epsilon( ) );
            }

        }
        twoWayDopplerObservationModel->setNormalizeWithSpeedOfLight( 0 );
        std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >(
                    uplinkDopplerObservationModel )->setNormalizeWithSpeedOfLight( 0 );
        std::dynamic_pointer_cast< OneWayDopplerObservationModel< double, double> >(
                    downlinkDopplerObservationModel )->setNormalizeWithSpeedOfLight( 0 );
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
        linkEndsStationSpacecraft[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "Station1"  );
        linkEndsStationSpacecraft[ reflector1 ] = std::make_pair< std::string, std::string >( "Spacecraft" , ""  );
        linkEndsStationSpacecraft[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth" , receivingStation  ) );

        LinkEnds uplinkLinkEndsStationSpacecraft;
        uplinkLinkEndsStationSpacecraft[ transmitter ] = std::make_pair< std::string, std::string >( "Earth" , "Station1"  );
        uplinkLinkEndsStationSpacecraft[ receiver ] = std::make_pair< std::string, std::string >( "Spacecraft" , ""  );

        LinkEnds downlinkLinkEndsStationSpacecraft;
        downlinkLinkEndsStationSpacecraft[ receiver ] = std::pair< std::string, std::string >( std::make_pair(  "Earth" , receivingStation ) );
        downlinkLinkEndsStationSpacecraft[ transmitter ] = std::make_pair< std::string, std::string >( "Spacecraft" , ""  );

        // Create observation settings
        std::shared_ptr< ObservationModelSettings > observableSettingsWithoutCorrections = std::make_shared< ObservationModelSettings >
                ( two_way_doppler, linkEndsStationSpacecraft );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithoutCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    observableSettingsWithoutCorrections, bodies );

        // Create observation settings
        std::shared_ptr< OneWayDopplerObservationSettings > oneWayObservableUplinkSettingsWithCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  uplinkLinkEndsStationSpacecraft, std::shared_ptr< LightTimeCorrectionSettings >( ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );
        std::shared_ptr< OneWayDopplerObservationSettings > oneWayObservableDownlinkSettingsWithCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  downlinkLinkEndsStationSpacecraft, std::shared_ptr< LightTimeCorrectionSettings >( ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );

        std::shared_ptr< OneWayDopplerObservationSettings > oneWayObservableUplinkSettingsWithoutCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  uplinkLinkEndsStationSpacecraft, nullptr );
        std::shared_ptr< OneWayDopplerObservationSettings > oneWayObservableDownlinkSettingsWithoutCorrections =
                std::make_shared< OneWayDopplerObservationSettings >
                (  downlinkLinkEndsStationSpacecraft, nullptr );

        std::shared_ptr< ObservationModelSettings > twoWayObservableSettingsWithCorrections =
                std::make_shared< TwoWayDopplerObservationSettings >
                ( oneWayObservableUplinkSettingsWithCorrections, oneWayObservableDownlinkSettingsWithCorrections );
        std::shared_ptr< ObservationModelSettings > twoWayObservableSettingsWithoutCorrections =
                std::make_shared< TwoWayDopplerObservationSettings >
                ( oneWayObservableUplinkSettingsWithoutCorrections, oneWayObservableDownlinkSettingsWithoutCorrections );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithCorrections =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    twoWayObservableSettingsWithCorrections, bodies );
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithoutCorrectionsDirect =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    twoWayObservableSettingsWithoutCorrections, bodies );


        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double observationWithoutCorrections = observationModelWithoutCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );
        double observationWithoutCorrectionsDirect = observationModelWithoutCorrectionsDirect->computeIdealObservations(
                    observationTime, receiver ).x( );
        BOOST_CHECK_SMALL( std::fabs( observationWithoutCorrections - observationWithoutCorrectionsDirect ),
                           static_cast< double >( std::numeric_limits< long double >::epsilon( ) ) );
        std::vector< double > linkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > linkEndStates;

        double observationWithCorrections = observationModelWithCorrections->computeIdealObservationsWithLinkEndData(
                    observationTime, receiver, linkEndTimes, linkEndStates ).x( );

        std::shared_ptr< RotationalEphemeris > earthRotationModel =
                bodies.at( "Earth" )->getRotationalEphemeris( );

        Eigen::Vector6d groundStationEarthFixedStateTransmission = Eigen::Vector6d::Zero( );
        groundStationEarthFixedStateTransmission.segment( 0, 3 ) = stationCartesianPosition;
        Eigen::Vector6d groundStationGeocentricStateTransmission = ephemerides::transformStateToInertialOrientation(
            groundStationEarthFixedStateTransmission, linkEndTimes.at( 0 ), earthRotationModel );
        Eigen::Vector6d groundStationStateTransmission = groundStationGeocentricStateTransmission;
        groundStationStateTransmission += bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 0 ) );


        Eigen::Vector6d groundStationEarthFixedStateReception = Eigen::Vector6d::Zero( );
        groundStationEarthFixedStateReception.segment( 0, 3 ) = receivingStationPosition;
        Eigen::Vector6d groundStationGeocentricStateReception = ephemerides::transformStateToInertialOrientation(
            groundStationEarthFixedStateReception, linkEndTimes.at( 2 ), earthRotationModel );
        Eigen::Vector6d groundStationStateReception = groundStationGeocentricStateReception;
        groundStationStateReception += bodies.at( "Earth" )->getStateInBaseFrameFromEphemeris( linkEndTimes.at( 2 ) );


        long double groundStationProperTimeRateAtTransmission = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( groundStationStateTransmission.segment( 3, 3 ).norm( ), 2 ) +
                      earthGravitationalParameter / groundStationGeocentricStateTransmission.segment( 0, 3 ).norm( ) ) );

        long double groundStationProperTimeRateAtReception = 1.0L - static_cast< long double >(
                    physical_constants::INVERSE_SQUARE_SPEED_OF_LIGHT *
                    ( 0.5 * std::pow( groundStationStateReception.segment( 3, 3 ).norm( ), 2 ) +
                      earthGravitationalParameter / groundStationGeocentricStateReception.segment( 0, 3 ).norm( ) ) );


        if( test == 0 )
        {
            BOOST_CHECK_SMALL( std::fabs( observationWithCorrections - observationWithoutCorrections ), 1.0E-6 );
            BOOST_CHECK_SMALL( std::fabs( static_cast< double >(
                                              groundStationProperTimeRateAtTransmission - groundStationProperTimeRateAtReception ) ),
                               10.0 * static_cast< double >( std::numeric_limits< double >::epsilon( ) ) );
        }
        else

        {
            long double properTimeRatioDeviation =
                    groundStationProperTimeRateAtTransmission / groundStationProperTimeRateAtReception - 1.0L;
            long double observableDifference =
                    observationWithCorrections / physical_constants::SPEED_OF_LIGHT -
                    ( observationWithoutCorrections / physical_constants::SPEED_OF_LIGHT + properTimeRatioDeviation +
                      observationWithoutCorrections / physical_constants::SPEED_OF_LIGHT * properTimeRatioDeviation );
            BOOST_CHECK_SMALL( std::fabs( static_cast< double >( observableDifference ) ), 1.0E-6 );
        }


    }
}


BOOST_AUTO_TEST_SUITE_END( )

}

}


