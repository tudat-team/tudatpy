/*    Copyright (c) 2010-2018, Delft University of Technology
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


BOOST_AUTO_TEST_SUITE( test_one_way_doppler_model )


BOOST_AUTO_TEST_CASE( testOneWayDoppplerModel )
{
    // Load Spice kernels
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
    SystemOfBodies bodies = createBodies( defaultBodySettings );

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

    bodies[ "Spacecraft" ] = std::make_shared< Body >( );
    bodies[ "Spacecraft" ]->setEphemeris(
                createBodyEphemeris( std::make_shared< KeplerEphemerisSettings >(
                                         spacecraftOrbitalElements, 0.0, earthGravitationalParameter, "Earth" ), "Spacecraft" ) );

    setGlobalFrameBodyEphemerides( bodies, "SSB", "ECLIPJ2000" );

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
                linkEnds, observableSettings, bodies );

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
                createLightTimeCalculator( linkEnds, transmitter, receiver, bodies );
        Eigen::Vector6d transmitterState, receiverState;        
        // Compute light time
        double lightTime = lightTimeCalculator->calculateLightTimeWithLinkEndsStates(
                    receiverState, transmitterState, observationTime, testCase, ancilliarySetings );

        // Compare light time calculator link end conditions with observation model
        {
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( receiverState, linkEndStates.at( 1 ), 1.0E-15 );
            TUDAT_CHECK_MATRIX_CLOSE_FRACTION( transmitterState, linkEndStates.at( 0 ), 1.0E-15 );

            if( testCase == 0 )
            {
                BOOST_CHECK_SMALL( std::fabs( observationTime  - linkEndTimes.at( 0 ) ), 1.0E-15 );
                BOOST_CHECK_SMALL( std::fabs( observationTime + lightTime - linkEndTimes.at( 1 ) ), 1.0E-10 );
            }
            else
            {
                BOOST_CHECK_SMALL( std::fabs( observationTime - linkEndTimes.at( 1 ) ), 1.0E-15 );
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
        biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d( 1.0E-6 ) ) );
        biasSettingsList.push_back( std::make_shared< ConstantRelativeObservationBiasSettings >( Eigen::Vector1d( 2.5E-4 ) ) );
        std::shared_ptr< ObservationBiasSettings > biasSettings = std::make_shared< MultipleObservationBiasSettings >(
                    biasSettingsList );

        std::shared_ptr< ObservationSettings > biasedObservableSettings = std::make_shared< ObservationSettings >
                ( one_way_doppler, std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings );

        // Create observation model
        std::shared_ptr< ObservationModel< 1, double, double> > biasedObservationModel =
                ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEnds, biasedObservableSettings, bodies );

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
                    linkEndsStationSpacecraft, observableSettingsWithoutCorrections, bodies );

        // Create observation settings
        std::shared_ptr< ObservationSettings > observableSettingsWithCorrections =
                std::make_shared< OneWayDopperObservationSettings >
                (  std::shared_ptr< LightTimeCorrectionSettings >( ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ),
                   std::make_shared< DirectFirstOrderDopplerProperTimeRateSettings >( "Earth" ) );

        // Create observation model.
        std::shared_ptr< ObservationModel< 1, double, double> > observationModelWithCorrections =
               ObservationModelCreator< 1, double, double>::createObservationModel(
                    linkEndsStationSpacecraft, observableSettingsWithCorrections, bodies );

        double observationTime = ( finalEphemerisTime + initialEphemerisTime ) / 2.0;

        double observationWithoutCorrections = observationModelWithoutCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );
        double observationWithCorrections = observationModelWithCorrections->computeIdealObservations(
                    observationTime, receiver ).x( );

        std::shared_ptr< RotationalEphemeris > earthRotationModel =
                bodies.at( "Earth" )->getRotationalEphemeris( );
        Eigen::Vector3d groundStationVelocityVector =
                earthRotationModel->getDerivativeOfRotationToTargetFrame( observationTime ) *
                ( earthRotationModel->getRotationToBaseFrame( observationTime ) * stationCartesianPosition );

        std::shared_ptr< Ephemeris > spacecraftEphemeris =
                bodies.at( "Spacecraft" )->getEphemeris( );
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

BOOST_AUTO_TEST_SUITE_END( )

}

}


