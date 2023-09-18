/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <limits>
#include <string>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/estimation.h"
#include "tudat/simulation/estimation_setup.h"

#include "tudat/io/readOdfFile.h"
#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/io/readTabulatedWeatherData.h"
#include "tudat/simulation/estimation_setup/processOdfFile.h"

#include <boost/date_time/gregorian/gregorian.hpp>

#include "tudat/astro/ground_stations/transmittingFrequencies.h"

using namespace tudat::propagators;
using namespace tudat::estimatable_parameters;
using namespace tudat::spice_interface;
using namespace tudat::ephemerides;
using namespace tudat::input_output;
using namespace tudat::observation_models;
using namespace tudat::simulation_setup;
using namespace tudat::numerical_integrators;
using namespace tudat::basic_astrodynamics;
using namespace tudat;

//template< typename StateScalarType, TimeType >
SystemOfBodies createEnvironment(
    const double initialTimeEnvironment,
    const double finalTimeEnvironment,
    const std::string globalFrameOrigin,
    const std::string moonEphemerisOrigin,
    const bool useInterpolatedEphemerides,
    const double planetEphemerisInterpolationStep,
    const double spacecraftEphemerisInterpolationStep )
{
    // Create settings for default bodies
    std::vector< std::string > bodiesToCreate = { "Earth", "Sun", "Moon" };
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings;
    if( useInterpolatedEphemerides )
    {
        bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialTimeEnvironment, finalTimeEnvironment, globalFrameOrigin, globalFrameOrientation, planetEphemerisInterpolationStep );
    }
    else
    {
        bodySettings = getDefaultBodySettings(
            bodiesToCreate, globalFrameOrigin, globalFrameOrientation );
    }

    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( moonEphemerisOrigin );

    // Add spacecraft settings
    std::string spacecraftName = "GRAIL-A";
    std::string spacecraftCentralBody = "Moon";
    bodySettings.addSettings( spacecraftName );
    if( useInterpolatedEphemerides )
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared<InterpolatedSpiceEphemerisSettings>(
                initialTimeEnvironment, finalTimeEnvironment, spacecraftEphemerisInterpolationStep, spacecraftCentralBody, globalFrameOrientation );
    }
    else
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared<DirectSpiceEphemerisSettings>(
                spacecraftCentralBody, globalFrameOrientation, false, false, false );
    }

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies< long double, Time >( bodySettings );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), coordinate_conversions::geodetic_position );
    return bodies;

}

std::shared_ptr< LightTimeCalculator< long double, Time > > getOneWayLightTimeCalculator(
    const SystemOfBodies& bodies,
    const std::vector< std::shared_ptr< ObservationModelSettings > > oneWayObservationSettingsList )
{
    // Create observation simulators
    std::shared_ptr< ObservationSimulator< 1, long double, Time > > observationSimulator =
        std::dynamic_pointer_cast< ObservationSimulator< 1, long double, Time > >(
            createObservationSimulators< long double, Time >( oneWayObservationSettingsList, bodies ).at( 0 ) );
    std::shared_ptr< OneWayRangeObservationModel< long double, Time > > oneWayRangeObservationModel =
        std::dynamic_pointer_cast< OneWayRangeObservationModel< long double, Time > >(
            observationSimulator->getObservationModels( ).begin( )->second );
    return oneWayRangeObservationModel->getLightTimeCalculator( );
}


int main( )
{
    double initialTimeEnvironment = Time(107561, 2262.19) - 2.0 * 3600.0;
    double finalTimeEnvironment = Time(108258, 2771.19) + 2.0 * 3600.0;

    // Load spice kernels
    spice_interface::loadStandardSpiceKernels( );
    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_120301_120529_sci_v02.bsp" );

    SystemOfBodies earthCenteredBodies = createEnvironment(
        initialTimeEnvironment, finalTimeEnvironment,
        "Earth", "Earth", false, TUDAT_NAN, TUDAT_NAN );
    SystemOfBodies earthCenteredInterpolatedBodies = createEnvironment(
        initialTimeEnvironment, finalTimeEnvironment,
        "Earth", "Earth", true, 120.0, 10.0 );
    SystemOfBodies barycentricBodies = createEnvironment(
        initialTimeEnvironment, finalTimeEnvironment,
        "SSB", "Earth", false, TUDAT_NAN, TUDAT_NAN );
    SystemOfBodies barycentricInterpolatedBodies = createEnvironment(
        initialTimeEnvironment, finalTimeEnvironment,
        "SSB", "Earth", true, 120.0, 10.0 );
    SystemOfBodies barycentricMoonInterpolatedBodies = createEnvironment(
        initialTimeEnvironment, finalTimeEnvironment,
        "SSB", "SSB", true, 120.0, 10.0 );

    // Define link ends.
    LinkEnds testLinkEnds;
    testLinkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", "" ) );
    testLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "GRAIL-A", "" );

    Time testTime = initialTimeEnvironment + 86400.1;

//    std::cout<<std::setprecision( 20 );
//    Eigen::Vector6ld barycentricMoonState = barycentricInterpolatedBodies.at( "Moon" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime );
//
//    Eigen::Vector6ld manualBarycentricMoonState =
//        barycentricInterpolatedBodies.at( "Moon" )->getEphemeris( )->getTemplatedStateFromEphemeris< long double, Time >( testTime ) +
//        barycentricInterpolatedBodies.at( "Earth" )->getEphemeris( )->getTemplatedStateFromEphemeris< long double, Time >( testTime );
//
//    Eigen::Vector6ld directEarthState = barycentricInterpolatedBodies.at( "Earth" )->getEphemeris( )->getTemplatedStateFromEphemeris< long double, Time >( testTime );
//
//    std::cout<<"Direct Earth: "<<directEarthState.transpose( )<<std::endl;
//
//    std::cout<<barycentricMoonState.transpose( )<<std::endl;
//    std::cout<<( manualBarycentricMoonState - barycentricMoonState ).transpose( )<<std::endl<<std::endl<<std::endl;


//    {
//        std::cout<<"Barycentric planets (Moon)"<<std::endl;
//        for ( int i = -20; i <= 0; i++ )
//        {
//            long double perturbationStep = std::pow( 10, i );
//            Eigen::Vector6ld upperturbedState =
//                barycentricInterpolatedBodies.at( "Moon" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime + perturbationStep );
//            Eigen::Vector6ld downperturbedState =
//                barycentricInterpolatedBodies.at( "Moon" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime - perturbationStep );
//
//            std::cout << i << " " << ( upperturbedState - downperturbedState ).transpose( )  / ( 2.0 * perturbationStep )<< std::endl;
//        }
//    }





    // Create observation settings
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, testLinkEnds ) );

    std::shared_ptr< LightTimeCalculator< long double, Time > > earthCenteredLightTimeCalculator =
        getOneWayLightTimeCalculator( earthCenteredBodies, observationSettingsList );
    std::shared_ptr< LightTimeCalculator< long double, Time > > barycentricLightTimeCalculator =
        getOneWayLightTimeCalculator( barycentricBodies, observationSettingsList );
    std::shared_ptr< LightTimeCalculator< long double, Time > > barycentricInterpolatedLightTimeCalculator =
        getOneWayLightTimeCalculator( barycentricInterpolatedBodies, observationSettingsList );
    std::shared_ptr< LightTimeCalculator< long double, Time > > earthCenteredInterpolatedLightTimeCalculator =
        getOneWayLightTimeCalculator( earthCenteredInterpolatedBodies, observationSettingsList );


    {
        long double originalLightTime = earthCenteredInterpolatedLightTimeCalculator->calculateLightTime( testTime );

        for ( int i = -20; i <= 0; i++ )
        {
            long double perturbationStep = std::pow( 10, i );

            Eigen::Vector6ld upperturbedState =
                barycentricInterpolatedLightTimeCalculator->getStateFunctionOfReceivingBody( )( testTime + perturbationStep );

            Eigen::Vector6ld downperturbedState =
                barycentricInterpolatedLightTimeCalculator->getStateFunctionOfReceivingBody( )(  testTime - perturbationStep );


            std::cout << i << " " << ( upperturbedState - downperturbedState ).transpose( )  / ( 2.0 * perturbationStep )<< std::endl;
        }
    }

    {
        long double originalLightTime = earthCenteredInterpolatedLightTimeCalculator->calculateLightTime( testTime );

        for ( int i = -20; i <= 0; i++ )
        {
            long double perturbationStep = std::pow( 10, i );

            Eigen::Vector6ld upperturbedState =
                barycentricInterpolatedLightTimeCalculator->getStateFunctionOfTransmittingBody( )( testTime + perturbationStep );

            Eigen::Vector6ld downperturbedState =
                barycentricInterpolatedLightTimeCalculator->getStateFunctionOfTransmittingBody( )(  testTime - perturbationStep );


            std::cout << i << " " << ( upperturbedState - downperturbedState ).transpose( )  / ( 2.0 * perturbationStep )<< std::endl;
        }
    }

    {
        long double originalLightTime = earthCenteredInterpolatedLightTimeCalculator->calculateLightTime( testTime );

        for ( int i = -20; i <= 0; i++ )
        {
            long double perturbationStep = std::pow( 10, i );

            long double upperturbedLightTime =
                barycentricInterpolatedLightTimeCalculator->calculateLightTime( testTime + perturbationStep );

            long double downperturbedLightTime =
                barycentricInterpolatedLightTimeCalculator->calculateLightTime( testTime - perturbationStep );


            std::cout << i << " " << ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * perturbationStep )<< std::endl;
        }
    }
//
//    {
//        for ( int i = -20; i <= 0; i++ )
//        {
//            long double perturbationStep = std::pow( 10, i );
////
////            Eigen::Vector6ld upperturbedState =
////                barycentricBodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime + perturbationStep );
////
////            Eigen::Vector6ld downperturbedState =
////                barycentricBodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime - perturbationStep );
//
//            Eigen::Vector6ld upperturbedState =
//                barycentricInterpolatedBodies.at( "Earth" )->getEphemeris( )->getCartesianLongStateFromExtendedTime( testTime + perturbationStep );
//
//            Eigen::Vector6ld downperturbedState =
//                barycentricInterpolatedBodies.at( "Earth" )->getEphemeris( )->getCartesianLongStateFromExtendedTime( testTime - perturbationStep );
//
//
//
//
//            std::cout << i << " " << ( upperturbedState - downperturbedState ).transpose( )  / ( 2.0 * perturbationStep )<< std::endl;
//        }
//    }
//
//    std::shared_ptr< ephemerides::TabulatedCartesianEphemeris< long double, Time > > earthEphemeris = std::dynamic_pointer_cast<
//        ephemerides::TabulatedCartesianEphemeris< long double, Time > >( barycentricInterpolatedBodies.at( "Moon" )->getEphemeris( ) );
//
//    std::cout<<earthEphemeris<<std::endl;
//
//    std::shared_ptr< interpolators::LagrangeInterpolator< Time, Eigen::Vector6ld > > stateInterpolator =
//        std::dynamic_pointer_cast< interpolators::LagrangeInterpolator< Time, Eigen::Vector6ld > >(
//            earthEphemeris->getInterpolator( ) );
//
//    std::cout<<stateInterpolator<<std::endl;
//
//    {
//        for ( int i = -20; i <= 0; i++ )
//        {
//            long double perturbationStep = std::pow( 10, i );
////
////            Eigen::Vector6ld upperturbedState =
////                barycentricBodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime + perturbationStep );
////
////            Eigen::Vector6ld downperturbedState =
////                barycentricBodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< long double, Time >( testTime - perturbationStep );
//
//            Eigen::Vector6ld upperturbedState = stateInterpolator->interpolate( testTime + perturbationStep );
//
//            Eigen::Vector6ld downperturbedState = stateInterpolator->interpolate( testTime - perturbationStep );
//
//
//
//
//            std::cout << i << " " << ( upperturbedState - downperturbedState ).transpose( ) / ( 2.0 * perturbationStep )<< std::endl;
//        }
//    }



}

//
//int main( )
//{
//    double initialTimeEnvironment = Time(107561, 2262.19) - 2.0 * 3600.0;
//    double finalTimeEnvironment = Time(108258, 2771.19) + 2.0 * 3600.0;
//
//    // Load spice kernels
//    spice_interface::loadStandardSpiceKernels( );
//    spice_interface::loadSpiceKernelInTudat( "/home/dominic/Tudat/Data/GRAIL_Spice/grail_120301_120529_sci_v02.bsp" );
//
//    SystemOfBodies earthCenteredBodies = createEnvironment(
//        initialTimeEnvironment, finalTimeEnvironment,
//        "Earth", "Earth", false, TUDAT_NAN, TUDAT_NAN );
////    SystemOfBodies earthCenteredInterpolatedBodies = createEnvironment(
////        initialTimeEnvironment, finalTimeEnvironment,
////        "Earth", "Earth", true, 120.0, 10.0 );
//    SystemOfBodies barycentricBodies = createEnvironment(
//        initialTimeEnvironment, finalTimeEnvironment,
//        "Sun", "Earth", false, TUDAT_NAN, TUDAT_NAN );
//    SystemOfBodies barycentricInterpolatedBodies = createEnvironment(
//        initialTimeEnvironment, finalTimeEnvironment,
//        "Sun", "Earth", true, 120.0, 10.0 );
//
//
//    // Define link ends.
//    LinkEnds testLinkEnds;
//    testLinkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", "" ) );
//    testLinkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "GRAIL-A", "" );
//
//    // Create observation settings
//    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
//    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >( one_way_range, testLinkEnds ) );
//
//    std::shared_ptr< LightTimeCalculator< long double, Time > > earthCenteredLightTimeCalculator =
//        getOneWayLightTimeCalculator( earthCenteredBodies, observationSettingsList );
//    std::shared_ptr< LightTimeCalculator< long double, Time > > barycentricLightTimeCalculator =
//        getOneWayLightTimeCalculator( barycentricBodies, observationSettingsList );
//    std::shared_ptr< LightTimeCalculator< long double, Time > > barycentricInterpolatedLightTimeCalculator =
//        getOneWayLightTimeCalculator( barycentricInterpolatedBodies, observationSettingsList );
//
////    std::shared_ptr< LightTimeCalculator< long double, Time > > earthCenteredInterpolatedLightTimeCalculator =
////        getOneWayLightTimeCalculator( earthCenteredInterpolatedBodies, observationSettingsList );
//
//    Time testTime = initialTimeEnvironment + 86400.0;
//
//
//    std::cout<<std::setprecision( 19 );
//    double lightTimeBarycenter = barycentricLightTimeCalculator->calculateLightTime( testTime );
//    double lightTimeEarth = earthCenteredLightTimeCalculator->calculateLightTime( testTime );
//
//    Time receptionTime = testTime;
//
//    typedef Eigen::Matrix< long double, 6, 1 > Vector6ld;
//
//    Vector6ld receiverStateEarth = earthCenteredLightTimeCalculator->getStateFunctionOfReceivingBody( )( receptionTime );
//    Vector6ld transmitterStateEarth = earthCenteredLightTimeCalculator->getStateFunctionOfTransmittingBody( )( receptionTime - lightTimeEarth );
//
//    Vector6ld receiverStateBarycenter = barycentricLightTimeCalculator->getStateFunctionOfReceivingBody( )( receptionTime ) -
//        barycentricBodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< long double, Time >( receptionTime );
//    Vector6ld transmitterStateBarycenter = barycentricLightTimeCalculator->getStateFunctionOfTransmittingBody( )( receptionTime - lightTimeBarycenter ) -
//        barycentricBodies.at( "Earth" )->getStateInBaseFrameFromEphemeris< long double, Time >( receptionTime - lightTimeBarycenter );
//
//    Vector6ld receiverStateBarycenterUncorrected = barycentricLightTimeCalculator->getStateFunctionOfReceivingBody( )( receptionTime );
//    Vector6ld transmitterStateBarycenterUncorrected = barycentricLightTimeCalculator->getStateFunctionOfTransmittingBody( )( receptionTime - lightTimeBarycenter );
//
//
//    std::cout<<receiverStateEarth.transpose( )<<std::endl;
//    std::cout<<receiverStateBarycenter.transpose( )<<std::endl;
//    std::cout<<receiverStateEarth.transpose( )-receiverStateBarycenter.transpose( )<<std::endl<<std::endl;
//
//    std::cout<<transmitterStateEarth.transpose( )<<std::endl;
//    std::cout<<transmitterStateBarycenter.transpose( )<<std::endl;
//    std::cout<<transmitterStateEarth.transpose( ) - transmitterStateBarycenter.transpose( )<<std::endl<<std::endl;
//
//    std::cout<<lightTimeBarycenter - lightTimeEarth<<std::endl;
//
//
//}