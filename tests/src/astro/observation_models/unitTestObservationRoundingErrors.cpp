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


namespace tudat
{

}
namespace unit_tests
{

BOOST_AUTO_TEST_SUITE( test_observation_model_rounding_error )


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
    std::vector<std::string> bodiesToCreate = { "Earth", "Sun", "Moon" };
    std::string globalFrameOrientation = "J2000";
    BodyListSettings bodySettings;
    if ( useInterpolatedEphemerides )
    {
        bodySettings = getDefaultBodySettings(
            bodiesToCreate, initialTimeEnvironment, finalTimeEnvironment, globalFrameOrigin, globalFrameOrientation,
            planetEphemerisInterpolationStep );
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
    if ( useInterpolatedEphemerides )
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared<InterpolatedSpiceEphemerisSettings>(
                initialTimeEnvironment, finalTimeEnvironment, spacecraftEphemerisInterpolationStep,
                spacecraftCentralBody, globalFrameOrientation );
    }
    else
    {
        bodySettings.at( spacecraftName )->ephemerisSettings =
            std::make_shared<DirectSpiceEphemerisSettings>(
                spacecraftCentralBody, globalFrameOrientation, false, false, false );
    }

    // Create bodies
    SystemOfBodies bodies = createSystemOfBodies<long double, Time>( bodySettings );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    return bodies;

}

std::shared_ptr<LightTimeCalculator<long double, Time> > getOneWayLightTimeCalculator(
    const SystemOfBodies &bodies,
    const std::vector<std::shared_ptr<ObservationModelSettings> > oneWayObservationSettingsList )
{
    // Create observation simulators
    std::shared_ptr<ObservationSimulator<1, long double, Time> > observationSimulator =
        std::dynamic_pointer_cast<ObservationSimulator<1, long double, Time> >(
            createObservationSimulators<long double, Time>( oneWayObservationSettingsList, bodies ).at( 0 ));
    std::shared_ptr<OneWayRangeObservationModel<long double, Time> > oneWayRangeObservationModel =
        std::dynamic_pointer_cast<OneWayRangeObservationModel<long double, Time> >(
            observationSimulator->getObservationModels( ).begin( )->second );
    return oneWayRangeObservationModel->getLightTimeCalculator( );
}

template< typename StateScalarType, typename TimeType >
void checkStateFunctionNumericalErrors(
    const std::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > stateFunction,
    const TimeType testTime,
    const std::vector< int > timeExponents = { 0, -6, -7, -8, -9, -10, -11, -12} )
{
    // Compute nominal state at test time
    Eigen::Vector6ld nominalState = stateFunction( testTime );

    // Time steps at which the relative numerical error will be around 1
    Eigen::Vector3ld limitTimes = std::numeric_limits< StateScalarType >::epsilon( ) *
                                  nominalState.segment( 0, 3 ).cwiseQuotient( nominalState.segment( 3, 3 ));

    std::cout<<"Limits: "<<limitTimes.transpose( )<<std::endl;
    // Compute numerical position partial, and compute error w.r.t. computation for Delta t = 1 s
    Eigen::Vector3ld nominalPartial = Eigen::Vector3ld::Zero( );
    for ( int i = 0; i < timeExponents.size( ); i++ )
    {
        // Define time step
        long double perturbationStep = std::pow( 10, timeExponents.at( i ) );

        // Calculate numerical partial
        Eigen::Vector3ld upperturbedState = stateFunction(
                testTime + perturbationStep ).segment( 0, 3 );
        Eigen::Vector3ld downperturbedState = stateFunction(
                testTime - perturbationStep ).segment( 0, 3 );
        Eigen::Vector3ld currentPartial = ( upperturbedState - downperturbedState ) / ( 2.0 * perturbationStep );

        if ( i == 0 )
        {
            // Set nominal partial
            nominalPartial = currentPartial;
        }
        else
        {
            // Test real (w.r.t. Delta t = 1 s) errors against theoretical limit (with safety factor of 2 to account for randomness in rounding error)
            Eigen::Vector3ld expectedRelativeErrorLevel = limitTimes / perturbationStep;
            Eigen::Vector3ld realRelativeErrorLevels = ( ( currentPartial - nominalPartial ).cwiseQuotient( nominalPartial ) ).segment( 0, 3 );
            std::cout<<timeExponents.at( i )<<" "<<expectedRelativeErrorLevel.transpose( )<<std::endl;
            std::cout<<timeExponents.at( i )<<" "<<realRelativeErrorLevels.transpose( )<<std::endl<<std::endl;

            for ( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( realRelativeErrorLevels( j ) ), 2.0 * std::fabs( expectedRelativeErrorLevel( j ) ) );
            }
        }
    }
}

BOOST_AUTO_TEST_CASE( test_ObservationModelContinuity )
{
    double initialTimeEnvironment = Time( 107561, 2262.19 ) - 2.0 * 3600.0;
    double finalTimeEnvironment = Time( 108258, 2771.19 ) + 2.0 * 3600.0;

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
    testLinkEnds[ receiver ] = std::pair<std::string, std::string>( std::make_pair( "Earth", "Station1" ));
    testLinkEnds[ transmitter ] = std::make_pair<std::string, std::string>( "GRAIL-A", "" );

    Time testTime = initialTimeEnvironment + 86400.1;
    // Create observation settings
    std::vector<std::shared_ptr<ObservationModelSettings> > observationSettingsList;
    observationSettingsList.push_back( std::make_shared<ObservationModelSettings>( one_way_range, testLinkEnds ));
    
    std::shared_ptr<LightTimeCalculator<long double, Time> > earthCenteredLightTimeCalculator =
        getOneWayLightTimeCalculator( earthCenteredBodies, observationSettingsList );
    std::shared_ptr<LightTimeCalculator<long double, Time> > barycentricLightTimeCalculator =
        getOneWayLightTimeCalculator( barycentricBodies, observationSettingsList );
    std::shared_ptr<LightTimeCalculator<long double, Time> > barycentricInterpolatedLightTimeCalculator =
        getOneWayLightTimeCalculator( barycentricInterpolatedBodies, observationSettingsList );
    std::shared_ptr<LightTimeCalculator<long double, Time> > earthCenteredInterpolatedLightTimeCalculator =
        getOneWayLightTimeCalculator( earthCenteredInterpolatedBodies, observationSettingsList );

    checkStateFunctionNumericalErrors< long double, Time >(
        barycentricInterpolatedLightTimeCalculator->getStateFunctionOfTransmittingBody( ), testTime );
    checkStateFunctionNumericalErrors< long double, Time >(
        barycentricInterpolatedLightTimeCalculator->getStateFunctionOfReceivingBody( ), testTime );

    checkStateFunctionNumericalErrors< long double, Time >(
        earthCenteredInterpolatedLightTimeCalculator->getStateFunctionOfTransmittingBody( ), testTime );
//    checkStateFunctionNumericalErrors< long double, Time >(
//        earthCenteredInterpolatedLightTimeCalculator->getStateFunctionOfReceivingBody( ), testTime );

//
//    {
//
//        Eigen::Vector6ld nominalPartial = Eigen::Vector6ld::Zero( );
//        for ( int i = 0; i >= -11; i-- )
//        {
//            long double perturbationStep = std::pow( 10, i );
//
//            Eigen::Vector6ld upperturbedState =
//                barycentricInterpolatedLightTimeCalculator->getStateFunctionOfTransmittingBody( )(
//                    testTime + perturbationStep );
//
//            Eigen::Vector6ld downperturbedState =
//                barycentricInterpolatedLightTimeCalculator->getStateFunctionOfTransmittingBody( )(
//                    testTime - perturbationStep );
//
//            Eigen::Vector6ld currentPartial =  ( upperturbedState - downperturbedState ) / ( 2.0 * perturbationStep );
//            if( i == 0 )
//            {
//                nominalPartial =  currentPartial;
//            }
//            else
//            {
//                TUDAT_CHECK_MATRIX_CLOSE_FRACTION( nominalPartial, currentPartial, 0.05 );
//            }
//        }
//    }
//
//    {
//        long double nominalPartial = 0.0L;
//
//        for ( int i = 0; i >= -9; i-- )
//        {
//            long double perturbationStep = std::pow( 10, i );
//
//            long double upperturbedLightTime =
//                barycentricInterpolatedLightTimeCalculator->calculateLightTime( testTime + perturbationStep );
//
//            long double downperturbedLightTime = barycentricInterpolatedLightTimeCalculator->calculateLightTime( testTime - perturbationStep );
//
//            long double currentPartial = ( upperturbedLightTime - downperturbedLightTime ) / ( 2.0 * perturbationStep );
//            if( i == 0 )
//            {
//                nominalPartial =  currentPartial;
//            }
//            else
//            {
//                BOOST_CHECK_CLOSE( nominalPartial, currentPartial, 0.05 );
//            }
//        }
//    }
}

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