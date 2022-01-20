/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_multi_arc_state_estimation )

//Using declarations.
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

template< typename ObservationScalarType = double , typename TimeType = double , typename StateScalarType  = double >
Eigen::VectorXd  executeParameterEstimation(
        const int linkArcs )
{
    //Load spice kernels.f
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    //Define setting for total number of bodies and those which need to be integrated numerically.
    //The first numberOfNumericalBodies from the bodyNames vector will be integrated numerically.

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = TimeType( 6.0E7 );
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings.at( "Earth" )->ephemerisSettings->resetMakeMultiArcEphemeris( true );
    bodySettings.at( "Moon" )->ephemerisSettings->resetFrameOrigin( "Sun" );
    bodySettings.at( "Mars" )->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Mars",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    

    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );

    createGroundStation( bodies.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    groundStations.push_back( mslStation );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Earth" ] = accelerationsOfEarth;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToEstimate;
    bodiesToEstimate.push_back( "Earth" );
    std::vector< std::string > bodiesToIntegrate;
    bodiesToIntegrate.push_back( "Earth" );
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );

    // Define propagator settings.
    std::vector< std::string > centralBodies;
    std::map< std::string, std::string > centralBodyMap;

    centralBodies.resize( numberOfNumericalBodies );
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies[ i ] = "SSB";
    }

    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );


    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;

    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime + 1.0E4;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double arcDuration = 1.6E7;
    double arcOverlap = 2.0E4;

    double currentStartTime = integrationStartTime, currentEndTime = integrationStartTime + arcDuration;

    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );
    integrationArcLimits.push_back( currentStartTime + arcOverlap );


    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsList;
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState =
                getInitialStateOfBody< TimeType, StateScalarType>(
                    bodiesToIntegrate.at( 0 ), centralBodies.at( 0 ), bodies, integrationArcStartTimes.at( i ) );
        propagatorSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                      currentInitialState,
                      integrationArcEndTimes.at( i ) ) );
    }
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettingsList, linkArcs );


    std::cout<<"************************************* RUNNING TEST *************************"<<std::endl;
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< double >( propagatorSettings, bodies, integrationArcStartTimes );

//    parameterNames = getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                              ( "Mars", constant_rotation_rate ) );
    parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >
                               ( "Mars", rotation_pole_position ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodies );


    // Define links in simulation.
    std::vector< LinkEnds > linkEnds2;
    linkEnds2.resize( 2 );
    linkEnds2[ 0 ][ transmitter ] = grazStation;
    linkEnds2[ 0 ][ receiver ] = mslStation;

    linkEnds2[ 1 ][ receiver ] = grazStation;
    linkEnds2[ 1 ][ transmitter ] = mslStation;

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                       one_way_range, linkEnds2[ 0 ] ) );
    observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                       one_way_range, linkEnds2[ 1 ] ) );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime - 4.0 * maximumTimeStep ), 3600.0 );

    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodies, parametersToEstimate,
                observationSettingsList, integratorSettings, propagatorSettings );

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


    TimeType observationTime;
    int numberOfObservationsPerArc = 5000;
    double timeBuffer = 9000.0;


    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservationsPerArc * integrationArcStartTimes.size( ) );

    for( unsigned int i = 0; i < integrationArcLimits.size( ) - 1; i++ )
    {
        double currentTimeStep = ( integrationArcLimits[ i + 1 ] - integrationArcLimits[ i ] - 2.0 * timeBuffer ) /
                static_cast< double >( numberOfObservationsPerArc - 1 );
        observationTime = integrationArcLimits[ i ] + timeBuffer;
        for( int j = 0; j < numberOfObservationsPerArc; j++ )
        {
            initialObservationTimes[ j + i * numberOfObservationsPerArc ] = observationTime;
            observationTime += currentTimeStep;
        }
    }

    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< > >(
                    one_way_range, linkEnds2[ 0 ], initialObservationTimes, receiver ) );

    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies  );


    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;

    for( unsigned int i = 0; i < numberOfNumericalBodies * integrationArcStartTimes.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
    }
    for( unsigned int i = numberOfNumericalBodies * integrationArcStartTimes.size( );
         i < static_cast< unsigned int >( initialParameterEstimate.rows( ) ); i++ )
    {
        initialParameterEstimate[ i ] *= ( 1.0 + 1.0E-6 );
    }

    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    std::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput =
            std::make_shared< PodInput< ObservationScalarType, TimeType > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ) );

    std::shared_ptr< PodOutput< StateScalarType, TimeType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput );

    return ( podOutput->parameterEstimate_ - truthParameters ).template cast< double >( );
}


BOOST_AUTO_TEST_CASE( test_MultiArcStateEstimation )
{
    // Execute test for linked arcs and separate arcs.
    for( unsigned int testCase = 0; testCase < 2; testCase++ )
    {
#if( TUDAT_BUILD_WITH_EXTENDED_PRECISION_PROPAGATION_TOOLS )
        Eigen::VectorXd parameterError = executeParameterEstimation< long double, tudat::Time, long double >(
                    testCase );
        int numberOfEstimatedArcs = ( parameterError.rows( ) - 3 ) / 6;

        std::cout <<"Estimation error: "<< parameterError.transpose( ) << std::endl;
        for( int i = 0; i < numberOfEstimatedArcs; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 1E-4 );
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 1.0E-10  );
            }
        }

        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 3 ) ), 1.0E-20 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 2 ) ), 1.0E-12 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 1 ) ), 1.0E-12 );
#else
        Eigen::VectorXd parameterError = executeParameterEstimation< double, double, double >(
                    testCase );
        int numberOfEstimatedArcs = ( parameterError.rows( ) - 3 ) / 6;

        std::cout << parameterError.transpose( ) << std::endl;
        for( int i = 0; i < numberOfEstimatedArcs; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 1E-1 );
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 1.0E-7  );
            }
        }

        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 3 ) ), 1.0E-17 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 2 ) ), 1.0E-9 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 1 ) ), 1.0E-9 );
#endif
    }

}

template< typename ObservationScalarType = double , typename TimeType = double , typename StateScalarType  = double >
Eigen::VectorXd  executeMultiBodyMultiArcParameterEstimation( )
{
    //Load spice kernels
    spice_interface::loadStandardSpiceKernels( );

    // Specify simulation time
    TimeType initialEphemerisTime = TimeType( 1.0E7 );
    TimeType finalEphemerisTime = initialEphemerisTime + 4.0 * 86400.0;

    //Define setting for celestial bodies
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - 86400.0, finalEphemerisTime + 86400.0,
                                    "Earth", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // CReate vehicles
    std::vector< std::string > vehicleNames = { "Borzi1", "Borzi2" };
    int numberOfVehicles = vehicleNames.size( );
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        bodies.createEmptyBody( vehicleNames.at( i ) );
        bodies.at( vehicleNames.at( i ) )->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                                           std::map< double, std::shared_ptr< Ephemeris > >( ), "Earth", "ECLIPJ2000" ) );
    }


    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
        accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
        accelerationMap[ vehicleNames.at( i ) ] = accelerationsOfVehicle;
    }

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate = vehicleNames;
    unsigned int numberOfNumericalBodies = bodiesToIntegrate.size( );
    std::vector< std::string > centralBodies;
    for( unsigned int i = 0; i < numberOfNumericalBodies; i++ )
    {
        centralBodies.push_back( "Earth" );
    }

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
                bodies, accelerationMap, bodiesToIntegrate, centralBodies );


    // Define integration arc limits
    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;
    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime;
    double integrationEndTime = finalEphemerisTime;
    double arcDuration = 1.01 * 86400.0;
    double arcOverlap = 300.0;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );

    integrationArcLimits.push_back( currentStartTime + arcOverlap );


    // Create initial state container per arc
    std::vector< Eigen::VectorXd > allBodiesPerArcInitialStates;
    allBodiesPerArcInitialStates.resize( integrationArcStartTimes.size( ) );
    for( unsigned int j = 0; j < integrationArcStartTimes.size( ); j++ )
    {
        allBodiesPerArcInitialStates[ j ] = Eigen::VectorXd::Zero( 6 * numberOfVehicles );
    }

    // Create initial state container per body
    std::vector< Eigen::VectorXd > singleBodyForAllArcsInitialStates;
    singleBodyForAllArcsInitialStates.resize( numberOfVehicles );
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        singleBodyForAllArcsInitialStates[ i ] = Eigen::VectorXd::Zero( 6 * integrationArcStartTimes.size( ) );
    }

    // Define (semi-arbitrary) initial states for vehicles
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        Eigen::Vector6d nominalInitialStateInKeplerianElements;
        nominalInitialStateInKeplerianElements( semiMajorAxisIndex ) = 7500.0E3;
        nominalInitialStateInKeplerianElements( eccentricityIndex ) = 0.1;
        nominalInitialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
        nominalInitialStateInKeplerianElements( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
        nominalInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
        nominalInitialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

        for( unsigned int j = 0; j < integrationArcStartTimes.size( ); j++ )
        {
            Eigen::Vector6d currentInitialStateInKeplerianElements = nominalInitialStateInKeplerianElements;
            currentInitialStateInKeplerianElements( inclinationIndex ) = 10.0 * static_cast< double >( i );
            currentInitialStateInKeplerianElements( trueAnomalyIndex ) = 10.0 * static_cast< double >( j );
            double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
            const Eigen::Vector6d currentInitialState = convertKeplerianToCartesianElements(
                        currentInitialStateInKeplerianElements, earthGravitationalParameter );
            allBodiesPerArcInitialStates[ j ].segment( i * 6, 6 ) = currentInitialState;
            singleBodyForAllArcsInitialStates[ i ].segment( j * 6, 6 ) = currentInitialState;
        }
    }


    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime  ), 30.0 );

    // Define propagator settings.
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsList;
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        propagatorSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                      allBodiesPerArcInitialStates.at( i ),
                      integrationArcEndTimes.at( i ), cowell, std::shared_ptr< DependentVariableSaveSettings >( ), 60.0 ) );
    }
    std::shared_ptr< MultiArcPropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettingsList );


    // Set parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
            getInitialMultiArcParameterSettings< >( propagatorSettings, bodies, integrationArcStartTimes );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodies );


    // Define links and observations in simulation.
    std::vector< LinkEnds > linkEndsList;
    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    linkEndsList.resize( numberOfVehicles );
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        linkEndsList[ i ][ observed_body ] = std::make_pair( vehicleNames.at( i ), "" );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                           position_observable, linkEndsList[ i ] ) );
    }


    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodies, parametersToEstimate,
                observationSettingsList, integratorSettings, propagatorSettings );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );


    // Define times at which observations are to be simulated
    TimeType observationTime;
    int numberOfObservationsPerArc = 100;
    double timeBuffer = 300.0;

    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservationsPerArc * integrationArcStartTimes.size( ) );
    for( unsigned int i = 0; i < integrationArcLimits.size( ) - 1; i++ )
    {
        double currentTimeStep = ( integrationArcLimits[ i + 1 ] - integrationArcLimits[ i ] - 2.0 * timeBuffer ) /
                static_cast< double >( numberOfObservationsPerArc - 1 );
        observationTime = integrationArcLimits[ i ] + timeBuffer;
        for( int j = 0; j < numberOfObservationsPerArc; j++ )
        {
            initialObservationTimes[ j + i * numberOfObservationsPerArc ] = observationTime;
            observationTime += currentTimeStep;
        }
    }


    // Define observation simulation settings
    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( int i = 0; i < numberOfVehicles; i++ )
    {
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        position_observable, linkEndsList[ i ], initialObservationTimes, observed_body ) );
    }

    // Simulate observations
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies  );

    // Perturb initial states
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    for( unsigned int i = 0; i < numberOfNumericalBodies * integrationArcStartTimes.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 1.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 4 + 6 * i ] += 1.0E-5;
        initialParameterEstimate[ 5 + 6 * i ] += 1.0E-5;
    }
    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define POD input
    std::shared_ptr< PodInput< ObservationScalarType, TimeType > > podInput =
            std::make_shared< PodInput< ObservationScalarType, TimeType > >(
                observationsAndTimes, ( initialParameterEstimate ).rows( ) );

    // Estimate parameters
    std::shared_ptr< PodOutput< StateScalarType, TimeType > > podOutput = orbitDeterminationManager.estimateParameters(
                podInput, std::make_shared< EstimationConvergenceChecker >( 3 ) );

    std::string outputFolder = "/home/dominic/Software/tudatBundleTest/tudatBundle/tudatApplications/master_thesis/Output/";

    return ( podOutput->parameterEstimate_ - truthParameters ).template cast< double >( );
}

BOOST_AUTO_TEST_CASE( test_MultiArcMultiBodyStateEstimation )
{

    Eigen::VectorXd parameterError = executeMultiBodyMultiArcParameterEstimation< double, double, double >( );
    int numberOfEstimatedArcs = ( parameterError.rows( ) ) / 6;

    std::cout <<"estimation error: "<< parameterError.transpose( ) << std::endl;
    for( int i = 0; i < numberOfEstimatedArcs; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 1E-4 );
            BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 1.0E-7  );
        }
    }
}

BOOST_AUTO_TEST_SUITE_END( )

}

}
