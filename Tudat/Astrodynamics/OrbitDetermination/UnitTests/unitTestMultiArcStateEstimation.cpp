/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#define BOOST_TEST_MAIN

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatSimulationHeader.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/simulateObservations.h"
#include "Tudat/SimulationSetup/EstimationSetup/orbitDeterminationManager.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createGroundStations.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_multi_arc_state_estimation )

//Using declarations.
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
    std::string kernelsPath = input_output::getSpiceKernelPath( );
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

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings[ "Earth" ]->ephemerisSettings-> resetMakeMultiArcEphemeris( true );
    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrigin( "Sun" );
    bodySettings[ "Mars" ]->rotationModelSettings = std::make_shared< SimpleRotationModelSettings >(
                "ECLIPJ2000", "IAU_Mars",
                spice_interface::computeRotationQuaternionBetweenFrames(
                    "ECLIPJ2000", "IAU_Mars", initialEphemerisTime ),
                initialEphemerisTime, 2.0 * mathematical_constants::PI /
                ( physical_constants::JULIAN_DAY + 40.0 * 60.0 ) );
    NamedBodyMap bodyMap = createBodies( bodySettings );

    setGlobalFrameBodyEphemerides( bodyMap, "SSB", "ECLIPJ2000" );

    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );

    createGroundStation( bodyMap.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );

    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    groundStations.push_back( mslStation );

    // Set accelerations between bodies that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfEarth;
    accelerationsOfEarth[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
    accelerationsOfEarth[ "Moon" ].push_back( std::make_shared< AccelerationSettings >( central_gravity ) );
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
                bodyMap, accelerationMap, bodiesToIntegrate, centralBodies );


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

    // Set parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< StateScalarType > >(
                                  "Earth", integrationArcStartTimes ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >
                              ( "Mars", constant_rotation_rate ) );
    parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >
                               ( "Mars", rotation_pole_position ) );

    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType >( parameterNames, bodyMap );


    // Define links in simulation.
    std::vector< LinkEnds > linkEnds2;
    linkEnds2.resize( 2 );
    linkEnds2[ 0 ][ transmitter ] = grazStation;
    linkEnds2[ 0 ][ receiver ] = mslStation;

    linkEnds2[ 1 ][ receiver ] = grazStation;
    linkEnds2[ 1 ][ transmitter ] = mslStation;

    observation_models::ObservationSettingsMap observationSettingsMap;
    observationSettingsMap.insert( std::make_pair( linkEnds2[ 0 ], std::make_shared< ObservationSettings >(
                                       one_way_range ) ) );
    observationSettingsMap.insert( std::make_pair( linkEnds2[ 1 ], std::make_shared< ObservationSettings >(
                                       one_way_range ) ) );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< TimeType > > integratorSettings =
            std::make_shared< IntegratorSettings< TimeType > >
            ( rungeKutta4, TimeType( initialEphemerisTime - 4.0 * maximumTimeStep ), 3600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< StateScalarType > > > propagatorSettingsList;
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > currentInitialState =
                getInitialStateOfBody< TimeType, StateScalarType>(
                    bodiesToIntegrate.at( 0 ), centralBodies.at( 0 ), bodyMap, integrationArcStartTimes.at( i ) );
        propagatorSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< StateScalarType > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate,
                      currentInitialState,
                      integrationArcEndTimes.at( i ) ) );
    }
    std::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< StateScalarType > >( propagatorSettingsList, linkArcs );

    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodyMap, parametersToEstimate,
                observationSettingsMap, integratorSettings, propagatorSettings );

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


    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > > measurementSimulationInput;
    std::map< LinkEnds, std::pair< std::vector< TimeType >, LinkEndType > > singleObservableSimulationInput;


    singleObservableSimulationInput[ linkEnds2[ 0 ] ] = std::make_pair( initialObservationTimes, receiver );
    measurementSimulationInput[ one_way_range ] = singleObservableSimulationInput;


    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;
    typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< TimeType >, LinkEndType > > > SingleObservablePodInputType;
    typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

    PodInputDataType observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( )  );


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
#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
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

BOOST_AUTO_TEST_SUITE_END( )

}

}
