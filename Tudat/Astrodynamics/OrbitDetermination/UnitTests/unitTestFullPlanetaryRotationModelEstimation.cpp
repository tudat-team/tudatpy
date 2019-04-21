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

#include <string>
#include <thread>

#include <limits>

#include <boost/make_shared.hpp>
#include <boost/test/unit_test.hpp>

#include "Tudat/Basics/testMacros.h"
#include "Tudat/SimulationSetup/tudatEstimationHeader.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/periodicSpinVariation.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_full_planetary_rotational_parameters_estimation )

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
using namespace tudat::coordinate_conversions;
using namespace tudat::ground_stations;
using namespace tudat::observation_models;


//! Unit test to check if periodic spin variation (for a full planetary rotational model) is estimated correctly
BOOST_AUTO_TEST_CASE( test_FullPlanetaryRotationalParameters )
{

    //Load spice kernels.f
    std::string kernelsPath = input_output::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    //Define environment settings

    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = 6.0E7;
    double maximumTimeStep = 3600.0;

    double buffer = 10.0 * maximumTimeStep;

    std::map< std::string, std::shared_ptr< BodySettings > > bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    bodySettings[ "Earth" ]->ephemerisSettings-> resetMakeMultiArcEphemeris( true );
    bodySettings[ "Moon" ]->ephemerisSettings->resetFrameOrigin( "Sun" );
    bodySettings[ "Mars" ]->rotationModelSettings = getHighAccuracyMarsRotationModel( initialEphemerisTime, finalEphemerisTime );

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

    // Define propagator settings.
    std::vector< std::string > centralBodies; centralBodies.push_back( "SSB" );

    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodyMap, accelerationMap, bodiesToIntegrate, centralBodies );


    std::vector< double > integrationArcStartTimes;
    std::vector< double > integrationArcEndTimes;

    std::vector< double > integrationArcLimits;

    double integrationStartTime = initialEphemerisTime + 1.0E4;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double arcDuration = 1.6E7;
    double arcOverlap = 2.0E4;
    int numberEstimationArcs = 0;

    double currentStartTime = integrationStartTime, currentEndTime = integrationStartTime + arcDuration;

    while ( currentEndTime < integrationEndTime )
    {
        integrationArcLimits.push_back( currentStartTime );
        integrationArcEndTimes.push_back( currentEndTime );
        integrationArcStartTimes.push_back( currentStartTime );
        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;

        numberEstimationArcs++;

    }
    integrationArcLimits.push_back( currentStartTime + arcOverlap );


    // Define links in simulation.
    std::vector< LinkEnds > linkEnds; linkEnds.resize( 2 );
    linkEnds[ 0 ][ transmitter ] = grazStation;
    linkEnds[ 0 ][ receiver ] = mslStation;

    linkEnds[ 1 ][ receiver ] = grazStation;
    linkEnds[ 1 ][ transmitter ] = mslStation;

    observation_models::ObservationSettingsMap observationSettingsMap;
    observationSettingsMap.insert( std::make_pair( linkEnds[ 0 ], std::make_shared< ObservationSettings >( one_way_range ) ) );
    observationSettingsMap.insert( std::make_pair( linkEnds[ 1 ], std::make_shared< ObservationSettings >( one_way_range ) ) );

    // Define integrator settings.
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< IntegratorSettings< double > >
            ( rungeKutta4, initialEphemerisTime - 4.0 * maximumTimeStep, 3600.0 );

    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > propagatorSettingsList;
    for( unsigned int i = 0; i < integrationArcStartTimes.size( ); i++ )
    {
        Eigen::Matrix< double, Eigen::Dynamic, 1 > currentInitialState = getInitialStateOfBody< double, double>(
                    bodiesToIntegrate.at( 0 ), centralBodies.at( 0 ), bodyMap, integrationArcStartTimes.at( i ) );

        propagatorSettingsList.push_back( std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( centralBodies, accelerationModelMap, bodiesToIntegrate, currentInitialState,
                      integrationArcEndTimes.at( i ) ) );
    }
    std::shared_ptr< PropagatorSettings< double > > propagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< double > >( propagatorSettingsList, 1 );


    // Define observation times.
    double observationTime;
    int numberOfObservationsPerArc = 5000;
    double timeBuffer = 9000.0;


    std::vector< double > initialObservationTimes;
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


    std::map< ObservableType, std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > > measurementSimulationInput;
    std::map< LinkEnds, std::pair< std::vector< double >, LinkEndType > > singleObservableSimulationInput;


    singleObservableSimulationInput[ linkEnds[ 0 ] ] = std::make_pair( initialObservationTimes, receiver );
    measurementSimulationInput[ one_way_range ] = singleObservableSimulationInput;



    // Set-up different cases with various parameters to estimate
    for ( int testCase = 0 ; testCase < 3 ; testCase++ ){

        // Set parameters that are to be estimated.
        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< ArcWiseInitialTranslationalStateEstimatableParameterSettings< double > >(
                                      "Earth", integrationArcStartTimes ) );

        // Estimate core factor and free core nutation rate
        if ( testCase == 0 ){
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", core_factor ) );
            parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", free_core_nutation_rate ) );
        }

        // Estimate periodic spin variation
        else if ( testCase == 1 ){
            parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >( "Mars", periodic_spin_variation ) );
        }

        // Estimate polar motion amplitude
        else if ( testCase == 2 ){
            parameterNames.push_back(  std::make_shared< EstimatableParameterSettings >( "Mars", polar_motion_amplitude ) );
        }


        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double >( parameterNames, bodyMap );


        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                    bodyMap, parametersToEstimate, observationSettingsMap, integratorSettings, propagatorSettings );


        // Define initial parameter estimate.
        Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
                parametersToEstimate->template getFullParameterValues< double >( );


        typedef Eigen::Matrix< double, Eigen::Dynamic, 1 > ObservationVectorType;
        typedef std::map< LinkEnds, std::pair< ObservationVectorType, std::pair< std::vector< double >, LinkEndType > > > SingleObservablePodInputType;
        typedef std::map< ObservableType, SingleObservablePodInputType > PodInputDataType;

        PodInputDataType observationsAndTimes = simulateObservations< double, double >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( )  );


        // Define perturbation of parameter estimate
        Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
        for (int i = 0 ; i < numberEstimationArcs ; i++)
        {
            initialParameterEstimate.segment( 0 + i * 6, 3 ) += 1.0E1 * Eigen::Vector3d::Identity();
            initialParameterEstimate.segment( 3 + i * 6, 3 ) += 1.0E-3 * Eigen::Vector3d::Identity();
        }


        if ( testCase == 0 ){
            initialParameterEstimate[ 6 * numberEstimationArcs ] += 1.0E-4;
            initialParameterEstimate[ 6 * numberEstimationArcs + 1 ] += 1.0E-8;
        }

        if ( testCase == 1 ){
            for( int i = 6 * numberEstimationArcs + 0 ; i < static_cast< int >( initialParameterEstimate.rows( ) ) ; i++ )
            {
                initialParameterEstimate[ i ] += 1.0E-8;
            }
        }

        if ( testCase == 2 ){
            for( int i = 6 * numberEstimationArcs ; i < static_cast< int >( initialParameterEstimate.rows( ) ) ; i++ )
            {
                initialParameterEstimate[ i ] += 1.0E-9;
            }
        }

        parametersToEstimate->resetParameterValues( initialParameterEstimate );

        std::shared_ptr< PodInput< double, double > > podInput = std::make_shared< PodInput< double, double > >(
                    observationsAndTimes, ( initialParameterEstimate ).rows( ) );

        std::shared_ptr< PodOutput< double, double > > podOutput = orbitDeterminationManager.estimateParameters( podInput );

        Eigen::VectorXd parameterError = podOutput->parameterEstimate_ - truthParameters;

//        std::cout << "TEST " << testCase << ": " << "\n\n";
//        std::cout << "truth parameter: " << truthParameters << "\n\n";
//        std::cout << "Initial parameter error: " << ( initialParameterEstimate - truthParameters ).transpose() << "\n\n";
//        std::cout << "Final parameter error: " << parameterError.transpose() << "\n\n";

        for( int i = 0; i < numberEstimationArcs; i++ )
        {
            for( unsigned int j = 0; j < 3; j++ )
            {
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j ) ), 3.0E-2 );
                BOOST_CHECK_SMALL( std::fabs( parameterError( i * 6 + j + 3 ) ), 1.0E-7  );
            }
        }


        if ( testCase == 0 ) {

            BOOST_CHECK_SMALL( std::fabs( parameterError( 6 * numberEstimationArcs ) ), 1.0E-5 );
            BOOST_CHECK_SMALL( std::fabs( parameterError( 6 * numberEstimationArcs + 1 ) ), 1.0E-8 );

        }
        else {

            for( int i = 6 * numberEstimationArcs ; i < static_cast< int >( initialParameterEstimate.rows( ) ); i++ )
            {
                if ( testCase == 1 ){
                    BOOST_CHECK_SMALL( std::fabs( parameterError( i ) ), 1.0E-12 );
                }
                if ( testCase == 2 ){
                    BOOST_CHECK_SMALL( std::fabs( parameterError( i ) ), 1.0E-6 );
                }
            }
        }

    }




}

BOOST_AUTO_TEST_SUITE_END( )

}

}
