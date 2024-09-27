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


#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/simulation/simulation.h"
#include "tudat/astro/observation_models/linkTypeDefs.h"
#include "tudat/simulation/estimation_setup/simulateObservations.h"
#include "tudat/simulation/estimation_setup/orbitDeterminationManager.h"
#include "tudat/simulation/environment_setup/createGroundStations.h"


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_hybrid_arc_state_estimation )

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
using namespace tudat::unit_conversions;

//! Perform parameter estimation for a hybrid-arc case (Mars single-arc w.r.t. SSB and orbiter multi-arc w.r.t Mars)
template< typename ObservationScalarType = double , typename TimeType = double , typename StateScalarType  = double >
Eigen::VectorXd  executeParameterEstimation(
        const Eigen::Matrix< StateScalarType, 12, 1 > initialStateDifference = Eigen::Matrix< StateScalarType, 12, 1 >::Zero( ),
        const bool patchMultiArcs = 0,
        const std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > forcedMultiArcInitialStates =
        std::vector< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >( ),
        const double arcDuration = 4.0 * 86400.0,
        const double arcOverlap  = 5.0E3 )
{
    //Load spice kernels.
    std::string kernelsPath = paths::getSpiceKernelPath( );
    spice_interface::loadStandardSpiceKernels( );

    // Specify initial time
    double initialEphemerisTime = 1.0E7;
    double finalEphemerisTime = initialEphemerisTime + 20.0 * 86400.0;
    double maximumTimeStep = 3600.0;
    double buffer = 5.0 * maximumTimeStep;

    // Create bodies needed in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Mars" );
    bodyNames.push_back( "Jupiter" );
    bodyNames.push_back( "Earth" );
    BodyListSettings bodySettings =
            getDefaultBodySettings( bodyNames, initialEphemerisTime - buffer, finalEphemerisTime + buffer );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    // Create orbiter
    bodies.createEmptyBody( "Orbiter" );
    bodies.at( "Orbiter" )->setConstantBodyMass( 5.0E3 );
    bodies.at( "Orbiter" )->setEphemeris( std::make_shared< MultiArcEphemeris >(
                                            std::map< double, std::shared_ptr< Ephemeris > >( ),
                                            "Mars", "ECLIPJ2000" ) );
    bodies.processBodyFrameDefinitions( );

    // Create and set radiation pressure settings
    double referenceAreaRadiation = 4.0;
    double radiationPressureCoefficient = 1.2;
    std::vector< std::string > occultingBodies;
    occultingBodies.push_back( "Earth" );
    std::shared_ptr< RadiationPressureInterfaceSettings > orbiterRadiationPressureSettings =
            std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                "Sun", referenceAreaRadiation, radiationPressureCoefficient, occultingBodies );
    bodies.at( "Orbiter" )->setRadiationPressureInterface(
                "Sun", createRadiationPressureInterface(
                    orbiterRadiationPressureSettings, "Orbiter", bodies ) );




    // Create ground stations
    std::pair< std::string, std::string > grazStation = std::pair< std::string, std::string >( "Earth", "" );
    std::pair< std::string, std::string > mslStation = std::pair< std::string, std::string >( "Mars", "MarsStation" );

    createGroundStation( bodies.at( "Mars" ), "MarsStation", ( Eigen::Vector3d( ) << 100.0, 0.5, 2.1 ).finished( ),
                         coordinate_conversions::geodetic_position );
    std::vector< std::pair< std::string, std::string > > groundStations;
    groundStations.push_back( grazStation );
    groundStations.push_back( mslStation );

    // Set accelerations to act on Mars
    SelectedAccelerationMap singleArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfMars;
    accelerationsOfMars[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfMars[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    singleArcAccelerationMap[ "Mars" ] = accelerationsOfMars;

    std::vector< std::string > singleArcBodiesToIntegrate, singleArcCentralBodies;
    singleArcBodiesToIntegrate.push_back( "Mars" );
    singleArcCentralBodies.push_back( "SSB" );
    AccelerationMap singleArcAccelerationModelMap = createAccelerationModelsMap(
                bodies, singleArcAccelerationMap, singleArcBodiesToIntegrate, singleArcCentralBodies );

    // Define and perturb Mars initial states
    Eigen::VectorXd singleArcInitialStates = getInitialStatesOfBodies(
                singleArcBodiesToIntegrate, singleArcCentralBodies, bodies, initialEphemerisTime );
    singleArcInitialStates += initialStateDifference.segment(
                0, singleArcInitialStates.rows( ) );

    // Create Mars propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< > > singleArcPropagatorSettings =
            std::make_shared< TranslationalStatePropagatorSettings< > >(
                singleArcCentralBodies, singleArcAccelerationModelMap, singleArcBodiesToIntegrate,
                singleArcInitialStates, finalEphemerisTime );

    // Set accelerations to act on orbiter
    SelectedAccelerationMap multiArcAccelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfOrbiter;
    accelerationsOfOrbiter[ "Mars" ].push_back( std::make_shared< SphericalHarmonicAccelerationSettings >( 2, 2 ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationsOfOrbiter[ "Sun" ].push_back( std::make_shared< AccelerationSettings >( cannon_ball_radiation_pressure ) );
    accelerationsOfOrbiter[ "Jupiter" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    multiArcAccelerationMap[ "Orbiter" ] = accelerationsOfOrbiter;

    std::vector< std::string > multiArcBodiesToIntegrate, multiArcCentralBodies;
    multiArcBodiesToIntegrate.push_back( "Orbiter" );
    multiArcCentralBodies.push_back( "Mars" );
    AccelerationMap multiArcAccelerationModelMap = createAccelerationModelsMap(
                bodies, multiArcAccelerationMap, multiArcBodiesToIntegrate, multiArcCentralBodies );


    // Creater orbiter arc times
    std::vector< double > integrationArcStarts, integrationArcEnds, integrationArcLimits;
    double integrationStartTime = initialEphemerisTime;
    double integrationEndTime = finalEphemerisTime - 1.0E4;
    double currentStartTime = integrationStartTime;
    double currentEndTime = integrationStartTime + arcDuration;
    do
    {
        integrationArcStarts.push_back( currentStartTime );
        integrationArcLimits.push_back( currentStartTime );

        integrationArcEnds.push_back( currentEndTime );

        currentStartTime = currentEndTime - arcOverlap;
        currentEndTime = currentStartTime + arcDuration;
    }
    while( currentEndTime < integrationEndTime );
    integrationArcLimits.push_back( currentStartTime + arcOverlap );

    // Create list of multi-arc initial states
    unsigned int numberOfIntegrationArcs = integrationArcStarts.size( );
    std::vector< Eigen::VectorXd > multiArcSystemInitialStates;
    multiArcSystemInitialStates.resize( numberOfIntegrationArcs );

    // Define and perturb (quasi-arbitrary) arc initial states
    double marsGravitationalParameter =  bodies.at( "Mars" )->getGravityFieldModel( )->getGravitationalParameter( );
    if( forcedMultiArcInitialStates.size( ) == 0 )
    {
        for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            Eigen::Vector6d orbiterInitialStateInKeplerianElements;
            orbiterInitialStateInKeplerianElements( semiMajorAxisIndex ) = 6000.0E3;
            orbiterInitialStateInKeplerianElements( eccentricityIndex ) = 0.05;
            orbiterInitialStateInKeplerianElements( inclinationIndex ) = convertDegreesToRadians( 85.3 );
            orbiterInitialStateInKeplerianElements( argumentOfPeriapsisIndex )
                    = convertDegreesToRadians( 235.7 - j );
            orbiterInitialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
                    = convertDegreesToRadians( 23.4 + j );
            orbiterInitialStateInKeplerianElements( trueAnomalyIndex ) = convertDegreesToRadians( 139.87 + j * 10.0 );

            // Convert state from Keplerian elements to Cartesian elements.
            multiArcSystemInitialStates[ j ]  = convertKeplerianToCartesianElements(
                        orbiterInitialStateInKeplerianElements,
                        marsGravitationalParameter ) + initialStateDifference.segment(
                        singleArcInitialStates.rows( ), 6 ).template cast< double >( );
        }
    }
    else
    {
        for( unsigned int j = 0; j < numberOfIntegrationArcs; j++ )
        {
            multiArcSystemInitialStates[ j ] = forcedMultiArcInitialStates.at( j ) + initialStateDifference.segment(
                        singleArcInitialStates.rows( ), 6 );
        }
    }

    // Create propagation settings for each arc
    std::vector< std::shared_ptr< SingleArcPropagatorSettings< double > > > arcPropagationSettingsList;
    for( unsigned int i = 0; i < numberOfIntegrationArcs; i++ )
    {
        arcPropagationSettingsList.push_back(
                    std::make_shared< TranslationalStatePropagatorSettings< double > >
                    ( multiArcCentralBodies, multiArcAccelerationModelMap, multiArcBodiesToIntegrate,
                      multiArcSystemInitialStates.at( i ), integrationArcEnds.at( i ) ) );
    }

    // Define propagator settings (multi- and hybrid-arc)
    std::shared_ptr< MultiArcPropagatorSettings< > > multiArcPropagatorSettings =
            std::make_shared< MultiArcPropagatorSettings< > >( arcPropagationSettingsList, patchMultiArcs );
    std::shared_ptr< HybridArcPropagatorSettings< > > hybridArcPropagatorSettings =
            std::make_shared< HybridArcPropagatorSettings< > >(
                singleArcPropagatorSettings, multiArcPropagatorSettings );

    // Define integrator settings
    std::shared_ptr< IntegratorSettings< > > integratorSettings =
            std::make_shared< IntegratorSettings< > >
            ( rungeKutta4, initialEphemerisTime, 60.0 );

    // Set parameters that are to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
    parameterNames =
            getInitialHybridArcParameterSettings< >( hybridArcPropagatorSettings, bodies, integrationArcStarts );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Sun", gravitational_parameter ) );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Mars", gravitational_parameter ) );
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< StateScalarType > > parametersToEstimate =
            createParametersToEstimate< StateScalarType, TimeType >( parameterNames, bodies );


    // Define links and observables in simulation.
    std::vector< LinkDefinition > linkEnds2;
    linkEnds2.resize( 2 );
    linkEnds2[ 0 ][ transmitter ] = grazStation;
    linkEnds2[ 0 ][ receiver ] = std::make_pair< std::string, std::string >( "Orbiter", "" );

    linkEnds2[ 1 ][ receiver ] = grazStation;
    linkEnds2[ 1 ][ transmitter ] = std::make_pair< std::string, std::string >( "Orbiter", "" );

//    linkEnds2[ 2 ][ transmitter ] = grazStation;
//    linkEnds2[ 2 ][ receiver ] = std::make_pair< std::string, std::string >( "Mars", "" );

//    linkEnds2[ 3 ][ receiver ] = grazStation;
//    linkEnds2[ 3 ][ transmitter ] = std::make_pair< std::string, std::string >( "Mars", "" );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( unsigned int i = 0; i  < linkEnds2.size( ); i++ )
    {
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                           one_way_range, linkEnds2[ i ] ) );
        observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                                                           angular_position, linkEnds2[ i ] ) );
    }


    // Create orbit determination object.
    OrbitDeterminationManager< ObservationScalarType, TimeType > orbitDeterminationManager =
            OrbitDeterminationManager< ObservationScalarType, TimeType >(
                bodies, parametersToEstimate,
                observationSettingsList, integratorSettings, hybridArcPropagatorSettings );
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > initialParameterEstimate =
            parametersToEstimate->template getFullParameterValues< StateScalarType >( );

    // Define observation simulation settings
    TimeType observationTime;
    int numberOfObservationsPerArc = 5000;
    double timeBuffer = 9000.0;

    std::vector< TimeType > initialObservationTimes;
    initialObservationTimes.resize( numberOfObservationsPerArc * integrationArcStarts.size( ) );
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
    for( unsigned int i = 0; i < linkEnds2.size( ); i++ )
    {
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        one_way_range, linkEnds2[ i ], initialObservationTimes, receiver ) );
        measurementSimulationInput.push_back(
                    std::make_shared< TabulatedObservationSimulationSettings< > >(
                        angular_position, linkEnds2[ i ], initialObservationTimes, receiver ) );
    }



    // Simulate observations
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< ObservationScalarType, TimeType >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Set weights
    std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightPerObservationParser;
    weightPerObservationParser[ observationParser( one_way_range ) ] = 1.0E-4;
    weightPerObservationParser[ observationParser( angular_position ) ] = 1.0E-20;
    observationsAndTimes->setConstantWeightPerObservable( weightPerObservationParser );

    // Perturb parameter vector
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    initialParameterEstimate[ 0 ] += 1.0E2;
    initialParameterEstimate[ 1 ] += 1.0E2;
    initialParameterEstimate[ 2 ] += 1.0E2;
    initialParameterEstimate[ 3 ] += 1E-3;
    initialParameterEstimate[ 4 ] += 1E-3;
    initialParameterEstimate[ 5 ] += 1E-3;

    for( unsigned int i = 1; i < multiArcBodiesToIntegrate.size( ) * integrationArcStarts.size( ); i++ )
    {
        initialParameterEstimate[ 0 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 1 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 2 + 6 * i ] += 100.0E0;
        initialParameterEstimate[ 3 + 6 * i ] += 1E-3;
        initialParameterEstimate[ 4 + 6 * i ] += 1E-3;
        initialParameterEstimate[ 5 + 6 * i ] += 1E-3;
    }

    for( unsigned int i = 6 * ( 1 + multiArcBodiesToIntegrate.size( ) * integrationArcStarts.size( ) );
         i < static_cast< unsigned int >( initialParameterEstimate.rows( ) ); i++ )
    {
        initialParameterEstimate[ i ] *= ( 1.0 + 1.0E-6 );
    }

    parametersToEstimate->resetParameterValues( initialParameterEstimate );

    // Define estimation settings
    std::shared_ptr< EstimationInput< ObservationScalarType, TimeType > > estimationInput =
            std::make_shared< EstimationInput< ObservationScalarType, TimeType > >(
                observationsAndTimes );

    estimationInput->setConvergenceChecker(
                std::make_shared< EstimationConvergenceChecker >( 6 ) );
    // Estimate parameters and return postfit error
    std::shared_ptr< EstimationOutput< StateScalarType > > estimationOutput = orbitDeterminationManager.estimateParameters(
                estimationInput  );

//    input_output::writeMatrixToFile( estimationOutput->normalizedDesignMatrix_, "hybridArcPartials.dat" );
//    input_output::writeMatrixToFile( estimationOutput->designMatrixTransformationDiagonal_, "hybridArcNormalization.dat" );
//    input_output::writeMatrixToFile( estimationOutput->inverseNormalizedCovarianceMatrix_, "hybridArcNormalizedInverseCovariance.dat" );
//    input_output::writeMatrixToFile( estimationOutput->getResidualHistoryMatrix( ), "hybridArcResidualHistory.dat" );

    return ( estimationOutput->parameterEstimate_ - truthParameters ).template cast< double >( );
}


BOOST_AUTO_TEST_CASE( test_HybridArcStateEstimation )
{
    // Perform estimation
    Eigen::VectorXd parameterError = executeParameterEstimation< double, double, double >( );
    int numberOfEstimatedArcs = ( parameterError.rows( ) - 8 ) / 6;

    std::cout<<std::endl<<std::endl<<"Final error: "<<parameterError.transpose( )<<std::endl;
    // Test error range: 5 m in-plane position and 2 micron/s in-plane velocity for Mars
    for( unsigned int j = 0; j < 2; j++ )
    {
        BOOST_CHECK_SMALL( std::fabs( parameterError( j ) ), 5.0 );
        BOOST_CHECK_SMALL( std::fabs( parameterError( j + 3 ) ), 2.5E-6  );
    }

    // Test error range: 1000 m in-plane position and 1.0 mm/s in-plane velocity for Mars (poor values due to short arc)
    BOOST_CHECK_SMALL( std::fabs( parameterError( 2 ) ), 1000.0 );
    BOOST_CHECK_SMALL( std::fabs( parameterError( 5 ) ), 1.0E-3  );

    // Test error range: 0.1 m position and 50 micron/s velocity for orbiter
    for( int i = 0; i < numberOfEstimatedArcs; i++ )
    {
        for( unsigned int j = 0; j < 3; j++ )
        {
            BOOST_CHECK_SMALL( std::fabs( parameterError( ( i + 1 )* 6 + j ) ), 1E-1 );
            BOOST_CHECK_SMALL( std::fabs( parameterError( ( i + 1 ) * 6 + j + 3 ) ), 5.0E-5  );
        }
    }

    // Test errors for gravitational parameters
    BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 2 ) ), 1.5E11 );
    BOOST_CHECK_SMALL( std::fabs( parameterError( parameterError.rows( ) - 1 ) ), 1.0E6 );

}

BOOST_AUTO_TEST_SUITE_END( )

}

}
