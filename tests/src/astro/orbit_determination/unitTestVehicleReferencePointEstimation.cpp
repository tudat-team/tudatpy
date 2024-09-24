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

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"

#include "tudat/simulation/estimation_setup/orbitDeterminationTestCases.h"
#include "tudat/simulation/estimation_setup/podProcessing.h"

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


namespace tudat
{
namespace unit_tests
{
BOOST_AUTO_TEST_SUITE( test_reference_point_estimation )

BOOST_AUTO_TEST_CASE( test_ReferencePointEstimation )
{

    // Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );
    bodyNames.push_back( "Sun" );
    bodyNames.push_back( "Moon" );
    bodyNames.push_back( "Mars" );

    // Specify number of observation days.
    int numberOfDaysOfData = 1;

    // Specify initial time
    double initialEphemerisTime = 1.0e7;
    double finalEphemerisTime = initialEphemerisTime + numberOfDaysOfData * 86400.0;

    // Create bodies needed in simulation
    BodyListSettings bodySettings =
        getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    bodySettings.addSettings( "Vehicle" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );

    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
        std::shared_ptr< interpolators::OneDimensionalInterpolator
            < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );
    bodies.at( "Vehicle" )->setRotationalEphemeris(
        createRotationModel( spiceRotationModelSettings( "ECLIPJ2000", "VehicleFixed", "IAU_Earth"), "Vehicle", bodies ) );
    bodies.at( "Vehicle" )->getVehicleSystems( )->setReferencePointPosition( "Antenna", Eigen::Vector3d::UnitX( ) );

    // Create ground stations.
    std::vector< std::string > groundStationNames;
    groundStationNames.push_back( "Station1" );
    groundStationNames.push_back( "Station2" );
    groundStationNames.push_back( "Station3" );

    createGroundStation( bodies.at( "Earth" ), "Station1", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station2", ( Eigen::Vector3d( ) << 0.0, -0.55, 2.0 ).finished( ),
                         coordinate_conversions::geodetic_position );
    createGroundStation( bodies.at( "Earth" ), "Station3", ( Eigen::Vector3d( ) << 0.0, 0.05, 4.0 ).finished( ),
                         coordinate_conversions::geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >(
        basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Sun" ].push_back( std::make_shared< AccelerationSettings >(
        basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Moon" ].push_back( std::make_shared< AccelerationSettings >(
        basic_astrodynamics::point_mass_gravity ) );
    accelerationsOfVehicle[ "Mars" ].push_back( std::make_shared< AccelerationSettings >(
        basic_astrodynamics::point_mass_gravity ) );

    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;


    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap(
        bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements for Vehicle.
    Eigen::Vector6d initialStateInKeplerianElements;
    initialStateInKeplerianElements( semiMajorAxisIndex ) = 7200.0E3;
    initialStateInKeplerianElements( eccentricityIndex ) = 0.05;
    initialStateInKeplerianElements( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialStateInKeplerianElements( argumentOfPeriapsisIndex )
        = unit_conversions::convertDegreesToRadians( 235.7 );
    initialStateInKeplerianElements( longitudeOfAscendingNodeIndex )
        = unit_conversions::convertDegreesToRadians( 23.4 );
    initialStateInKeplerianElements( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );

    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );

    // Set (perturbed) initial state.
    Eigen::Matrix< double, 6, 1 > systemInitialState = convertKeplerianToCartesianElements(
        initialStateInKeplerianElements, earthGravitationalParameter );


    // Create integrator settings.
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = rungeKuttaFixedStepSettings(
        40.0, CoefficientSets::rungeKuttaFehlberg78 );

    // Create propagator settings.
    std::shared_ptr< TranslationalStatePropagatorSettings< double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double > >
            ( centralBodies, accelerationModelMap, bodiesToIntegrate, systemInitialState, initialEphemerisTime,
              integratorSettings, std::make_shared< PropagationTimeTerminationSettings >( finalEphemerisTime ) );

    // Define link ends.
    std::vector< LinkDefinition > stationReceiverLinkEnds;
    std::vector< LinkDefinition > stationTransmitterLinkEnds;

    for( unsigned int i = 0; i < groundStationNames.size( ); i++ )
    {
        LinkDefinition linkEnds;
        linkEnds[ transmitter ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ receiver ] = std::make_pair< std::string, std::string >( "Vehicle", "Antenna" );
        stationTransmitterLinkEnds.push_back( linkEnds );

        linkEnds[ receiver ] = std::pair< std::string, std::string >( std::make_pair( "Earth", groundStationNames.at( i ) ) );
        linkEnds[ transmitter ] = std::make_pair< std::string, std::string >( "Vehicle", "Antenna" );
        stationReceiverLinkEnds.push_back( linkEnds );
    }

    std::map< ObservableType, std::vector< LinkDefinition > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationTransmitterLinkEnds[ 0 ] );
    linkEndsPerObservable[ one_way_range ].push_back( stationReceiverLinkEnds[ 1 ] );

    // Define parameters to be estimated.
    std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames =
        getInitialStateParameterSettings< double >( propagatorSettings, bodies );
    parameterNames.push_back( std::make_shared< EstimatableParameterSettings >( "Vehicle", reference_point_position, "Antenna" ) );


    // Create parameters
    std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
        createParametersToEstimate< double >( parameterNames, bodies , propagatorSettings );

    printEstimatableParameterEntries( parametersToEstimate );

    std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
    for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            observationSettingsList.push_back(
                std::make_shared< ObservationModelSettings >(
                    currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ) ) );
        }
    }

    // Create orbit determination object.
    OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
        bodies, parametersToEstimate, observationSettingsList, propagatorSettings );

    // Compute list of observation times.
    std::vector< double > baseTimeList;
    double observationTimeStart = initialEphemerisTime + 1000.0;
    double  observationInterval = 60.0;
    for( int i = 0; i < numberOfDaysOfData; i++ )
    {
        for( unsigned int j = 0; j < 500; j++ )
        {
            baseTimeList.push_back( observationTimeStart + static_cast< double >( i ) * 86400.0 +
                                    static_cast< double >( j ) * observationInterval );
        }
    }


    std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
    for( std::map< ObservableType, std::vector< LinkDefinition > >::iterator linkEndIterator = linkEndsPerObservable.begin( );
         linkEndIterator != linkEndsPerObservable.end( ); linkEndIterator++ )
    {
        ObservableType currentObservable = linkEndIterator->first;
        std::vector< LinkDefinition > currentLinkEndsList = linkEndIterator->second;
        for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
        {
            measurementSimulationInput.push_back(
                std::make_shared< TabulatedObservationSimulationSettings< > >(
                    currentObservable, currentLinkEndsList[ i ], baseTimeList, receiver ) );
        }
    }

    // Simulate observations.
    std::shared_ptr< ObservationCollection< > > observationsAndTimes = simulateObservations< double, double >(
        measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

    // Set observation weights
    std::map< std::shared_ptr< observation_models::ObservationCollectionParser >, double > weightsPerObservationParser;
    weightsPerObservationParser[ observationParser( one_way_range ) ] = 1.0 / ( 1.0 * 1.0 );
    weightsPerObservationParser[ observationParser( angular_position ) ] = 1.0 / ( 1.0E-5 * 1.0E-5 );
    weightsPerObservationParser[ observationParser( one_way_doppler ) ] = 1.0 / ( 1.0E-11 * 1.0E-11 * physical_constants::SPEED_OF_LIGHT * physical_constants::SPEED_OF_LIGHT  );
    observationsAndTimes->setConstantWeightPerObservable( weightsPerObservationParser );

    // Perturb parameter estimate.
    Eigen::Matrix< double, Eigen::Dynamic, 1 > initialParameterEstimate =
        parametersToEstimate->template getFullParameterValues< double >( );
    Eigen::Matrix< double, Eigen::Dynamic, 1 > truthParameters = initialParameterEstimate;
    Eigen::Matrix< double, Eigen::Dynamic, 1 > parameterPerturbation =
        Eigen::Matrix< double, Eigen::Dynamic, 1 >::Zero( truthParameters.rows( ) );

    // Perturbe initial state estimate.
    parameterPerturbation.segment( 0, 3 ) = Eigen::Vector3d::Constant( 1.0 );
    parameterPerturbation.segment( 3, 3 ) = Eigen::Vector3d::Constant( 1.E-3 );
    parameterPerturbation.segment( 6, 3 ) = Eigen::Vector3d::Constant( 10.0 );
    initialParameterEstimate += parameterPerturbation;
    parametersToEstimate->resetParameterValues( initialParameterEstimate );


    // Define estimation input
    std::shared_ptr< EstimationInput< double, double  > > estimationInput = std::make_shared< EstimationInput< double, double > >(
        observationsAndTimes );


    estimationInput->defineEstimationSettings( true, true, true, true, false );
    estimationInput->setConvergenceChecker(
        std::make_shared< EstimationConvergenceChecker >( 4 ) );

    // Perform estimation
    std::shared_ptr< EstimationOutput< double > > estimationOutput = orbitDeterminationManager.estimateParameters(
        estimationInput );

    Eigen::VectorXd estimationError = estimationOutput->parameterEstimate_ - truthParameters;
    std::cout <<"estimation error: "<< ( estimationError ).transpose( ) << std::endl;


    // Check if parameters are correctly estimated
    Eigen::VectorXd estimatedParametervalues = estimationOutput->parameterEstimate_;

    for( unsigned int i = 0; i < 3; i++ )
    {
        BOOST_CHECK_SMALL( std::fabs( estimationError( i ) ), 1.0E-6 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 3 ) ), 1.0E-9 );
        BOOST_CHECK_SMALL( std::fabs( estimationError( i + 6 ) ), 1.0E-6 );
    }

}

BOOST_AUTO_TEST_SUITE_END( )

}

}



