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


namespace tudat
{
namespace unit_tests
{


BOOST_AUTO_TEST_SUITE( test_time_bias_partials )

//! Test partial derivatives of angular position observable, using general test suite of observation partials.
BOOST_AUTO_TEST_CASE( testTimeBiasPartials )
{
    const int numberOfDaysOfData = 1;
    int numberOfIterations = 10;

    //Load spice kernels.
    spice_interface::loadStandardSpiceKernels( );

    // Define bodies in simulation
    std::vector< std::string > bodyNames;
    bodyNames.push_back( "Earth" );

    // Specify initial time
    double initialTime = 1.0e7;
    double finalTime = initialTime + numberOfDaysOfData * 86400.0;

    std::vector< double > arcs;
    arcs.push_back( initialTime );
    arcs.push_back( initialTime + 8.0 * 3600.0 );
    arcs.push_back( initialTime + 16.0 * 3600.0 );

    std::vector< double > timeBiasesPerArc = { 0.0, 0.0, 0.0 };

    std::vector< Eigen::VectorXd > biasesPerArc;
    biasesPerArc.push_back( Eigen::Vector1d::Zero( ) );
    biasesPerArc.push_back( Eigen::Vector1d::Zero( ) );
    biasesPerArc.push_back( Eigen::Vector1d::Zero( ) );

    // Create bodies needed in simulation
    BodyListSettings bodySettings = getDefaultBodySettings( bodyNames, "Earth", "ECLIPJ2000" );
    SystemOfBodies bodies = createSystemOfBodies( bodySettings );
    bodies.createEmptyBody( "Vehicle" );
    bodies.at( "Vehicle" )->setEphemeris( std::make_shared< TabulatedCartesianEphemeris< > >(
        std::shared_ptr< interpolators::OneDimensionalInterpolator < double, Eigen::Vector6d > >( ), "Earth", "ECLIPJ2000" ) );


    // Create ground station
    createGroundStation( bodies.at( "Earth" ), "Station", ( Eigen::Vector3d( ) << 0.0, 0.35, 0.0 ).finished( ), geodetic_position );

    // Set accelerations on Vehicle that are to be taken into account.
    SelectedAccelerationMap accelerationMap;
    std::map< std::string, std::vector< std::shared_ptr< AccelerationSettings > > > accelerationsOfVehicle;
    accelerationsOfVehicle[ "Earth" ].push_back( std::make_shared< AccelerationSettings >( point_mass_gravity ) );
    accelerationMap[ "Vehicle" ] = accelerationsOfVehicle;

    // Set bodies for which initial state is to be estimated and integrated.
    std::vector< std::string > bodiesToIntegrate;
    std::vector< std::string > centralBodies;
    bodiesToIntegrate.push_back( "Vehicle" );
    centralBodies.push_back( "Earth" );

    // Create acceleration models
    AccelerationMap accelerationModelMap = createAccelerationModelsMap( bodies, accelerationMap, bodiesToIntegrate, centralBodies );

    // Set Keplerian elements
    Eigen::Vector6d initialKeplerianState;
    initialKeplerianState( semiMajorAxisIndex ) = 7200.0E3;
    initialKeplerianState( eccentricityIndex ) = 0.05;
    initialKeplerianState( inclinationIndex ) = unit_conversions::convertDegreesToRadians( 85.3 );
    initialKeplerianState( argumentOfPeriapsisIndex ) = unit_conversions::convertDegreesToRadians( 235.7 );
    initialKeplerianState( longitudeOfAscendingNodeIndex ) = unit_conversions::convertDegreesToRadians( 23.4 );
    initialKeplerianState( trueAnomalyIndex ) = unit_conversions::convertDegreesToRadians( 139.87 );
    double earthGravitationalParameter = bodies.at( "Earth" )->getGravityFieldModel( )->getGravitationalParameter( );
    Eigen::Vector6d initialState = convertKeplerianToCartesianElements( initialKeplerianState, earthGravitationalParameter );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double, double > >
        ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialState, finalTime  );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
        initialTime, 120.0, CoefficientSets::rungeKuttaFehlberg78, 120.0, 120.0, 1.0, 1.0 );

    // Define link ends.
    LinkEnds testLinkEnds;
    testLinkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
    testLinkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );

    for( unsigned int observableCase = 1; observableCase < 2; observableCase++ )
    {
        ObservableType currentObservable = undefined_observation_model;
        switch( observableCase )
        {
        case 0:
            currentObservable = one_way_range;
            break;
        case 1:
            currentObservable = angular_position;
            break;
        default:
            throw std::runtime_error( "Error when testing time bias partials; observable not implemented" );
        }
        for ( unsigned int testCase = 0 ; testCase < 1 ; testCase++ )
        {
            std::cout<<observableCase<<" "<<testCase<<std::endl;
            bool multiArcBiases = testCase;

            // Add bias parameters to estimation
            std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
            parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", initialState, "Earth" ) );
            if ( !multiArcBiases )
            {
                parameterNames.push_back( std::make_shared< ConstantTimeBiasEstimatableParameterSettings >(
                    testLinkEnds, currentObservable, receiver ) );
            }
            else
            {
                parameterNames.push_back( std::make_shared< ArcWiseTimeBiasEstimatableParameterSettings >(
                    testLinkEnds, currentObservable, arcs, receiver ) );
            }

            // Create parameters
            std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
                createParametersToEstimate< double, double >( parameterNames, bodies );
            printEstimatableParameterEntries( parametersToEstimate );

            // Define bias settings for observation model
            std::shared_ptr< ObservationBiasSettings > biasSettings;
            if ( !multiArcBiases )
            {
                std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                biasSettingsList.push_back( std::make_shared< ConstantTimeBiasSettings >( 0.0, receiver ) );
                biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
            }
            else
            {
                std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                biasSettingsList.push_back( std::make_shared< ArcWiseTimeBiasSettings >( arcs, timeBiasesPerArc, receiver ) );
                biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
            }
            // Define observation model settings
            std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
            observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                currentObservable, testLinkEnds, std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ) );

            // Create orbit determination object.
            OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
                bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

            // Define list of observation times
            std::vector< double > baseTimeList;
            double observationInterval = 600.0;
            int numberOfObservations = 10;
            for( unsigned int j = 0; j < numberOfObservations; j++ )
            {
                baseTimeList.push_back( ( initialTime + 600.0 ) + static_cast< double >( j ) * observationInterval );
            }

            // Define observation simulation times
            std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
            measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                currentObservable, testLinkEnds, baseTimeList, receiver ) );

            // Simulate observations
            std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
                measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );
            Eigen::VectorXd observations = simulatedObservations->getObservationVector( );

            // Extract observation model
            std::cout<<"Number of simulators "<<orbitDeterminationManager.getObservationSimulators( ).size( )<<std::endl;
            std::shared_ptr< ObservationModel< 2, double, double > > observationModel =
                std::dynamic_pointer_cast< ObservationSimulator< 2, double, double > >( orbitDeterminationManager.getObservationSimulators( ).at( 0 ) )->getObservationModels( ).at(
                testLinkEnds );

            // Compute analytical partials
            auto observationManager = std::dynamic_pointer_cast< ObservationManager< 2, double, double > >( orbitDeterminationManager.getObservationManager( currentObservable ) );
            Eigen::MatrixXd partials = Eigen::MatrixXd::Zero( baseTimeList.size( ), 7 );
            observationManager->computeObservationsWithPartials(
                baseTimeList, testLinkEnds, receiver, nullptr, observations, partials, false, true );

            // Get nominal parameter values
            Eigen::VectorXd originalParameters = parametersToEstimate->getFullParameterValues< double >( );

            for( unsigned int parameterIndex = 6; parameterIndex < parametersToEstimate->getParameterSetSize( ); parameterIndex++ )
            {
                // Define perturbation for numerical partial
                double timeBiasPerturbation = 1.0;

                // Perturb bias value up and recompute observations
                Eigen::VectorXd perturbedParameters = originalParameters;
                perturbedParameters( parameterIndex ) += timeBiasPerturbation;
                parametersToEstimate->resetParameterValues( perturbedParameters );
                std::shared_ptr< ObservationCollection< double, double > > upperturbedSimulatedObservations = simulateObservations< double, double >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

                // Perturb bias value down and recompute observations
                perturbedParameters = originalParameters;
                perturbedParameters( parameterIndex ) -= timeBiasPerturbation;
                parametersToEstimate->resetParameterValues( perturbedParameters );
                std::shared_ptr< ObservationCollection< double, double > > downperturbedSimulatedObservations = simulateObservations< double, double >(
                    measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

                // Reset to nominal values
                parametersToEstimate->resetParameterValues( perturbedParameters );

                // Compute numerical partials
                Eigen::MatrixXd numericalPartials = ( upperturbedSimulatedObservations->getObservationVector( ) - downperturbedSimulatedObservations->getObservationVector( ) ) /
                    ( 2.0 * timeBiasPerturbation );

                std::cout<<numericalPartials<<std::endl<<std::endl<<partials<<std::endl<<std::endl;
                std::cout<<std::endl<<std::endl<<numericalPartials.cwiseQuotient( partials.block( 0, parameterIndex, partials.rows( ), 1 ) )-Eigen::VectorXd::Constant( numericalPartials.rows( ), 1 )<<std::endl;
            }
        }
    }
}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





