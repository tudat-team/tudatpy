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

    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariables;
    dependentVariables.push_back( std::make_shared< SingleDependentVariableSaveSettings >(
        total_acceleration_dependent_variable, "Vehicle" ) );

    // Create propagator settings
    std::shared_ptr< TranslationalStatePropagatorSettings< double, double > > propagatorSettings =
        std::make_shared< TranslationalStatePropagatorSettings< double, double > >
        ( centralBodies, accelerationModelMap, bodiesToIntegrate, initialState, finalTime, cowell, dependentVariables );

    // Create integrator settings
    std::shared_ptr< IntegratorSettings< double > > integratorSettings = std::make_shared< RungeKuttaVariableStepSizeSettings< double > >(
        initialTime, 120.0, CoefficientSets::rungeKuttaFehlberg78, 120.0, 120.0, 1.0, 1.0 );

//    // Define link ends.
//    LinkEnds stationTransmitterLinkEnds;
//    stationTransmitterLinkEnds[ transmitter ] = LinkEndId( "Earth", "Station" );
//    stationTransmitterLinkEnds[ receiver ] = LinkEndId( "Vehicle", "" );

    LinkEnds stationReceiverLinkEnds;
    stationReceiverLinkEnds[ receiver ] = LinkEndId( "Earth", "Station" );
    stationReceiverLinkEnds[ transmitter ] = LinkEndId( "Vehicle", "" );

    std::map< ObservableType, std::vector< LinkEnds > > linkEndsPerObservable;
    linkEndsPerObservable[ one_way_range ] = std::vector< LinkEnds >( { stationReceiverLinkEnds } );//, stationTransmitterLinkEnds } );

    for ( unsigned int testCase = 0 ; testCase < 1 ; testCase++ )
    {
        bool multiArcBiases = testCase;

        std::vector< std::shared_ptr< EstimatableParameterSettings > > parameterNames;
        parameterNames.push_back( std::make_shared< InitialTranslationalStateEstimatableParameterSettings< double > >( "Vehicle", initialState, "Earth" ) );
        if ( !multiArcBiases )
        {
            parameterNames.push_back( std::make_shared< ConstantTimeBiasEstimatableParameterSettings >(
                linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, receiver ) );
//            parameterNames.push_back( std::make_shared< ConstantTimeBiasEstimatableParameterSettings >(
//                linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, receiver ) );
//            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
//                linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, true ) );
//            parameterNames.push_back( std::make_shared< ConstantObservationBiasEstimatableParameterSettings >(
//                linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, true ) );
        }
        else
        {
            parameterNames.push_back( std::make_shared< ArcWiseTimeBiasEstimatableParameterSettings >(
                linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, arcs, receiver ) );
            parameterNames.push_back( std::make_shared< ArcWiseTimeBiasEstimatableParameterSettings >(
                linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, arcs, receiver ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                linkEndsPerObservable.at( one_way_range ).at( 0 ), one_way_range, arcs, receiver, true ) );
            parameterNames.push_back( std::make_shared< ArcWiseConstantObservationBiasEstimatableParameterSettings >(
                linkEndsPerObservable.at( one_way_range ).at( 1 ), one_way_range, arcs, receiver, true ) );
        }

        // Create parameters
        std::shared_ptr< estimatable_parameters::EstimatableParameterSet< double > > parametersToEstimate =
            createParametersToEstimate< double, double >( parameterNames, bodies );
        printEstimatableParameterEntries( parametersToEstimate );

        std::vector< std::shared_ptr< ObservationModelSettings > > observationSettingsList;
        for( auto linkEndIterator : linkEndsPerObservable )
        {
            ObservableType currentObservable = linkEndIterator.first;

            std::vector< LinkEnds > currentLinkEndsList = linkEndIterator.second;
            for( unsigned int i = 0; i < currentLinkEndsList.size( ); i++ )
            {
                std::shared_ptr< ObservationBiasSettings > biasSettings;
                if ( !multiArcBiases )
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                    biasSettingsList.push_back( std::make_shared< ConstantTimeBiasSettings >( 0.0, receiver ) );
                    biasSettingsList.push_back( std::make_shared< ConstantObservationBiasSettings >( Eigen::Vector1d::Zero( ), true ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
                }
                else
                {
                    std::vector< std::shared_ptr< ObservationBiasSettings > > biasSettingsList;
                    biasSettingsList.push_back( std::make_shared< ArcWiseTimeBiasSettings >( arcs, timeBiasesPerArc, receiver ) );
                    biasSettingsList.push_back( std::make_shared< ArcWiseConstantObservationBiasSettings >( arcs, biasesPerArc, receiver, true ) );
                    biasSettings = std::make_shared< MultipleObservationBiasSettings >( biasSettingsList );
                }
                observationSettingsList.push_back( std::make_shared< ObservationModelSettings >(
                    currentObservable, currentLinkEndsList.at( i ), std::shared_ptr< LightTimeCorrectionSettings >( ), biasSettings ) );
            }
        }

        // Create orbit determination object.
        OrbitDeterminationManager< double, double > orbitDeterminationManager = OrbitDeterminationManager< double, double >(
            bodies, parametersToEstimate, observationSettingsList, integratorSettings, propagatorSettings );

        std::vector< double > baseTimeList;
        double observationInterval = 600.0;
        for( unsigned int j = 0; j < 10; j++ )
        {
            baseTimeList.push_back( ( initialTime + 600.0 ) + static_cast< double >( j ) * observationInterval );
        }
        std::vector< std::shared_ptr< ObservationSimulationSettings< double > > > measurementSimulationInput;
        LinkEnds currentLinkEnds = linkEndsPerObservable.at( one_way_range ).at( 0 );
                measurementSimulationInput.push_back( std::make_shared< TabulatedObservationSimulationSettings< double > >(
                    one_way_range, currentLinkEnds, baseTimeList, receiver ) );

        // Simulate observations
        std::shared_ptr< ObservationCollection< double, double > > simulatedObservations = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );



        std::shared_ptr< ObservationModel< 1, double, double > > observationModel =
            std::dynamic_pointer_cast< ObservationSimulator< 1, double, double > >( orbitDeterminationManager.getObservationSimulators( ).at( 0 ) )->getObservationModels( ).at(
            linkEndsPerObservable.at( one_way_range ).at( 0 ) );
        std::shared_ptr< EstimatableParameter< Eigen::VectorXd > > timeBiasParameter = parametersToEstimate->getEstimatedVectorParameters( ).at( 0 );
        auto observationManager = std::dynamic_pointer_cast< ObservationManager< 1, double, double > >( orbitDeterminationManager.getObservationManager( one_way_range ) );
        Eigen::MatrixXd partials = Eigen::MatrixXd::Zero( baseTimeList.size( ), 7 );
        Eigen::VectorXd observations = simulatedObservations->getObservationVector( );
        observationManager->computeObservationsWithPartials(
            baseTimeList, currentLinkEnds, receiver, nullptr, observations, partials, false, true );
        std::cout<<partials<<std::endl;

        double timeBiasPerturbation = 0.2;
        Eigen::VectorXd originalParameters = parametersToEstimate->getFullParameterValues< double >( );

        Eigen::VectorXd perturbedParameters = originalParameters;
        perturbedParameters( 6 ) += timeBiasPerturbation;
        parametersToEstimate->resetParameterValues( perturbedParameters );
        std::shared_ptr< ObservationCollection< double, double > > upperturbedSimulatedObservations = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

        perturbedParameters = originalParameters;
        perturbedParameters( 6 ) -= timeBiasPerturbation;
        parametersToEstimate->resetParameterValues( perturbedParameters );
        std::shared_ptr< ObservationCollection< double, double > > downperturbedSimulatedObservations = simulateObservations< double, double >(
            measurementSimulationInput, orbitDeterminationManager.getObservationSimulators( ), bodies );

        Eigen::MatrixXd numericalPartials = ( upperturbedSimulatedObservations->getObservationVector( ) - downperturbedSimulatedObservations->getObservationVector( ) ) /
            ( 2.0 * timeBiasPerturbation );
        std::cout<<std::endl<<std::endl<<numericalPartials<<std::endl;
        std::cout<<std::endl<<std::endl<<numericalPartials.cwiseQuotient( partials.block( 0, 6, 10, 1 ) )-Eigen::VectorXd::Constant( numericalPartials.rows( ), 1 )<<std::endl;

    }

}


BOOST_AUTO_TEST_SUITE_END( )

} // namespace unit_tests

} // namespace tudat





